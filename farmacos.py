import networkx as nx
import plotly.graph_objects as go
import plotly.express as px
from collections import Counter
import numpy as np
import pandas as pd
import streamlit as st

# Configuración de Streamlit
st.set_page_config(layout="wide")

# Cargar datos
@st.cache_data
def load_data():
    return pd.read_csv(r"DDIBUENO.csv")


df = load_data()

@st.cache_data
def load_data2():
    return pd.read_csv(r"drugbank_ddi_label_map.csv")
meeaning_y=load_data2()

@st.cache_data
def load_data():
    essentials_df = pd.read_excel(r"GlobalEssentialMedicinesDatabase.xlsx")  # tu nuevo archivo
    return essentials_df

# DICCIONARIO CON NOMBRES COMPLETOS DE ATC
ATC_CATEGORIES = {
    'A': 'ALIMENTARY TRACT AND METABOLISM',
    'B': 'BLOOD AND BLOOD FORMING ORGANS',
    'C': 'CARDIOVASCULAR SYSTEM',
    'D': 'DERMATOLOGICALS',
    'G': 'GENITO URINARY SYSTEM AND SEX HORMONES',
    'H': 'SYSTEMIC HORMONAL PREPARATIONS, EXCL. SEX HORMONES AND INSULINS',
    'J': 'ANTIINFECTIVES FOR SYSTEMIC USE',
    'L': 'ANTINEOPLASTIC AND IMMUNOMODULATING AGENTS',
    'M': 'MUSCULO-SKELETAL SYSTEM',
    'N': 'NERVOUS SYSTEM',
    'P': 'ANTIPARASITIC PRODUCTS, INSECTICIDES AND REPELLENTS',
    'R': 'RESPIRATORY SYSTEM',
    'S': 'SENSORY ORGANS',
    'V': 'VARIOUS',
    'Sin ATC': 'No ATC',
    'Multi ATC': 'Multi ATC'
}

# Colores para cada categoría ATC
ATC_COLORS = {
    'A': '#FF6B35',
    'B': '#004E89',
    'C': '#FF0000',
    'D': '#FFA500',
    'G': '#9370DB',
    'H': '#FF69B4',
    'J': '#32CD32',
    'L': '#8B0000',
    'M': '#D2691E',
    'N': '#4B0082',
    'P': '#00CED1',
    'R': '#1E90FF',
    'S': '#FFD700',
    'V': '#808080',
    'Sin ATC': '#CCCCCC',
    'Multi ATC': '#800080'
}

def get_atc_color_and_category(atc_code, num_atc):
    if num_atc == 0:
        return ATC_COLORS['Sin ATC'], 'Sin ATC'
    
    if num_atc > 1:
        if pd.notna(atc_code) and '|' in str(atc_code):
            return ATC_COLORS['Multi ATC'], 'Multi ATC'
    
    if pd.isna(atc_code) or atc_code == '' or atc_code == 'No ATC':
        return ATC_COLORS['Sin ATC'], 'Sin ATC'
    
    first_char = str(atc_code)[0].upper()
    
    if '|' in str(atc_code):
        return ATC_COLORS['Multi ATC'], 'Multi ATC'
    
    return ATC_COLORS.get(first_char, '#CCCCCC'), ATC_CATEGORIES.get(first_char, 'Unknown')

# FUNCIÓN PARA OBTENER LA LETRA ATC DEL NOMBRE DE CATEGORÍA
def get_atc_letter_from_category(category_name):
    """Obtener la letra ATC a partir del nombre completo de la categoría"""
    for letter, name in ATC_CATEGORIES.items():
        if name == category_name:
            return letter
    if category_name == 'Sin ATC':
        return 'Sin ATC'
    elif category_name == 'Multi ATC':
        return 'Multi ATC'
    return None

@st.cache_data
def crear_grafo_completo(df):
    """Crear grafo completo con todos los fármacos"""
    G = nx.DiGraph()
    
    for _, row in df.iterrows():
        drug1 = row['Common_name_x']
        drug2 = row['Common_name_y']
        interaction_type = row['Y'] 
        
        atc1 = row['atc_code_x'] if pd.notna(row['atc_code_x']) else "No ATC"
        num_atc1 = row['num_atc_x']
        color1, category1 = get_atc_color_and_category(atc1, num_atc1)
        
        atc2 = row['atc_code_y'] if pd.notna(row['atc_code_y']) else "No ATC"
        num_atc2 = row['num_atc_y']
        color2, category2 = get_atc_color_and_category(atc2, num_atc2)
        
        if not G.has_node(drug1):
            tooltip_info = f"<b>Drug:</b> {drug1}<br><b>ATC Code:</b> {atc1}<br><b>ATC Category:</b> {category1}"
            G.add_node(drug1, 
                      atc_code=atc1,
                      atc_category=category1,
                      color=color1,
                      num_atc=num_atc1,
                      tooltip=tooltip_info)
        
        if not G.has_node(drug2):
            tooltip_info = f"<b>Drug:</b> {drug2}<br><b>ATC Code:</b> {atc2}<br><b>ATC Category:</b> {category2}"
            G.add_node(drug2, 
                      atc_code=atc2,
                      atc_category=category2,
                      color=color2,
                      num_atc=num_atc2,
                      tooltip=tooltip_info)
        
        if not G.has_edge(drug1, drug2):
            edge_tooltip = f"<b>From:</b> {drug1} → <b>To:</b> {drug2}"
            G.add_edge(drug1, drug2, 
                      interaction_type=interaction_type,
                      tooltip=edge_tooltip)
    
    return G

def crear_subgrafo_centrado(G, farmaco_objetivo):
    """Crear subgrafo centrado en un fármaco específico"""
    farmaco_encontrado = None
    for node in G.nodes():
        if farmaco_objetivo.lower() in node.lower():
            farmaco_encontrado = node
            break
    
    if farmaco_encontrado:
        predecessors = list(G.predecessors(farmaco_encontrado))
        successors = list(G.successors(farmaco_encontrado))
        subgraph_nodes = [farmaco_encontrado] + predecessors + successors
        G_sub = G.subgraph(subgraph_nodes).copy()
        return G_sub, farmaco_encontrado
    else:
        return None, None

def crear_grafo_plotly(G, farmaco_principal=None, active_categories=None):
    """Crear visualización del grafo con Plotly"""
    
    if active_categories is None or not any(active_categories.values()):
        fig = go.Figure()
        fig.update_layout(
            title=dict(
                text="No drugs to display<br>Select categories on the right",
                font=dict(size=16),
                x=0.5,
                y=0.5,
                xanchor='center',
                yanchor='middle'
            ),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor='black',
            height=600,
            width=800
        )
        return fig, None
    
    nodes_to_keep = []
    
    if farmaco_principal:
        nodes_to_keep.append(farmaco_principal)
    
    for node in G.nodes():
        if node == farmaco_principal:
            continue
        category = G.nodes[node]['atc_category']
        if active_categories.get(category, False):
            nodes_to_keep.append(node)
    
    if not nodes_to_keep:
        fig = go.Figure()
        fig.update_layout(
            title=dict(
                text="No drugs to display<br>Select categories on the right",
                font=dict(size=16),
                x=0.5,
                y=0.5,
                xanchor='center',
                yanchor='middle'
            ),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor='black',
            height=600,
            width=800
        )
        return fig, None
    
    edges_to_keep = []
    for u, v in G.edges():
        if u in nodes_to_keep and v in nodes_to_keep:
            edges_to_keep.append((u, v))
    
    G_filtered = nx.DiGraph()
    for node in nodes_to_keep:
        G_filtered.add_node(node, **G.nodes[node])
    
    for u, v in edges_to_keep:
        G_filtered.add_edge(u, v, **G[u][v])
    
    if farmaco_principal and farmaco_principal in G_filtered.nodes():
        pos = {}
        pos[farmaco_principal] = np.array([0, 0])
        
        neighbors = list(G_filtered.predecessors(farmaco_principal)) + \
                   list(G_filtered.successors(farmaco_principal))
        neighbors = list(set(neighbors))
        
        n_neighbors = len(neighbors)
        if n_neighbors > 0:
            radius = 1.5 + 0.2 * min(n_neighbors, 20)
            angle_step = 2 * np.pi / n_neighbors
            
            for i, neighbor in enumerate(neighbors):
                angle = i * angle_step
                x = radius * np.cos(angle)
                y = radius * np.sin(angle)
                pos[neighbor] = np.array([x, y])
    else:
        pos = nx.spring_layout(G_filtered, k=2/np.sqrt(len(G_filtered.nodes())), 
                              iterations=50, seed=42)
    
    node_x, node_y, node_colors, node_sizes, node_texts, node_names = [], [], [], [], [], []
    
    for node in G_filtered.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        node_colors.append(G_filtered.nodes[node]['color'])
        
        if farmaco_principal and node == farmaco_principal:
            node_sizes.append(30)
        else:
            node_sizes.append(15)
        
        node_texts.append(G_filtered.nodes[node]['tooltip'])
        
        if farmaco_principal and node == farmaco_principal:
            node_names.append(node)
        else:
            if len(node) > 20:
                node_names.append(node[:17] + "...")
            else:
                node_names.append(node)
    
    edge_x, edge_y, edge_texts = [], [], []
    
    for u, v, data in G_filtered.edges(data=True):
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)
        
        edge_texts.append(data['tooltip'])
    
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=1.5, color='gray'),
        hoverinfo='text',
        text=edge_texts,
        mode='lines',
        name='Interactions',
        hoverlabel=dict(bgcolor='black', font_size=12)
    )
    
    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text',
        text=node_names,
        textposition="top center",
        textfont=dict(size=10),
        hoverinfo='text',
        hovertext=node_texts,
        marker=dict(
            color=node_colors,
            size=node_sizes,
            line=dict(width=1, color='black')
        ),
        name='Drugs',
        hoverlabel=dict(bgcolor='black', font_size=12)
    )
    
    if farmaco_principal:
        title_text = f"Drug: {farmaco_principal}<br>Drugs: {len(G_filtered.nodes())} | Interactions: {len(G_filtered.edges())}"
    else:
        title_text = f"Complete Network<br>Drugs: {len(G_filtered.nodes())} | Interactions: {len(G_filtered.edges())}"
    
    fig = go.Figure(data=[edge_trace, node_trace],
                   layout=go.Layout(
                       title=dict(
                           text=title_text,
                           font=dict(size=14),
                           x=0.5,
                           xanchor='center'
                       ),
                       showlegend=False,
                       hovermode='closest',
                       margin=dict(b=20, l=5, r=5, t=60),
                       xaxis=dict(
                           showgrid=False, 
                           zeroline=False, 
                           showticklabels=False
                       ),
                       yaxis=dict(
                           showgrid=False, 
                           zeroline=False, 
                           showticklabels=False
                       ),
                       plot_bgcolor='black',
                       height=600,
                       width=800,
                       dragmode='pan'
                   ))
    
    return fig, G_filtered

def mostrar_analisis_interacciones(G_filtrado, farmaco_principal):
    """Mostrar análisis detallado de interacciones con valores de 'Y'"""
    
    if not G_filtrado or len(G_filtrado.edges()) == 0:
        st.info("No interactions to analyze.")
        return
    
    st.subheader("Interaction Analysis")
    
    # Preparar datos para análisis
    interacciones_data = []
    for u, v, data in G_filtrado.edges(data=True):
        y_value = data.get('interaction_type', 'N/A')
        interacciones_data.append({
            'From': u,
            'To': v,
            'Y Value': y_value,
            'Direction': f"{u} → {v}"
        })
    
    interacciones_df = pd.DataFrame(interacciones_data)
    
    # Mostrar estadísticas de valores Y
    st.write("**Y Value Distribution:**")
    y_counts = interacciones_df['Y Value'].value_counts().sort_index()
    
    col1, col2 = st.columns(2)
    
    with col1:
        for y_val, count in y_counts.items():
            percentage = (count / len(interacciones_df)) * 100
            st.write(f"**Y = {y_val}:** {count} interactions ({percentage:.1f}%)")
    
    with col2:
        if len(y_counts) > 0:
            fig_pie = px.pie(
                values=y_counts.values,
                names=[f"Y = {k}" for k in y_counts.index],
                title="Distribution of Y Values",
                color_discrete_sequence=px.colors.qualitative.Set3
            )
            st.plotly_chart(fig_pie, use_container_width=True)
    
    # Mostrar tabla de interacciones
    with st.expander("View All Interactions with Y value", expanded=False):
        # Ordenar por valor Y
        interacciones_df = interacciones_df.sort_values('Y Value')
        
        # Agregar colores según valor Y
        def color_y_value(val):
            
         return ''
        
        styled_df = interacciones_df.style.applymap(
            color_y_value, subset=['Y Value']
        )
        
        st.dataframe(
            styled_df,
            use_container_width=True,
            column_config={
                "From": st.column_config.TextColumn("From (Drug A)", width="medium"),
                "To": st.column_config.TextColumn("To (Drug B)", width="medium"),
                "Y Value": st.column_config.NumberColumn("Y Value", width="small"),
                "Direction": st.column_config.TextColumn("Direction", width="large")
            },
            hide_index=True
        )
    
    # Análisis por fármaco principal
    if farmaco_principal:
        st.write(f"**Interactions involving {farmaco_principal}:**")
        
        # Interacciones salientes (de farmaco_principal a otros)
        outgoing = []
        for _, row in interacciones_df.iterrows():
            if row['From'] == farmaco_principal:
                outgoing.append({
                    'Target': row['To'],
                    'Y Value': row['Y Value'],
                    'Type': 'Outgoing'
                })
        
        # Interacciones entrantes (de otros a farmaco_principal)
        incoming = []
        for _, row in interacciones_df.iterrows():
            if row['To'] == farmaco_principal:
                incoming.append({
                    'Source': row['From'],
                    'Y Value': row['Y Value'],
                    'Type': 'Incoming'
                })
        
        if outgoing:
            st.write(f"** {farmaco_principal} to drug 2 ({len(outgoing)}):**")
            outgoing_df = pd.DataFrame(outgoing)
            st.dataframe(
                outgoing_df.sort_values('Y Value'),
                use_container_width=True,
                hide_index=True
            )
        else:
            outgoing_df = None
            st.write("No outgoing interactions.")
        
        if incoming:
            st.write(f"**Drug 1 to {farmaco_principal} ({len(incoming)}):**")
            incoming_df = pd.DataFrame(incoming)
            st.dataframe(
                incoming_df.sort_values('Y Value'),
                use_container_width=True,
                hide_index=True
            )
        else:
            incoming_df = None
            st.write("No incoming interactions.")
        
        

def pestaña_visualizacion():
    """Pestaña principal de visualización de interacciones"""
    st.title("Drug Interaction Network Visualization")
    
    # Inicializar estado de sesión
    if 'G_completo' not in st.session_state:
        with st.spinner("Loading drug interaction data..."):
            st.session_state.G_completo = crear_grafo_completo(df)
    
    G_completo = st.session_state.G_completo
    
    # Sidebar para controles
    with st.sidebar:
        st.header("Controls")
        
        # Opción para seleccionar fármaco principal
        st.subheader("Select Main Drug (Optional)")
        all_drugs = sorted(set(df['Common_name_x'].tolist() + df['Common_name_y'].tolist()))
        
        farmaco_principal_seleccionado = st.selectbox(
            "Choose a drug to focus on:",
            ["None"] + all_drugs,
            index=0,
            help="Select 'None' to see the complete network without focus"
        )
        
        if farmaco_principal_seleccionado != "None":
            farmaco_principal = farmaco_principal_seleccionado
            G_actual, farmaco_real = crear_subgrafo_centrado(G_completo, farmaco_principal)
            if G_actual is None:
                st.error(f"Drug '{farmaco_principal}' not found!")
                st.stop()
            farmaco_principal = farmaco_real
            st.success(f"Showing network for: {farmaco_principal}")
        else:
            farmaco_principal = None
            G_actual = G_completo
        
        
        
        # Obtener todas las categorías presentes en el grafo actual
        all_categories = set()
        for node in G_actual.nodes():
            if farmaco_principal and node == farmaco_principal:
                continue
            all_categories.add(G_actual.nodes[node]['atc_category'])
        
        category_list = sorted(list(all_categories))
        
        # Inicializar checkboxes desmarcados
        if 'active_categories' not in st.session_state:
            st.session_state.active_categories = {cat: False for cat in category_list}
        
        # Mostrar checkboxes con contadores
        for category in category_list:
            count = sum(1 for node in G_actual.nodes() 
                       if G_actual.nodes[node]['atc_category'] == category and 
                       node != farmaco_principal)
            
            col1, col2 = st.columns([1, 3])
            with col1:
                checked = st.checkbox(
                    "",
                    value=st.session_state.active_categories.get(category, False),
                    key=f"checkbox_{category.replace(' ', '_').replace(',', '')}",
                    label_visibility="collapsed"
                )
                if checked != st.session_state.active_categories.get(category, False):
                    st.session_state.active_categories[category] = checked
            
            with col2:
                atc_letter = get_atc_letter_from_category(category)
                if atc_letter:
                    color = ATC_COLORS.get(atc_letter, '#CCCCCC')
                else:
                    color = '#CCCCCC'
                
                st.markdown(
                    f"<span style='color:{color}; font-weight:bold;'>■</span> "
                    f"{category} ({count})",
                    unsafe_allow_html=True
                )
        
        # Botones de selección rápida
        col1, col2 = st.columns(2)
        
        with col1:
            select_all_clicked = st.button("Select All", use_container_width=True, key="select_all_btn")
        
        with col2:
            deselect_all_clicked = st.button("Deselect All", use_container_width=True, key="deselect_all_btn")
        
        if select_all_clicked:
            for category in category_list:
                st.session_state.active_categories[category] = True
            st.rerun()
        
        if deselect_all_clicked:
            for category in category_list:
                st.session_state.active_categories[category] = False
            st.rerun()
    
    # Área principal
    st.subheader("Network Visualization")
    
    # Verificar si hay algo para mostrar
    if not any(st.session_state.active_categories.values()) and not farmaco_principal:
        st.info("**Please select at least one ATC category from the sidebar to display drugs.**")
        
        with st.expander("Dataset Statistics", expanded=True):
            col1, col2 = st.columns(2)
            with col1:
                st.metric("Total Drugs in Dataset", len(G_completo.nodes()))
            with col2:
                st.metric("Total Interactions", len(G_completo.edges()))
            
            category_counts = Counter()
            for node in G_completo.nodes():
                category_counts[G_completo.nodes[node]['atc_category']] += 1
            
            st.write("**Drugs by ATC Category:**")
            for category, count in sorted(category_counts.items()):
                atc_letter = get_atc_letter_from_category(category)
                if atc_letter:
                    color = ATC_COLORS.get(atc_letter, '#CCCCCC')
                else:
                    color = '#CCCCCC'
                
                st.markdown(
                    f"<span style='color:{color}; font-weight:bold;'>■</span> "
                    f"{category}: {count} drugs",
                    unsafe_allow_html=True
                )
    else:
        # Crear y mostrar gráfico
        fig, G_filtrado = crear_grafo_plotly(
            G_actual, 
            farmaco_principal, 
            st.session_state.active_categories
        )
        
        st.plotly_chart(fig, use_container_width=True, 
                       config={'displayModeBar': True, 'scrollZoom': True})
        
        # Mostrar estadísticas
        st.subheader("Statistics")
        
        if G_filtrado:
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Displayed Drugs", len(G_filtrado.nodes()))
            with col2:
                st.metric("Displayed Interactions", len(G_filtrado.edges()))
            with col3:
                if farmaco_principal:
                    degree = G_filtrado.degree(farmaco_principal)
                    st.metric(f"Connections of {farmaco_principal}", degree)
            
            # Mostrar análisis de interacciones
            mostrar_analisis_interacciones(G_filtrado, farmaco_principal)
        
        # Mostrar lista de fármacos visibles
        if G_filtrado and len(G_filtrado.nodes()) > 0:
            with st.expander("View list of displayed drugs", expanded=False):
                drugs_list = sorted(G_filtrado.nodes())
                for drug in drugs_list:
                    if drug == farmaco_principal:
                        st.markdown(f"**🔵 {drug}** (Main drug)")
                    else:
                        category = G_filtrado.nodes[drug]['atc_category']
                        color = G_filtrado.nodes[drug]['color']
                        st.markdown(
                            f"<span style='color:{color};'>●</span> {drug} ({category})",
                            unsafe_allow_html=True
                        )



def procesar_pais(dd_df, essentials_df, pais):
    """
    Filtrar interacciones para un país específico basado en fármacos esenciales
    """
    # 1. Obtener fármacos esenciales del país
    if pais in essentials_df.columns:
        farmacos_esenciales = essentials_df[essentials_df[pais] == 1]['Medicine'].tolist()
        atc_codes_esenciales = essentials_df[essentials_df[pais] == 1]['ATC code primary'].dropna().tolist()
        
        # También considerar códigos secundarios si existen
        if 'ATC code secondary' in essentials_df.columns:
            atc_secundarios = essentials_df[essentials_df[pais] == 1]['ATC code secondary'].dropna().tolist()
            atc_codes_esenciales.extend(atc_secundarios)
    
    else:
        st.error(f"País '{pais}' no encontrado en el dataset")
        return pd.DataFrame()
    
    # 2. Crear lista de ATC codes únicos
    atc_uniques = []
    for atc in atc_codes_esenciales:
        if isinstance(atc, str):
            # Separar por | si hay múltiples códigos
            atc_uniques.extend([a.strip() for a in atc.split('|')])
    
    # 3. Filtrar DDIBUENO por fármacos esenciales
    # Opción A: Filtrar por nombre de fármaco
    filtered_by_name = dd_df[
        dd_df['Common_name_x'].isin(farmacos_esenciales) | 
        dd_df['Common_name_y'].isin(farmacos_esenciales)
    ].copy()
    
    # Opción B: Filtrar por código ATC (más robusto)
    filtered_by_atc = dd_df[
        dd_df['atc_code_x'].isin(atc_uniques) | 
        dd_df['atc_code_y'].isin(atc_uniques)
    ].copy()
    
    # Combinar ambos filtros
    filtered_df = pd.concat([filtered_by_name, filtered_by_atc]).drop_duplicates()
    
    # 4. Agregar información del país
    filtered_df['Country'] = pais
    filtered_df['Drug_A_Essential'] = filtered_df['Common_name_x'].isin(farmacos_esenciales)
    filtered_df['Drug_B_Essential'] = filtered_df['Common_name_y'].isin(farmacos_esenciales)
    
    return filtered_df, farmacos_esenciales, atc_uniques






def pestaña_esenciales(dd_df, essentials_df):
    """Pestaña para análisis de fármacos esenciales por país"""
    
    st.title("🌍 Essential Drugs Analysis")
    
    # Sidebar para seleccionar país
    with st.sidebar:
        st.header("Country Selection")
        
        # Obtener lista de países del dataset
        # Las primeras columnas son: Medicine, ATC code primary, ATC code secondary, WHO Model List
        paises = [col for col in essentials_df.columns if col not in 
                 ['Medicine', 'ATC code primary', 'ATC code secondary', 'WHO Model List of Essential Medications']]
        
        # Seleccionar país
        pais_seleccionado = st.selectbox(
            "Select Country:",
            paises,
            index=paises.index('Mexico') if 'Mexico' in paises else 0
        )
        
        # Opción para múltiples países
        paises_seleccionados = st.multiselect(
            "Or select multiple countries:",
            paises,
            default=[pais_seleccionado] if pais_seleccionado else []
        )
        
        # Opciones de análisis
        st.header("Analysis Options")
        show_stats = st.checkbox("Show Statistics", True)
        show_network = st.checkbox("Show Network Visualization", True)
        show_interactions = st.checkbox("Show Drug Interactions", True)
        compare_countries = st.checkbox("Compare Countries", False)
    
    # Contenido principal
    if compare_countries and len(paises_seleccionados) > 1:
        mostrar_comparacion_paises(dd_df, essentials_df, paises_seleccionados)
    elif len(paises_seleccionados) == 1:
        mostrar_analisis_pais(dd_df, essentials_df, paises_seleccionados[0], 
                             show_stats, show_network, show_interactions)
    else:
        st.info("👈 Please select a country from the sidebar to start analysis")

def mostrar_analisis_pais(dd_df, essentials_df, pais, show_stats, show_network, show_interactions):
    """Mostrar análisis para un país específico"""
    
    st.header(f"🇺🇳 Analysis for {pais}")
    
    # Procesar datos del país
    with st.spinner(f"Processing data for {pais}..."):
        filtered_df, farmacos_esenciales, atc_codes = procesar_pais(dd_df, essentials_df, pais)
    
    if filtered_df.empty:
        st.warning(f"No interactions found for essential drugs in {pais}")
        return
    
    # 1. Estadísticas principales
    if show_stats:
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Essential Drugs", len(farmacos_esenciales))
        with col2:
            st.metric("Unique ATC Codes", len(atc_codes))
        with col3:
            st.metric("Total Interactions", len(filtered_df))
        with col4:
            essential_interactions = len(filtered_df[
                filtered_df['Drug_A_Essential'] & filtered_df['Drug_B_Essential']
            ])
            st.metric("Both Essential", essential_interactions)
    
    # 2. Top fármacos más interactivos
    st.subheader("📊 Top Interactive Essential Drugs")
    
    # Contar interacciones por fármaco
    interacciones_por_farmaco = {}
    for _, row in filtered_df.iterrows():
        if row['Drug_A_Essential']:
            interacciones_por_farmaco[row['Common_name_x']] = \
                interacciones_por_farmaco.get(row['Common_name_x'], 0) + 1
        if row['Drug_B_Essential']:
            interacciones_por_farmaco[row['Common_name_y']] = \
                interacciones_por_farmaco.get(row['Common_name_y'], 0) + 1
    
    if interacciones_por_farmaco:
        top_drugs = pd.DataFrame(
            interacciones_por_farmaco.items(),
            columns=['Drug', 'Interaction Count']
        ).sort_values('Interaction Count', ascending=False).head(10)
        
        st.dataframe(top_drugs, use_container_width=True)
    
    # 3. Visualización de red (similar a tu código original)
    if show_network:
        st.subheader("🌐 Interaction Network")
        
        # Crear grafo solo con fármacos esenciales
        G_pais = crear_grafo_completo(filtered_df)
        
        # Filtrar solo nodos esenciales
        nodes_to_keep = [
            node for node in G_pais.nodes() 
            if node in farmacos_esenciales
        ]
        
        G_filtered = G_pais.subgraph(nodes_to_keep).copy()
        
        # Crear visualización
        fig, _ = crear_grafo_plotly(G_filtered)
        st.plotly_chart(fig, use_container_width=True)
    
    # 4. Tabla de interacciones
    if show_interactions:
        st.subheader("💊 Drug Interactions Details")
        
        with st.expander("View All Interactions", expanded=False):
            # Agrupar por severidad
            filtered_df['Severity'] = filtered_df['Y'].map({
                0: 'None', 1: 'Mild', 2: 'Moderate', 3: 'Severe', 4: 'Contraindicated'
            })
            
            st.dataframe(
                filtered_df[['Common_name_x', 'Common_name_y', 'Y', 'Severity', 
                           'Drug_A_Essential', 'Drug_B_Essential']].rename(
                    columns={
                        'Common_name_x': 'Drug A',
                        'Common_name_y': 'Drug B',
                        'Y': 'Y Value',
                        'Drug_A_Essential': 'A is Essential',
                        'Drug_B_Essential': 'B is Essential'
                    }
                ),
                use_container_width=True,
                hide_index=True
            )
    
    # 5. Análisis por categoría ATC
    st.subheader("🏥 Analysis by ATC Category")
    
    # Mapear códigos ATC a categorías
    ATC_CATEGORIES = {
        'A': 'Alimentary Tract',
        'B': 'Blood',
        'C': 'Cardiovascular',
        'D': 'Dermatologicals',
        'G': 'Genito Urinary',
        'H': 'Hormones',
        'J': 'Antiinfectives',
        'L': 'Antineoplastic',
        'M': 'Musculoskeletal',
        'N': 'Nervous System',
        'P': 'Antiparasitic',
        'R': 'Respiratory',
        'S': 'Sensory Organs',
        'V': 'Various'
    }
    
    # Contar fármacos por categoría
    categorias_count = {}
    for atc in atc_codes:
        if atc and len(atc) > 0:
            cat = atc[0]
            categorias_count[ATC_CATEGORIES.get(cat, 'Other')] = \
                categorias_count.get(ATC_CATEGORIES.get(cat, 'Other'), 0) + 1
    
    if categorias_count:
        import plotly.express as px
        cats_df = pd.DataFrame(
            categorias_count.items(),
            columns=['Category', 'Count']
        ).sort_values('Count', ascending=False)
        
        fig = px.bar(
            cats_df,
            x='Category',
            y='Count',
            title=f"Essential Drugs by Category in {pais}",
            color='Count',
            color_continuous_scale='Blues'
        )
        st.plotly_chart(fig, use_container_width=True)

def mostrar_comparacion_paises(dd_df, essentials_df, paises):
    """Comparar múltiples países"""
    
    st.header("📊 Country Comparison")
    
    # Recopilar datos para cada país
    comparacion_data = []
    
    for pais in paises:
        filtered_df, farmacos_esenciales, atc_codes = procesar_pais(dd_df, essentials_df, pais)
        
        comparacion_data.append({
            'Country': pais,
            'Essential Drugs': len(farmacos_esenciales),
            'ATC Codes': len(atc_codes),
            'Total Interactions': len(filtered_df),
            'High Risk Interactions': len(filtered_df[filtered_df['Y'] >= 3]) if not filtered_df.empty else 0
        })
    
    comparacion_df = pd.DataFrame(comparacion_data)
    
    # Mostrar tabla comparativa
    st.dataframe(
        comparacion_df,
        use_container_width=True,
        hide_index=True
    )
    
    # Gráfico comparativo
    import plotly.express as px
    
    fig = px.bar(
        comparacion_df,
        x='Country',
        y=['Essential Drugs', 'Total Interactions', 'High Risk Interactions'],
        title="Country Comparison",
        barmode='group'
    )
    st.plotly_chart(fig, use_container_width=True)
    
    # Matriz de similitud
    st.subheader("Similarity between Countries")
    
    # Calcular similitud de Jaccard
    similarity_matrix = []
    for i, pais1 in enumerate(paises):
        row = []
        farmacos1 = set(essentials_df[essentials_df[pais1] == 1]['Medicine'].tolist())
        
        for j, pais2 in enumerate(paises):
            farmacos2 = set(essentials_df[essentials_df[pais2] == 1]['Medicine'].tolist())
            
            if len(farmacos1.union(farmacos2)) > 0:
                similitud = len(farmacos1.intersection(farmacos2)) / len(farmacos1.union(farmacos2))
            else:
                similitud = 0
            
            row.append(similitud)
        similarity_matrix.append(row)
    
    # Mostrar matriz de similitud
    import plotly.graph_objects as go
    
    fig = go.Figure(data=go.Heatmap(
        z=similarity_matrix,
        x=paises,
        y=paises,
        colorscale='Blues',
        text=[[f"{val:.2f}" for val in row] for row in similarity_matrix],
        texttemplate="%{text}",
        textfont={"size": 10}
    ))
    
    fig.update_layout(
        title="Drug List Similarity (Jaccard Index)",
        xaxis_title="Country",
        yaxis_title="Country"
    )
    
    st.plotly_chart(fig, use_container_width=True)










   

def pestaña_significado_y():
    st.title("Dataset, menaing of Y values")
    st.header('Dataset')
    st.dataframe(meeaning_y,use_container_width=True)
    




































































def main():
    # Crear pestañas de navegación
    tab2, tab3, tab4 = st.tabs([

        "Network interactions Visualization", 
        "Essential Drugs by Country",
        "Y dataset"
    ])
    
    
    
    with tab2:
        pestaña_visualizacion()
    
    with tab3:
        pestaña_esenciales()
    with tab4:
        pestaña_significado_y()
if __name__ == "__main__":
    main()
