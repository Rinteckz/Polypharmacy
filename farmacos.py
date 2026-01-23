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
def load_essential_drugs():
    
    essentials_df = pd.read_csv(r"farmacos_esenciales.csv")
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
st.title("**Drug-Drug interactions**")
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
    
    # Mostrar codigo de efecto (este se busca en el csv que pondré en la mismsa pag)
    st.write("**Effect by DDI code (Y):**")
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
    
    with st.expander("All Interactions in the network", expanded=False):
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
        st.write(f"**Interactions in the network involving {farmaco_principal}:**")
        
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
            st.write(f"**{farmaco_principal} to drug B ({len(outgoing)}):**")
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
            st.write(f"**Drug A to {farmaco_principal} ({len(incoming)}):**")
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
    st.title("Drug-Drug Interaction Network Visualization")
    
    # Inicializar estado de sesión
    if 'G_completo' not in st.session_state:
        with st.spinner("Loading drug interaction data..."):
            st.session_state.G_completo = crear_grafo_completo(df)
    
    G_completo = st.session_state.G_completo
    
    # Sidebar para controles
    with st.sidebar:
        st.header("Controls")
        
        
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
        
        
        
        all_categories = set()
        for node in G_actual.nodes():
            if farmaco_principal and node == farmaco_principal:
                continue
            all_categories.add(G_actual.nodes[node]['atc_category'])
        
        category_list = sorted(list(all_categories))
        
        if 'active_categories' not in st.session_state:
            st.session_state.active_categories = {cat: False for cat in category_list}
        
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
        
        # Botones de selección r
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
        st.info("<--- **Open the sidebar to choose a drug**")
        
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
            
            
            # Mostrar análisis de interacciones
            mostrar_analisis_interacciones(G_filtrado, farmaco_principal)
        
        if G_filtrado and len(G_filtrado.nodes()) > 0:
            with st.expander("View list of displayed drugs", expanded=False):
                drugs_list = sorted(G_filtrado.nodes())
                for drug in drugs_list:
                    if drug == farmaco_principal:
                        st.markdown(f"**{drug}** (Main drug)")
                    else:
                        category = G_filtrado.nodes[drug]['atc_category']
                        color = G_filtrado.nodes[drug]['color']
                        st.markdown(
                            f"<span style='color:{color};'>●</span> {drug} ({category})",
                            unsafe_allow_html=True
                        )


   

def pestaña_significado_y():
    st.title("Dataset, menaing of Y values")
    st.header('Dataset')
    st.dataframe(meeaning_y,use_container_width=True)
    





























































def pestaña_esenciales():
    """Pestaña simple para fármacos esenciales"""
    
    st.title("Essential Drugs By WHO")
    
    # Cargar datos
    essentials_df = load_essential_drugs()
    
    if essentials_df.empty:
        st.info("Please install openpyxl to load essential drugs data.")
        return
    
    # Obtener lista de países
    non_country_cols = ['Medicine', 'ATC code primary', 'ATC code secondary', 
                       'WHO Model List of Essential Medications']
    paises = [col for col in essentials_df.columns if col not in non_country_cols]
    
    if not paises:
        st.warning("No countries found in the dataset.")
        return
    
    # Seleccionar país
    st.subheader("Select Country")
    pais_seleccionado = st.selectbox(
        "Choose a country:",
        paises,
        index=paises.index('Mexico') if 'Mexico' in paises else 0
    )
    
    # Filtrar fármacos esenciales del país
    farmacos_pais = essentials_df[essentials_df[pais_seleccionado] == 1]
    
    if farmacos_pais.empty:
        st.info(f"No essential drugs found for {pais_seleccionado}")
        return
    
    st.subheader(f"Statistics for {pais_seleccionado}")
    
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Total Essential Drugs", len(farmacos_pais))
    
    
    codigos_atc_esenciales = set()
    
    codigos_atc_esenciales.update(farmacos_pais['ATC code primary'].dropna().unique())
    
    for codigos in farmacos_pais['ATC code secondary'].dropna():
        if isinstance(codigos, str):
            if '|' in codigos:
                for codigo in codigos.split('|'):
                    codigos_atc_esenciales.add(codigo.strip())
            else:
                codigos_atc_esenciales.add(codigos.strip())
    
    codigos_atc_esenciales = {codigo.strip() for codigo in codigos_atc_esenciales 
                              if pd.notna(codigo) and codigo.strip() != ''}
    
    #with col2:
        #st.metric("Unique ATC Codes", len(codigos_atc_esenciales))
    
    
    def verificar_farmaco(atc_code, num_atc):
        if pd.isna(atc_code) or atc_code == 'No ATC' or num_atc == 0:
            return False
        
        if '|' in str(atc_code):
            codigos = [c.strip() for c in str(atc_code).split('|')]
            return any(codigo in codigos_atc_esenciales for codigo in codigos)
        else:
            return str(atc_code).strip() in codigos_atc_esenciales
    
    # Aplicar verificación a ambas columnas
    df['esencial_x'] = df.apply(lambda row: verificar_farmaco(row['atc_code_x'], row['num_atc_x']), axis=1)
    df['esencial_y'] = df.apply(lambda row: verificar_farmaco(row['atc_code_y'], row['num_atc_y']), axis=1)
    
    # Contar cuántos son esenciales
    esencial_x_count = df['esencial_x'].sum()
    esencial_y_count = df['esencial_y'].sum()
    
    #with col3:
        #st.metric("Drugs in DDIBUENO", esencial_x_count + esencial_y_count)
    
    
    farmacos_esenciales_lista = []
    
    # Procesar fármacos de la columna X
    for idx, row in df[df['esencial_x']].iterrows():
        farmacos_esenciales_lista.append({
            'Common_Name': row['Common_name_x'],
            'ATC_Code': row['atc_code_x'],
            'Num_ATC': row['num_atc_x'],
            'Source': 'Column X',
            'Original_Index': idx,
            'Drug_Pair': f"{row['Common_name_x']} - {row['Common_name_y']}",
            'Y_Value': row['Y']
        })
    
    for idx, row in df[df['esencial_y']].iterrows():
        farmacos_esenciales_lista.append({
            'Common_Name': row['Common_name_y'],
            'ATC_Code': row['atc_code_y'],
            'Num_ATC': row['num_atc_y'],
            'Source': 'Column Y',
            'Original_Index': idx,
            'Drug_Pair': f"{row['Common_name_x']} - {row['Common_name_y']}",
            'Y_Value': row['Y']
        })
    
    # Crear DataFrame
    if not farmacos_esenciales_lista:
        st.info(f"No essential drugs from {pais_seleccionado} found in DDI database")
        return
    
    df_farmacos_esenciales = pd.DataFrame(farmacos_esenciales_lista)
    
    # Eliminar duplicados (mismo nombre y ATC)
    df_farmacos_esenciales = df_farmacos_esenciales.drop_duplicates(subset=['Common_Name', 'ATC_Code'])
    
    # Ordenar por nombre
    df_farmacos_esenciales = df_farmacos_esenciales.sort_values('Common_Name')
    
    
    st.subheader(f"Essential Drugs from {pais_seleccionado} found in DDIBUENO")
    
    # Mostrar tabla principal
    st.dataframe(
        df_farmacos_esenciales[['Common_Name', 'ATC_Code']],
        use_container_width=True,
        column_config={
            "Common_Name": st.column_config.TextColumn("Drug Name", width="medium"),
            "ATC_Code": st.column_config.TextColumn("ATC Code", width="small")
        },
        hide_index=True
    )
    
    # Mostrar estadísticas detalladas
    st.subheader("Detailed Statistics")
    
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Unique Essential Drugs", len(df_farmacos_esenciales))
    with col2:
        st.metric("From Column X", len(df_farmacos_esenciales[df_farmacos_esenciales['Source'] == 'Column X']))
    with col3:
        st.metric("From Column Y", len(df_farmacos_esenciales[df_farmacos_esenciales['Source'] == 'Column Y']))
    #with col4:
       # unique_drugs = set(df_farmacos_esenciales['Common_Name'])
        #st.metric("Unique Drug Names", len(unique_drugs))
    
    
    st.subheader("Distribution by ATC Group")
    
    def extraer_grupo_atc(atc_code):
        if pd.isna(atc_code) or atc_code == 'No ATC':
            return 'Sin ATC'
        if isinstance(atc_code, str) and len(atc_code) > 0:
            return atc_code[0]
        return 'Unknown'
    
    df_farmacos_esenciales['ATC_Group'] = df_farmacos_esenciales['ATC_Code'].apply(extraer_grupo_atc)
    
    ATC_GROUP_NAMES = {
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
    
    # Crear DataFrame para distribución
    distribucion = df_farmacos_esenciales['ATC_Group'].value_counts().reset_index()
    distribucion.columns = ['ATC_Group', 'Count']
    distribucion['Group_Name'] = distribucion['ATC_Group'].map(ATC_GROUP_NAMES)
    
    # Mostrar gráfico de barras
    fig = px.bar(
        distribucion,
        x='Group_Name',
        y='Count',
        title=f"Essential Drugs by ATC Group - {pais_seleccionado}",
        color='Count',
        color_continuous_scale='Blues',
        labels={'Group_Name': 'ATC Group', 'Count': 'Number of Drugs'}
    )
    st.plotly_chart(fig, use_container_width=True)
    
    # Mostrar tabla de distribución
    with st.expander("View Distribution Table", expanded=False):
        st.dataframe(
            distribucion,
            use_container_width=True,
            hide_index=True
        )
    
    
    st.subheader("Interactions Involving Essential Drugs")
    
    essential_drugs_names = set(df_farmacos_esenciales['Common_Name'])
    
    essential_interactions = df[
        df['Common_name_x'].isin(essential_drugs_names) | 
        df['Common_name_y'].isin(essential_drugs_names)
    ]
    
    if not essential_interactions.empty:
        st.write(f"**Found {len(essential_interactions)} interactions involving essential drugs**")
        
        with st.expander("View All Interactions", expanded=False):
            st.dataframe(
                essential_interactions[['Common_name_x', 'Common_name_y', 'Y', 'atc_code_x', 'atc_code_y']].rename(
                    columns={
                        'Common_name_x': 'Drug A',
                        'Common_name_y': 'Drug B',
                        'Y': 'Y Value',
                        'atc_code_x': 'ATC Code A',
                        'atc_code_y': 'ATC Code B'
                    }
                ),
                use_container_width=True,
                hide_index=True
            )
        
        st.write("**Severity Analysis (Y Values):**")
        
        y_counts = essential_interactions['Y'].value_counts().sort_index()
        
        col1, col2 = st.columns(2)
        
        with col1:
            for y_val, count in y_counts.items():
                percentage = (count / len(essential_interactions)) * 100
                severity = {
                    0: "None", 1: "Mild", 2: "Moderate", 
                    3: "Severe", 4: "Contraindicated"
                }.get(y_val, "Unknown")
                
                st.write(f"**Y={y_val} ({severity}):** {count} ({percentage:.1f}%)")
        
        with col2:
            fig_pie = px.pie(
                values=y_counts.values,
                names=[f"Y={k}" for k in y_counts.index],
                title="Severity Distribution",
                color_discrete_sequence=px.colors.sequential.Reds
            )
            st.plotly_chart(fig_pie, use_container_width=True)
        
        # Top fármacos más interactivos
        st.write("**Most Interactive Essential Drugs:**")
        
        interaction_counts = {}
        for _, row in essential_interactions.iterrows():
            if row['Common_name_x'] in essential_drugs_names:
                interaction_counts[row['Common_name_x']] = interaction_counts.get(row['Common_name_x'], 0) + 1
            if row['Common_name_y'] in essential_drugs_names:
                interaction_counts[row['Common_name_y']] = interaction_counts.get(row['Common_name_y'], 0) + 1
        
        if interaction_counts:
            top_drugs = pd.DataFrame(
                interaction_counts.items(),
                columns=['Drug', 'Interaction_Count']
            ).sort_values('Interaction_Count', ascending=False).head(10)
            
            st.dataframe(
                top_drugs,
                use_container_width=True,
                hide_index=True
            )
    
    else:
        st.info("No interactions found involving these essential drugs")
    

def about():
    st.title("About this web page")
    st.write("This application has been developed to visualize drug-drug interactions using network graphs. It allows user to navigate trhough a large dataset of drugs and their interactions.")
    st.write("We showcase essential drugs as defined by the World Healt Organization (WHO) and compare them with the dataset DDI, to identify which essential drugs have kwon interactions.")
    st.markdown(
    "The DDI dataset was obtained from <a href='https://tdcommons.ai/'>Therapeutic Data Commons</a>",
    unsafe_allow_html=True
    )

    st.markdown(
    "The essential drug list was obtained from <a href='https://pmc.ncbi.nlm.nih.gov/articles/PMC6560372/pdf/BLT.18.222448.pdf'>Comparison of essential medicines lists in 137 countries</a>",
    unsafe_allow_html=True
    )

    st.markdown("This web-page has been developed by Adolfo Guzmán Arenas, Jorge Luis Chavez Perez, Gilberto Lorenzo Martínez Luna y Edgardo ALberto Marin Colli")





def manual():
    """
    Función que contiene todo el contenido del manual de usuario
    """
    st.title("📖 Manual de Usuario")
    
    # Navegación interna del manual
    manual_tabs = st.tabs(["Guía Rápida", "Funcionalidades", "Tutorial", "FAQ", "Contacto"])
    
    with manual_tabs[0]:
        st.header("🚀 Guía Rápida de Inicio")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("Primeros Pasos")
            st.markdown("""
            1. **Carga tus datos** en la pestaña 'Y dataset'
            2. **Configura parámetros** en 'Essentials'
            3. **Visualiza resultados** en 'Network interactions Visualization'
            4. **Exporta** tus análisis
            """)
            
        with col2:
            st.subheader("Atajos Importantes")
            st.markdown("""
            - `Ctrl+R`: Recargar aplicación
            - `Ctrl+S`: Guardar configuración
            - **Hover**: Sobre iconos para info
            - **Click derecho**: En gráficos para opciones
            """)
    
    with manual_tabs[1]:
        st.header("🔧 Funcionalidades Detalladas")
        
        # Funcionalidades por pestaña
        with st.expander("📊 Network interactions Visualization", expanded=True):
            st.markdown("""
            ### Visualización de Redes
            - **Gráficos interactivos**: Zoom, pan, selección
            - **Filtros dinámicos**: Por tipo de interacción, peso, dirección
            - **Exportación**: PNG, SVG, PDF
            - **Análisis**: Métricas de centralidad, clustering
            """)
        
        with st.expander("⚙️ Essentials"):
            st.markdown("""
            ### Configuración Esencial
            - **Parámetros de red**: Umbrales, algoritmos
            - **Preprocesamiento**: Normalización, filtrado
            - **Opciones de visualización**: Colores, etiquetas, layout
            """)
        
        with st.expander("📁 Y dataset"):
            st.markdown("""
            ### Gestión de Datos
            - **Formatos soportados**: CSV, Excel, JSON, TXT
            - **Límites**: Hasta 1GB de datos
            - **Previsualización**: Tabla interactiva
            - **Limpieza automática**: Valores nulos, duplicados
            """)
    
    with manual_tabs[2]:
        st.header("🎬 Tutorial Paso a Paso")
        
        step = st.select_slider(
            "Selecciona el paso:",
            options=["Paso 1", "Paso 2", "Paso 3", "Paso 4", "Paso 5"]
        )
        
        if step == "Paso 1":
            st.subheader("Carga de Datos")
            st.markdown("""
            1. Ve a la pestaña **'Y dataset'**
            2. Haz clic en **'Browse files'**
            3. Selecciona tu archivo (CSV, Excel, etc.)
            4. Revisa la previsualización
            5. Confirma con **'Load Data'**
            """)
            st.image("tutorial_step1.png", caption="Interfaz de carga de datos")
            
        elif step == "Paso 2":
            st.subheader("Configuración de Parámetros")
            st.markdown("""
            1. Navega a **'Essentials'**
            2. Ajusta los umbrales según necesites
            3. Selecciona el algoritmo de análisis
            4. Configura opciones de visualización
            5. Guarda la configuración
            """)
    
    with manual_tabs[3]:
        st.header("❓ Preguntas Frecuentes")
        
        faq_items = {
            "¿Cómo interpreto los resultados de la red?": """
            **Métricas clave a considerar:**
            - **Degree centrality**: Nodos más conectados
            - **Betweenness**: Nodos que actúan como puentes
            - **Clustering coefficient**: Nivel de agrupamiento
            
            Los nodos más grandes/coloreados son más importantes.
            """,
            
            "¿Qué formatos de archivo acepta?": """
            **Formatos soportados:**
            - CSV (separado por comas)
            - Excel (.xlsx, .xls)
            - JSON
            - TXT (formato tabular)
            
            **Requisitos:**
            - Máximo 1GB por archivo
            - UTF-8 encoding recomendado
            """,
            
            "¿Puedo guardar mi sesión?": """
            **Sí, de dos formas:**
            
            1. **Configuración**: Exporta settings como JSON
            2. **Resultados**: Descarga gráficos y tablas
            3. **Sesión completa**: Disponible en versión premium
            
            *Nota: Los datos subidos no se almacenan permanentemente*
            """
        }
        
        for question, answer in faq_items.items():
            with st.expander(question):
                st.markdown(answer)
    
    with manual_tabs[4]:
        st.header("📞 Contacto y Soporte")
        
        st.subheader("¿Necesitas ayuda adicional?")
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.markdown("""
            ### 📧 Email
            **Soporte técnico:**  
            support@networkapp.com
            
            **Consultas generales:**  
            info@networkapp.com
            """)
        
        with col2:
            st.markdown("""
            ### 📚 Recursos
            - [Documentación completa](https://docs.example.com)
            - [Video tutoriales](https://youtube.com/example)
            - [Foro de la comunidad](https://forum.example.com)
            """)
        
        with col3:
            st.markdown("""
            ### 🐛 Reportar Problemas
            **GitHub Issues:**  
            [github.com/example/issues](https://github.com/example/issues)
            
            **Prioridad:**
            - Urgente: < 24 horas
            - Normal: 2-3 días
            """)
        
        # Formulario de contacto simple
        st.divider()
        st.subheader("Envíanos tu consulta")
        
        with st.form("contact_form"):
            name = st.text_input("Nombre")
            email = st.text_input("Email")
            issue_type = st.selectbox("Tipo de consulta", 
                                     ["Problema técnico", "Sugerencia", "Pregunta general"])
            message = st.text_area("Mensaje", height=150)
            
            submitted = st.form_submit_button("Enviar consulta")
            if submitted:
                st.success("¡Consulta enviada! Te responderemos en 24-48 horas.")
def main():
    # Crear pestañas de navegación
    tab1, tab2, tab3, tab4, tab5 = st.tabs([

        "Network interactions Visualization", 
        "Essentials",
        "Y dataset",
        "About"
        "User Manual"
    ])
    
    
    
    with tab1:
        pestaña_visualizacion()
    with tab2:
        pestaña_esenciales()
    with tab3:
        pestaña_significado_y()
    with tab4:
        about()
    with tab5:
        manual()

if __name__ == "__main__":
    main()
