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
    df = pd.read_csv(r"DDIBUENO.csv")
    
    # Mostrar información del dataset
    st.sidebar.info(f"Dataset loaded: {len(df)} rows")
    
    # Verificar duplicados en las interacciones
    unique_interactions = df[['Common_name_x', 'Common_name_y']].drop_duplicates()
    st.sidebar.info(f"Unique drug pairs: {len(unique_interactions)}")
    
    # Contar interacciones por tipo
    interaction_counts = df['Y'].value_counts()
    st.sidebar.info(f"Interaction types: Type 1: {interaction_counts.get(1, 0)}, Type 2: {interaction_counts.get(2, 0)}")
    
    return df

df = load_data()

ATC_CATEGORIES = {
    'A': 'Alimentary tract',
    'B': 'Blood organs',
    'C': 'Cardiovascular',
    'D': 'Dermatologicals',
    'G': 'Genito-urinary',
    'H': 'Hormonal',
    'J': 'Anti-infectives',
    'L': 'Antineoplastic',
    'M': 'Musculo-skeletal',
    'N': 'Nervous system',
    'P': 'Antiparasitic',
    'R': 'Respiratory',
    'S': 'Sensory organs',
    'V': 'Various',
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

@st.cache_data
def crear_grafo_completo(df):
    """Crear grafo completo con todos los fármacos - EXACTO como en Matplotlib"""
    G = nx.DiGraph()
    
    # Estadísticas para debug
    total_rows = len(df)
    edges_added = 0
    edges_skipped = 0
    
    for idx, row in df.iterrows():
        drug1 = row['Common_name_x']
        drug2 = row['Common_name_y']
        interaction_type = row['Y'] 
        
        atc1 = row['atc_code_x'] if pd.notna(row['atc_code_x']) else "No ATC"
        num_atc1 = row['num_atc_x']
        color1, category1 = get_atc_color_and_category(atc1, num_atc1)
        
        atc2 = row['atc_code_y'] if pd.notna(row['atc_code_y']) else "No ATC"
        num_atc2 = row['num_atc_y']
        color2, category2 = get_atc_color_and_category(atc2, num_atc2)
        
        # Agregar nodo 1 si no existe
        if not G.has_node(drug1):
            tooltip_info = f"<b>Drug:</b> {drug1}<br><b>ATC Code:</b> {atc1}<br><b>ATC Category:</b> {category1}"
            G.add_node(drug1, 
                      atc_code=atc1,
                      atc_category=category1,
                      color=color1,
                      num_atc=num_atc1,
                      tooltip=tooltip_info)
        
        # Agregar nodo 2 si no existe
        if not G.has_node(drug2):
            tooltip_info = f"<b>Drug:</b> {drug2}<br><b>ATC Code:</b> {atc2}<br><b>ATC Category:</b> {category2}"
            G.add_node(drug2, 
                      atc_code=atc2,
                      atc_category=category2,
                      color=color2,
                      num_atc=num_atc2,
                      tooltip=tooltip_info)
        
        # Agregar arista (EXACTAMENTE como en Matplotlib)
        if not G.has_edge(drug1, drug2):
            interaction_desc = "Type 1" if interaction_type == 1 else "Type 2" if interaction_type == 2 else "Unknown"
            edge_tooltip = f"<b>Interaction Type:</b> {interaction_desc}<br><b>From:</b> {drug1} → <b>To:</b> {drug2}"
            G.add_edge(drug1, drug2, 
                      interaction_type=interaction_type,
                      tooltip=edge_tooltip)
            edges_added += 1
        else:
            edges_skipped += 1
    
    # Guardar estadísticas
    st.session_state.graph_stats = {
        'total_rows': total_rows,
        'edges_added': edges_added,
        'edges_skipped': edges_skipped,
        'unique_drugs': len(G.nodes()),
        'unique_interactions': len(G.edges())
    }
    
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
    
    # Si no hay categorías activas, no mostrar nada
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
            plot_bgcolor='white',
            height=600,
            width=800
        )
        return fig, None
    
    # Filtrar nodos por categorías activas
    nodes_to_keep = []
    
    # El fármaco principal SIEMPRE se muestra si existe
    if farmaco_principal:
        nodes_to_keep.append(farmaco_principal)
    
    # Filtrar otros nodos por categorías activas
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
            plot_bgcolor='white',
            height=600,
            width=800
        )
        return fig, None
    
    # Crear subgrafo filtrado
    edges_to_keep = []
    for u, v in G.edges():
        if u in nodes_to_keep and v in nodes_to_keep:
            edges_to_keep.append((u, v))
    
    G_filtered = nx.DiGraph()
    for node in nodes_to_keep:
        G_filtered.add_node(node, **G.nodes[node])
    
    for u, v in edges_to_keep:
        G_filtered.add_edge(u, v, **G[u][v])
    
    # Calcular posiciones de los nodos
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
    
    # Preparar datos para los nodos
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
    
    # Preparar datos para las aristas
    edge_x, edge_y, edge_texts = [], [], []
    
    for u, v, data in G_filtered.edges(data=True):
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
        edge_texts.append(data['tooltip'])
    
    # Crear trazas de Plotly
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=1.5, color='gray'),
        hoverinfo='text',
        text=edge_texts,
        mode='lines',
        hoverlabel=dict(bgcolor='white', font_size=12)
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
        hoverlabel=dict(bgcolor='white', font_size=12)
    )
    
    # Título
    if farmaco_principal:
        title_text = f"Drug: {farmaco_principal}<br>Drugs: {len(G_filtered.nodes())} | Interactions: {len(G_filtered.edges())}"
    else:
        title_text = f"Complete Network<br>Drugs: {len(G_filtered.nodes())} | Interactions: {len(G_filtered.edges())}"
    
    # Crear figura
    fig = go.Figure(data=[edge_trace, node_trace],
                   layout=go.Layout(
                       title=dict(text=title_text, font=dict(size=14), x=0.5, xanchor='center'),
                       showlegend=False,
                       hovermode='closest',
                       margin=dict(b=20, l=5, r=5, t=60),
                       xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                       yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                       plot_bgcolor='white',
                       height=600,
                       width=800,
                       dragmode='pan'
                   ))
    
    return fig, G_filtered

def main():
    st.title("Drug Interaction Network Visualization")
    
    # Inicializar estados
    if 'G_completo' not in st.session_state:
        with st.spinner("Loading drug interaction data..."):
            st.session_state.G_completo = crear_grafo_completo(df)
    
    if 'selected_drug' not in st.session_state:
        st.session_state.selected_drug = "None"
    
    if 'select_all_clicked' not in st.session_state:
        st.session_state.select_all_clicked = False
    
    if 'deselect_all_clicked' not in st.session_state:
        st.session_state.deselect_all_clicked = False
    
    if 'active_categories' not in st.session_state:
        st.session_state.active_categories = {}
    
    G_completo = st.session_state.G_completo
    
    # Sidebar
    with st.sidebar:
        st.header("Controls")
        
        # Mostrar estadísticas del dataset
        with st.expander("📊 Dataset Info", expanded=False):
            st.write(f"**CSV Rows:** {len(df):,}")
            st.write(f"**Unique Drugs:** {len(G_completo.nodes()):,}")
            st.write(f"**Unique Interactions:** {len(G_completo.edges()):,}")
            
            if 'graph_stats' in st.session_state:
                stats = st.session_state.graph_stats
                st.write(f"**Edges added:** {stats['edges_added']:,}")
                st.write(f"**Duplicates skipped:** {stats['edges_skipped']:,}")
                if stats['edges_skipped'] > 0:
                    st.info(f"Note: {stats['edges_skipped']:,} duplicate interactions were skipped")
        
        # Selección de fármaco principal
        st.subheader("Select Main Drug (Optional)")
        all_drugs = sorted(set(df['Common_name_x'].tolist() + df['Common_name_y'].tolist()))
        
        farmaco_principal_seleccionado = st.selectbox(
            "Choose a drug to focus on:",
            ["None"] + all_drugs,
            index=0,
            key="drug_selector"
        )
        
        st.session_state.selected_drug = farmaco_principal_seleccionado
        
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
        
        # Filtrado por categoría ATC
        st.subheader("Filter by ATC Category")
        st.write("Select categories to display (none selected by default):")
        
        # Obtener categorías
        all_categories = set()
        for node in G_actual.nodes():
            if farmaco_principal and node == farmaco_principal:
                continue
            all_categories.add(G_actual.nodes[node]['atc_category'])
        
        category_list = sorted(list(all_categories))
        
        # Inicializar estados de categorías
        for category in category_list:
            if category not in st.session_state.active_categories:
                st.session_state.active_categories[category] = False
        
        # Botones de selección rápida
        col1, col2 = st.columns(2)
        with col1:
            if st.button("Select All", use_container_width=True, key="select_all_btn"):
                st.session_state.select_all_clicked = True
                st.session_state.deselect_all_clicked = False
                st.rerun()
        
        with col2:
            if st.button("Deselect All", use_container_width=True, key="deselect_all_btn"):
                st.session_state.deselect_all_clicked = True
                st.session_state.select_all_clicked = False
                st.rerun()
        
        # Aplicar cambios de botones
        if st.session_state.select_all_clicked:
            for category in category_list:
                st.session_state.active_categories[category] = True
            st.session_state.select_all_clicked = False
        
        if st.session_state.deselect_all_clicked:
            for category in category_list:
                st.session_state.active_categories[category] = False
            st.session_state.deselect_all_clicked = False
        
        # Mostrar checkboxes
        st.write("---")
        st.write("**Categories:**")
        
        for category in category_list:
            count = sum(1 for node in G_actual.nodes() 
                       if G_actual.nodes[node]['atc_category'] == category and 
                       node != farmaco_principal)
            
            checked = st.checkbox(
                f"{category} ({count:,})",
                value=st.session_state.active_categories[category],
                key=f"cat_checkbox_{category}"
            )
            
            if checked != st.session_state.active_categories[category]:
                st.session_state.active_categories[category] = checked
    
    # Área principal
    st.subheader("Network Visualization")
    
    # Determinar grafo actual
    if st.session_state.selected_drug != "None":
        farmaco_principal = st.session_state.selected_drug
        G_actual, farmaco_real = crear_subgrafo_centrado(G_completo, farmaco_principal)
        if G_actual is not None:
            farmaco_principal = farmaco_real
    else:
        farmaco_principal = None
        G_actual = G_completo
    
    # Verificar si hay algo para mostrar
    if not any(st.session_state.active_categories.values()) and not farmaco_principal:
        st.info("👈 **Please select at least one ATC category from the sidebar to display drugs.**")
        
        # Mostrar estadísticas DETALLADAS
        with st.expander("📊 Detailed Dataset Statistics", expanded=True):
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Total CSV Rows", f"{len(df):,}")
            with col2:
                st.metric("Unique Drugs", f"{len(G_completo.nodes()):,}")
            with col3:
                st.metric("Unique Interactions", f"{len(G_completo.edges()):,}")
            
            # Información adicional
            st.write("**Additional Statistics:**")
            
            # Contar fármacos únicos
            all_unique_drugs = set(df['Common_name_x'].tolist() + df['Common_name_y'].tolist())
            st.write(f"- Total unique drug names in CSV: {len(all_unique_drugs):,}")
            
            # Contar pares únicos (sin dirección)
            df['drug_pair'] = df.apply(lambda row: tuple(sorted([row['Common_name_x'], row['Common_name_y']])), axis=1)
            unique_pairs = df['drug_pair'].nunique()
            st.write(f"- Unique drug pairs (undirected): {unique_pairs:,}")
            
            # Distribución por categoría ATC
            st.write("**Drugs by ATC Category:**")
            category_counts = Counter()
            for node in G_completo.nodes():
                category_counts[G_completo.nodes[node]['atc_category']] += 1
            
            for category, count in sorted(category_counts.items()):
                color = ATC_COLORS.get(category[:1], ATC_COLORS.get(category, '#CCCCCC'))
                percentage = (count / len(G_completo.nodes())) * 100
                st.markdown(
                    f"<span style='color:{color}; font-weight:bold;'>■</span> "
                    f"{category}: {count:,} drugs ({percentage:.1f}%)",
                    unsafe_allow_html=True
                )
            
            # Estadísticas de duplicados si existen
            if 'graph_stats' in st.session_state:
                stats = st.session_state.graph_stats
                if stats['edges_skipped'] > 0:
                    st.warning(f"⚠️ **Note:** {stats['edges_skipped']:,} duplicate interactions were removed from the graph")
                    st.write(f"- Graph shows unique interactions only")
                    st.write(f"- Raw CSV contains {stats['total_rows']:,} rows")
                    st.write(f"- After removing duplicates: {stats['edges_added']:,} unique interactions")
    else:
        # Crear y mostrar gráfico
        fig, G_filtrado = crear_grafo_plotly(
            G_actual, 
            farmaco_principal, 
            st.session_state.active_categories
        )
        
        st.plotly_chart(fig, use_container_width=True, 
                       config={'displayModeBar': True, 'scrollZoom': True})
        
        # Mostrar estadísticas del gráfico filtrado
        st.subheader("Current View Statistics")
        
        if G_filtrado:
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Displayed Drugs", f"{len(G_filtrado.nodes()):,}")
            with col2:
                st.metric("Displayed Interactions", f"{len(G_filtrado.edges()):,}")
            with col3:
                if farmaco_principal:
                    degree = G_filtrado.degree(farmaco_principal)
                    st.metric(f"Connections of {farmaco_principal}", f"{degree:,}")
            
            # Detalles por categoría
            with st.expander("View category details", expanded=False):
                category_counts = Counter()
                for node in G_filtrado.nodes():
                    if node != farmaco_principal:
                        category_counts[G_filtrado.nodes[node]['atc_category']] += 1
                
                if category_counts:
                    st.write("**Drugs in view by category:**")
                    for category, count in sorted(category_counts.items()):
                        color = ATC_COLORS.get(category[:1], ATC_COLORS.get(category, '#CCCCCC'))
                        st.markdown(
                            f"<span style='color:{color}; font-weight:bold;'>■</span> "
                            f"{category}: {count:,} drugs",
                            unsafe_allow_html=True
                        )
                else:
                    st.write("No drugs from selected categories (only main drug shown)")
        
        # Mostrar lista de fármacos
        if G_filtrado and len(G_filtrado.nodes()) > 0:
            with st.expander("View list of displayed drugs", expanded=False):
                drugs_list = sorted(G_filtrado.nodes())
                st.write(f"**Total: {len(drugs_list):,} drugs**")
                
                # Mostrar en columnas
                cols = 3
                for i in range(0, len(drugs_list), cols):
                    col_items = drugs_list[i:i+cols]
                    col_objects = st.columns(cols)
                    for j, drug in enumerate(col_items):
                        with col_objects[j]:
                            if drug == farmaco_principal:
                                st.markdown(f"**🔵 {drug}**")
                                st.caption("(Main drug)")
                            else:
                                category = G_filtrado.nodes[drug]['atc_category']
                                color = G_filtrado.nodes[drug]['color']
                                st.markdown(
                                    f"<span style='color:{color};'>●</span> {drug}",
                                    unsafe_allow_html=True
                                )
                                st.caption(f"({category})")
    
    # Información
    st.markdown("---")
    st.markdown("**About the data:**")
    st.markdown("""
    - **Interactions are unique**: If the same drug pair appears multiple times in the CSV, 
      only one interaction is shown in the graph
    - **Directed graph**: Arrows show the direction of interactions (Drug A → Drug B)
    - **ATC Categories**: Colors represent Anatomical Therapeutic Chemical classification
    - **Interaction Types**: Type 1 and Type 2 based on the 'Y' column in the CSV
    """)

if __name__ == "__main__":
    main()




