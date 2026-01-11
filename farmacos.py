# app.py
import streamlit as st
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import time
import plotly.graph_objects as go
import plotly.express as px

# Configuración de la página
st.set_page_config(
    page_title="Drug Interaction Network",
    page_icon="💊",
    layout="wide",
    initial_sidebar_state="expanded"
)

# CSS personalizado
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        color: #1E3A8A;
        text-align: center;
        margin-bottom: 1rem;
    }
    .stButton > button {
        width: 100%;
        background-color: #3B82F6;
        color: white;
    }
    .metric-card {
        background-color: #F8FAFC;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #3B82F6;
    }
    .drug-card {
        background-color: #F0F9FF;
        padding: 0.5rem;
        margin: 0.5rem 0;
        border-radius: 0.5rem;
        border: 1px solid #E0F2FE;
    }
    .dataframe tbody tr:hover {
        background-color: #f0f9ff;
    }
</style>
""", unsafe_allow_html=True)

# Título principal
st.markdown('<h1 class="main-header">💊 Drug Interaction Network Analyzer</h1>', unsafe_allow_html=True)

# Definiciones ATC
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

# Función para colores ATC
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

# Función para crear grafo optimizada - MODIFICADA para usar TODO el dataset
@st.cache_data(ttl=3600)
def crear_grafo_completo(_df, farmaco_objetivo=None, max_interactions=None):
    """Versión optimizada para Streamlit que usa TODO el dataset"""
    G = nx.DiGraph()
    
    # Determinar qué datos usar
    if farmaco_objetivo:
        # Filtrar solo interacciones del fármaco objetivo
        mask = (_df['Common_name_x'].str.contains(farmaco_objetivo, case=False, na=False) | 
                _df['Common_name_y'].str.contains(farmaco_objetivo, case=False, na=False))
        filtered_df = _df[mask]
        
        # Si hay límite, aplicarlo solo para visualización
        if max_interactions and len(filtered_df) > max_interactions:
            filtered_df = filtered_df.head(max_interactions)
    else:
        # Para grafo completo, usar TODO el dataset
        filtered_df = _df
    
    st.info(f"Processing {len(filtered_df):,} interactions...")
    
    # Progreso
    progress_bar = st.progress(0)
    total_rows = len(filtered_df)
    
    for idx, row in filtered_df.iterrows():
        drug1 = row['Common_name_x']
        drug2 = row['Common_name_y']
        interaction_type = row['Y']
        
        atc1 = row['atc_code_x'] if pd.notna(row['atc_code_x']) else "No ATC"
        num_atc1 = row['num_atc_x'] if 'num_atc_x' in row else 0
        color1, category1 = get_atc_color_and_category(atc1, num_atc1)
        
        atc2 = row['atc_code_y'] if pd.notna(row['atc_code_y']) else "No ATC"
        num_atc2 = row['num_atc_y'] if 'num_atc_y' in row else 0
        color2, category2 = get_atc_color_and_category(atc2, num_atc2)
        
        if not G.has_node(drug1):
            G.add_node(drug1, 
                      atc_code=atc1,
                      atc_category=category1,
                      color=color1,
                      num_atc=num_atc1)
        
        if not G.has_node(drug2):
            G.add_node(drug2, 
                      atc_code=atc2,
                      atc_category=category2,
                      color=color2,
                      num_atc=num_atc2)
        
        if not G.has_edge(drug1, drug2):
            G.add_edge(drug1, drug2, 
                      interaction_type=interaction_type)
        
        # Actualizar progreso cada 1000 filas
        if idx % 1000 == 0:
            progress_bar.progress(min(idx / total_rows, 1.0))
    
    progress_bar.progress(1.0)
    time.sleep(0.5)
    progress_bar.empty()
    
    return G

# Función para visualización con Plotly optimizada para grandes datasets
def visualizar_con_plotly(G, farmaco_principal=None):
    """Visualización interactiva con Plotly optimizada"""
    
    if len(G.nodes()) == 0:
        fig = go.Figure()
        fig.add_annotation(
            text="No drugs to display",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=20)
        )
        return fig
    
    # Para grafos grandes, usar layout más eficiente
    if len(G.nodes()) > 100:
        # Layout más eficiente para grafos grandes
        pos = nx.spring_layout(G, k=3/np.sqrt(len(G.nodes())), iterations=30, seed=42)
        
        # Simplificar visualización para muchos nodos
        edge_trace = None  # Omitir aristas para mejor rendimiento
        
        # Solo mostrar nodos importantes
        if farmaco_principal and farmaco_principal in G.nodes():
            # Incluir fármaco principal y sus vecinos directos
            important_nodes = [farmaco_principal] + list(G.neighbors(farmaco_principal))
            G_filtered = G.subgraph(important_nodes)
            pos_filtered = {node: pos[node] for node in important_nodes if node in pos}
        else:
            # Para grafo completo, usar todos los nodos
            G_filtered = G
            pos_filtered = pos
    else:
        # Para grafos pequeños, visualización completa
        if farmaco_principal and farmaco_principal in G.nodes():
            pos = nx.spring_layout(G, k=1, iterations=50, seed=42)
        else:
            pos = nx.spring_layout(G, seed=42)
        
        G_filtered = G
        pos_filtered = pos
        
        # Preparar aristas para grafos pequeños
        edge_x = []
        edge_y = []
        
        for edge in G_filtered.edges():
            if edge[0] in pos_filtered and edge[1] in pos_filtered:
                x0, y0 = pos_filtered[edge[0]]
                x1, y1 = pos_filtered[edge[1]]
                edge_x.extend([x0, x1, None])
                edge_y.extend([y0, y1, None])
        
        if edge_x and edge_y:
            edge_trace = go.Scatter(
                x=edge_x, y=edge_y,
                mode='lines',
                line=dict(width=0.5, color='#888'),
                hoverinfo='none'
            )
        else:
            edge_trace = None
    
    # Preparar datos de nodos
    node_x = []
    node_y = []
    node_texts = []
    node_colors = []
    node_sizes = []
    
    for node in G_filtered.nodes():
        if node in pos_filtered:
            x, y = pos_filtered[node]
            node_x.append(x)
            node_y.append(y)
            
            node_data = G_filtered.nodes[node]
            node_texts.append(
                f"<b>{node}</b><br>"
                f"ATC: {node_data.get('atc_code', 'N/A')}<br>"
                f"Category: {node_data.get('atc_category', 'Unknown')}<br>"
                f"Connections: {G_filtered.degree(node)}"
            )
            
            node_colors.append(node_data.get('color', '#CCCCCC'))
            
            # Tamaño basado en grado de conexión
            if node == farmaco_principal:
                node_sizes.append(30)
            else:
                degree = G_filtered.degree(node)
                node_sizes.append(10 + min(degree * 2, 20))
    
    # Crear traza de nodos
    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text' if len(G_filtered.nodes()) < 50 else 'markers',
        text=[node[:12] + ('...' if len(node) > 12 else '') for node in G_filtered.nodes()] 
              if len(G_filtered.nodes()) < 50 else [],
        textposition="top center",
        hovertext=node_texts,
        hoverinfo='text',
        marker=dict(
            color=node_colors,
            size=node_sizes,
            line=dict(width=1, color='darkgray')
        )
    )
    
    # Crear figura
    if edge_trace:
        fig = go.Figure(data=[edge_trace, node_trace])
    else:
        fig = go.Figure(data=[node_trace])
    
    # Configurar layout
    title = f"Drug Interaction Network"
    if farmaco_principal:
        title = f"Interactions for: {farmaco_principal}"
    
    fig.update_layout(
        title=f"{title}<br>{len(G.nodes())} drugs, {len(G.edges())} interactions",
        showlegend=False,
        hovermode='closest',
        margin=dict(b=20, l=5, r=5, t=60),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        height=700,
        plot_bgcolor='white',
        hoverlabel=dict(
            bgcolor="white",
            font_size=12,
            font_family="Arial"
        )
    )
    
    return fig

# Función principal de carga de datos - CARGA COMPLETA
@st.cache_data(ttl=3600)
def load_complete_data(file_path=None, use_cols=None):
    """Carga completa de datos sin límites"""
    try:
        st.info("Loading complete dataset... This may take a moment for large files.")
        
        if hasattr(file_path, 'name'):  # Si es archivo subido
            if use_cols:
                df = pd.read_csv(file_path, usecols=use_cols)
            else:
                df = pd.read_csv(file_path)
        else:  # Si es archivo local
            if use_cols:
                df = pd.read_csv(file_path, usecols=use_cols)
            else:
                df = pd.read_csv(file_path)
        
        st.success(f"✓ Dataset loaded successfully: {len(df):,} rows")
        return df
    except Exception as e:
        st.error(f"Error loading data: {e}")
        return None

# Sidebar
st.sidebar.header("⚙️ Configuration")
uploaded_file = st.sidebar.file_uploader("Upload CSV file", type=['csv'])

if uploaded_file is not None:
    file_path = uploaded_file
    st.sidebar.success("Using uploaded file")
else:
    file_path = "DDIBUENO.csv"
    st.sidebar.info("Using default file path")

st.sidebar.header("🔍 Filters")

# Opción para limitar visualización (no datos)
max_interactions = st.sidebar.slider(
    "Max interactions to display (0 = all)",
    min_value=0,
    max_value=50000,
    value=0,
    step=1000,
    help="Limit display for better performance. 0 shows all interactions."
)

visualization_method = st.sidebar.radio(
    "Visualization method",
    ["Plotly (Interactive)", "Matplotlib (Static)"],
    index=0
)

show_stats = st.sidebar.checkbox("Show detailed statistics", value=True)

# Opción para procesamiento completo
process_all = st.sidebar.checkbox("Process entire dataset", value=True, 
                                 help="Uncheck to limit processing for testing")

# Contenido principal
tab1, tab2, tab3 = st.tabs(["Network Visualization", "Data Explorer", "About"])

with tab1:
    col1, col2 = st.columns([3, 1])
    
    with col2:
        st.subheader("🔎 Drug Search")
        
        @st.cache_data
        def get_drug_list(_df):
            all_drugs = set(_df['Common_name_x'].tolist() + _df['Common_name_y'].tolist())
            return sorted(list(all_drugs))
        
        with st.spinner("Loading complete dataset..."):
            use_cols = ['Common_name_x', 'Common_name_y', 'Y', 'atc_code_x', 'atc_code_y']
            try:
                if hasattr(file_path, 'name'):
                    file_content = file_path.getvalue().decode('utf-8')
                    first_lines = file_content.split('\n')[:2]
                    first_line = first_lines[0] if first_lines else ''
                    if 'num_atc_x' in first_line:
                        use_cols.extend(['num_atc_x', 'num_atc_y'])
                else:
                    df_sample = pd.read_csv(file_path, nrows=1)
                    if 'num_atc_x' in df_sample.columns:
                        use_cols.extend(['num_atc_x', 'num_atc_y'])
            except:
                pass
            
            # Cargar TODO el dataset
            df = load_complete_data(file_path, use_cols)
        
        if df is not None:
            drug_list = get_drug_list(df)
            
            st.success(f"✅ Loaded {len(df):,} total interactions")
            st.info(f"📊 {len(drug_list):,} unique drugs available")
            
            # Buscador mejorado
            search_term = st.text_input("Search drug:", 
                                       placeholder="Type drug name...",
                                       help="Search across all drugs")
            
            # Filtrar lista basada en búsqueda
            if search_term:
                filtered_drugs = [drug for drug in drug_list 
                                if search_term.lower() in drug.lower()]
            else:
                filtered_drugs = drug_list[:1000]  # Mostrar primeros 1000
            
            target_drug = st.selectbox(
                "Select target drug:",
                options=["Complete Network (All Drugs)"] + filtered_drugs,
                index=0,
                help="Select 'Complete Network' to see all interactions"
            )
            
            # Opciones adicionales
            col_opt1, col_opt2 = st.columns(2)
            with col_opt1:
                show_direct_only = st.checkbox("Direct connections only", value=False,
                                             help="Show only direct interactions")
            with col_opt2:
                min_degree = st.number_input("Min connections", min_value=1, value=1,
                                           help="Minimum number of connections to display")
            
            if st.button("🚀 Generate Complete Network", type="primary", use_container_width=True):
                with st.spinner(f"Building network from {len(df):,} interactions..."):
                    selected_drug = None
                    if target_drug != "Complete Network (All Drugs)":
                        selected_drug = target_drug
                    
                    # Determinar límite
                    limit = max_interactions if max_interactions > 0 else None
                    
                    # Crear grafo COMPLETO o filtrado
                    if selected_drug:
                        G = crear_grafo_completo(df, selected_drug, limit)
                        farmaco_principal = selected_drug
                    else:
                        G = crear_grafo_completo(df, None, limit)
                        farmaco_principal = None
                    
                    # Filtrar por grado mínimo si se especifica
                    if min_degree > 1:
                        nodes_to_remove = [node for node in G.nodes() 
                                         if G.degree(node) < min_degree]
                        G.remove_nodes_from(nodes_to_remove)
                    
                    # Guardar en estado de sesión
                    st.session_state['graph'] = G
                    st.session_state['main_drug'] = farmaco_principal
                    st.session_state['df'] = df
                    
                    # Estadísticas DETALLADAS
                    if show_stats:
                        st.subheader("📈 Network Statistics")
                        
                        col1_stat, col2_stat, col3_stat, col4_stat = st.columns(4)
                        
                        with col1_stat:
                            st.metric("Total Drugs", f"{len(G.nodes()):,}")
                        with col2_stat:
                            st.metric("Total Interactions", f"{len(G.edges()):,}")
                        with col3_stat:
                            if farmaco_principal:
                                degree = G.degree(farmaco_principal) if farmaco_principal in G else 0
                                st.metric("Target Drug Connections", degree)
                        with col4_stat:
                            density = nx.density(G) if len(G.nodes()) > 1 else 0
                            st.metric("Network Density", f"{density:.4f}")
                        
                        # Más estadísticas
                        if len(G.nodes()) > 0:
                            degrees = [deg for _, deg in G.degree()]
                            avg_degree = np.mean(degrees) if degrees else 0
                            max_degree = max(degrees) if degrees else 0
                            
                            col5_stat, col6_stat = st.columns(2)
                            with col5_stat:
                                st.metric("Average Connections", f"{avg_degree:.1f}")
                            with col6_stat:
                                st.metric("Max Connections", max_degree)
                    
                    # Visualización
                    if visualization_method == "Plotly (Interactive)":
                        fig = visualizar_con_plotly(G, farmaco_principal)
                        st.plotly_chart(fig, use_container_width=True)
                    else:
                        # Versión Matplotlib simplificada para grandes grafos
                        fig, ax = plt.subplots(figsize=(14, 10))
                        
                        if len(G.nodes()) > 100:
                            # Para grafos grandes, mostrar solo subgrafo
                            if farmaco_principal and farmaco_principal in G.nodes():
                                # Mostrar fármaco principal y sus vecinos
                                important_nodes = [farmaco_principal] + list(G.neighbors(farmaco_principal))[:50]
                                G_display = G.subgraph(important_nodes)
                            else:
                                # Mostrar nodos con mayor grado
                                degrees = dict(G.degree())
                                top_nodes = sorted(degrees.items(), key=lambda x: x[1], reverse=True)[:100]
                                G_display = G.subgraph([node for node, _ in top_nodes])
                        else:
                            G_display = G
                        
                        pos = nx.spring_layout(G_display, seed=42)
                        
                        node_colors = [G_display.nodes[node]['color'] for node in G_display.nodes()]
                        node_sizes = []
                        for node in G_display.nodes():
                            if node == farmaco_principal:
                                node_sizes.append(1200)
                            else:
                                node_sizes.append(300 + min(G_display.degree(node) * 50, 1000))
                        
                        nx.draw(G_display, pos, 
                               node_color=node_colors,
                               node_size=node_sizes,
                               edge_color='gray',
                               width=0.5,
                               with_labels=True,
                               font_size=8,
                               font_weight='bold',
                               ax=ax)
                        
                        title = f"Drug Interaction Network"
                        if farmaco_principal:
                            title = f"Interactions for: {farmaco_principal}"
                        
                        ax.set_title(f"{title}\n{len(G_display.nodes())} of {len(G.nodes())} drugs shown", 
                                   fontsize=14, pad=20)
                        ax.axis('off')
                        
                        # Añadir leyenda de tamaño
                        ax.text(1.05, 0.5, "Node size ≈ Connections", 
                               transform=ax.transAxes, fontsize=10,
                               verticalalignment='center')
                        
                        st.pyplot(fig)
                    
                    # Información detallada en expansores
                    with st.expander("📊 Detailed Network Analysis", expanded=False):
                        tab_analysis, tab_nodes, tab_edges = st.tabs(["Analysis", "Nodes", "Edges"])
                        
                        with tab_analysis:
                            # Distribución de grados
                            if len(G.nodes()) > 0:
                                degrees = [deg for _, deg in G.degree()]
                                
                                fig_hist = px.histogram(x=degrees, 
                                                       title="Distribution of Connections",
                                                       labels={'x': 'Number of Connections', 'y': 'Count'},
                                                       nbins=30)
                                fig_hist.update_layout(showlegend=False)
                                st.plotly_chart(fig_hist, use_container_width=True)
                            
                            # Top fármacos más conectados
                            if len(G.nodes()) > 0:
                                degrees = dict(G.degree())
                                top_drugs = sorted(degrees.items(), key=lambda x: x[1], reverse=True)[:20]
                                
                                top_df = pd.DataFrame(top_drugs, columns=['Drug', 'Connections'])
                                fig_bar = px.bar(top_df, x='Drug', y='Connections',
                                               title="Top 20 Most Connected Drugs",
                                               color='Connections')
                                fig_bar.update_layout(xaxis_tickangle=-45)
                                st.plotly_chart(fig_bar, use_container_width=True)
                        
                        with tab_nodes:
                            st.subheader(f"All Drugs in Network ({len(G.nodes()):,})")
                            
                            # Búsqueda en nodos
                            node_search = st.text_input("Search in network drugs:", 
                                                       placeholder="Filter drugs...")
                            
                            # Mostrar nodos paginados
                            nodes_list = list(G.nodes())
                            if node_search:
                                nodes_list = [node for node in nodes_list 
                                            if node_search.lower() in node.lower()]
                            
                            items_per_page = 50
                            total_pages = (len(nodes_list) // items_per_page) + 1
                            
                            if total_pages > 1:
                                page = st.number_input("Page", min_value=1, max_value=total_pages, value=1)
                                start_idx = (page - 1) * items_per_page
                                end_idx = min(start_idx + items_per_page, len(nodes_list))
                            else:
                                start_idx, end_idx = 0, len(nodes_list)
                            
                            for node in nodes_list[start_idx:end_idx]:
                                node_data = G.nodes[node]
                                degree = G.degree(node)
                                
                                col_node1, col_node2 = st.columns([3, 1])
                                with col_node1:
                                    st.markdown(f"""
                                    <div class="drug-card">
                                    <b>{node}</b><br>
                                    <small>ATC: {node_data.get('atc_code', 'N/A')} | 
                                    Category: {node_data.get('atc_category', 'Unknown')} | 
                                    Connections: {degree}</small>
                                    </div>
                                    """, unsafe_allow_html=True)
                                with col_node2:
                                    if st.button("📊", key=f"btn_{node}", help=f"Show {node} network"):
                                        st.session_state['selected_drug'] = node
                                        st.rerun()
                        
                        with tab_edges:
                            st.subheader(f"All Interactions ({len(G.edges()):,})")
                            
                            # Mostrar algunas interacciones
                            edges_list = list(G.edges(data=True))[:100]  # Mostrar primeras 100
                            for u, v, data in edges_list:
                                interaction_type = data.get('interaction_type', 'Unknown')
                                st.markdown(f"**{u}** → **{v}** (Type: {interaction_type})")

with tab2:
    st.subheader("📁 Complete Data Explorer")
    
    if 'df' in st.session_state:
        df = st.session_state['df']
        
        st.success(f"Exploring dataset with {len(df):,} rows")
        
        # Estadísticas rápidas
        col_quick1, col_quick2, col_quick3, col_quick4 = st.columns(4)
        with col_quick1:
            st.metric("Total Rows", f"{len(df):,}")
        with col_quick2:
            unique_drugs = len(set(df['Common_name_x'].tolist() + df['Common_name_y'].tolist()))
            st.metric("Unique Drugs", f"{unique_drugs:,}")
        with col_quick3:
            type1_count = len(df[df['Y'] == 1])
            st.metric("Type 1 Interactions", f"{type1_count:,}")
        with col_quick4:
            type2_count = len(df[df['Y'] == 2])
            st.metric("Type 2 Interactions", f"{type2_count:,}")
        
        # Filtros avanzados
        st.subheader("🔍 Advanced Filters")
        
        col_filt1, col_filt2, col_filt3 = st.columns(3)
        
        with col_filt1:
            selected_atc = st.selectbox(
                "ATC Category:",
                ["All Categories"] + list(ATC_CATEGORIES.values()),
                key="filter_atc_main"
            )
        
        with col_filt2:
            interaction_type = st.selectbox(
                "Interaction Type:",
                ["All Types", "Type 1", "Type 2"],
                key="filter_type_main"
            )
        
        with col_filt3:
            drug_filter = st.text_input("Filter by drug name:", 
                                       placeholder="Any drug name...")
        
        # Aplicar filtros
        if st.button("Apply Filters to Entire Dataset", type="primary"):
            filtered_df = df.copy()
            
            if selected_atc != "All Categories":
                atc_code = None
                for code, category in ATC_CATEGORIES.items():
                    if category == selected_atc:
                        atc_code = code
                        break
                
                if atc_code:
                    mask = (filtered_df['atc_code_x'].astype(str).str.startswith(atc_code)) | \
                           (filtered_df['atc_code_y'].astype(str).str.startswith(atc_code))
                    filtered_df = filtered_df[mask]
            
            if interaction_type != "All Types":
                type_value = 1 if interaction_type == "Type 1" else 2
                filtered_df = filtered_df[filtered_df['Y'] == type_value]
            
            if drug_filter:
                mask = (filtered_df['Common_name_x'].str.contains(drug_filter, case=False, na=False) | 
                       filtered_df['Common_name_y'].str.contains(drug_filter, case=False, na=False))
                filtered_df = filtered_df[mask]
            
            st.success(f"Filtered to {len(filtered_df):,} rows")
            
            # Mostrar datos filtrados con paginación
            st.subheader("Filtered Data Preview")
            
            total_filtered = len(filtered_df)
            items_per_page = 100
            total_pages = (total_filtered // items_per_page) + 1
            
            if total_pages > 1:
                page = st.number_input("Page", min_value=1, max_value=total_pages, value=1, 
                                     key="page_filtered")
                start_idx = (page - 1) * items_per_page
                end_idx = min(start_idx + items_per_page, total_filtered)
                
                st.caption(f"Showing rows {start_idx+1:,} to {end_idx:,} of {total_filtered:,}")
            else:
                start_idx, end_idx = 0, min(100, total_filtered)
            
            st.dataframe(filtered_df.iloc[start_idx:end_idx], use_container_width=True)
            
            # Opciones de exportación
            st.subheader("📥 Export Options")
            
            col_exp1, col_exp2 = st.columns(2)
            
            with col_exp1:
                csv = filtered_df.to_csv(index=False).encode('utf-8')
                st.download_button(
                    label=f"Download Filtered Data ({len(filtered_df):,} rows)",
                    data=csv,
                    file_name=f"filtered_interactions_{len(filtered_df)}_rows.csv",
                    mime="text/csv",
                    use_container_width=True
                )
            
            with col_exp2:
                csv_all = df.to_csv(index=False).encode('utf-8')
                st.download_button(
                    label=f"Download Complete Dataset ({len(df):,} rows)",
                    data=csv_all,
                    file_name=f"complete_dataset_{len(df)}_rows.csv",
                    mime="text/csv",
                    use_container_width=True
                )
    else:
        st.info("Please load and generate a network in the 'Network Visualization' tab first.")

with tab3:
    st.subheader("ℹ️ About This Tool")
    
    st.markdown("""
    ## Drug Interaction Network Analyzer
    
    This tool visualizes drug-drug interactions from your dataset.
    
    ### ATC Categories:
    """)
    
    atc_table_data = []
    for code, category in ATC_CATEGORIES.items():
        if code in ATC_COLORS:
            atc_table_data.append({
                "Code": code,
                "Category": category,
                "Color": ATC_COLORS[code]
            })
    
    atc_table = pd.DataFrame(atc_table_data)
    
    def color_cell(color):
        return f'<div style="background-color:{color}; width:20px; height:20px; border-radius:3px;"></div>'
    
    atc_table['Color Sample'] = atc_table['Color'].apply(color_cell)
    
    st.markdown(atc_table[['Code', 'Category', 'Color Sample']].to_html(
        escape=False, 
        index=False
    ), unsafe_allow_html=True)
    
    st.markdown("""
    ### How to Use:
    1. Upload your CSV file or use the default
    2. Select a target drug or view all drugs
    3. Adjust filters in the sidebar
    4. Explore the interactive network
    5. Use the Data Explorer tab for detailed analysis
    """)


# Footer
st.markdown("---")

st.markdown(
    "<div style='text-align: center; color: gray;'>"
    "💊 Drug Interaction Network Analyzer | Complete Dataset Analysis | Built with Streamlit"
    "</div>",
    unsafe_allow_html=True
)

# Inicialización
if 'graph' not in st.session_state:
    st.session_state['graph'] = None
if 'main_drug' not in st.session_state:
    st.session_state['main_drug'] = None
if 'df' not in st.session_state:
    st.session_state['df'] = None
if 'selected_drug' not in st.session_state:
    st.session_state['selected_drug'] = None

# Manejar selección de fármaco desde la lista
if 'selected_drug' in st.session_state and st.session_state['selected_drug']:
    st.info(f"Selected drug: {st.session_state['selected_drug']}")
    # Aquí podrías agregar lógica para regenerar el grafo con este fármaco
