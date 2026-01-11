# app_streamlit_completo.py
import streamlit as st
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import io
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
</style>
""", unsafe_allow_html=True)

# Título principal
st.markdown('<h1 class="main-header">💊 Drug Interaction Network Analyzer</h1>', unsafe_allow_html=True)

# Definiciones ATC (mantener tu código original)
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

# Función para crear grafo optimizada
@st.cache_data(ttl=3600)
def crear_grafo_optimizado(_df, farmaco_objetivo=None, max_interactions=100):
    """Versión optimizada para Streamlit"""
    G = nx.DiGraph()
    
    # Filtrar primero si hay un fármaco objetivo
    if farmaco_objetivo:
        mask = (_df['Common_name_x'].str.contains(farmaco_objetivo, case=False, na=False) | 
                _df['Common_name_y'].str.contains(farmaco_objetivo, case=False, na=False))
        filtered_df = _df[mask].head(max_interactions)
    else:
        filtered_df = _df.head(max_interactions)
    
    for _, row in filtered_df.iterrows():
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
    
    return G

# Función para visualización con Plotly (más rápida)
def visualizar_con_plotly(G, farmaco_principal=None):
    """Visualización interactiva con Plotly"""
    
    if len(G.nodes()) == 0:
        return go.Figure()
    
    # Posiciones de los nodos
    if farmaco_principal and farmaco_principal in G.nodes():
        pos = nx.spring_layout(G, k=1, iterations=50, seed=42)
    else:
        pos = nx.spring_layout(G, seed=42)
    
    # Preparar datos de aristas
    edge_x = []
    edge_y = []
    edge_texts = []
    
    for edge in G.edges(data=True):
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
        
        interaction_type = edge[2].get('interaction_type', 'Unknown')
        edge_texts.append(f"{edge[0]} → {edge[1]}<br>Type: {interaction_type}")
    
    # Preparar datos de nodos
    node_x = []
    node_y = []
    node_texts = []
    node_colors = []
    node_sizes = []
    
    for node in G.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        
        node_data = G.nodes[node]
        node_texts.append(
            f"<b>{node}</b><br>"
            f"ATC: {node_data.get('atc_code', 'N/A')}<br>"
            f"Category: {node_data.get('atc_category', 'Unknown')}"
        )
        
        node_colors.append(node_data.get('color', '#CCCCCC'))
        
        # Tamaño más grande para fármaco principal
        if node == farmaco_principal:
            node_sizes.append(25)
        else:
            node_sizes.append(15)
    
    # Crear trazas
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        mode='lines',
        line=dict(width=1, color='#888'),
        hoverinfo='none'
    )
    
    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text',
        text=[node[:15] + ('...' if len(node) > 15 else '') for node in G.nodes()],
        textposition="top center",
        hovertext=node_texts,
        hoverinfo='text',
        marker=dict(
            color=node_colors,
            size=node_sizes,
            line=dict(width=2, color='darkgray')
        )
    )
    
    # Crear figura
    fig = go.Figure(data=[edge_trace, node_trace])
    
    # Configurar layout
    fig.update_layout(
        title=f"Drug Interaction Network ({len(G.nodes())} drugs, {len(G.edges())} interactions)",
        showlegend=False,
        hovermode='closest',
        margin=dict(b=20, l=5, r=5, t=40),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        height=600,
        plot_bgcolor='white'
    )
    
    return fig

# Función para visualización con Matplotlib (para opción alternativa)
def visualizar_con_matplotlib(G, farmaco_principal=None):
    """Visualización con Matplotlib (opcional)"""
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Posiciones
    if farmaco_principal and farmaco_principal in G.nodes():
        pos = nx.spring_layout(G, k=1, iterations=50, seed=42)
    else:
        pos = nx.spring_layout(G, seed=42)
    
    # Colores y tamaños
    node_colors = []
    node_sizes = []
    
    for node in G.nodes():
        node_colors.append(G.nodes[node]['color'])
        if node == farmaco_principal:
            node_sizes.append(800)
        else:
            node_sizes.append(400)
    
    # Dibujar
    nx.draw_networkx_nodes(G, pos, 
                          node_color=node_colors,
                          node_size=node_sizes,
                          alpha=0.9,
                          ax=ax)
    
    nx.draw_networkx_edges(G, pos,
                          width=1,
                          alpha=0.6,
                          edge_color='gray',
                          arrows=True,
                          arrowsize=10,
                          ax=ax)
    
    # Etiquetas
    labels = {}
    for node in G.nodes():
        if len(node) > 20:
            labels[node] = node[:17] + "..."
        else:
            labels[node] = node
    
    nx.draw_networkx_labels(G, pos, labels, font_size=8, ax=ax)
    
    ax.set_title(f"Drug Network: {len(G.nodes())} nodes, {len(G.edges())} edges")
    ax.axis('off')
    
    return fig

# Función principal de carga de datos
@st.cache_data(ttl=3600)
def load_data(file_path=None, sample_size=None, use_cols=None):
    """Carga optimizada de datos"""
    try:
        if sample_size:
            # Leer solo una muestra
            df = pd.read_csv(file_path, nrows=sample_size, usecols=use_cols)
        else:
            # Leer todo pero solo columnas necesarias
            if use_cols:
                df = pd.read_csv(file_path, usecols=use_cols)
            else:
                df = pd.read_csv(file_path)
        
        return df
    except Exception as e:
        st.error(f"Error loading data: {e}")
        return None

# Sidebar - Configuración
st.sidebar.header("⚙️ Configuration")

# Opción de carga de archivo
uploaded_file = st.sidebar.file_uploader("Upload CSV file", type=['csv'])

if uploaded_file is not None:
    # Usar archivo subido
    file_path = uploaded_file
    st.sidebar.success("Using uploaded file")
else:
    # Usar ruta local
    file_path = "DDIBUENO.csv"
    st.sidebar.info("Using default file path")

# Filtros en sidebar
st.sidebar.header("🔍 Filters")

# Limitar tamaño de datos para prueba
sample_size = st.sidebar.slider(
    "Sample size (rows)",
    min_value=1000,
    max_value=50000,
    value=10000,
    step=1000,
    help="Reduce for faster testing"
)

max_interactions = st.sidebar.slider(
    "Max interactions per drug",
    min_value=10,
    max_value=200,
    value=50,
    step=10
)

# Opciones de visualización
visualization_method = st.sidebar.radio(
    "Visualization method",
    ["Plotly (Interactive)", "Matplotlib (Static)"],
    index=0
)

show_stats = st.sidebar.checkbox("Show detailed statistics", value=True)

# Contenido principal
tab1, tab2, tab3 = st.tabs(["Network Visualization", "Data Explorer", "About"])

with tab1:
    col1, col2 = st.columns([3, 1])
    
    with col2:
        st.subheader("🔎 Drug Search")
        
        # Cargar lista de fármacos
        @st.cache_data
        def get_drug_list(_df):
            all_drugs = set(_df['Common_name_x'].tolist() + _df['Common_name_y'].tolist())
            return sorted(list(all_drugs))
        
        # Cargar datos
        with st.spinner("Loading data..."):
            # Columnas necesarias
            use_cols = ['Common_name_x', 'Common_name_y', 'Y', 'atc_code_x', 'atc_code_y']
            if 'num_atc_x' in pd.read_csv(file_path, nrows=1).columns:
                use_cols.extend(['num_atc_x', 'num_atc_y'])
            
            df = load_data(file_path, sample_size, use_cols)
        
        if df is not None:
            drug_list = get_drug_list(df)
            
            # Selector de fármaco
            target_drug = st.selectbox(
                "Select target drug:",
                options=["All Drugs"] + drug_list[:500],  # Limitar a 500 para rendimiento
                index=0
            )
            
            # Búsqueda manual
            search_term = st.text_input("Or search by name:")
            
            # Botón para generar
            if st.button("Generate Network", type="primary"):
                with st.spinner("Building network..."):
                    # Determinar fármaco objetivo
                    selected_drug = None
                    if target_drug != "All Drugs":
                        selected_drug = target_drug
                    elif search_term:
                        # Buscar fármaco que contenga el término
                        for drug in drug_list:
                            if search_term.lower() in drug.lower():
                                selected_drug = drug
                                break
                    
                    # Crear grafo
                    if selected_drug:
                        G = crear_grafo_optimizado(df, selected_drug, max_interactions)
                        farmaco_principal = selected_drug
                        st.session_state['graph'] = G
                        st.session_state['main_drug'] = farmaco_principal
                    else:
                        G = crear_grafo_optimizado(df, None, max_interactions)
                        farmaco_principal = None
                        st.session_state['graph'] = G
                        st.session_state['main_drug'] = None
                    
                    # Estadísticas
                    if show_stats:
                        col_metrics = st.columns(3)
                        with col_metrics[0]:
                            st.metric("Total Drugs", len(G.nodes()))
                        with col_metrics[1]:
                            st.metric("Total Interactions", len(G.edges()))
                        with col_metrics[2]:
                            if farmaco_principal:
                                st.metric("Target Drug", farmaco_principal)
                    
                    # Visualización
                    if visualization_method == "Plotly (Interactive)":
                        fig = visualizar_con_plotly(G, farmaco_principal)
                        st.plotly_chart(fig, use_container_width=True)
                    else:
                        fig = visualizar_con_matplotlib(G, farmaco_principal)
                        st.pyplot(fig)
                    
                    # Información detallada
                    with st.expander("📊 Network Details"):
                        col_info1, col_info2 = st.columns(2)
                        
                        with col_info1:
                            st.subheader("Drugs in Network")
                            for node in list(G.nodes())[:20]:  # Mostrar primeros 20
                                with st.container():
                                    node_data = G.nodes[node]
                                    st.markdown(f"""
                                    <div class="drug-card">
                                    <b>{node}</b><br>
                                    ATC: {node_data.get('atc_code', 'N/A')}<br>
                                    Category: {node_data.get('atc_category', 'Unknown')}
                                    </div>
                                    """, unsafe_allow_html=True)
                        
                        with col_info2:
                            st.subheader("ATC Category Distribution")
                            categories = Counter()
                            for node in G.nodes():
                                categories[G.nodes[node]['atc_category']] += 1
                            
                            # Gráfico de torta
                            if categories:
                                fig_pie = px.pie(
                                    values=list(categories.values()),
                                    names=list(categories.keys()),
                                    title="Drugs by ATC Category"
                                )
                                st.plotly_chart(fig_pie, use_container_width=True)

with tab2:
    st.subheader("📁 Data Explorer")
    
    if 'df' in locals():
        # Mostrar vista previa de datos
        st.dataframe(df.head(100), use_container_width=True)
        
        # Estadísticas básicas
        col_stats1, col_stats2, col_stats3 = st.columns(3)
        
        with col_stats1:
            st.metric("Total Rows", len(df))
        
        with col_stats2:
            unique_drugs = len(set(df['Common_name_x'].tolist() + df['Common_name_y'].tolist()))
            st.metric("Unique Drugs", unique_drugs)
        
        with col_stats3:
            st.metric("Total Interactions", len(df))
        
        # Filtros interactivos
        st.subheader("Filter Data")
        
        col_filter1, col_filter2 = st.columns(2)
        
        with col_filter1:
            selected_atc = st.selectbox(
                "Filter by ATC Category:",
                ["All"] + list(ATC_CATEGORIES.values())
            )
        
        with col_filter2:
            interaction_type = st.selectbox(
                "Filter by Interaction Type:",
                ["All", "Type 1", "Type 2"]
            )
        
        if st.button("Apply Filters"):
            filtered_df = df.copy()
            
            if selected_atc != "All":
                # Convertir categoría a código ATC
                atc_code = None
                for code, category in ATC_CATEGORIES.items():
                    if category == selected_atc:
                        atc_code = code
                        break
                
                if atc_code:
                    mask = (filtered_df['atc_code_x'].astype(str).str.startswith(atc_code)) | \
                           (filtered_df['atc_code_y'].astype(str).str.startswith(atc_code))
                    filtered_df = filtered_df[mask]
            
            if interaction_type != "All":
                type_value = 1 if interaction_type == "Type 1" else 2
                filtered_df = filtered_df[filtered_df['Y'] == type_value]
            
            st.dataframe(filtered_df.head(50), use_container_width=True)

with tab3:
    st.subheader("ℹ️ About This Tool")
    
    st.markdown("""
    ## Drug Interaction Network Analyzer
    
    This tool visualizes drug-drug interactions from your dataset.
    
    ### Features:
    - **Interactive Network Visualization**: Explore drug interactions
    - **ATC Category Filtering**: Filter by therapeutic categories
    - **Search Functionality**: Find specific drugs and their interactions
    - **Statistical Analysis**: View distribution and metrics
    
    ### ATC Categories:
    """)
    
    # Mostrar tabla de categorías ATC
    atc_table = pd.DataFrame([
        {"Code": code, "Category": category, "Color": color}
        for code, category in ATC_CATEGORIES.items()
        if code in ATC_COLORS
    ])
    
    # Añadir muestra de color
    def color_cell(color):
        return f'<div style="background-color:{color}; width:20px; height:20px; border-radius:3px;"></div>'
    
    atc_table['Color Sample'] = atc_table['Color'].apply(color_cell)
    
    st.markdown(atc_table.to_html(escape=False, index=False), unsafe_allow_html=True)
    
    st.markdown("""
    ### How to Use:
    1. Upload your CSV file or use the default
    2. Select a target drug or view all drugs
    3. Adjust filters in the sidebar
    4. Explore the interactive network
    5. Use the Data Explorer tab for detailed analysis
    
    ### Tips for Large Datasets:
    - Start with a small sample size
    - Use the 'Max interactions' slider to limit results
    - Plotly visualization is faster for large networks
    """)

# Footer
st.markdown("---")
st.markdown(
    "<div style='text-align: center; color: gray;'>"
    "Drug Interaction Network Analyzer | Built with Streamlit"
    "</div>",
    unsafe_allow_html=True
)

# Inicialización
if 'graph' not in st.session_state:
    st.session_state['graph'] = None
if 'main_drug' not in st.session_state:
    st.session_state['main_drug'] = None

