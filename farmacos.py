# app.py - VERSIÓN SIMPLIFICADA
import streamlit as st
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import plotly.graph_objects as go

# Configuración de la página
st.set_page_config(
    page_title="Drug Interaction Network",
    page_icon="💊",
    layout="wide"
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
</style>
""", unsafe_allow_html=True)

# Título principal
st.markdown('<h1 class="main-header">💊 Drug Interaction Network</h1>', unsafe_allow_html=True)

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

# Función para crear grafo COMPLETO
@st.cache_data(ttl=3600)
def crear_grafo_completo(_df, farmaco_objetivo=None):
    """Crea grafo con TODO el dataset"""
    G = nx.DiGraph()
    
    # Determinar qué datos usar
    if farmaco_objetivo:
        # Filtrar solo interacciones del fármaco objetivo
        mask = (_df['Common_name_x'].str.contains(farmaco_objetivo, case=False, na=False) | 
                _df['Common_name_y'].str.contains(farmaco_objetivo, case=False, na=False))
        filtered_df = _df[mask]
    else:
        # Usar TODO el dataset
        filtered_df = _df
    
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

# Función para visualización con Plotly
def visualizar_con_plotly(G, farmaco_principal=None):
    """Visualización interactiva con Plotly"""
    
    if len(G.nodes()) == 0:
        fig = go.Figure()
        fig.add_annotation(
            text="No drugs to display",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=20)
        )
        return fig
    
    # Posiciones de los nodos
    pos = nx.spring_layout(G, seed=42)
    
    # Preparar aristas
    edge_x = []
    edge_y = []
    
    for edge in G.edges():
        if edge[0] in pos and edge[1] in pos:
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])
    
    # Preparar nodos
    node_x = []
    node_y = []
    node_texts = []
    node_colors = []
    node_sizes = []
    
    for node in G.nodes():
        if node in pos:
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
            
            # Tamaño basado en si es fármaco principal
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
    title = "Drug Interaction Network"
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
        plot_bgcolor='white'
    )
    
    return fig

# Función para carga de datos
@st.cache_data(ttl=3600)
def load_data(file_path=None):
    """Carga completa de datos"""
    try:
        if hasattr(file_path, 'name'):  # Si es archivo subido
            df = pd.read_csv(file_path)
        else:  # Si es archivo local
            df = pd.read_csv(file_path)
        
        return df
    except Exception as e:
        st.error(f"Error loading data: {e}")
        return None

# Sidebar - Configuración simple
st.sidebar.header("⚙️ Configuration")

uploaded_file = st.sidebar.file_uploader("Upload CSV file", type=['csv'])

if uploaded_file is not None:
    file_path = uploaded_file
    st.sidebar.success("Using uploaded file")
else:
    file_path = "DDIBUENO.csv"
    st.sidebar.info("Using default file path")

# Contenido principal - SOLO VISUALIZACIÓN
st.subheader("🔍 Drug Search and Visualization")

# Cargar datos
@st.cache_data
def get_drug_list(_df):
    all_drugs = set(_df['Common_name_x'].tolist() + _df['Common_name_y'].tolist())
    return sorted(list(all_drugs))

with st.spinner("Loading dataset..."):
    df = load_data(file_path)

if df is not None:
    drug_list = get_drug_list(df)
    
    col1, col2 = st.columns([3, 1])
    
    with col2:
        st.info(f"📊 {len(df):,} interactions")
        st.info(f"💊 {len(drug_list):,} unique drugs")
        
        # Buscador
        search_term = st.text_input("Search drug:", placeholder="Type drug name...")
        
        # Filtrar lista basada en búsqueda
        if search_term:
            filtered_drugs = [drug for drug in drug_list 
                            if search_term.lower() in drug.lower()]
        else:
            filtered_drugs = drug_list[:500]
        
        target_drug = st.selectbox(
            "Select drug:",
            options=["Show All Drugs"] + filtered_drugs,
            index=0
        )
        
        if st.button("Generate Network", type="primary", use_container_width=True):
            with st.spinner("Building network..."):
                selected_drug = None
                if target_drug != "Show All Drugs":
                    selected_drug = target_drug
                
                # Crear grafo COMPLETO
                G = crear_grafo_completo(df, selected_drug)
                farmaco_principal = selected_drug
                
                # Guardar en estado
                st.session_state['graph'] = G
                st.session_state['main_drug'] = farmaco_principal
                st.session_state['df'] = df
    
    with col1:
        # Mostrar estadísticas básicas si hay grafo
        if 'graph' in st.session_state and st.session_state['graph']:
            G = st.session_state['graph']
            
            # Estadísticas rápidas
            st.metric("Drugs in Network", len(G.nodes()))
            st.metric("Interactions", len(G.edges()))
            
            if st.session_state['main_drug']:
                degree = G.degree(st.session_state['main_drug']) if st.session_state['main_drug'] in G else 0
                st.metric("Connections for selected drug", degree)
            
            # Visualización
            fig = visualizar_con_plotly(G, st.session_state['main_drug'])
            st.plotly_chart(fig, use_container_width=True)
            
            # Información básica de categorías
            if len(G.nodes()) > 0:
                with st.expander("View ATC Categories in Network"):
                    categories = Counter()
                    for node in G.nodes():
                        categories[G.nodes[node]['atc_category']] += 1
                    
                    for category, count in sorted(categories.items()):
                        st.write(f"**{category}**: {count} drugs")
        else:
            st.info("Select a drug and click 'Generate Network' to visualize")

# Footer simple
st.markdown("---")
st.markdown(
    "<div style='text-align: center; color: gray;'>"
    "Drug Interaction Network Visualizer | Streamlit"
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
