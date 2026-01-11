# app_streamlit_completo.py
import streamlit as st
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
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
    .warning-box {
        background-color: #FFF3CD;
        border: 1px solid #FFEAA7;
        border-radius: 0.5rem;
        padding: 1rem;
        margin: 1rem 0;
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
        # Para grafo completo, usar muestra limitada
        filtered_df = _df.head(1000)  # Límite para evitar sobrecarga
    
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
        fig = go.Figure()
        fig.add_annotation(
            text="No drugs to display",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=20)
        )
        return fig
    
    # Posiciones de los nodos
    if farmaco_principal and farmaco_principal in G.nodes():
        pos = nx.spring_layout(G, k=1, iterations=50, seed=42)
    else:
        pos = nx.spring_layout(G, k=1.5/np.sqrt(len(G.nodes())), iterations=50, seed=42)
    
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
    node_names = []
    
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
        node_names.append(node)
        
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
        text=[node[:15] + ('...' if len(node) > 15 else '') for node in node_names],
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
    network_title = "Drug Interaction Network"
    if farmaco_principal:
        network_title = f"Interactions for: {farmaco_principal}"
    
    fig.update_layout(
        title=f"{network_title}<br>{len(G.nodes())} drugs, {len(G.edges())} interactions",
        showlegend=False,
        hovermode='closest',
        margin=dict(b=20, l=5, r=5, t=60),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        height=600,
        plot_bgcolor='white'
    )
    
    return fig

# Función principal de carga de datos - CARGA COMPLETA
@st.cache_data(ttl=3600)
def load_complete_data(file_path=None, use_cols=None):
    """Carga completa de datos sin muestreo"""
    try:
        # Mostrar progreso
        progress_text = st.empty()
        progress_bar = st.progress(0)
        
        progress_text.text("Loading data...")
        
        if use_cols:
            df = pd.read_csv(file_path, usecols=use_cols)
        else:
            df = pd.read_csv(file_path)
        
        progress_bar.progress(100)
        progress_text.text("Data loaded successfully!")
        time.sleep(0.5)
        progress_text.empty()
        progress_bar.empty()
        
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

# REMOVIDO: Slider de sample_size
# En su lugar, mantener solo max_interactions
max_interactions = st.sidebar.slider(
    "Max interactions per drug",
    min_value=10,
    max_value=500,
    value=100,
    step=10,
    help="Limit interactions for better performance"
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
        
        # Cargar datos COMPLETOS
        with st.spinner("Loading complete dataset..."):
            # Columnas necesarias
            use_cols = ['Common_name_x', 'Common_name_y', 'Y', 'atc_code_x', 'atc_code_y']
            
            # Verificar si las columnas num_atc existen
            try:
                df_sample = pd.read_csv(file_path, nrows=1)
                if 'num_atc_x' in df_sample.columns:
                    use_cols.extend(['num_atc_x', 'num_atc_y'])
            except:
                pass
            
            # CARGAR TODO EL DATASET
            df = load_complete_data(file_path, use_cols)
        
        if df is not None:
            # Mostrar advertencia si el dataset es muy grande
            if len(df) > 100000:
                st.warning(f"⚠️ Large dataset detected: {len(df):,} rows. Visualization may be slower.")
            
            # Obtener lista de fármacos
            drug_list = get_drug_list(df)
            
            st.info(f"📊 Dataset loaded: {len(df):,} interactions, {len(drug_list)} unique drugs")
            
            # Selector de fármaco con búsqueda
            search_term = st.text_input("Search drug by name:", 
                                       placeholder="Type drug name...")
            
            # Filtrar fármacos basado en búsqueda
            filtered_drugs = []
            if search_term:
                filtered_drugs = [drug for drug in drug_list 
                                 if search_term.lower() in drug.lower()]
            else:
                filtered_drugs = drug_list[:200]  # Mostrar primeros 200 por defecto
            
            target_drug = st.selectbox(
                "Select target drug:",
                options=["All Drugs"] + filtered_drugs,
                index=0,
                help="Select 'All Drugs' for complete network overview"
            )
            
            # Botón para generar
            if st.button("Generate Network", type="primary", use_container_width=True):
                with st.spinner("Building network..."):
                    # Determinar fármaco objetivo
                    selected_drug = None
                    if target_drug != "All Drugs" and target_drug:
                        selected_drug = target_drug
                    
                    # Crear grafo
                    if selected_drug:
                        G = crear_grafo_optimizado(df, selected_drug, max_interactions)
                        farmaco_principal = selected_drug
                        st.session_state['graph'] = G
                        st.session_state['main_drug'] = farmaco_principal
                        st.session_state['df'] = df
                    else:
                        # Para grafo completo, usar muestra limitada
                        st.info("Showing overview of first 1000 interactions")
                        G = crear_grafo_optimizado(df.head(1000), None, max_interactions)
                        farmaco_principal = None
                        st.session_state['graph'] = G
                        st.session_state['main_drug'] = None
                        st.session_state['df'] = df
                    
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
                            else:
                                st.metric("Network Type", "Complete")
                    
                    # Visualización
                    if visualization_method == "Plotly (Interactive)":
                        fig = visualizar_con_plotly(G, farmaco_principal)
                        st.plotly_chart(fig, use_container_width=True)
                    else:
                        # Versión Matplotlib simplificada
                        fig, ax = plt.subplots(figsize=(12, 8))
                        pos = nx.spring_layout(G, seed=42)
                        
                        node_colors = [G.nodes[node]['color'] for node in G.nodes()]
                        node_sizes = [800 if node == farmaco_principal else 400 
                                    for node in G.nodes()]
                        
                        nx.draw(G, pos, 
                               node_color=node_colors,
                               node_size=node_sizes,
                               edge_color='gray',
                               with_labels=True,
                               font_size=8,
                               ax=ax)
                        ax.set_title(f"Drug Network: {len(G.nodes())} drugs")
                        ax.axis('off')
                        st.pyplot(fig)
                    
                    # Información detallada
                    with st.expander("📊 Network Details", expanded=False):
                        col_info1, col_info2 = st.columns(2)
                        
                        with col_info1:
                            st.subheader("Drugs in Network")
                            nodes_to_show = list(G.nodes())[:20]  # Mostrar primeros 20
                            for node in nodes_to_show:
                                node_data = G.nodes[node]
                                st.markdown(f"""
                                <div class="drug-card">
                                <b>{node}</b><br>
                                <small>ATC: {node_data.get('atc_code', 'N/A')}<br>
                                Category: {node_data.get('atc_category', 'Unknown')}</small>
                                </div>
                                """, unsafe_allow_html=True)
                            
                            if len(G.nodes()) > 20:
                                st.caption(f"... and {len(G.nodes()) - 20} more drugs")
                        
                        with col_info2:
                            st.subheader("ATC Category Distribution")
                            categories = Counter()
                            for node in G.nodes():
                                categories[G.nodes[node]['atc_category']] += 1
                            
                            if categories:
                                # Crear gráfico de barras en lugar de torta para muchas categorías
                                cat_df = pd.DataFrame({
                                    'Category': list(categories.keys()),
                                    'Count': list(categories.values())
                                }).sort_values('Count', ascending=False)
                                
                                fig_bar = px.bar(cat_df, 
                                               x='Category', 
                                               y='Count',
                                               title="Drugs by ATC Category",
                                               color='Count',
                                               color_continuous_scale='Blues')
                                fig_bar.update_layout(xaxis_tickangle=-45)
                                st.plotly_chart(fig_bar, use_container_width=True)
                                
                                # También mostrar tabla
                                with st.expander("View as table"):
                                    st.dataframe(cat_df, use_container_width=True)

with tab2:
    st.subheader("📁 Data Explorer")
    
    if 'df' in st.session_state:
        df = st.session_state['df']
        
        st.info(f"📈 Dataset statistics: {len(df):,} total interactions")
        
        # Mostrar vista previa de datos
        with st.expander("View Raw Data (First 100 rows)", expanded=False):
            st.dataframe(df.head(100), use_container_width=True)
        
        # Estadísticas básicas
        col_stats1, col_stats2, col_stats3 = st.columns(3)
        
        with col_stats1:
            st.metric("Total Rows", f"{len(df):,}")
        
        with col_stats2:
            unique_drugs = len(set(df['Common_name_x'].tolist() + df['Common_name_y'].tolist()))
            st.metric("Unique Drugs", f"{unique_drugs:,}")
        
        with col_stats3:
            st.metric("Total Interactions", f"{len(df):,}")
        
        # Filtros interactivos
        st.subheader("Filter Data")
        
        col_filter1, col_filter2, col_filter3 = st.columns(3)
        
        with col_filter1:
            interaction_type = st.selectbox(
                "Interaction Type:",
                ["All", "Type 1", "Type 2"],
                key="filter_type"
            )
        
        with col_filter2:
            min_interactions = st.number_input(
                "Min interactions per drug:",
                min_value=1,
                max_value=100,
                value=1,
                step=1
            )
        
        with col_filter3:
            atc_options = ["All"] + list(ATC_CATEGORIES.values())
            selected_atc = st.selectbox(
                "ATC Category:",
                atc_options,
                key="filter_atc"
            )
        
        if st.button("Apply Filters", type="primary"):
            with st.spinner("Applying filters..."):
                filtered_df = df.copy()
                
                # Aplicar filtro de tipo de interacción
                if interaction_type != "All":
                    type_value = 1 if interaction_type == "Type 1" else 2
                    filtered_df = filtered_df[filtered_df['Y'] == type_value]
                
                # Aplicar filtro ATC
                if selected_atc != "All":
                    atc_code = None
                    for code, category in ATC_CATEGORIES.items():
                        if category == selected_atc:
                            atc_code = code
                            break
                    
                    if atc_code:
                        mask = (filtered_df['atc_code_x'].astype(str).str.startswith(atc_code)) | \
                               (filtered_df['atc_code_y'].astype(str).str.startswith(atc_code))
                        filtered_df = filtered_df[mask]
                
                st.success(f"Filtered to {len(filtered_df):,} interactions")
                
                # Mostrar datos filtrados
                st.dataframe(filtered_df.head(50), use_container_width=True)
                
                # Exportar opción
                csv = filtered_df.to_csv(index=False).encode('utf-8')
                st.download_button(
                    label="Download filtered data as CSV",
                    data=csv,
                    file_name="filtered_interactions.csv",
                    mime="text/csv"
                )

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
    
    # Mostrar tabla de categorías ATC - CORREGIDO
    atc_table_data = []
    for code, category in ATC_CATEGORIES.items():
        if code in ATC_COLORS:
            atc_table_data.append({
                "Code": code,
                "Category": category,
                "Color": ATC_COLORS[code]
            })
    
    atc_table = pd.DataFrame(atc_table_data)
    
    # Añadir muestra de color
    def color_cell(color_hex):
        return f'<div style="background-color:{color_hex}; width:20px; height:20px; border-radius:3px;"></div>'
    
    atc_table['Color Sample'] = atc_table['Color'].apply(color_cell)
    
    # Mostrar tabla
    st.markdown(atc_table[['Code', 'Category', 'Color Sample']].to_html(
        escape=False, 
        index=False
    ), unsafe_allow_html=True)
    
    st.markdown("""
    
    ### How to Use:
    1. **Upload** your CSV file or use the default
    2. **Search** for a specific drug or select from the list
    3. **Generate** the network visualization
    4. **Explore** interactions and statistics
    5. **Use filters** in Data Explorer for detailed analysis
    
    ### Dataset Requirements:
    - CSV format with columns: `Common_name_x`, `Common_name_y`, `Y`, `atc_code_x`, `atc_code_y`
    - Optional: `num_atc_x`, `num_atc_y` for ATC multiplicity
    
    ### Performance Notes:
    - Complete dataset is loaded for searching
    - Network visualization shows limited interactions for performance
    - Use filters to focus on specific subsets
    """)
    
    # Información técnica
    with st.expander("Technical Details"):
        st.markdown("""
        **Technologies Used:**
        - Streamlit for web interface
        - NetworkX for graph operations
        - Plotly for interactive visualizations
        - Pandas for data manipulation
        
        **Graph Features:**
        - Directed graph showing drug interactions
        - Color-coded by ATC therapeutic category
        - Interactive hover information
        - Filtering by ATC categories
        """)

# Footer
st.markdown("---")
st.markdown(
    "<div style='text-align: center; color: gray;'>"
    "💊 Drug Interaction Network Analyzer | Built with Streamlit"
    "</div>",
    unsafe_allow_html=True
)

# Inicialización de estado
if 'graph' not in st.session_state:
    st.session_state['graph'] = None
if 'main_drug' not in st.session_state:
    st.session_state['main_drug'] = None
if 'df' not in st.session_state:
    st.session_state['df'] = None

# Nota importante sobre rendimiento
if uploaded_file:
    file_size = uploaded_file.size / (1024 * 1024)  # Convertir a MB
    if file_size > 50:  # Si el archivo es mayor a 50MB
        st.sidebar.warning(f"⚠️ Large file: {file_size:.1f} MB. Loading may take time.")

