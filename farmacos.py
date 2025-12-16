import streamlit as st
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from collections import Counter
import numpy as np
import pandas as pd
import io
import base64

# Configurar la página
st.set_page_config(
    page_title="Drug Interaction Network",
    page_icon="💊",
    layout="wide",
    initial_sidebar_state="expanded"
)

# CSS personalizado
st.markdown("""
<style>
    .main {
        padding: 2rem;
    }
    .stButton>button {
        width: 100%;
    }
    .drug-title {
        color: #1E88E5;
        font-size: 2rem;
        margin-bottom: 1rem;
    }
    .info-box {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 0.5rem;
        margin: 1rem 0;
    }
    .category-tag {
        display: inline-block;
        padding: 3px 10px;
        margin: 2px;
        border-radius: 15px;
        font-size: 12px;
        font-weight: 500;
    }
    .stCheckbox > div {
        margin-bottom: 0.5rem;
    }
</style>
""", unsafe_allow_html=True)

# Título principal
st.title("💊 Drug Interaction Network Explorer")
st.markdown("---")

# ATC Categories with descriptions (igual que tu código)
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

# Colors for each ATC category (igual que tu código)
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
    """Exactamente igual que tu código original"""
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

def crear_grafo_con_informacion(df, farmaco_objetivo=None):
    """Versión adaptada para Streamlit"""
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
    
    if farmaco_objetivo:
        farmaco_encontrado = None
        for node in G.nodes():
            if farmaco_objetivo.lower() in node.lower():
                farmaco_encontrado = node
                break
        
        if farmaco_encontrado:
            predecessors = list(G.predecessors(farmaco_encontrado))
            successors = list(G.successors(farmaco_encontrado))
            subgraph_nodes = [farmaco_encontrado] + predecessors + successors
            G = G.subgraph(subgraph_nodes).copy()
        else:
            st.warning(f"Drug '{farmaco_objetivo}' not found in the data")
            return None
    
    return G

def generar_visualizacion_streamlit(df, farmaco_objetivo=None, selected_categories=None):
    """Versión Streamlit de tu visualización"""
    
    # Crear grafo
    G = crear_grafo_con_informacion(df, farmaco_objetivo)
    
    if G is None:
        return None, None
    
    # Filtrar por categorías si hay selección
    if selected_categories:
        nodes_to_keep = []
        for node in G.nodes():
            if node == farmaco_objetivo:
                nodes_to_keep.append(node)
                continue
            
            node_category = G.nodes[node]['atc_category']
            if node_category in selected_categories:
                nodes_to_keep.append(node)
        
        G = G.subgraph(nodes_to_keep).copy()
    
    if len(G.nodes()) == 0:
        return None, None
    
    # Crear figura
    fig, ax = plt.subplots(figsize=(14, 10))
    
    # Layout
    farmaco_principal = farmaco_objetivo
    if farmaco_principal and farmaco_principal in G.nodes():
        pos = {}
        pos[farmaco_principal] = np.array([0, 0])
        
        neighbors = list(G.predecessors(farmaco_principal)) + \
                   list(G.successors(farmaco_principal))
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
        pos = nx.spring_layout(G, k=2/np.sqrt(len(G.nodes())), 
                              iterations=50, seed=42)
    
    # Preparar datos para visualización
    node_colors = []
    node_sizes = []
    
    for node in G.nodes():
        node_colors.append(G.nodes[node]['color'])
        if farmaco_principal and node == farmaco_principal:
            node_sizes.append(1200)
        else:
            node_sizes.append(600)
    
    # Dibujar nodos
    nx.draw_networkx_nodes(G, pos, 
                          node_color=node_colors,
                          node_size=node_sizes,
                          alpha=0.9,
                          edgecolors='black',
                          linewidths=1,
                          ax=ax)
    
    # Dibujar aristas
    edge_colors = []
    for u, v, data in G.edges(data=True):
        edge_colors.append('gray')
    
    nx.draw_networkx_edges(G, pos,
                          width=1.5,
                          alpha=0.6,
                          edge_color=edge_colors,
                          arrows=True,
                          arrowsize=10,
                          arrowstyle='-|>',
                          node_size=node_sizes,
                          ax=ax)
    
    # Etiquetas
    labels = {}
    for node in G.nodes():
        if farmaco_principal and node == farmaco_principal:
            labels[node] = node
        else:
            if len(node) > 20:
                labels[node] = node[:17] + "..."
            else:
                labels[node] = node
    
    nx.draw_networkx_labels(
        G,
        pos,
        labels,
        font_size=8,
        font_weight='normal',
        bbox=dict(
            facecolor='white',
            edgecolor='none',
            alpha=0.6
        ),
        ax=ax
    )
    
    # Título y estadísticas
    stats_text = (f"Drugs: {len(G.nodes())} | Interactions: {len(G.edges())}")
    
    if farmaco_principal:
        title = f"Drug: {farmaco_principal}\n{stats_text}"
    else:
        title = f"Complete Network\n{stats_text}"
    
    ax.set_title(title, fontsize=14, pad=20)
    ax.axis('off')
    
    plt.tight_layout()
    
    # Información para tooltips (simulada con Streamlit)
    info_data = {
        'nodes': {},
        'edges': []
    }
    
    for node in G.nodes():
        info_data['nodes'][node] = {
            'atc_code': G.nodes[node]['atc_code'],
            'atc_category': G.nodes[node]['atc_category'],
            'color': G.nodes[node]['color']
        }
    
    for u, v, data in G.edges(data=True):
        info_data['edges'].append({
            'from': u,
            'to': v,
            'y_value': data.get('interaction_type', 'Unknown')
        })
    
    return fig, info_data

# Función principal de la app
def main():
    # Sidebar - Controles (lado izquierdo)
    with st.sidebar:
        st.header("⚙️ Configuration")
        
        # Cargar datos
        @st.cache_data
        def load_data():
            try:
                df = pd.read_csv("DDIBUENO.csv")
                return df
            except Exception as e:
                st.error(f"Error loading CSV: {e}")
                return None
        
        df = load_data()
        
        if df is None:
            st.stop()
        
        # Obtener lista de fármacos
        all_drugs = sorted(set(list(df['Common_name_x'].unique()) + 
                              list(df['Common_name_y'].unique())))
        
        # Selectbox para fármaco objetivo
        target_drug = st.selectbox(
            "Select Target Drug",
            options=["(Show all)"] + all_drugs,
            index=0,
            help="Choose a drug to analyze its interactions"
        )
        
        st.markdown("---")
        st.header("🎨 Filter by ATC Category")
        
        # Obtener todas las categorías presentes
        all_categories = set()
        for _, row in df.iterrows():
            for suffix in ['x', 'y']:
                atc_code = row[f'atc_code_{suffix}']
                num_atc = row[f'num_atc_{suffix}']
                color, category = get_atc_color_and_category(atc_code, num_atc)
                all_categories.add(category)
        
        category_list = sorted(list(all_categories))
        
        # Checkboxes para categorías
        selected_categories = []
        for category in category_list:
            if st.checkbox(f"{category}", value=True, key=f"cat_{category}"):
                selected_categories.append(category)
        
        st.markdown("---")
        
        # Botones de acción
        col1, col2 = st.columns(2)
        with col1:
            if st.button("Select All", use_container_width=True):
                for category in category_list:
                    st.session_state[f"cat_{category}"] = True
                st.rerun()
        
        with col2:
            if st.button("Deselect All", use_container_width=True):
                for category in category_list:
                    st.session_state[f"cat_{category}"] = False
                st.rerun()
        
        st.markdown("---")
        
        # Leyenda de colores
        st.header("🌈 Color Legend")
        for category, color in ATC_COLORS.items():
            if category in ['Sin ATC', 'Multi ATC'] or category in ATC_CATEGORIES:
                st.markdown(f"""
                <div style='display: flex; align-items: center; margin-bottom: 5px;'>
                    <div style='width: 20px; height: 20px; background-color: {color}; 
                         margin-right: 10px; border: 1px solid #ddd; border-radius: 4px;'></div>
                    <span>{category}</span>
                </div>
                """, unsafe_allow_html=True)
    
    # Contenido principal
    if df is not None:
        # Preparar parámetros
        farmaco_param = None if target_drug == "(Show all)" else target_drug
        
        # Generar visualización
        with st.spinner("Generating network visualization..."):
            fig, info_data = generar_visualizacion_streamlit(
                df, 
                farmaco_param, 
                selected_categories
            )
        
        if fig is not None:
            # Mostrar gráfico
            st.pyplot(fig)
            
            # Estadísticas
            col1, col2, col3 = st.columns(3)
            
            if info_data and 'nodes' in info_data:
                with col1:
                    st.metric("Total Drugs", len(info_data['nodes']))
                
                with col2:
                    st.metric("Total Interactions", 
                             len(info_data['edges']) if 'edges' in info_data else 0)
                
                with col3:
                    if farmaco_param:
                        st.metric("Target Drug", farmaco_param)
                    else:
                        st.metric("Network Type", "Complete")
            
            # Información detallada
            st.markdown("---")
            st.subheader("🔍 Detailed Information")
            
            # Pestañas para nodos y aristas
            tab1, tab2 = st.tabs(["📋 Drugs Information", "🔗 Interactions"])
            
            with tab1:
                if info_data and 'nodes' in info_data:
                    st.dataframe(
                        pd.DataFrame.from_dict(info_data['nodes'], orient='index')
                        .reset_index()
                        .rename(columns={'index': 'Drug Name'}),
                        use_container_width=True,
                        height=300
                    )
            
            with tab2:
                if info_data and 'edges' in info_data:
                    edges_df = pd.DataFrame(info_data['edges'])
                    st.dataframe(
                        edges_df,
                        use_container_width=True,
                        height=300
                    )
                    
                    # Mostrar valor Y de las interacciones
                    if 'y_value' in edges_df.columns:
                        st.subheader("Interaction Types (Y values)")
                        y_counts = edges_df['y_value'].value_counts()
                        st.bar_chart(y_counts)
        else:
            st.warning("No data to display with current filters. Try selecting different categories.")
        
        # Footer
        st.markdown("---")
        st.caption("""
        **Interactive Features:**
        - Hover over the graph to see details
        - Click checkboxes to filter by ATC category
        - Select a target drug to focus on specific interactions
        """)
        
        # Información del dataset
        with st.expander("📊 Dataset Information"):
            st.write(f"**Total records:** {len(df)}")
            st.write(f"**Unique drugs:** {len(all_drugs)}")
            
            # Contar categorías
            category_counts = Counter()
            for _, row in df.iterrows():
                for suffix in ['x', 'y']:
                    atc_code = row[f'atc_code_{suffix}']
                    num_atc = row[f'num_atc_{suffix}']
                    color, category = get_atc_color_and_category(atc_code, num_atc)
                    category_counts[category] += 1
            
            st.write("**Drugs by ATC Category:**")
            for category, count in sorted(category_counts.items()):
                st.write(f"  - {category}: {count}")

if __name__ == "__main__":
    main()

