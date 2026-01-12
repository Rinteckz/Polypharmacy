import networkx as nx
import plotly.graph_objects as go
import plotly.express as px
from collections import Counter
import numpy as np
import pandas as pd
import streamlit as st

# Configuración de Streamlit
st.set_page_config(layout="wide")

# Cargar datos (cambia la ruta según sea necesario)
@st.cache_data
def load_data():
    return pd.read_csv(r"C:\Users\edjom\OneDrive\Escritorio\TherapeuticDataCommons\DDIBUENO.csv")

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

def crear_grafo_con_informacion(df, farmaco_objetivo=None):
    """Versión con información adicional para tooltips"""
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
            interaction_desc = "Type 1" if interaction_type == 1 else "Type 2" if interaction_type == 2 else "Unknown"
            edge_tooltip = f"<b>Interaction Type:</b> {interaction_desc}<br><b>From:</b> {drug1} → <b>To:</b> {drug2}"
            G.add_edge(drug1, drug2, 
                      interaction_type=interaction_type,
                      tooltip=edge_tooltip)
    
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

def crear_grafo_plotly(G, farmaco_principal=None, active_categories=None):
    """Crear visualización del grafo con Plotly"""
    
    if active_categories is None:
        active_categories = {}
    
    # Filtrar nodos por categorías activas
    nodes_to_keep = []
    for node in G.nodes():
        if node == farmaco_principal:
            nodes_to_keep.append(node)
        else:
            category = G.nodes[node]['atc_category']
            if active_categories.get(category, False):
                nodes_to_keep.append(node)
    
    if not nodes_to_keep:
        return None
    
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
    node_x = []
    node_y = []
    node_colors = []
    node_sizes = []
    node_texts = []
    node_names = []
    
    for node in G_filtered.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        node_colors.append(G_filtered.nodes[node]['color'])
        
        # Tamaño diferente para fármaco principal
        if farmaco_principal and node == farmaco_principal:
            node_sizes.append(25)
        else:
            node_sizes.append(15)
        
        node_texts.append(G_filtered.nodes[node]['tooltip'])
        node_names.append(node)
    
    # Preparar datos para las aristas
    edge_x = []
    edge_y = []
    edge_texts = []
    
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
    
    # Crear trazas de Plotly
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=1.5, color='gray'),
        hoverinfo='text',
        text=edge_texts,
        mode='lines',
        name='Interactions'
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
        name='Drugs'
    )
    
    # Crear figura
    fig = go.Figure(data=[edge_trace, node_trace],
                   layout=go.Layout(
                       title=dict(
                           text=f'Drug Interaction Network<br>Drugs: {len(G_filtered.nodes())} | Interactions: {len(G_filtered.edges())}',
                           font=dict(size=16)
                       ),
                       showlegend=False,
                       hovermode='closest',
                       margin=dict(b=20, l=5, r=5, t=40),
                       xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                       yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                       plot_bgcolor='white',
                       height=700,
                       width=1000
                   ))
    
    return fig, G_filtered

def main():
    st.title("Drug Interaction Network Visualization")
    
    # Sidebar para controles
    with st.sidebar:
        st.header("Controls")
        
        # Buscar fármaco
        all_drugs = sorted(set(df['Common_name_x'].tolist() + df['Common_name_y'].tolist()))
        target_drug = st.selectbox(
            "Select a drug to focus on:",
            ["All drugs"] + all_drugs,
            index=0
        )
        
        if target_drug != "All drugs":
            farmaco_objetivo = target_drug
        else:
            farmaco_objetivo = None
        
        st.subheader("Filter by ATC Category")
        
        # Obtener todas las categorías únicas
        G_full = crear_grafo_con_informacion(df, farmaco_objetivo)
        if G_full is None and farmaco_objetivo:
            st.stop()
        
        if G_full:
            all_categories = set()
            farmaco_principal = None
            
            if farmaco_objetivo:
                for node in G_full.nodes():
                    if farmaco_objetivo.lower() in node.lower():
                        farmaco_principal = node
                        break
            
            for node in G_full.nodes():
                if node != farmaco_principal:
                    all_categories.add(G_full.nodes[node]['atc_category'])
            
            category_list = sorted(list(all_categories))
            
            # Checkboxes para categorías
            active_categories = {}
            for category in category_list:
                active_categories[category] = st.checkbox(
                    category,
                    value=True,
                    key=f"cat_{category}"
                )
            
            # Botones de selección rápida
            col1, col2 = st.columns(2)
            with col1:
                if st.button("Select All"):
                    for category in category_list:
                        st.session_state[f"cat_{category}"] = True
            
            with col2:
                if st.button("Deselect All"):
                    for category in category_list:
                        st.session_state[f"cat_{category}"] = False
        else:
            active_categories = {}
    
    # Área principal para el gráfico
    if G_full:
        # Actualizar estados de checkboxes desde session_state
        if 'active_categories' in locals():
            for category in category_list:
                if f"cat_{category}" in st.session_state:
                    active_categories[category] = st.session_state[f"cat_{category}"]
        
        # Crear gráfico
        fig, G_filtered = crear_grafo_plotly(G_full, farmaco_principal, active_categories)
        
        if fig:
            st.plotly_chart(fig, use_container_width=True)
            
            # Mostrar estadísticas
            st.subheader("Network Statistics")
            col1, col2, col3 = st.columns(3)
            
            with col1:
                st.metric("Total Drugs", len(G_full.nodes()))
            with col2:
                st.metric("Filtered Drugs", len(G_filtered.nodes()))
            with col3:
                st.metric("Interactions", len(G_filtered.edges()))
            
            # Contar por categoría
            category_counts = Counter()
            for node in G_filtered.nodes():
                if node != farmaco_principal:
                    category_counts[G_filtered.nodes[node]['atc_category']] += 1
            
            if category_counts:
                st.subheader("Drugs by ATC Category")
                category_df = pd.DataFrame({
                    'Category': list(category_counts.keys()),
                    'Count': list(category_counts.values())
                })
                category_df = category_df.sort_values('Count', ascending=False)
                
                # Gráfico de barras
                fig_bar = px.bar(
                    category_df,
                    x='Category',
                    y='Count',
                    color='Category',
                    color_discrete_map=ATC_COLORS,
                    title="Drug Distribution by ATC Category"
                )
                st.plotly_chart(fig_bar, use_container_width=True)
        else:
            st.warning("No drugs to display. Please select at least one category.")
    else:
        st.info("Select a drug or use 'All drugs' to visualize the network.")

if __name__ == "__main__":
    main()
