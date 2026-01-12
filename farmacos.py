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
            interaction_desc = "Type 1" if interaction_type == 1 else "Type 2" if interaction_type == 2 else "Unknown"
            edge_tooltip = f"<b>Interaction Type:</b> {interaction_desc}<br><b>From:</b> {drug1} → <b>To:</b> {drug2}"
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
    """Crear visualización del grafo con Plotly - FIEL AL COMPORTAMIENTO ORIGINAL"""
    
    # Si no hay categorías activas, no mostrar nada (como en el original)
    if active_categories is None or not any(active_categories.values()):
        # Crear figura vacía con mensaje
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
    
    # Filtrar nodos por categorías activas (EXACTAMENTE como en el original)
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
        # Crear figura vacía
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
    
    # Calcular posiciones de los nodos (EXACTAMENTE como en el original)
    if farmaco_principal and farmaco_principal in G_filtered.nodes():
        # Layout centrado en fármaco principal
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
        # Layout spring para cuando no hay fármaco principal
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
        
        # Tamaño diferente para fármaco principal (como en el original)
        if farmaco_principal and node == farmaco_principal:
            node_sizes.append(30)  # Equivalente a 1200 en matplotlib
        else:
            node_sizes.append(15)  # Equivalente a 600 en matplotlib
        
        node_texts.append(G_filtered.nodes[node]['tooltip'])
        
        # Nombre abreviado para etiqueta (como en el original)
        if farmaco_principal and node == farmaco_principal:
            node_names.append(node)  # Nombre completo para fármaco principal
        else:
            if len(node) > 20:
                node_names.append(node[:17] + "...")
            else:
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
    
    # Título basado en si hay fármaco principal o no
    if farmaco_principal:
        title_text = f"Drug: {farmaco_principal}<br>Drugs: {len(G_filtered.nodes())} | Interactions: {len(G_filtered.edges())}"
    else:
        title_text = f"Complete Network<br>Drugs: {len(G_filtered.nodes())} | Interactions: {len(G_filtered.edges())}"
    
    # Crear figura
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

def main():
    st.title("Drug Interaction Network Visualization")
    
    # Inicializar estado de sesión
    if 'G_completo' not in st.session_state:
        with st.spinner("Loading drug interaction data..."):
            st.session_state.G_completo = crear_grafo_completo(df)
    
    G_completo = st.session_state.G_completo
    
    # Sidebar para controles
    with st.sidebar:
        st.header("Controls")
        
        # Opción para seleccionar fármaco principal (NUEVO para Streamlit)
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
            # Crear subgrafo centrado en el fármaco
            G_actual, farmaco_real = crear_subgrafo_centrado(G_completo, farmaco_principal)
            if G_actual is None:
                st.error(f"Drug '{farmaco_principal}' not found!")
                st.stop()
            farmaco_principal = farmaco_real
            st.success(f"Showing network for: {farmaco_principal}")
        else:
            farmaco_principal = None
            G_actual = G_completo
        
        st.subheader("Filter by ATC Category")
        st.write("Select categories to display (none selected by default):")
        
        # Obtener todas las categorías presentes en el grafo actual
        all_categories = set()
        for node in G_actual.nodes():
            # Si hay fármaco principal, no incluirlo en las categorías (como en el original)
            if farmaco_principal and node == farmaco_principal:
                continue
            all_categories.add(G_actual.nodes[node]['atc_category'])
        
        category_list = sorted(list(all_categories))
        
        # Inicializar checkboxes desmarcados (como en el original)
        if 'active_categories' not in st.session_state:
            st.session_state.active_categories = {cat: False for cat in category_list}
        
        # Mostrar checkboxes con contadores
        for category in category_list:
            count = sum(1 for node in G_actual.nodes() 
                       if G_actual.nodes[node]['atc_category'] == category and 
                       node != farmaco_principal)
            
            col1, col2 = st.columns([1, 3])
            with col1:
                # Checkbox con estado guardado
                checked = st.checkbox(
                    "",
                    value=st.session_state.active_categories.get(category, False),
                    key=f"cat_{category}",
                    label_visibility="collapsed"
                )
                st.session_state.active_categories[category] = checked
            
            with col2:
                # Mostrar categoría con color y contador
                color = ATC_COLORS.get(
                    category if category in ATC_COLORS else category[0],
                    '#CCCCCC'
                )

                st.markdown(
                    f"<span style='color:{color}; font-weight:bold;'>■</span> "
                    f"{category} ({count})",
                    unsafe_allow_html=True
                )
        
        # Botones de selección rápida
        col1, col2 = st.columns(2)
        with col1:
            if st.button("Select All", use_container_width=True):
                for category in category_list:
                    st.session_state.active_categories[category] = True
                    st.session_state[f"cat_{category}"] = True
        
        with col2:
            if st.button("Deselect All", use_container_width=True):
                for category in category_list:
                    st.session_state.active_categories[category] = False
                    st.session_state[f"cat_{category}"] = False
    
    # Área principal
    st.subheader("Network Visualization")
    
    # Verificar si hay algo para mostrar
    if not any(st.session_state.active_categories.values()) and not farmaco_principal:
        # Mostrar mensaje como en el original
        st.info("👈 **Please select at least one ATC category from the sidebar to display drugs.**")
        
        # Mostrar estadísticas generales
        with st.expander("Dataset Statistics", expanded=True):
            col1, col2 = st.columns(2)
            with col1:
                st.metric("Total Drugs in Dataset", len(G_completo.nodes()))
            with col2:
                st.metric("Total Interactions", len(G_completo.edges()))
            
            # Mostrar distribución de categorías
            category_counts = Counter()
            for node in G_completo.nodes():
                category_counts[G_completo.nodes[node]['atc_category']] += 1
            
            st.write("**Drugs by ATC Category:**")
            for category, count in sorted(category_counts.items()):
                color = ATC_COLORS.get(category[:1], ATC_COLORS.get(category, '#CCCCCC'))
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
            
            # Mostrar detalles por categoría
            with st.expander("View details by category", expanded=False):
                category_counts = Counter()
                for node in G_filtrado.nodes():
                    if node != farmaco_principal:
                        category_counts[G_filtrado.nodes[node]['atc_category']] += 1
                
                if category_counts:
                    for category, count in sorted(category_counts.items()):
                        color = ATC_COLORS.get(category[:1], ATC_COLORS.get(category, '#CCCCCC'))
                        st.markdown(
                            f"<span style='color:{color}; font-weight:bold;'>■</span> "
                            f"{category}: {count} drugs",
                            unsafe_allow_html=True
                        )
                else:
                    st.write("No drugs from selected categories (only main drug shown)")
        
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
    
    # Información adicional
    st.markdown("---")
    st.markdown("**How to use:**")
    st.markdown("""
    1. **Select a main drug** (optional) from the dropdown to focus on its interactions
    2. **Select ATC categories** you want to visualize (all are deselected by default)
    3. **Hover over nodes** to see drug details
    4. **Hover over edges** to see interaction details
    5. **Use the mouse** to pan and zoom the graph
    """)

if __name__ == "__main__":
    main()


