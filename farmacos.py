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

def crear_grafo_plotly(G, farmaco_principal=None, active_categories=None, show_all=True):
    """Crear visualización del grafo con Plotly"""
    
    if active_categories is None:
        active_categories = {}
    
    # Si show_all es True, mostrar todos los nodos (excepto si hay un fármaco principal)
    if show_all:
        nodes_to_keep = list(G.nodes())
    else:
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
        return None, None
    
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
        
        # Añadir nodos no conectados directamente
        remaining_nodes = set(G_filtered.nodes()) - set([farmaco_principal] + neighbors)
        if remaining_nodes:
            n_remaining = len(remaining_nodes)
            outer_radius = radius * 1.8 if radius > 0 else 2.0
            outer_angle_step = 2 * np.pi / n_remaining
            
            for i, node in enumerate(remaining_nodes):
                angle = i * outer_angle_step
                x = outer_radius * np.cos(angle)
                y = outer_radius * np.sin(angle)
                pos[node] = np.array([x, y])
    else:
        # Layout spring para grafo completo
        pos = nx.spring_layout(G_filtered, k=1.5/np.sqrt(len(G_filtered.nodes())), 
                              iterations=100, seed=42)
    
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
            node_sizes.append(30)
        else:
            # Tamaño basado en grado (número de conexiones)
            degree = G_filtered.degree(node)
            node_sizes.append(10 + min(degree, 10) * 2)
        
        node_texts.append(G_filtered.nodes[node]['tooltip'])
        
        # Nombre abreviado para etiqueta
        if len(node) > 20:
            node_names.append(node[:17] + "...")
        else:
            node_names.append(node)
    
    # Preparar datos para las aristas
    edge_x = []
    edge_y = []
    edge_texts = []
    edge_colors = []
    
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
        
        # Color basado en tipo de interacción
        interaction_type = data.get('interaction_type', 0)
        if interaction_type == 1:
            edge_colors.append('#FF0000')  # Rojo para tipo 1
        elif interaction_type == 2:
            edge_colors.append('#0000FF')  # Azul para tipo 2
        else:
            edge_colors.append('#808080')  # Gris para desconocido
    
    # Crear trazas de Plotly para aristas
    edge_traces = []
    for i in range(0, len(edge_x), 3):
        edge_trace = go.Scatter(
            x=edge_x[i:i+3], y=edge_y[i:i+3],
            line=dict(width=1.5, color=edge_colors[i//3]),
            hoverinfo='text',
            text=[edge_texts[i//3]],
            mode='lines',
            hoverlabel=dict(bgcolor='white', font_size=12),
            showlegend=False
        )
        edge_traces.append(edge_trace)
    
    # Crear traza para nodos
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
            line=dict(width=1.5, color='black')
        ),
        name='Drugs',
        hoverlabel=dict(bgcolor='white', font_size=12)
    )
    
    # Crear figura
    fig = go.Figure(data=edge_traces + [node_trace],
                   layout=go.Layout(
                       title=dict(
                           text=f'Drug Interaction Network<br>Drugs: {len(G_filtered.nodes())} | Interactions: {len(G_filtered.edges())}',
                           font=dict(size=16),
                           x=0.5,
                           xanchor='center'
                       ),
                       showlegend=False,
                       hovermode='closest',
                       margin=dict(b=20, l=5, r=5, t=60),
                       xaxis=dict(
                           showgrid=False, 
                           zeroline=False, 
                           showticklabels=False,
                           range=[min(node_x) - 0.5, max(node_x) + 0.5]
                       ),
                       yaxis=dict(
                           showgrid=False, 
                           zeroline=False, 
                           showticklabels=False,
                           range=[min(node_y) - 0.5, max(node_y) + 0.5]
                       ),
                       plot_bgcolor='white',
                       height=800,
                       width=1200,
                       dragmode='pan'
                   ))
    
    # Añadir cuadro de leyenda para colores ATC
    legend_items = []
    y_pos = 1.0
    for category, color in ATC_COLORS.items():
        if category in ATC_CATEGORIES:
            legend_items.append(
                dict(x=1.02, y=y_pos, xref="paper", yref="paper",
                     text=f"<span style='color:{color}'><b>■</b></span> {ATC_CATEGORIES[category]}",
                     showarrow=False, font=dict(size=10))
            )
            y_pos -= 0.04
    
    # Añadir leyenda para tipos de interacción
    legend_items.append(
        dict(x=1.02, y=y_pos, xref="paper", yref="paper",
             text="<b>Interaction Types:</b>",
             showarrow=False, font=dict(size=11, weight='bold'))
    )
    y_pos -= 0.03
    legend_items.append(
        dict(x=1.02, y=y_pos, xref="paper", yref="paper",
             text="<span style='color:#FF0000'><b>—</b></span> Type 1",
             showarrow=False, font=dict(size=10))
    )
    y_pos -= 0.02
    legend_items.append(
        dict(x=1.02, y=y_pos, xref="paper", yref="paper",
             text="<span style='color:#0000FF'><b>—</b></span> Type 2",
             showarrow=False, font=dict(size=10))
    )
    
    fig.update_layout(annotations=legend_items)
    
    return fig, G_filtered

def main():
    st.title("📊 Drug Interaction Network Visualization")
    
    # Cargar grafo completo una vez
    if 'G_completo' not in st.session_state:
        with st.spinner("Loading drug interaction data..."):
            st.session_state.G_completo = crear_grafo_completo(df)
            st.success(f"Loaded {len(st.session_state.G_completo.nodes())} drugs and {len(st.session_state.G_completo.edges())} interactions")
    
    G_completo = st.session_state.G_completo
    
    # Sidebar para controles
    with st.sidebar:
        st.header("⚙️ Controls")
        
        # Modo de visualización
        visualization_mode = st.radio(
            "Visualization Mode:",
            ["Complete Network", "Focused on Specific Drug"]
        )
        
        # Buscar fármaco
        all_drugs = sorted(set(df['Common_name_x'].tolist() + df['Common_name_y'].tolist()))
        farmaco_principal = None
        
        if visualization_mode == "Focused on Specific Drug":
            target_drug = st.selectbox(
                "Select drug to focus on:",
                all_drugs,
                index=0
            )
            
            # Crear subgrafo centrado
            G_filtrado, farmaco_principal = crear_subgrafo_centrado(G_completo, target_drug)
            
            if G_filtrado is None:
                st.error(f"Drug '{target_drug}' not found!")
                st.stop()
            
            st.success(f"Found {len(G_filtrado.nodes())} drugs connected to {farmaco_principal}")
        else:
            G_filtrado = G_completo
            farmaco_principal = None
        
        st.subheader("🔍 Filter by ATC Category")
        
        # Obtener todas las categorías presentes
        if G_filtrado:
            all_categories = set()
            for node in G_filtrado.nodes():
                all_categories.add(G_filtrado.nodes[node]['atc_category'])
            
            category_list = sorted(list(all_categories))
            
            # Inicializar estados de checkboxes
            if 'active_categories' not in st.session_state:
                st.session_state.active_categories = {cat: True for cat in category_list}
            
            # Checkboxes para categorías
            st.session_state.active_categories = {}
            for category in category_list:
                st.session_state.active_categories[category] = st.checkbox(
                    f"{category} ({sum(1 for node in G_filtrado.nodes() if G_filtrado.nodes[node]['atc_category'] == category)})",
                    value=True,
                    key=f"cat_{category}"
                )
            
            # Botones de selección rápida
            col1, col2 = st.columns(2)
            with col1:
                if st.button("✅ Select All", use_container_width=True):
                    for category in category_list:
                        st.session_state[f"cat_{category}"] = True
            
            with col2:
                if st.button("❌ Deselect All", use_container_width=True):
                    for category in category_list:
                        st.session_state[f"cat_{category}"] = False
            
            # Opción para mostrar todos sin filtrar
            show_all_categories = st.checkbox("Show all drugs (ignore filters)", value=True)
    
    # Área principal para el gráfico
    st.subheader("📈 Network Visualization")
    
    if G_filtrado:
        # Crear gráfico
        fig, G_visualizado = crear_grafo_plotly(
            G_filtrado, 
            farmaco_principal, 
            st.session_state.active_categories,
            show_all=show_all_categories
        )
        
        if fig:
            st.plotly_chart(fig, use_container_width=True, config={'displayModeBar': True})
            
            # Mostrar estadísticas
            st.subheader("📊 Statistics")
            
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.metric("Total Drugs in Dataset", len(G_completo.nodes()))
            with col2:
                st.metric("Filtered Drugs", len(G_visualizado.nodes()))
            with col3:
                st.metric("Total Interactions", len(G_visualizado.edges()))
            with col4:
                if farmaco_principal:
                    degree = G_completo.degree(farmaco_principal)
                    st.metric(f"Connections of {farmaco_principal[:15]}...", degree)
            
            # Estadísticas detalladas
            with st.expander("📋 Detailed Statistics", expanded=True):
                # Contar por categoría
                category_counts = Counter()
                for node in G_visualizado.nodes():
                    category_counts[G_visualizado.nodes[node]['atc_category']] += 1
                
                # Preparar datos para tabla
                stats_data = []
                for category, count in sorted(category_counts.items()):
                    percentage = (count / len(G_visualizado.nodes())) * 100
                    stats_data.append({
                        'Category': category,
                        'Drug Count': count,
                        'Percentage': f"{percentage:.1f}%",
                        'Color': ATC_COLORS.get(category[:1], ATC_COLORS.get(category, '#CCCCCC'))
                    })
                
                # Mostrar tabla
                stats_df = pd.DataFrame(stats_data)
                st.dataframe(stats_df, use_container_width=True, hide_index=True)
                
                # Gráfico de distribución
                if len(stats_data) > 0:
                    fig_bar = px.bar(
                        stats_df,
                        x='Category',
                        y='Drug Count',
                        color='Category',
                        color_discrete_map=ATC_COLORS,
                        title="Drug Distribution by ATC Category",
                        text='Drug Count'
                    )
                    fig_bar.update_traces(textposition='outside')
                    fig_bar.update_layout(showlegend=False)
                    st.plotly_chart(fig_bar, use_container_width=True)
            
            # Mostrar lista de fármacos
            with st.expander("🧪 List of Drugs in View"):
                drugs_in_view = sorted(G_visualizado.nodes())
                cols = 4
                rows = (len(drugs_in_view) + cols - 1) // cols
                
                for i in range(rows):
                    cols_list = st.columns(cols)
                    for j in range(cols):
                        idx = i * cols + j
                        if idx < len(drugs_in_view):
                            drug = drugs_in_view[idx]
                            color = G_visualizado.nodes[drug]['color']
                            category = G_visualizado.nodes[drug]['atc_category']
                            cols_list[j].markdown(
                                f"<div style='background-color:{color}; padding:5px; border-radius:5px; margin:2px;'>"
                                f"<b>{drug}</b><br><small>{category}</small>"
                                f"</div>",
                                unsafe_allow_html=True
                            )
        else:
            st.warning("No drugs to display. Please select at least one category.")
    else:
        st.info("Select a visualization mode and drug to begin.")

    # Footer
    st.markdown("---")
    st.markdown("**Drug Interaction Network** • Powered by NetworkX and Plotly")

if __name__ == "__main__":
    main()
