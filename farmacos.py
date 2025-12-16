import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from collections import Counter
import numpy as np
import pandas as pd
import threading
def auto_remove_annotation(annotation, canvas, delay=3):
    def _remove():
        try:
            annotation.remove()
            canvas.draw_idle()
        except Exception:
            pass
    threading.Timer(delay, _remove).start()

# Intentar importar mplcursors, si no está instalado, usar alternativa
try:
    import mplcursors
    MPLCURSORS_AVAILABLE = True
except ImportError:
    print("mplcursors no está instalado. Usando tooltips alternativos...")
    MPLCURSORS_AVAILABLE = False

df = pd.read_csv(r"DDIBUENO.csv")

# ATC Categories with descriptions
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

# Colors for each ATC category
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
            # Información detallada para tooltips
            tooltip_info = f"Drug: {drug1}\nATC Code: {atc1}\nATC Category: {category1}"
            G.add_node(drug1, 
                      atc_code=atc1,
                      atc_category=category1,
                      color=color1,
                      num_atc=num_atc1,
                      tooltip=tooltip_info)
        
        if not G.has_node(drug2):
            tooltip_info = f"Drug: {drug2}\nATC Code: {atc2}\nATC Category: {category2}"
            G.add_node(drug2, 
                      atc_code=atc2,
                      atc_category=category2,
                      color=color2,
                      num_atc=num_atc2,
                      tooltip=tooltip_info)
        
        # Solo agregar arista si no existe, con información de tooltip
        if not G.has_edge(drug1, drug2):
            edge_tooltip = f"Interaction Type: {interaction_type}\nFrom: {drug1} → To: {drug2}"
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
            # Solo incluir conexiones directas para mejor rendimiento
            predecessors = list(G.predecessors(farmaco_encontrado))
            successors = list(G.successors(farmaco_encontrado))
            subgraph_nodes = [farmaco_encontrado] + predecessors + successors
            G = G.subgraph(subgraph_nodes).copy()
        else:
            print(f"Drug '{farmaco_objetivo}' not found in the data")
            return None
    
    return G

def visualizar_grafo_con_tooltips(df, farmaco_objetivo=None, figsize=(16, 10)):
    """Visualización con tooltips interactivos para nodos y aristas"""
    
    # Crear grafo con información para tooltips
    G = crear_grafo_con_informacion(df, farmaco_objetivo)
    
    if G is None:
        return
    
    # Encontrar fármaco principal
    farmaco_principal = None
    if farmaco_objetivo:
        for node in G.nodes():
            if farmaco_objetivo.lower() in node.lower():
                farmaco_principal = node
                break
    
    # Obtener todas las categorías ATC presentes
    all_categories = set()
    for node in G.nodes():
        if node != farmaco_principal:  # No filtrar el fármaco principal
            all_categories.add(G.nodes[node]['atc_category'])
    
    # Ordenar categorías
    category_list = sorted(list(all_categories))
    
    # Crear figura con espacio para checkboxes
    fig = plt.figure(figsize=figsize)
    
    # Grid para organización
    gs = plt.GridSpec(1, 5, figure=fig, width_ratios=[4, 0.1, 1, 0.1, 1])
    
    # Área principal del gráfico
    ax_graph = fig.add_subplot(gs[0, 0])
    
    # Área para checkboxes (lado derecho)
    ax_checkboxes = fig.add_subplot(gs[0, 2])
    ax_checkboxes.axis('off')
    
    # Área para controles adicionales
    ax_controls = fig.add_subplot(gs[0, 4])
    ax_controls.axis('off')
    
    # Estado inicial: todas las categorías desmarcadas (mejor rendimiento)
    active_categories = {cat: False for cat in category_list}
    
    # Variables globales para tooltips
    global node_positions, edge_info_list
    
    # Función para filtrar nodos por categorías seleccionadas
    def filtrar_por_categorias():
        # Siempre incluir el fármaco principal
        nodes_to_keep = [farmaco_principal] if farmaco_principal else []
        
        # Agregar nodos de categorías seleccionadas
        for node in G.nodes():
            if node == farmaco_principal:
                continue
            
            category = G.nodes[node]['atc_category']
            if active_categories.get(category, False):  # Solo si está seleccionada
                nodes_to_keep.append(node)
        
        # Crear subgrafo
        G_filtered = G.subgraph(nodes_to_keep).copy()
        
        # También filtrar aristas donde ambos nodos estén presentes
        edges_to_keep = []
        for u, v in G.edges():
            if u in nodes_to_keep and v in nodes_to_keep:
                edges_to_keep.append((u, v))
        
        # Crear nuevo grafo con solo las aristas válidas
        G_final = nx.DiGraph()
        for node in nodes_to_keep:
            G_final.add_node(node, **G.nodes[node])
        
        for u, v in edges_to_keep:
            G_final.add_edge(u, v, **G[u][v])
        
        return G_final
    
    # Función para dibujar el grafo
    def dibujar_grafo():
        ax_graph.clear()
        
        # Filtrar grafo
        G_filtered = filtrar_por_categorias()
        
        if len(G_filtered.nodes()) == 0:
            ax_graph.text(0.5, 0.5, "No drugs to display\nSelect categories on the right", 
                         ha='center', va='center', fontsize=12)
            ax_graph.axis('off')
            return
        
        # Crear layout
        if farmaco_principal and farmaco_principal in G_filtered.nodes():
            # Layout centrado en fármaco principal
            pos = {}
            pos[farmaco_principal] = np.array([0, 0])
            
            neighbors = list(G_filtered.predecessors(farmaco_principal)) + \
                       list(G_filtered.successors(farmaco_principal))
            neighbors = list(set(neighbors))  # Eliminar duplicados
            
            n_neighbors = len(neighbors)
            if n_neighbors > 0:
                radius = 1.5 + 0.2 * min(n_neighbors, 20)  # Limitar radio
                angle_step = 2 * np.pi / n_neighbors
                
                for i, neighbor in enumerate(neighbors):
                    angle = i * angle_step
                    x = radius * np.cos(angle)
                    y = radius * np.sin(angle)
                    pos[neighbor] = np.array([x, y])
        else:
            # Layout spring para grafo completo
            pos = nx.spring_layout(G_filtered, k=2/np.sqrt(len(G_filtered.nodes())), 
                                  iterations=50, seed=42)
        
        # Guardar posiciones para tooltips
        global node_positions
        node_positions = pos
        
        # Preparar datos para visualización
        node_colors = []
        node_sizes = []
        
        for node in G_filtered.nodes():
            node_colors.append(G_filtered.nodes[node]['color'])
            
            # Tamaño diferente para fármaco principal
            if farmaco_principal and node == farmaco_principal:
                node_sizes.append(1200)
            else:
                node_sizes.append(600)
        
        # Dibujar nodos
        nodes_collection = nx.draw_networkx_nodes(G_filtered, pos, 
                                      node_color=node_colors,
                                      node_size=node_sizes,
                                      alpha=0.9,
                                      edgecolors='black',
                                      linewidths=1,
                                      ax=ax_graph)
        
        # Preparar información de aristas para tooltips
        global edge_info_list
        edge_info_list = []
        edge_colors = []
        
        for u, v, data in G_filtered.edges(data=True):
            
            edge_colors.append('gray')
            
            edge_info_list.append({
                'u': u, 
                'v': v, 
                'data': data,
                'color': edge_colors[-1]
            })
        
        # Dibujar aristas - IMPORTANTE: dibujar todas juntas
        edges_collection = nx.draw_networkx_edges(G_filtered, pos,
                                      width=1.5,
                                      alpha=0.6,
                                      edge_color=edge_colors,
                                      arrows=True,
                                      arrowsize=10,
                                      arrowstyle='-|>',
                                      node_size=node_sizes,
                                      ax=ax_graph)
        
        # Dibujar etiquetas de nodos (nombres siempre visibles)
        labels = {}
        for node in G_filtered.nodes():
            if farmaco_principal and node == farmaco_principal:
                labels[node] = node  # nombre completo
            else:
                # acortar nombres largos
                if len(node) > 20:
                    labels[node] = node[:17] + "..."
                else:
                    labels[node] = node

        nx.draw_networkx_labels(
            G_filtered,
            pos,
            labels,
            font_size=6,              # más pequeño
            font_weight='normal',
            bbox=dict(                
                facecolor='white',
                edgecolor='none',
                alpha=0.6
            ),
            ax=ax_graph
        )
        
        # Información de estadísticas
        stats_text = (f"Drugs: {len(G_filtered.nodes())} | "
                     f"Interactions: {len(G_filtered.edges())}")
        
        if farmaco_principal:
            title = f"Drug: {farmaco_principal}\n{stats_text}"

        else:
            title = f"Complete Network\n{stats_text}"
        
        ax_graph.set_title(title, fontsize=12, pad=20,y=0.9)
        ax_graph.axis('off')
        
        # Configurar tooltips
        if MPLCURSORS_AVAILABLE:
            configurar_tooltips_mplcursors(G_filtered, nodes_collection, edges_collection)
        else:
            configurar_tooltips_alternativos(G_filtered, pos, nodes_collection, edges_collection)
        
        return G_filtered
    
    # Función para configurar tooltips con mplcursors (si está disponible)
    def configurar_tooltips_mplcursors(G_filtered, nodes_collection, edges_collection):
        try:
            
            # Configurar tooltips para nodos
            cursor_nodes = mplcursors.cursor(nodes_collection, hover=True)
            

            @cursor_nodes.connect("add")
            def on_add_node(sel):
                node_idx = sel.index
                nodes_list = list(G_filtered.nodes())
                if node_idx < len(nodes_list):
                    node_name = nodes_list[node_idx]
                    node_data = G_filtered.nodes[node_name]
                    
                    tooltip_text = (f"DRUG INFORMATION\n"
                                  f"Name: {node_name}\n"
                                  f"ATC Code: {node_data['atc_code']}\n"
                                  f"ATC Category: {node_data['atc_category']}")
                    
                    sel.annotation.set_text(tooltip_text)
                    sel.annotation.set_bbox(dict(
                        boxstyle="round,pad=0.5",
                        facecolor="lightyellow",
                        alpha=0.95,
                        edgecolor="orange"
                    ))
                    sel.annotation.set_fontsize(9)
                    auto_remove_annotation(sel.annotation, sel.annotation.figure.canvas, delay=5)

            
            # Configurar tooltips para aristas
            cursor_edges = mplcursors.cursor(edges_collection, hover=True)
            
            @cursor_edges.connect("add")
            def on_add_edge(sel):
                edge_idx = sel.index
                if edge_idx < len(edge_info_list):
                    edge_info = edge_info_list[edge_idx]
                    u = edge_info['u']
                    v = edge_info['v']
                    data = edge_info['data']
                    
                    interaction_type = data.get('interaction_type', 'Unknown')
                    interaction_desc = "Type 1" if interaction_type == 1 else "Type 2" if interaction_type == 2 else "Unknown"
                    
                    tooltip_text = (f"INTERACTION DETAILS\n"
                                  f"From: {u}\n"
                                  f"To: {v}\n"
                                  f"Interaction Type (Y): {interaction_desc}")
                    
                    sel.annotation.set_text(tooltip_text)
                    sel.annotation.set_bbox(dict(
                        boxstyle="round,pad=0.5",
                        facecolor="lightblue",
                        alpha=0.95,
                        edgecolor="blue"
                    ))
                    sel.annotation.set_fontsize(9)
                    
        except Exception as e:
            print(f"Error configurando tooltips: {e}")
            print("Usando tooltips alternativos...")
    
    # Función alternativa para tooltips (sin mplcursors)
    def configurar_tooltips_alternativos(G_filtered, pos, nodes_collection, edges_collection):
        """Tooltips alternativos usando eventos de matplotlib"""
        
        # Crear área para mostrar información
        info_box = ax_controls.text(0.05, 0.05, 
                                   "Hover over nodes/edges for info...",
                                   fontsize=8,
                                   transform=ax_controls.transAxes,
                                   bbox=dict(boxstyle="round,pad=0.5", 
                                            facecolor="lightyellow", 
                                            alpha=0.8))
        
        def on_motion(event):
            if event.inaxes != ax_graph:
                info_box.set_text("Hover over nodes/edges for info...")
                fig.canvas.draw_idle()
                return
            
            # Verificar si está sobre un nodo
            node_hit = None
            min_dist = 0.05
            
            for node, (x, y) in node_positions.items():
                dist = np.sqrt((event.xdata - x)**2 + (event.ydata - y)**2)
                if dist < min_dist:
                    node_hit = node
                    min_dist = dist
                    break
            
            if node_hit:
                node_data = G_filtered.nodes[node_hit]
                info_text = (f"DRUG: {node_hit}\n"
                           f"ATC Code: {node_data['atc_code']}\n"
                           f"ATC Category: {node_data['atc_category']}")
                info_box.set_text(info_text)
                fig.canvas.draw_idle()
                return
            
            # Verificar si está sobre una arista (aproximación)
            edge_hit = None
            for edge_info in edge_info_list:
                u, v = edge_info['u'], edge_info['v']
                if u in node_positions and v in node_positions:
                    x1, y1 = node_positions[u]
                    x2, y2 = node_positions[v]
                    
                    # Distancia punto-recta simplificada
                    line_length = np.sqrt((x2-x1)**2 + (y2-y1)**2)
                    if line_length == 0:
                        continue
                    
                    t = ((event.xdata - x1)*(x2 - x1) + (event.ydata - y1)*(y2 - y1)) / (line_length**2)
                    t = max(0, min(1, t))
                    
                    proj_x = x1 + t * (x2 - x1)
                    proj_y = y1 + t * (y2 - y1)
                    
                    dist = np.sqrt((event.xdata - proj_x)**2 + (event.ydata - proj_y)**2)
                    
                    if dist < 0.03:
                        edge_hit = edge_info
                        break
            
            if edge_hit:
                interaction_type = edge_hit['data'].get('interaction_type', 'Unknown')
                interaction_desc = "Type 1" if interaction_type == 1 else "Type 2" if interaction_type == 2 else "Unknown"
                
                info_text = (f"INTERACTION:\n"
                           f"From: {edge_hit['u']}\n"
                           f"To: {edge_hit['v']}\n"
                           f"Type (Y): {interaction_desc}")
                info_box.set_text(info_text)
                fig.canvas.draw_idle()
                return
            
            # Si no está sobre nada
            info_box.set_text("Hover over nodes/edges for info...")
            fig.canvas.draw_idle()
        
        # Conectar evento de movimiento del mouse
        fig.canvas.mpl_connect('motion_notify_event', on_motion)
    
    # Crear checkboxes
    checkbox_labels = [f"{cat}" for cat in category_list]
    checkbox_status = [active_categories[cat] for cat in category_list]
    
    # Ajustar posición de checkboxes
    y_positions = np.linspace(0.9, 0.1, len(checkbox_labels))
    
    # Dibujar título de checkboxes
    ax_checkboxes.text(0.1, 0.95, "Filter by ATC Category:", 
                      fontsize=11, fontweight='bold',
                      transform=ax_checkboxes.transAxes)
    
    # Crear checkboxes manualmente
    checkboxes = []
    for i, (label, status, y_pos) in enumerate(zip(checkbox_labels, checkbox_status, y_positions)):
        # Color del checkbox según categoría
        color = ATC_COLORS.get(label[:1], '#CCCCCC') if label != 'Sin ATC' and label != 'Multi ATC' else ATC_COLORS.get(label, '#CCCCCC')
        
        # Dibujar checkbox (blanco porque todos están desmarcados inicialmente)
        checkbox_rect = plt.Rectangle((0.1, y_pos-0.03), 0.05, 0.05, 
                                     facecolor='white',  # Blanco = desmarcado
                                     edgecolor='black',
                                     transform=ax_checkboxes.transAxes)
        ax_checkboxes.add_patch(checkbox_rect)
        
        # Texto de la categoría
        ax_checkboxes.text(0.2, y_pos, label, 
                          fontsize=9,
                          transform=ax_checkboxes.transAxes,
                          verticalalignment='center')
        
        # Guardar referencia
        checkboxes.append({
            'rect': checkbox_rect,
            'label': label,
            'y_pos': y_pos,
            'category': category_list[i],
            'color': color
        })
    
    # Función para manejar clics en checkboxes
    def onclick(event):
        if event.inaxes != ax_checkboxes:
            return
        
        # Verificar si el clic fue en algún checkbox
        for checkbox in checkboxes:
            rect = checkbox['rect']
            
            # Coordenadas del rectángulo
            x_min, y_min = rect.get_xy()
            x_max = x_min + rect.get_width()
            y_max = y_min + rect.get_height()
            
            # Transformar coordenadas
            x_min_t, y_min_t = ax_checkboxes.transAxes.transform((x_min, y_min))
            x_max_t, y_max_t = ax_checkboxes.transAxes.transform((x_max, y_max))
            
            # Verificar si el clic está dentro del rectángulo
            if (x_min_t <= event.x <= x_max_t and y_min_t <= event.y <= y_max_t):
                # Cambiar estado
                category = checkbox['category']
                new_state = not active_categories[category]
                active_categories[category] = new_state
                
                # Actualizar color del checkbox
                if new_state:  # Marcado
                    checkbox['rect'].set_facecolor(checkbox['color'])
                else:  # Desmarcado
                    checkbox['rect'].set_facecolor('white')
                
                # Redibujar gráfico
                dibujar_grafo()
                fig.canvas.draw_idle()
                break
    
    # Conectar evento de clic
    fig.canvas.mpl_connect('button_press_event', onclick)
    
   
    # Función para seleccionar/deseleccionar todo
    def select_all(event):
        for checkbox in checkboxes:
            category = checkbox['category']
            active_categories[category] = True
            checkbox['rect'].set_facecolor(checkbox['color'])
        
        dibujar_grafo()
        fig.canvas.draw_idle()
    
    def deselect_all(event):
        for checkbox in checkboxes:
            category = checkbox['category']
            active_categories[category] = False
            checkbox['rect'].set_facecolor('white')
        
        dibujar_grafo()
        fig.canvas.draw_idle()
    
    # Agregar botones (simulados con texto clickeable)
    select_all_text = ax_controls.text(0.1, 0.25, "✓ Select All", 
                                      fontsize=10, color='blue',
                                      bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue"),
                                      transform=ax_controls.transAxes)
    
    deselect_all_text = ax_controls.text(0.1, 0.15, "✗ Deselect All", 
                                        fontsize=10, color='red',
                                        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightcoral"),
                                        transform=ax_controls.transAxes)
    
    # Conectar eventos a los "botones" de texto
    select_all_text.set_picker(True)
    deselect_all_text.set_picker(True)
    
    def on_pick(event):
        if event.artist == select_all_text:
            select_all(event)
        elif event.artist == deselect_all_text:
            deselect_all(event)
    
    fig.canvas.mpl_connect('pick_event', on_pick)
    
    # Dibujar gráfico inicial (vacío ya que todas las categorías están desmarcadas)
    dibujar_grafo()
    
    plt.tight_layout()
    plt.show()
    
    # Mostrar estadísticas iniciales
    print("="*60)
    print(f"DRUG INTERACTION NETWORK WITH TOOLTIPS")
    print("="*60)
    
    if farmaco_principal:
        print(f"\nMain Drug: {farmaco_principal}")
    
    print(f"\nTotal drugs in dataset: {len(G.nodes())}")
    print(f"Total interactions in dataset: {len(G.edges())}")
    
    # Contar por categoría
    category_counts = Counter()
    for node in G.nodes():
        if node != farmaco_principal:
            category_counts[G.nodes[node]['atc_category']] += 1
    
    print(f"\nDrugs by ATC Category:")
    for category, count in sorted(category_counts.items()):
        print(f"  {category}: {count} drugs")
    
    print(f"\n🔍 INTERACTIVE FEATURES:")
    if MPLCURSORS_AVAILABLE:
        print("  • Hover over NODES to see drug details (name, ATC code, category)")
        print("  • Hover over EDGES to see interaction type (Y value)")
    else:
        print("  • Move mouse over NODES/EDGES to see information at bottom")
    print("  • Click checkboxes to filter by ATC category")
    print("  • Start with all categories deselected for better performance")
    print("  • Use 'Select All' to show everything")
    print("="*60)

# Función principal para probar
if __name__ == "__main__":
    print("Loading drug interaction data...")
    
    # Probar con phenobarbital
    target_drug = "verteporfin"
    
    print(f"\nAnalyzing interactions for: {target_drug}")
    
    if not MPLCURSORS_AVAILABLE:
        print("Note: mplcursors not installed. Using alternative tooltips.")
        print("To install: pip install mplcursors")
    
    # Usar la versión con tooltips

    visualizar_grafo_con_tooltips(df, target_drug)
