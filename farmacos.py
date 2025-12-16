# app.py - VERSIÓN STREAMLIT FUNCIONAL
import streamlit as st
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

# Configurar la página PRIMERO
st.set_page_config(
    page_title="Drug Interaction Network",
    page_icon="💊",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Título principal
st.title("💊 Drug Interaction Network Explorer")
st.markdown("---")

# Cargar datos CON CACHÉ (IMPORTANTE)
@st.cache_data
def load_data():
    try:
        df = pd.read_csv("DDIBUENO.csv")
        st.success(f"✅ Data loaded: {len(df)} interactions")
        return df
    except Exception as e:
        st.error(f"❌ Error loading CSV: {e}")
        return None

# Cargar datos
df = load_data()

if df is None:
    st.stop()

# Sidebar para controles
with st.sidebar:
    st.header("⚙️ Configuration")
    
    # Obtener lista de fármacos
    all_drugs = sorted(set(list(df['Common_name_x'].unique()) + 
                          list(df['Common_name_y'].unique())))
    
    target_drug = st.selectbox(
        "Select Target Drug",
        options=["(Show all)"] + all_drugs,
        index=0
    )
    
    st.markdown("---")
    st.header("🎨 Filter by ATC Category")
    
    # Checkbox simple para empezar
    show_filter = st.checkbox("Enable ATC Filtering", value=False)

# Contenido principal
col1, col2 = st.columns([3, 1])

with col1:
    st.subheader("📊 Network Visualization")
    
    # Crear grafo simple para probar
    if target_drug == "(Show all)":
        st.info("Showing complete network")
        # Crear grafo completo (limitado para demo)
        G = nx.DiGraph()
        sample_df = df.head(50)  # Limitar para prueba
        
        for _, row in sample_df.iterrows():
            G.add_edge(row['Common_name_x'], row['Common_name_y'])
        
        # Dibujar
        fig, ax = plt.subplots(figsize=(10, 8))
        pos = nx.spring_layout(G, seed=42)
        nx.draw(G, pos, with_labels=True, ax=ax, node_color='lightblue',
                node_size=500, font_size=8, arrows=True)
        ax.set_title("Drug Interaction Network")
        
        st.pyplot(fig)
        
    else:
        st.info(f"Showing interactions for: **{target_drug}**")
        # Filtrar interacciones del fármaco seleccionado
        filtered = df[(df['Common_name_x'] == target_drug) | 
                     (df['Common_name_y'] == target_drug)]
        
        if len(filtered) > 0:
            G = nx.DiGraph()
            for _, row in filtered.iterrows():
                G.add_edge(row['Common_name_x'], row['Common_name_y'])
            
            fig, ax = plt.subplots(figsize=(10, 8))
            pos = nx.spring_layout(G, seed=42)
            nx.draw(G, pos, with_labels=True, ax=ax, node_color='lightblue',
                    node_size=800, font_size=10, arrows=True)
            ax.set_title(f"Interactions for {target_drug}")
            
            st.pyplot(fig)
            st.metric("Interactions found", len(filtered))
        else:
            st.warning(f"No interactions found for {target_drug}")

with col2:
    st.subheader("ℹ️ Information")
    st.metric("Total drugs in dataset", len(all_drugs))
    st.metric("Total interactions", len(df))
    
    if target_drug != "(Show all)":
        # Mostrar información básica del fármaco
        st.info(f"""
        **Selected Drug:** {target_drug}
        
        **Interactions:** {
            len(df[(df['Common_name_x'] == target_drug) | 
                  (df['Common_name_y'] == target_drug)])
        }
        """)

# Footer
st.markdown("---")
st.caption("Drug Interaction Network Explorer | Powered by Streamlit & NetworkX")

