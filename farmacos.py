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

@st.cache_data
def load_data2():
    return pd.read_csv(r"drugbank_ddi_label_map.csv")
meeaning_y=load_data2()

# Cargar datos de fármacos esenciales
@st.cache_data
def load_essential_drugs():
    
        
    essentials_df = pd.read_csv(r"farmacos_esenciales.csv")
    return essentials_df
    

# Resto de tu código original (sin cambios hasta las funciones de pestañas)
# ... [todo tu código original de visualización se mantiene igual] ...

def pestaña_visualizacion():
    # ... [tu código original de visualización] ...
    pass

def pestaña_significado_y():
    st.title("Dataset, meaning of Y values")
    st.header('Dataset')
    st.dataframe(meeaning_y, use_container_width=True)

# =================================================================
# PESTAÑA DE FÁRMACOS ESENCIALES - VERSIÓN CORREGIDA
# =================================================================

def pestaña_esenciales():
    """Pestaña simple para fármacos esenciales"""
    
    st.title("🌍 Essential Drugs")
    
    # Cargar datos
    essentials_df = load_essential_drugs()
    
    if essentials_df.empty:
        st.info("Please install openpyxl to load essential drugs data.")
        return
    
    # Obtener lista de países
    non_country_cols = ['Medicine', 'ATC code primary', 'ATC code secondary', 
                       'WHO Model List of Essential Medications']
    paises = [col for col in essentials_df.columns if col not in non_country_cols]
    
    if not paises:
        st.warning("No countries found in the dataset.")
        return
    
    # Seleccionar país
    st.subheader("Select Country")
    pais_seleccionado = st.selectbox(
        "Choose a country:",
        paises,
        index=paises.index('Mexico') if 'Mexico' in paises else 0
    )
    
    # Filtrar fármacos esenciales del país
    farmacos_pais = essentials_df[essentials_df[pais_seleccionado] == 1]
    
    if farmacos_pais.empty:
        st.info(f"No essential drugs found for {pais_seleccionado}")
        return
    
    # Mostrar estadísticas iniciales
    st.subheader(f"Statistics for {pais_seleccionado}")
    
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Total Essential Drugs", len(farmacos_pais))
    
    # =============================================================
    # PASO 1: Obtener todos los códigos ATC esenciales (como tu código)
    # =============================================================
    codigos_atc_esenciales = set()
    
    # Agregar códigos ATC primarios
    codigos_atc_esenciales.update(farmacos_pais['ATC code primary'].dropna().unique())
    
    # Agregar códigos ATC secundarios
    for codigos in farmacos_pais['ATC code secondary'].dropna():
        if isinstance(codigos, str):
            if '|' in codigos:
                for codigo in codigos.split('|'):
                    codigos_atc_esenciales.add(codigo.strip())
            else:
                codigos_atc_esenciales.add(codigos.strip())
    
    # Limpiar espacios en blanco
    codigos_atc_esenciales = {codigo.strip() for codigo in codigos_atc_esenciales 
                              if pd.notna(codigo) and codigo.strip() != ''}
    
    with col2:
        st.metric("Unique ATC Codes", len(codigos_atc_esenciales))
    
    # =============================================================
    # PASO 2: Verificar fármacos en DDIBUENO (como tu código)
    # =============================================================
    def verificar_farmaco(atc_code, num_atc):
        if pd.isna(atc_code) or atc_code == 'No ATC' or num_atc == 0:
            return False
        
        if '|' in str(atc_code):
            codigos = [c.strip() for c in str(atc_code).split('|')]
            return any(codigo in codigos_atc_esenciales for codigo in codigos)
        else:
            return str(atc_code).strip() in codigos_atc_esenciales
    
    # Aplicar verificación a ambas columnas
    df['esencial_x'] = df.apply(lambda row: verificar_farmaco(row['atc_code_x'], row['num_atc_x']), axis=1)
    df['esencial_y'] = df.apply(lambda row: verificar_farmaco(row['atc_code_y'], row['num_atc_y']), axis=1)
    
    # Contar cuántos son esenciales
    esencial_x_count = df['esencial_x'].sum()
    esencial_y_count = df['esencial_y'].sum()
    
    with col3:
        st.metric("Drugs in DDIBUENO", esencial_x_count + esencial_y_count)
    
    # =============================================================
    # PASO 3: Crear lista de fármacos esenciales encontrados
    # =============================================================
    farmacos_esenciales_lista = []
    
    # Procesar fármacos de la columna X
    for idx, row in df[df['esencial_x']].iterrows():
        farmacos_esenciales_lista.append({
            'Common_Name': row['Common_name_x'],
            'ATC_Code': row['atc_code_x'],
            'Num_ATC': row['num_atc_x'],
            'Source': 'Column X',
            'Original_Index': idx,
            'Drug_Pair': f"{row['Common_name_x']} - {row['Common_name_y']}",
            'Y_Value': row['Y']
        })
    
    # Procesar fármacos de la columna Y
    for idx, row in df[df['esencial_y']].iterrows():
        farmacos_esenciales_lista.append({
            'Common_Name': row['Common_name_y'],
            'ATC_Code': row['atc_code_y'],
            'Num_ATC': row['num_atc_y'],
            'Source': 'Column Y',
            'Original_Index': idx,
            'Drug_Pair': f"{row['Common_name_x']} - {row['Common_name_y']}",
            'Y_Value': row['Y']
        })
    
    # Crear DataFrame
    if not farmacos_esenciales_lista:
        st.info(f"No essential drugs from {pais_seleccionado} found in DDIBUENO database")
        return
    
    df_farmacos_esenciales = pd.DataFrame(farmacos_esenciales_lista)
    
    # Eliminar duplicados (mismo nombre y ATC)
    df_farmacos_esenciales = df_farmacos_esenciales.drop_duplicates(subset=['Common_Name', 'ATC_Code'])
    
    # Ordenar por nombre
    df_farmacos_esenciales = df_farmacos_esenciales.sort_values('Common_Name')
    
    # =============================================================
    # PASO 4: Mostrar resultados principales
    # =============================================================
    st.subheader(f"Essential Drugs from {pais_seleccionado} found in DDIBUENO")
    
    # Mostrar tabla principal
    st.dataframe(
        df_farmacos_esenciales[['Common_Name', 'ATC_Code', 'Source', 'Y_Value']],
        use_container_width=True,
        column_config={
            "Common_Name": st.column_config.TextColumn("Drug Name", width="medium"),
            "ATC_Code": st.column_config.TextColumn("ATC Code", width="small"),
            "Source": st.column_config.TextColumn("Found in", width="small"),
            "Y_Value": st.column_config.NumberColumn("Y Value", width="small")
        },
        hide_index=True
    )
    
    # Mostrar estadísticas detalladas
    st.subheader("Detailed Statistics")
    
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Unique Essential Drugs", len(df_farmacos_esenciales))
    with col2:
        st.metric("From Column X", len(df_farmacos_esenciales[df_farmacos_esenciales['Source'] == 'Column X']))
    with col3:
        st.metric("From Column Y", len(df_farmacos_esenciales[df_farmacos_esenciales['Source'] == 'Column Y']))
    with col4:
        unique_drugs = set(df_farmacos_esenciales['Common_Name'])
        st.metric("Unique Drug Names", len(unique_drugs))
    
    # =============================================================
    # PASO 5: Mostrar distribución por grupo ATC
    # =============================================================
    st.subheader("Distribution by ATC Group")
    
    def extraer_grupo_atc(atc_code):
        if pd.isna(atc_code) or atc_code == 'No ATC':
            return 'Sin ATC'
        if isinstance(atc_code, str) and len(atc_code) > 0:
            return atc_code[0]
        return 'Unknown'
    
    df_farmacos_esenciales['ATC_Group'] = df_farmacos_esenciales['ATC_Code'].apply(extraer_grupo_atc)
    
    # Diccionario de nombres de grupos ATC
    ATC_GROUP_NAMES = {
        'A': 'Alimentary Tract',
        'B': 'Blood',
        'C': 'Cardiovascular',
        'D': 'Dermatologicals',
        'G': 'Genito Urinary',
        'H': 'Hormones',
        'J': 'Antiinfectives',
        'L': 'Antineoplastic',
        'M': 'Musculoskeletal',
        'N': 'Nervous System',
        'P': 'Antiparasitic',
        'R': 'Respiratory',
        'S': 'Sensory Organs',
        'V': 'Various',
        'Sin ATC': 'No ATC',
        'Unknown': 'Unknown'
    }
    
    # Crear DataFrame para distribución
    distribucion = df_farmacos_esenciales['ATC_Group'].value_counts().reset_index()
    distribucion.columns = ['ATC_Group', 'Count']
    distribucion['Group_Name'] = distribucion['ATC_Group'].map(ATC_GROUP_NAMES)
    
    # Mostrar gráfico de barras
    fig = px.bar(
        distribucion,
        x='Group_Name',
        y='Count',
        title=f"Essential Drugs by ATC Group - {pais_seleccionado}",
        color='Count',
        color_continuous_scale='Blues',
        labels={'Group_Name': 'ATC Group', 'Count': 'Number of Drugs'}
    )
    st.plotly_chart(fig, use_container_width=True)
    
    # Mostrar tabla de distribución
    with st.expander("View Distribution Table", expanded=False):
        st.dataframe(
            distribucion,
            use_container_width=True,
            hide_index=True
        )
    
    # =============================================================
    # PASO 6: Mostrar interacciones de estos fármacos
    # =============================================================
    st.subheader("Interactions Involving Essential Drugs")
    
    # Encontrar todas las interacciones que involucran fármacos esenciales
    essential_drugs_names = set(df_farmacos_esenciales['Common_Name'])
    
    # Filtrar interacciones donde al menos uno de los fármacos es esencial
    essential_interactions = df[
        df['Common_name_x'].isin(essential_drugs_names) | 
        df['Common_name_y'].isin(essential_drugs_names)
    ]
    
    if not essential_interactions.empty:
        st.write(f"**Found {len(essential_interactions)} interactions involving essential drugs**")
        
        # Mostrar tabla de interacciones
        with st.expander("View All Interactions", expanded=False):
            st.dataframe(
                essential_interactions[['Common_name_x', 'Common_name_y', 'Y', 'atc_code_x', 'atc_code_y']].rename(
                    columns={
                        'Common_name_x': 'Drug A',
                        'Common_name_y': 'Drug B',
                        'Y': 'Y Value',
                        'atc_code_x': 'ATC Code A',
                        'atc_code_y': 'ATC Code B'
                    }
                ),
                use_container_width=True,
                hide_index=True
            )
        
        # Análisis de severidad
        st.write("**Severity Analysis (Y Values):**")
        
        y_counts = essential_interactions['Y'].value_counts().sort_index()
        
        col1, col2 = st.columns(2)
        
        with col1:
            for y_val, count in y_counts.items():
                percentage = (count / len(essential_interactions)) * 100
                severity = {
                    0: "None", 1: "Mild", 2: "Moderate", 
                    3: "Severe", 4: "Contraindicated"
                }.get(y_val, "Unknown")
                
                st.write(f"**Y={y_val} ({severity}):** {count} ({percentage:.1f}%)")
        
        with col2:
            fig_pie = px.pie(
                values=y_counts.values,
                names=[f"Y={k}" for k in y_counts.index],
                title="Severity Distribution",
                color_discrete_sequence=px.colors.sequential.Reds
            )
            st.plotly_chart(fig_pie, use_container_width=True)
        
        # Top fármacos más interactivos
        st.write("**Most Interactive Essential Drugs:**")
        
        interaction_counts = {}
        for _, row in essential_interactions.iterrows():
            if row['Common_name_x'] in essential_drugs_names:
                interaction_counts[row['Common_name_x']] = interaction_counts.get(row['Common_name_x'], 0) + 1
            if row['Common_name_y'] in essential_drugs_names:
                interaction_counts[row['Common_name_y']] = interaction_counts.get(row['Common_name_y'], 0) + 1
        
        if interaction_counts:
            top_drugs = pd.DataFrame(
                interaction_counts.items(),
                columns=['Drug', 'Interaction_Count']
            ).sort_values('Interaction_Count', ascending=False).head(10)
            
            st.dataframe(
                top_drugs,
                use_container_width=True,
                hide_index=True
            )
    
    else:
        st.info("No interactions found involving these essential drugs")
    
    # =============================================================
    # PASO 7: Opción para descargar datos
    # =============================================================
    st.subheader("Download Data")
    
    col1, col2 = st.columns(2)
    
    with col1:
        # Descargar lista de fármacos esenciales
        csv_farmacos = df_farmacos_esenciales.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="📥 Download Essential Drugs List",
            data=csv_farmacos,
            file_name=f"essential_drugs_{pais_seleccionado}.csv",
            mime="text/csv"
        )
    
    with col2:
        if not essential_interactions.empty:
            # Descargar interacciones
            csv_interacciones = essential_interactions.to_csv(index=False).encode('utf-8')
            st.download_button(
                label="📥 Download Interactions",
                data=csv_interacciones,
                file_name=f"essential_interactions_{pais_seleccionado}.csv",
                mime="text/csv"
            )

# =================================================================
# MAIN FUNCTION
# =================================================================

def main():
    # Crear pestañas
    tab1, tab2, tab3 = st.tabs([
        "🔍 Drug Interactions", 
        "📊 Y Values", 
        "🌍 Essential Drugs"
    ])
    
    with tab1:
        pestaña_visualizacion()
    
    with tab2:
        pestaña_significado_y()
    
    with tab3:
        pestaña_esenciales()

if __name__ == "__main__":
    main()
