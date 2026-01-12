# pages/3_📈_Estadisticas.py
import streamlit as st
import pandas as pd
import plotly.express as px

# Configuración de la página
st.set_page_config(layout="wide")

# Título de la página
st.title("📈 Estadísticas de Interacciones de Fármacos")

# 1. CARGAR DATOS
@st.cache_data
def cargar_datos():
    return pd.read_csv(r"DDIBUENO.csv")

df = cargar_datos()

# 2. MOSTRAR ESTADÍSTICAS BÁSICAS
st.header("📊 Resumen General del Dataset")

col1, col2, col3, col4 = st.columns(4)

with col1:
    total_interacciones = len(df)
    st.metric("Total Interacciones", total_interacciones)

with col2:
    farmacos_unicos = len(set(df['Common_name_x'].tolist() + df['Common_name_y'].tolist()))
    st.metric("Fármacos Únicos", farmacos_unicos)

with col3:
    st.metric("Columnas en Dataset", len(df.columns))

with col4:
    st.metric("Valores Y Diferentes", df['Y'].nunique())

# 3. DISTRIBUCIÓN DE VALORES Y
st.header("📈 Distribución de Valores Y")

# Gráfico de barras
fig = px.histogram(df, x='Y', 
                   title="Frecuencia de Valores Y",
                   labels={'Y': 'Valor Y', 'count': 'Frecuencia'},
                   color='Y',
                   color_discrete_sequence=px.colors.qualitative.Set3)

st.plotly_chart(fig, use_container_width=True)

# 4. TOP FÁRMACOS MÁS INTERACTIVOS
st.header("🏆 Top 10 Fármacos con Más Interacciones")

# Contar interacciones por fármaco
interacciones_salientes = df['Common_name_x'].value_counts()
interacciones_entrantes = df['Common_name_y'].value_counts()

# Combinar
total_interacciones_farmaco = (interacciones_salientes + interacciones_entrantes).sort_values(ascending=False)

# Crear dataframe
top_farmacos = pd.DataFrame({
    'Fármaco': total_interacciones_farmaco.head(10).index,
    'Total Interacciones': total_interacciones_farmaco.head(10).values,
    'Como Drug A': [interacciones_salientes.get(f, 0) for f in total_interacciones_farmaco.head(10).index],
    'Como Drug B': [interacciones_entrantes.get(f, 0) for f in total_interacciones_farmaco.head(10).index]
})

# Mostrar tabla
st.dataframe(top_farmacos, use_container_width=True, hide_index=True)

# 5. GRÁFICO DE TOP FÁRMACOS
fig2 = px.bar(top_farmacos, x='Fármaco', y='Total Interacciones',
              title="Top 10 Fármacos más Interactivos",
              color='Total Interacciones',
              color_continuous_scale='Viridis')
st.plotly_chart(fig2, use_container_width=True)

# 6. ANÁLISIS POR ATC (si tu dataset tiene ATC)
st.header("🧪 Análisis por Categoría ATC")

if 'atc_code_x' in df.columns:
    # Extraer primera letra ATC
    df['ATC_Group'] = df['atc_code_x'].astype(str).str[0]
    
    # Contar por grupo ATC
    atc_counts = df['ATC_Group'].value_counts().reset_index()
    atc_counts.columns = ['Grupo ATC', 'Cantidad']
    
    # Gráfico
    fig3 = px.pie(atc_counts, names='Grupo ATC', values='Cantidad',
                  title="Distribución por Grupo ATC")
    st.plotly_chart(fig3, use_container_width=True)
else:
    st.info("El dataset no contiene información ATC para este análisis")

# 7. BÚSQUEDA ESPECÍFICA
st.header("🔍 Buscar Estadísticas por Fármaco")

farmaco_buscar = st.selectbox(
    "Selecciona un fármaco para ver sus estadísticas:",
    sorted(set(df['Common_name_x'].tolist() + df['Common_name_y'].tolist()))
)

if farmaco_buscar:
    # Calcular estadísticas para ese fármaco
    como_a = df[df['Common_name_x'] == farmaco_buscar]
    como_b = df[df['Common_name_y'] == farmaco_buscar]
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.metric("Como Drug A (afecta a otros)", len(como_a))
        st.metric("Valor Y promedio como A", como_a['Y'].mean().round(2) if len(como_a) > 0 else 0)
    
    with col2:
        st.metric("Como Drug B (es afectado)", len(como_b))
        st.metric("Valor Y promedio como B", como_b['Y'].mean().round(2) if len(como_b) > 0 else 0)
    
    st.metric("Total Interacciones", len(como_a) + len(como_b))

# 8. EXPORTAR DATOS
st.header("💾 Exportar Datos")

with st.expander("Exportar estadísticas"):
    st.download_button(
        label="📥 Descargar estadísticas de top fármacos (CSV)",
        data=top_farmacos.to_csv(index=False).encode('utf-8'),
        file_name="top_farmacos_estadisticas.csv",
        mime="text/csv"
    )
    
    st.download_button(
        label="📥 Descargar distribución de valores Y (CSV)",
        data=df['Y'].value_counts().reset_index().to_csv(index=False).encode('utf-8'),
        file_name="distribucion_valores_y.csv",
        mime="text/csv"
    )

# 9. LINK DE REGRESO
st.markdown("---")
st.page_link("pages/1_🏠_Home.py", label="🏠 Volver a la Visualización Principal", icon="🏠")