# %% Librerias
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams["figure.dpi"] = 150     # resolución en pantalla
plt.rcParams["savefig.dpi"] = 300    # resolución al guardar con savefig

# %% Cargar las tablas limpias
cinetica = pd.read_csv("../datos/procesados/cinetica.csv", dtype={"pubmed": str})
inhibicion = pd.read_csv("../datos/procesados/inhibicion.csv", dtype={"pubmed": str})

print("Cinética:", cinetica.shape)
print("Inhibición:", inhibicion.shape)

# %% Variables en escala logarítmica
cinetica["log_km"] = np.log10(cinetica["km"])
cinetica["log_kcat"] = np.log10(cinetica["kcat"])
cinetica["log_kcat_km"] = np.log10(cinetica["kcat_km"])
inhibicion["log_ki"] = np.log10(inhibicion["ki"])

# %% Estadísticos descriptivos de las variables numéricas
logs = pd.concat([cinetica[["log_km", "log_kcat", "log_kcat_km"]],
                  inhibicion[["log_ki"]]], axis=1)

resumen = logs.describe().T
resumen["asimetria"] = logs.skew()
resumen["mediana_original"] = 10 ** resumen["50%"]
resumen.round(2)
# %% Estadísticos del Km
km = cinetica["log_km"].dropna()

media = km.mean()
desv = km.std()
q1 = km.quantile(0.25)
mediana = km.quantile(0.50)
q3 = km.quantile(0.75)

print(f"Media:    {media:.2f}  ->  {10**media:.0f} µM")
print(f"Mediana:  {mediana:.2f}  ->  {10**mediana:.0f} µM")
print(f"Q1:       {q1:.2f}  ->  {10**q1:.0f} µM")
print(f"Q3:       {q3:.2f}  ->  {10**q3:.0f} µM")
print(f"Desv. estándar: {desv:.2f}")
# %% Km: cada medición como un punto
fig, eje = plt.subplots(figsize=(11, 3))

sns.stripplot(x=km, color="steelblue", alpha=0.45, size=4, jitter=0.3,
              edgecolor="white", linewidth=0.3, ax=eje)

eje.axvline(mediana, color="crimson", linewidth=2, label=f"Mediana = {mediana:.2f} ({10**mediana:.0f} µM)")
eje.axvline(media, color="goldenrod", linestyle="--", linewidth=2, label=f"Media = {media:.2f} ({10**media:.0f} µM)")
eje.axvline(q1, color="darkorange", linestyle=":", linewidth=2, label=f"Q1 = {q1:.2f} ({10**q1:.0f} µM)")
eje.axvline(q3, color="darkorange", linestyle=":", linewidth=2, label=f"Q3 = {q3:.2f} ({10**q3:.0f} µM)")

eje.set_xlabel("log₁₀ Km (µM)")
eje.set_title(f"Km: cada punto es una medición (n = {len(km)})")
eje.legend(fontsize=8, loc="upper left")
plt.show()
# %% Km: la desviación estándar
fig, eje = plt.subplots(figsize=(11, 3.5))

# Franjas de ±2 y ±1 desviación estándar alrededor de la media
eje.axvspan(media - 2*desv, media + 2*desv, color="gold", alpha=0.12,
            label=f"Media ± 2 desv. ({10**(media-2*desv):.1f} a {10**(media+2*desv):.0f} µM)")
eje.axvspan(media - desv, media + desv, color="gold", alpha=0.35,
            label=f"Media ± 1 desv. ({10**(media-desv):.0f} a {10**(media+desv):.0f} µM)")

sns.stripplot(x=km, color="steelblue", alpha=0.3, size=3, jitter=0.3, ax=eje)
eje.axvline(media, color="goldenrod", linestyle="--", linewidth=2, label=f"Media = {media:.2f}")

eje.set_yticks([])
eje.set_xlabel("log₁₀ Km (µM)")
eje.set_title("Km: la desviación estándar alrededor de la media")
eje.legend(fontsize=8, loc="upper left")
plt.show()

# ¿Qué porcentaje de los puntos cae dentro de cada franja?
dentro_1 = ((km >= media - desv) & (km <= media + desv)).mean() * 100
dentro_2 = ((km >= media - 2*desv) & (km <= media + 2*desv)).mean() * 100
print(f"Dentro de ±1 desviación: {dentro_1:.1f} % de los puntos")
print(f"Dentro de ±2 desviaciones: {dentro_2:.1f} % de los puntos")
# %% Km: valores atípicos (regla del IQR)
iqr = q3 - q1
limite_bajo = q1 - 1.5 * iqr
limite_alto = q3 + 1.5 * iqr

print(f"Límites: {limite_bajo:.2f} a {limite_alto:.2f}  "
      f"({10**limite_bajo:.3f} a {10**limite_alto:.0f} µM)")

atipicos_km = cinetica[(cinetica["log_km"] < limite_bajo) | (cinetica["log_km"] > limite_alto)]
print("Atípicos bajos:", (cinetica["log_km"] < limite_bajo).sum())
print("Atípicos altos:", (cinetica["log_km"] > limite_alto).sum())

# ¿Qué tienen en común?
print("\nClase de antibiótico:")
print(atipicos_km["clase_antibiotico"].value_counts())
print("\nClase de enzima:")
print(atipicos_km["clase_enzima"].value_counts())
print("\nSilvestre o mutante:")
print(atipicos_km["tipo_enzima"].value_counts())
print("\nArtículos distintos:", atipicos_km["pubmed"].nunique())

atipicos_km.sort_values("log_km")[["log_km", "km", "antibiotico", "clase_enzima",
                                   "tipo_enzima", "mutacion", "pubmed"]].head(15)

# %% Regla del IQR (rango intercuartílico) para el Km
q1 = cinetica["log_km"].quantile(0.25)
q3 = cinetica["log_km"].quantile(0.75)
iqr = q3 - q1

limite_bajo = q1 - 1.5 * iqr
limite_alto = q3 + 1.5 * iqr

atipicos_bajos = cinetica[cinetica["log_km"] < limite_bajo]
atipicos_altos = cinetica[cinetica["log_km"] > limite_alto]

print(f"Q1 = {q1:.2f} | Q3 = {q3:.2f} | IQR = {iqr:.2f}")
print(f"Límites: {limite_bajo:.2f} a {limite_alto:.2f}")
print(f"Atípicos bajos: {len(atipicos_bajos)} | Atípicos altos: {len(atipicos_altos)}")


# %% Exploración automática: gráfico, estadísticos y atípicos de cada variable
from pathlib import Path

carpeta = Path("../figuras/exploracion")
carpeta.mkdir(parents=True, exist_ok=True)

def revisar_variable(tabla, columna, nombre, unidad, agrupar):
    datos = tabla[columna].dropna()
    media, desv = datos.mean(), datos.std()
    q1, mediana, q3 = datos.quantile([0.25, 0.5, 0.75])
    iqr = q3 - q1
    bajo, alto = q1 - 1.5 * iqr, q3 + 1.5 * iqr

    # --- Gráfico ---
    fig, eje = plt.subplots(figsize=(11, 3.5))
    eje.axvspan(media - desv, media + desv, color="gold", alpha=0.3, label="Media ± 1 desv.")
    sns.stripplot(x=datos, color="steelblue", alpha=0.45, size=4, jitter=0.3,
                  edgecolor="white", linewidth=0.3, ax=eje)
    eje.axvline(mediana, color="crimson", linewidth=2,
                label=f"Mediana = {mediana:.2f} ({10**mediana:.3g} {unidad})")
    eje.axvline(media, color="goldenrod", linestyle="--", linewidth=2, label=f"Media = {media:.2f}")
    eje.axvline(q1, color="darkorange", linestyle=":", linewidth=2, label=f"Q1 = {q1:.2f}")
    eje.axvline(q3, color="darkorange", linestyle=":", linewidth=2, label=f"Q3 = {q3:.2f}")
    eje.axvline(bajo, color="gray", linestyle="-.", linewidth=1, label="Límites de atípicos (IQR)")
    eje.axvline(alto, color="gray", linestyle="-.", linewidth=1)
    eje.set_yticks([])
    eje.set_xlabel(f"log₁₀ {nombre} ({unidad})")
    eje.set_title(f"{nombre} (n = {len(datos)})")
    eje.legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(carpeta / f"{columna}.png")
    plt.show()

    # --- Estadísticos ---
    print(f"===== {nombre} =====")
    print(f"n = {len(datos)} | media = {media:.2f} | desv = {desv:.2f} | asimetría = {datos.skew():.2f}")
    print(f"Q1 = {q1:.2f} ({10**q1:.3g}) | mediana = {mediana:.2f} ({10**mediana:.3g}) "
          f"| Q3 = {q3:.2f} ({10**q3:.3g}) {unidad}")
    dentro_1 = ((datos >= media - desv) & (datos <= media + desv)).mean() * 100
    dentro_2 = ((datos >= media - 2*desv) & (datos <= media + 2*desv)).mean() * 100
    print(f"Dentro de ±1 desv: {dentro_1:.1f} % | dentro de ±2 desv: {dentro_2:.1f} %")

    # --- Atípicos ---
    bajos = tabla[tabla[columna] < bajo]
    altos = tabla[tabla[columna] > alto]
    print(f"Límites IQR: {bajo:.2f} a {alto:.2f} ({10**bajo:.3g} a {10**alto:.3g} {unidad})")
    print(f"Atípicos: {len(bajos)} bajos | {len(altos)} altos "
          f"({(len(bajos) + len(altos)) / len(datos) * 100:.1f} % del total)")
    for etiqueta, grupo in [("BAJOS", bajos), ("ALTOS", altos)]:
        if len(grupo) > 0:
            print(f"  {etiqueta} — {grupo['pubmed'].nunique()} artículos distintos")
            for col in agrupar:
                print(f"    {col}: {grupo[col].value_counts().head(4).to_dict()}")
    print()


cinetica_cols = ["clase_enzima", "clase_antibiotico", "tipo_enzima"]
variables = [(cinetica, "log_km", "Km", "µM", cinetica_cols),
             (cinetica, "log_kcat", "kcat", "s⁻¹", cinetica_cols),
             (cinetica, "log_kcat_km", "kcat/Km", "M⁻¹ s⁻¹", cinetica_cols),
             (inhibicion, "log_ki", "Ki", "µM", ["tipo_inhibidor", "clase_enzima", "inhibidores"])]

for tabla, columna, nombre, unidad, agrupar in variables:
    revisar_variable(tabla, columna, nombre, unidad, agrupar)
# %% Correlaciones entre las variables numéricas (solo enzimas silvestres)
silvestres = cinetica[cinetica["tipo_enzima"] == "Wildtype"]
numericas = ["log_km", "log_kcat", "log_kcat_km", "ph", "temperatura", "año"]

pearson = silvestres[numericas].corr(method="pearson")
spearman = silvestres[numericas].corr(method="spearman")

print("Entradas silvestres:", len(silvestres))
print("\nPearson:")
print(pearson.round(2))
print("\nSpearman:")
print(spearman.round(2))
# %% Eficiencia (kcat/Km) según la clase de enzima, solo silvestres
clases_validas = ["A", "B", "C", "D"]
comparacion = silvestres[silvestres["clase_enzima"].isin(clases_validas)]

resumen_clases = (comparacion.groupby("clase_enzima")["log_kcat_km"]
                  .describe()[["count", "25%", "50%", "75%"]])
resumen_clases["mediana_real"] = 10 ** resumen_clases["50%"]
resumen_clases.round(2)

# %% ¿De dónde vienen los datos de la clase D?
clase_d = comparacion[comparacion["clase_enzima"] == "D"]
print("Artículos distintos:", clase_d["pubmed"].nunique())
print("Proteínas distintas:", clase_d["uniprot"].nunique())
print(clase_d["antibiotico"].value_counts())

# %% Clase de enzima × clase de antibiótico (silvestres, antibióticos clínicos)
clinicos_ab = ["penicilina", "cefalosporina", "carbapenemico", "monobactamico"]
cruce = comparacion[comparacion["clase_antibiotico"].isin(clinicos_ab)]

# 1. Conteos: cuántas mediciones hay de cada combinación
conteos = pd.crosstab(cruce["clase_enzima"], cruce["clase_antibiotico"])
print("Número de mediciones:")
print(conteos)

# 2. Porcentajes por fila: de lo que se midió en cada clase de enzima, qué % fue con cada antibiótico
porcentajes = pd.crosstab(cruce["clase_enzima"], cruce["clase_antibiotico"],
                          normalize="index") * 100
print("\nPorcentaje de cada clase de enzima:")
print(porcentajes.round(1))

# 3. Eficiencia mediana de cada combinación (lo que responde la pregunta 1)
eficiencia = cruce.pivot_table(index="clase_enzima", columns="clase_antibiotico",
                               values="log_kcat_km", aggfunc="median")
print("\nMediana de log10 kcat/Km:")
print(eficiencia.round(2))

# %%# %% Cuántas mediciones de kcat/Km hay en cada combinación
n_eficiencia = cruce.pivot_table(index="clase_enzima", columns="clase_antibiotico",
                                 values="log_kcat_km", aggfunc="count")
n_eficiencia
# %% Figura 1: distribución de los parámetros cinéticos
variables = [(cinetica, "log_km", r"$\log_{10}$ Km (µM)", "Afinidad (Km)"),
             (cinetica, "log_kcat", r"$\log_{10}$ kcat (s$^{-1}$)", "Velocidad (kcat)"),
             (cinetica, "log_kcat_km", r"$\log_{10}$ kcat/Km (M$^{-1}$ s$^{-1}$)", "Eficiencia (kcat/Km)"),
             (inhibicion, "log_ki", r"$\log_{10}$ Ki (µM)", "Potencia de inhibidores (Ki)")]

sns.set_style("whitegrid")
fig, ejes = plt.subplots(2, 2, figsize=(11, 8))

for eje, (tabla, columna, etiqueta, titulo) in zip(ejes.flat, variables):
    datos = tabla[columna].dropna()
    sns.histplot(datos, bins=30, kde=True, color="steelblue", ax=eje)

    q1, mediana, q3 = datos.quantile([0.25, 0.5, 0.75])
    eje.axvline(mediana, color="crimson", linewidth=2, label=f"Mediana ({mediana:.2f})")
    eje.axvline(q1, color="darkorange", linestyle="--", label=f"Q1 ({q1:.2f})")
    eje.axvline(q3, color="darkorange", linestyle="--", label=f"Q3 ({q3:.2f})")

    eje.set_title(f"{titulo}  (n = {len(datos)})")
    eje.set_xlabel(etiqueta)
    eje.set_ylabel("Número de mediciones")
    eje.legend(fontsize=8)

fig.suptitle("Distribución de los parámetros cinéticos de β-lactamasas (SABIO-RK)", fontsize=13)
fig.tight_layout()
fig.savefig("../figuras/fig1_distribuciones.png", dpi=300)
plt.show()

# %% Figura 2: eficiencia por clase de enzima (silvestres)
fig, eje = plt.subplots(figsize=(8, 5))

orden = ["A", "B", "C", "D"]
sns.boxplot(data=comparacion, x="clase_enzima", y="log_kcat_km", order=orden,
            color="lightsteelblue", showfliers=False, ax=eje)
sns.stripplot(data=comparacion, x="clase_enzima", y="log_kcat_km", order=orden,
              color="steelblue", alpha=0.35, size=3, jitter=0.25, ax=eje)

# Número de mediciones debajo de cada clase
n_por_clase = comparacion.groupby("clase_enzima")["log_kcat_km"].count()
eje.set_xticks(range(4))
eje.set_xticklabels([f"Clase {c}\n(n = {n_por_clase[c]})" for c in orden])

eje.set_xlabel("Clase molecular de β-lactamasa")
eje.set_ylabel(r"$\log_{10}$ kcat/Km (M$^{-1}$ s$^{-1}$)")
eje.set_title("Eficiencia catalítica según la clase de β-lactamasa (enzimas silvestres)")
sns.despine()
fig.tight_layout()
fig.savefig("../figuras/fig2_eficiencia_por_clase.png")
plt.show()

# %% Figura 3: relación entre Km y kcat (silvestres)
datos_rel = silvestres.dropna(subset=["log_km", "log_kcat"])
r_pearson = datos_rel["log_km"].corr(datos_rel["log_kcat"], method="pearson")
r_spearman = datos_rel["log_km"].rank().corr(datos_rel["log_kcat"].rank())

fig, eje = plt.subplots(figsize=(8, 6))
sns.regplot(data=datos_rel, x="log_km", y="log_kcat", ax=eje,
            scatter_kws={"alpha": 0.3, "s": 15, "color": "steelblue", "edgecolor": "none"},
            line_kws={"color": "crimson", "linewidth": 2})

eje.text(0.03, 0.97,
         f"n = {len(datos_rel)}\nPearson r = {r_pearson:.2f}\nSpearman ρ = {r_spearman:.2f}",
         transform=eje.transAxes, va="top", fontsize=10,
         bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "gray"})

eje.set_xlabel(r"$\log_{10}$ Km (µM)")
eje.set_ylabel(r"$\log_{10}$ kcat (s$^{-1}$)")
eje.set_title("Relación entre afinidad (Km) y velocidad (kcat) en enzimas silvestres")
sns.despine()
fig.tight_layout()
fig.savefig("../figuras/fig3_km_vs_kcat.png")
plt.show()

# %%# %% Figura 4: evolución en el tiempo de lo que se estudia
cinetica["periodo"] = (cinetica["año"] // 5) * 5
clinicos_ab = ["penicilina", "cefalosporina", "carbapenemico", "monobactamico"]
por_periodo = (cinetica[cinetica["clase_antibiotico"].isin(clinicos_ab)]
               .groupby(["periodo", "clase_antibiotico"]).size()
               .unstack(fill_value=0))

pct_mutantes = (cinetica.groupby("periodo")["tipo_enzima"]
                .apply(lambda x: (x == "Mutant").mean() * 100))
articulos = cinetica.groupby("periodo")["pubmed"].nunique()
# Incluir todos los quinquenios, aunque no tengan datos
todos = range(1980, 2020, 5)
por_periodo = por_periodo.reindex(todos, fill_value=0)
pct_mutantes = pct_mutantes.reindex(todos)          # sin datos -> vacío: la línea se corta
articulos = articulos.reindex(todos, fill_value=0)
fig, (eje1, eje2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
# Panel A: mediciones por tipo de antibiótico
for clase in clinicos_ab:
    eje1.plot(por_periodo.index, por_periodo[clase], marker="o", linewidth=2, label=clase)
eje1.set_ylabel("Número de mediciones")
eje1.set_title("A. Mediciones por tipo de antibiótico")
eje1.legend(title="Antibiótico")
# Panel B: porcentaje de enzimas mutantes
eje2.plot(pct_mutantes.index, pct_mutantes.values, marker="o", color="purple", linewidth=2)
for x, y, n in zip(pct_mutantes.index, pct_mutantes.values, articulos.values):
    if n > 0:
        eje2.annotate(f"{n} art.", (x, y), textcoords="offset points", xytext=(0, 8),
                      ha="center", fontsize=8, color="gray")
eje2.set_ylabel("Mediciones con enzimas mutantes (%)")
eje2.set_title("B. Proporción de estudios con enzimas mutantes")
eje2.set_xlabel("Periodo de publicación (inicio del quinquenio)")
eje2.set_xticks(list(todos))
fig.suptitle("Evolución de las mediciones cinéticas de β-lactamasas (1984–2019)", fontsize=13)
fig.tight_layout()
fig.savefig("../figuras/fig4_evolucion_temporal.png")
plt.show()
print(por_periodo)
print(articulos)

# %% Figura 5: heatmap de eficiencia, clase de enzima × tipo de antibiótico
from matplotlib.colors import ListedColormap

orden_ab = ["penicilina", "cefalosporina", "carbapenemico", "monobactamico"]
orden_enz = ["A", "B", "C", "D"]

mediana = (cruce.pivot_table(index="clase_enzima", columns="clase_antibiotico",
                             values="log_kcat_km", aggfunc="median")
           .reindex(index=orden_enz, columns=orden_ab))
n = (cruce.pivot_table(index="clase_enzima", columns="clase_antibiotico",
                       values="log_kcat_km", aggfunc="count")
     .reindex(index=orden_enz, columns=orden_ab).fillna(0).astype(int))

texto = mediana.round(1).astype(str) + "\n(n = " + n.astype(str) + ")"
texto = texto.where(n > 0, "sin datos")
confiable = mediana.where(n >= 5)

fig, eje = plt.subplots(figsize=(8, 5.5))

# Capa de fondo: blanco = sin datos, gris = menos de 5 mediciones (con su texto)
fondo = (n > 0).astype(float)
sns.heatmap(fondo, cmap=ListedColormap(["white", "lightgray"]), cbar=False,
            annot=texto, fmt="", annot_kws={"color": "dimgray"},
            linewidths=1, linecolor="white", ax=eje)

# Capa de color: solo las celdas con 5 o más mediciones
sns.heatmap(confiable, annot=texto, fmt="", cmap="viridis", vmin=3, vmax=7,
            linewidths=1, linecolor="white", ax=eje,
            cbar_kws={"label": r"Mediana de $\log_{10}$ kcat/Km (M$^{-1}$ s$^{-1}$)"})

eje.grid(False)
eje.set_xticklabels(["Penicilinas", "Cefalosporinas", "Carbapenémicos", "Monobactámicos"])
eje.set_yticklabels([f"Clase {c}" for c in orden_enz], rotation=0)
eje.set_xlabel("Tipo de antibiótico")
eje.set_ylabel("Clase de β-lactamasa")
eje.set_title("Eficiencia catalítica por clase de enzima y tipo de antibiótico\n"
              "(enzimas silvestres; celdas grises: menos de 5 mediciones)")
fig.tight_layout()
fig.savefig("../figuras/fig5_heatmap_eficiencia.png")
plt.show()
# %% Pregunta 2: potencia de los inhibidores según la clase de enzima
inh = inhibicion[inhibicion["clase_enzima"].isin(["A", "B", "C", "D"])]

tabla_ki = inh.pivot_table(index="tipo_inhibidor", columns="clase_enzima",
                           values="log_ki", aggfunc="median").round(2)
n_ki = inh.pivot_table(index="tipo_inhibidor", columns="clase_enzima",
                       values="log_ki", aggfunc="count").fillna(0).astype(int)
print("Mediana de log10 Ki (µM):")
print(tabla_ki)
print("\nNúmero de mediciones:")
print(n_ki)

# Inhibidores clínicos: qué clases se midieron con ellos
clinicos_df = inh[inh["tipo_inhibidor"] == "inhibidor clinico"]
print("\nInhibidores clínicos por clase de enzima:")
print(pd.crosstab(clinicos_df["inhibidores"], clinicos_df["clase_enzima"]))

# %%
