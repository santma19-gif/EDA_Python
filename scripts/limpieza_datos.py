# %% Importar librerias
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# %% Cargar datos
df = pd.read_csv('../datos/crudos/sabio_betalactamasas_2026-11-19.csv',
                 low_memory=False,
                 )
print(df.shape)
df.head()
# %%
diagnostico = pd.DataFrame({
    "tipo": df.dtypes,
    "faltantes": df.isnull().sum(),
    "porcentaje_faltantes": (df.isna().mean()*100).round(1),
    "distintos": df.nunique(),
})

diagnostico = diagnostico.sort_values("porcentaje_faltantes", ascending=False)
diagnostico


# %% Selección de columnas (43 a 19)
seleccion = ["EntryID",
             "PubMedID",
             "Organism",
             "Comment",
             "UniprotIDs",
             "Cofactors",
             "Substrates",
             "Inhibitors",
             "KinlawType",
             "EnzymeType",
             "Specification",
             "pH",
             "Temperature",
             "Year",
             "ParameterType",
             "ParameterSpecies",
             "ParameterStartValue",
             "ParameterStandDev",
             "ParameterUnit"]

print("Antes de la selección:", df.shape)
df = df[seleccion]
print("Después de la selección:", df.shape)

# %% Normalizar nombres de columnas
df = df.rename(columns={"EntryID": "entrada",
                        "PubMedID": "pubmed",
                        "Organism": "organismo",
                        "Comment": "comentario",
                        "UniprotIDs": "uniprot",
                        "Cofactors": "cofactor",
                        "Substrates": "sustratos",
                        "Inhibitors": "inhibidores",
                        "KinlawType": "tipo_cinetica",
                        "EnzymeType": "tipo_enzima",
                        "Specification": "mutacion",
                        "pH": "ph",
                        "Temperature": "temperatura",
                        "Year": "año",
                        "ParameterType": "parametro",
                        "ParameterSpecies": "especie_parametro",
                        "ParameterStartValue": "valor",
                        "ParameterStandDev": "desviacion",
                        "ParameterUnit": "unidad",
                        })

df.columns.tolist()

# %% Tipos de parametros
df["parametro"].value_counts()
# %% Parametros seleccionados
parametros_seleccionados = ["Km",
                            "kcat", 
                            "kcat/Km", 
                            "Ki",
                            ]

filas_antes = len(df)
entradas_antes = df["entrada"].nunique()

seleccion = df["parametro"].isin(parametros_seleccionados)
df = df[seleccion]

print(f"Filas: {filas_antes} -> {len(df)}")
print(f"Entradas: {entradas_antes} -> {df['entrada'].nunique()}")

# %% Unidades por tipo de parámetro
df.groupby("parametro")["unidad"].value_counts()

# %% Revisión filas sin unidad ("-") y sin valor
sin_unidad = df["unidad"] == "-"
sin_valor = df["valor"].isna()

print("Filas sin unidad:", sin_unidad.sum())
print("Filas sin valor:", sin_valor.sum())
print("Sin unidad Y sin valor:", (sin_unidad & sin_valor).sum())
print("Sin unidad pero CON valor:", (sin_unidad & ~sin_valor).sum())
# %% Filas de kcat/Km con unidad µM/s
df[df["unidad"] == "µM/s"][["entrada", "pubmed", "sustratos", "parametro", "valor", "unidad"]]

# %%
# %% ¿Coincide kcat/Km calculado con los valores en "µM/s"?
entradas_sospechosas = df[df["unidad"] == "µM/s"]["entrada"]


revision = (df[df["entrada"].isin(entradas_sospechosas)]
            .pivot_table(index=["entrada", "sustratos"],
                         columns="parametro",
                         values="valor"))

revision["kcat/Km calculado (µM-1 s-1)"] = revision["kcat"] / revision["Km"]
revision
# %% Corregir la unidad mal escrita (PubMed 10390216)
df.loc[df["unidad"] == "µM/s", "unidad"] = "µM^(-1)*s^(-1)"

# %% Excluir filas sin unidad ("-": 364 sin valor y 2 con valor)
filas_antes = len(df)
df = df[df["unidad"] != "-"]
print(f"Filas: {filas_antes} -> {len(df)}")

# %%
# %% Unificar unidades
# Factor por el que se multiplica cada valor para llevarlo a la unidad final
factores = {"µM": 1,                       # Km y Ki ya en µM
            "mM": 1_000,                   # mM -> µM
            "nM": 0.001,                   # nM -> µM
            "M": 1_000_000,                # M -> µM
            "s^(-1)": 1,                   # kcat ya en s^-1
            "min^(-1)": 1 / 60,            # min^-1 -> s^-1
            "M^(-1)*s^(-1)": 1,            # kcat/Km ya en M^-1 s^-1
            "µM^(-1)*s^(-1)": 1_000_000,   # µM^-1 s^-1 -> M^-1 s^-1
            "mM^(-1)*s^(-1)": 1_000,       # mM^-1 s^-1 -> M^-1 s^-1
            "M^(-1)*min^(-1)": 1 / 60,     # M^-1 min^-1 -> M^-1 s^-1
            "µM^(-1)*min^(-1)": 1_000_000 / 60,
            }

unidad_final = {"Km": "µM",
                "Ki": "µM",
                "kcat": "s^(-1)",
                "kcat/Km": "M^(-1)*s^(-1)",
                }

df["factor"] = df["unidad"].map(factores)
print("Filas sin factor (unidad no reconocida):", df["factor"].isna().sum())

df["valor"] = df["valor"] * df["factor"]
df["desviacion"] = df["desviacion"] * df["factor"]
df["unidad"] = df["parametro"].map(unidad_final)
df = df.drop(columns="factor")

df.groupby("parametro")["unidad"].value_counts()

# %% ¿Cuántos valores de cada parámetro tiene cada entrada?
conteo = (df.groupby(["entrada", "parametro"])
          .size()
          .unstack(fill_value=0))

print("Entradas con más de un valor del mismo tipo:")
print((conteo > 1).sum())

# %% Entradas con más de un Ki
entradas_varios_ki = conteo[conteo["Ki"] > 1].index

(df[df["entrada"].isin(entradas_varios_ki) & (df["parametro"] == "Ki")]
 [["entrada", "pubmed", "inhibidores", "especie_parametro", "valor", "unidad"]])
# %% ¿Qué tipo de cinética tienen esas entradas?
df[df["entrada"].isin(entradas_varios_ki)][["entrada", "tipo_cinetica"]].drop_duplicates()
# %%
# %% Excluir entradas con dos Ki (inhibición mixta: no se sabe cuál es cuál)
filas_antes = len(df)
df = df[~(df["entrada"].isin(entradas_varios_ki) & (df["parametro"] == "Ki"))]
print(f"Filas: {filas_antes} -> {len(df)}")
# %%
# %% Formato ancho: tabla de cinética (una fila por entrada)
cinetica_largo = df[df["parametro"].isin(["Km", "kcat", "kcat/Km"])]

# 1. Datos descriptivos: una fila por entrada
columnas_entrada = ["entrada", "pubmed", "organismo", "comentario", "uniprot",
                    "cofactor", "sustratos", "tipo_cinetica", "tipo_enzima",
                    "mutacion", "ph", "temperatura", "año"]

atributos = cinetica_largo[columnas_entrada].drop_duplicates("entrada")

# 2. Valores: Km, kcat y kcat/Km como columnas
valores = cinetica_largo.pivot(index="entrada",
                               columns="parametro",
                               values="valor")

valores = valores.rename(columns={"Km": "km",
                                  "kcat": "kcat",
                                  "kcat/Km": "kcat_km"})

# 3. Unir las dos partes
cinetica = atributos.merge(valores, on="entrada")

print("Filas en formato largo:", len(cinetica_largo))
print("Entradas en formato ancho:", cinetica.shape)
cinetica.head()
# %% Formato ancho: tabla de inhibición (una fila por entrada)
columnas_inhibicion = ["entrada", "pubmed", "organismo", "comentario", "uniprot",
                       "cofactor", "sustratos", "inhibidores", "tipo_cinetica",
                       "tipo_enzima", "mutacion", "ph", "temperatura", "año",
                       "valor"]

inhibicion = (df[df["parametro"] == "Ki"][columnas_inhibicion]
              .rename(columns={"valor": "ki"}))

print("Entradas con Ki:", inhibicion.shape)
print("¿Alguna entrada repetida?", inhibicion["entrada"].duplicated().any())
inhibicion.head()
# %% Consistencia de kcat/Km: reportado vs calculado
# km está en µM; para que kcat/km quede en M^-1 s^-1, se pasa km a M (× 1e-6)
cinetica["kcat_km_calculado"] = cinetica["kcat"] / (cinetica["km"] * 1e-6)

# Cociente reportado / calculado, en escala logarítmica:
#   0 = coinciden, 1 = el reportado es 10 veces mayor, -10 = 10^10 veces menor
cinetica["desfase_log10"] = np.log10(cinetica["kcat_km"] / cinetica["kcat_km_calculado"])

con_los_tres = cinetica["desfase_log10"].notna()
print("Entradas con Km, kcat y kcat/Km:", con_los_tres.sum())

cinetica["desfase_log10"].round().value_counts().sort_index()
# %%
# %% ¿De qué artículos vienen los errores grandes?
errores = cinetica[cinetica["desfase_log10"].abs() >= 4]
print("Entradas con desfase >= 10^4:", len(errores))
print("Artículos distintos:", errores["pubmed"].nunique())
errores.groupby(errores["desfase_log10"].round())["pubmed"].unique()

# %% ¿Cuántos ceros hay?
(cinetica[["km", "kcat", "kcat_km"]] == 0).sum()
# %% Corregir kcat/Km
# 1. Los ceros no son mediciones válidas: pasan a vacío
for columna in ["km", "kcat", "kcat_km"]:
    ceros = cinetica[columna] == 0
    print(f"Ceros en {columna}: {ceros.sum()}")
    cinetica.loc[ceros, columna] = np.nan

# 2. Recalcular con los ceros ya quitados
cinetica["kcat_km_calculado"] = cinetica["kcat"] / (cinetica["km"] * 1e-6)

# 3. Si hay Km y kcat, se usa el calculado; si no, el reportado
hay_calculado = cinetica["kcat_km_calculado"].notna()

cinetica["origen_kcat_km"] = np.where(hay_calculado, "calculado", "reportado")
cinetica["kcat_km"] = np.where(hay_calculado,
                               cinetica["kcat_km_calculado"],
                               cinetica["kcat_km"])

print(cinetica["origen_kcat_km"].value_counts())

# 4. Quitar las columnas auxiliares
cinetica = cinetica.drop(columns=["kcat_km_calculado", "desfase_log10"])
cinetica.head()
# %%
# %% Etiquetar las entradas sin kcat/Km
cinetica.loc[cinetica["kcat_km"].isna(), "origen_kcat_km"] = "sin dato"
cinetica["origen_kcat_km"].value_counts()
# %% Confirmar los problemas de comentario y uniprot con pandas
# 1. Comentarios con "-"
print("Comentarios '-' en cinética:", (cinetica["comentario"] == "-").sum())
print("Comentarios '-' en inhibición:", (inhibicion["comentario"] == "-").sum())

# 2. UniProt con el ID repetido (ej. 'P29723|P29723')
partes = cinetica["uniprot"].str.split("|")
repetidos = partes.apply(lambda ids: isinstance(ids, list) and len(ids) != len(set(ids)))
print("\nUniProt con ID repetido:", repetidos.sum())

# 3. UniProt con notación de oligómero (ej. '(P52700)*4')
oligomero = cinetica["uniprot"].str.contains(r"\(", na=False)
print("UniProt con notación de oligómero:", oligomero.sum())

# Ejemplos de cada caso
print("\nEjemplos:")
print(cinetica.loc[repetidos, "uniprot"].value_counts().head(3))
print(cinetica.loc[oligomero, "uniprot"].value_counts().head(3))

# %%
# %% Paso A: detalles menores (comentario, uniprot, pubmed)
import re

def limpiar_uniprot(texto):
    """'(P52700)*4|P52700' -> 'P52700'   |   'P29723|P29723' -> 'P29723'"""
    if pd.isna(texto):
        return texto
    ids = []
    for parte in texto.split("|"):
        parte = re.sub(r"[()]|\*\d+", "", parte).strip()   # quita paréntesis y '*4'
        if parte and parte not in ids:
            ids.append(parte)
    return "|".join(ids)


def limpiar_detalles(tabla):
    tabla = tabla.copy()
    tabla["comentario"] = tabla["comentario"].replace("-", np.nan)
    tabla["uniprot"] = tabla["uniprot"].apply(limpiar_uniprot)
    tabla["pubmed"] = tabla["pubmed"].astype(str)
    return tabla


print("UniProt antes:", cinetica["uniprot"].nunique())

cinetica = limpiar_detalles(cinetica)
inhibicion = limpiar_detalles(inhibicion)

print("UniProt después:", cinetica["uniprot"].nunique())
print("Tipo de pubmed:", cinetica["pubmed"].dtype)
print("Comentarios '-' que quedan:", (cinetica["comentario"] == "-").sum())

cinetica["uniprot"].value_counts().head(10)

# %% Paso B: extraer el antibiótico (quitar el H2O de sustratos)
def quitar_agua(texto):
    if pd.isna(texto):
        return texto
    partes = [p.strip() for p in texto.split("|") if p.strip() != "H2O"]
    return "|".join(partes)

cinetica["antibiotico"] = cinetica["sustratos"].apply(quitar_agua)

print("Antibióticos distintos:", cinetica["antibiotico"].nunique())
print("Entradas con más de un sustrato (además del agua):",
      cinetica["antibiotico"].str.contains(r"\|", na=False).sum())

print(cinetica["antibiotico"].value_counts().to_string())

# %%
