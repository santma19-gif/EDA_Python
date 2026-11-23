# %% Importar librerias
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# %% Cargar datos (toma la descarga más reciente)
from pathlib import Path

archivos = sorted(Path("../datos/crudos").glob("sabio_betalactamasas_*.csv"))
print("Usando:", archivos[-1].name)

df = pd.read_csv(archivos[-1],
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

# %% Paso B: unificar sinónimos y errores de escritura
sinonimos = {"Benzylpenicillin": "Penicillin G",
             "Phenoxymethylpenicillin": "Penicillin V",
             "Cephaloridin": "Cephaloridine",
             "Cefaloridine": "Cephaloridine",
             "Ticaracillin": "Ticarcillin",
             "Amoxycillin": "Amoxicillin",
             "Cephazoline": "Cefazolin",
             "6-Aminopenicillanate": "6-Aminopenicillanic acid",
             }

genericos = ["Cephalosporin", "beta-Lactam"]   # no dicen qué compuesto se midió

antes = cinetica["antibiotico"].nunique()
cinetica["antibiotico"] = cinetica["antibiotico"].replace(sinonimos)
cinetica.loc[cinetica["antibiotico"].isin(genericos), "antibiotico"] = np.nan

print(f"Nombres distintos: {antes} -> {cinetica['antibiotico'].nunique()}")
cinetica["antibiotico"].value_counts().head(10)

# %% Paso C: clase de antibiótico
clases = {
    "penicilina": ["Penicillin G", "Penicillin V", "Ampicillin", "Amoxicillin",
                   "Piperacillin", "Carbenicillin", "Ticarcillin", "Oxacillin",
                   "Cloxacillin", "Azlocillin", "Mezlocillin", "Methicillin",
                   "Mecillinam", "Sulbenicillin", "Epicillin", "Temocillin",
                   "Talampicillin"],
    "cefalosporina": ["Cephalothin", "Cefotaxime", "Cephaloridine", "Cefuroxime",
                      "Ceftazidime", "Cefoxitin", "Cefazolin", "Cephalexin",
                      "Cefaclor", "Cefepime", "Cefdinir", "Cefamandole",
                      "Cefpodoxime", "Cephapirin", "Ceftriaxone", "Cefsulodin",
                      "Cefpirome", "Cefoperazone", "Cephradine", "Cephaloglycin",
                      "Ceftizoxime", "Cefteram", "Cefixime", "Cefcapene",
                      "Cefotiam", "Cefotetan", "Cefprozil", "Cefmetazole",
                      "Cefadroxil", "Desacetyl cefotaxime",
                      "Moxalactam",      # oxacefem: se agrupa con las cefalosporinas
                      "Loracarbef"],     # carbacefem: se agrupa con las cefalosporinas
    "carbapenemico": ["Imipenem", "Meropenem", "Biapenem", "Ertapenem",
                      "Doripenem", "Panipenem", "Asparenomycin A"],
    "monobactamico": ["Aztreonam", "Carumonam"],
    "sustrato de laboratorio": ["Nitrocefin", "CENTA", "Bocillin FL",
                                "Pyridine-2-azo-4-dimethylaniline cephalosporin"],
    "inhibidor": ["Sulbactam", "Tazobactam", "Clavulanic acid", "BRL42715"],
}

# Invertir el diccionario: de {clase: [antibióticos]} a {antibiótico: clase}
mapa_clases = {antibiotico: clase
               for clase, lista in clases.items()
               for antibiotico in lista}

cinetica["clase_antibiotico"] = cinetica["antibiotico"].map(mapa_clases)

# Lo que no está en ninguna lista (y no está vacío) pasa a "otros"
sin_clase = cinetica["antibiotico"].notna() & cinetica["clase_antibiotico"].isna()
print("Compuestos clasificados como 'otros':")
print(sorted(cinetica.loc[sin_clase, "antibiotico"].unique()))
cinetica.loc[sin_clase, "clase_antibiotico"] = "otros"

cinetica["clase_antibiotico"].value_counts(dropna=False)

# %% Paso D: lista de inhibidores
print("Inhibidores distintos:", inhibicion["inhibidores"].nunique())
print("Entradas sin inhibidor:", inhibicion["inhibidores"].isna().sum())
print("Entradas con más de un inhibidor:",
      inhibicion["inhibidores"].str.contains(r"\|", na=False).sum())

print(inhibicion["inhibidores"].value_counts().to_string())

# %% Paso D: tipo de inhibidor
# 1. Recuperar el inhibidor desde especie_parametro (formato "1 | Nombre | Rol")
especie_ki = (df[df["parametro"] == "Ki"]
              .drop_duplicates("entrada")
              .set_index("entrada")["especie_parametro"])

sin_inhibidor = inhibicion["inhibidores"].isna()
print("Sin inhibidor, antes:", sin_inhibidor.sum())
print(especie_ki[inhibicion.loc[sin_inhibidor, "entrada"]].value_counts().to_string())

# %% Paso D: tipo de inhibidor
# 1. Excluir los Ki que no son de un inhibidor (del sustrato o del producto)
filas_antes = len(inhibicion)
inhibicion = inhibicion[inhibicion["inhibidores"].notna()]
print(f"Entradas: {filas_antes} -> {len(inhibicion)}")

# 2. Sinónimos
sinonimos_inhibidores = {
    "Clavulanate": "Clavulanic acid",
    "Benzylpenicillin": "Penicillin G",
    "Latamoxef": "Moxalactam",
    "Sodium benzyl 2-benzoxazolone-5-acyl phosphonate":
        "Sodium benzyl 2-benzoxazolone-5-acylphosphonate",
}
antes = inhibicion["inhibidores"].nunique()
inhibicion["inhibidores"] = inhibicion["inhibidores"].replace(sinonimos_inhibidores)
print(f"Inhibidores distintos: {antes} -> {inhibicion['inhibidores'].nunique()}")

# 3. Clasificar
clinicos = ["Clavulanic acid", "Sulbactam", "Tazobactam"]
clases_clinicas = ["penicilina", "cefalosporina", "carbapenemico", "monobactamico"]
betalactamicos = [ab for ab, clase in mapa_clases.items() if clase in clases_clinicas]
betalactamicos.append("Cefmenoxime")   # cefalosporina que no estaba en el paso C

def tipo_inhibidor(nombre):
    if nombre in clinicos:
        return "inhibidor clinico"
    if nombre in betalactamicos:
        return "antibiotico betalactamico"
    return "experimental"

inhibicion["tipo_inhibidor"] = inhibicion["inhibidores"].apply(tipo_inhibidor)
inhibicion["tipo_inhibidor"].value_counts()

# %% Paso E: tipo de enzima (metalo vs serina), a partir del cofactor y el UniProt
# 1. Proteínas (UniProt) que tienen Zn en al menos una entrada, en cualquiera de las dos tablas
todas = pd.concat([cinetica[["uniprot", "cofactor"]],
                   inhibicion[["uniprot", "cofactor"]]])
con_zinc = todas["cofactor"].str.contains("Zn", na=False)
uniprot_metalo = set(todas.loc[con_zinc, "uniprot"].dropna())
print("Proteínas con zinc:", sorted(uniprot_metalo))

# 2. Clasificar cada entrada
def tipo_betalactamasa(fila):
    if "Zn" in str(fila["cofactor"]):
        return "metalo (clase B)"
    if fila["uniprot"] in uniprot_metalo:
        return "metalo (clase B)"
    if pd.notna(fila["uniprot"]):
        return "serina (clases A, C, D)"
    return "sin clasificar"

cinetica["tipo_betalactamasa"] = cinetica.apply(tipo_betalactamasa, axis=1)
inhibicion["tipo_betalactamasa"] = inhibicion.apply(tipo_betalactamasa, axis=1)

print("\nCinética:")
print(cinetica["tipo_betalactamasa"].value_counts())
print("\nInhibición:")
print(inhibicion["tipo_betalactamasa"].value_counts())

# %% Exportar los UniProt distintos para consultarlos en UniProt
ids = (pd.concat([cinetica["uniprot"], inhibicion["uniprot"]])
       .dropna()
       .str.split("|")
       .explode()
       .unique())

print("Proteínas distintas:", len(ids))
pd.Series(ids).to_csv("../datos/crudos/uniprot_ids.txt", index=False, header=False)
# %%
# %% Paso E: clase de enzima (A, B, C, D) según la familia en UniProt
familias = pd.read_csv("../datos/crudos/idmapping_2026_11_22.tsv.gz", sep="\t")

def clase_desde_familia(familia):
    if pd.isna(familia):
        return np.nan
    for letra in ["A", "B", "C", "D"]:
        if f"Class-{letra} beta-lactamase family" in familia:
            return letra
    return "otra enzima"

familias["clase"] = familias["Protein families"].apply(clase_desde_familia)
mapa_enzimas = dict(zip(familias["From"], familias["clase"]))

# Correcciones a mano (ver justificación en la bitácora)
mapa_enzimas.update({"Q9K578": "B",            # FEZ-1: metalo-β-lactamasa sin familia anotada
                     "Q79F59": "B",            # eliminada de UniProt; tenía Zn en SABIO-RK
                     "P29723": "otra enzima"}) # Tp47: carboxipeptidasa con actividad β-lactamasa

def asignar_clase(tabla):
    tabla = tabla.copy()
    primer_id = tabla["uniprot"].str.split("|").str[0]
    tabla["clase_enzima"] = primer_id.map(mapa_enzimas).fillna("sin clasificar")
    return tabla.drop(columns="tipo_betalactamasa", errors="ignore")

cinetica = asignar_clase(cinetica)
inhibicion = asignar_clase(inhibicion)

print("Cinética:")
print(cinetica["clase_enzima"].value_counts())
print("\nInhibición:")
print(inhibicion["clase_enzima"].value_counts())
# %%
# %% Paso F: cierre
# 1. Excluir "otra enzima": no son β-lactamasas de las clases A-D (incluye las humanas)
for nombre, tabla in [("cinetica", cinetica), ("inhibicion", inhibicion)]:
    print(f"{nombre}: {len(tabla)} -> {(tabla['clase_enzima'] != 'otra enzima').sum()}")

cinetica = cinetica[cinetica["clase_enzima"] != "otra enzima"].copy()
inhibicion = inhibicion[inhibicion["clase_enzima"] != "otra enzima"].copy()

# 2. Valores imposibles: negativos y rangos de pH y temperatura
print("\nValores negativos:")
print((cinetica[["km", "kcat", "kcat_km"]] < 0).sum())
print("Ki negativos:", (inhibicion["ki"] < 0).sum())
print("\nRangos de condiciones:")
print(cinetica[["ph", "temperatura"]].describe().loc[["min", "max"]])
# Ceros en Ki: no son mediciones válidas
print("Ceros en Ki:", (inhibicion["ki"] == 0).sum())
inhibicion.loc[inhibicion["ki"] == 0, "ki"] = np.nan

# 3. Guardar las tablas limpias
cinetica.to_csv("../datos/procesados/cinetica.csv", index=False)
inhibicion.to_csv("../datos/procesados/inhibicion.csv", index=False)

print("\nTamaño final:")
print("Cinética:", cinetica.shape)
print("Inhibición:", inhibicion.shape)
# %%
