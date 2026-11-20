"""
Descarga los datos de cinética de β-lactamasas (EC 3.5.2.6) desde SABIO-RK.

Fuente: SABIO-RK, Export API (https://sabio.h-its.org/ui/export-api)
Licencia de los datos: uso no comercial (términos de SABIO-RK). Por eso los
archivos descargados NO se suben al repositorio; este script los regenera.

Qué hace:
1. Pide a la API todas las entradas con ECNumber:3.5.2.6, en páginas de 1000.
2. Guarda la respuesta tal cual en datos/crudos/ como JSON (copia exacta).
3. Aplana el JSON en un CSV en formato largo (una fila por parámetro),
   con las mismas columnas que la exportación manual de la página.

El CSV NO está limpio: no se corrigen unidades, valores ni nombres.
Toda la limpieza se hace después, en otro script.

Uso (desde la carpeta principal del repositorio):
    python scripts/descargar_datos.py
"""

import json
import time
from datetime import date
from pathlib import Path

import pandas as pd
import requests

# ---------------------------------------------------------------------------
# Configuración
# ---------------------------------------------------------------------------
URL_API = "https://sabio.h-its.org/export-api/sabio/kinlaw-entry/json"
CONSULTA = "ECNumber:3.5.2.6"   # todas las entradas de β-lactamasa
TAMANO_PAGINA = 1000            # máximo permitido por la API para JSON
PAUSA_SEGUNDOS = 2              # la API permite 60 peticiones por minuto

CARPETA_CRUDOS = Path(__file__).resolve().parent.parent / "datos" / "crudos"
HOY = date.today().isoformat()
ARCHIVO_JSON = CARPETA_CRUDOS / f"sabio_betalactamasas_{HOY}.json"
ARCHIVO_CSV = CARPETA_CRUDOS / f"sabio_betalactamasas_{HOY}.csv"


# ---------------------------------------------------------------------------
# 1. Descarga
# ---------------------------------------------------------------------------
def pedir_pagina(pagina):
    """Pide una página a la API. Si se pasa del límite (error 429), espera y reintenta."""
    parametros = {"q": CONSULTA, "page": pagina, "pageSize": TAMANO_PAGINA}
    for intento in range(3):
        respuesta = requests.get(URL_API, params=parametros, timeout=120)
        if respuesta.status_code == 429:
            print("  Límite de peticiones alcanzado; esperando 60 s...")
            time.sleep(60)
            continue
        respuesta.raise_for_status()
        return respuesta.json()
    raise RuntimeError(f"No se pudo descargar la página {pagina} después de 3 intentos.")


def descargar_todo():
    """Recorre todas las páginas y devuelve la lista completa de entradas."""
    primera = pedir_pagina(1)
    total = primera["meta"]["total_count"]
    paginas = primera["meta"]["total_pages"]
    print(f"La consulta '{CONSULTA}' tiene {total} entradas en {paginas} página(s).")

    entradas = list(primera["data"])
    for pagina in range(2, paginas + 1):
        time.sleep(PAUSA_SEGUNDOS)
        print(f"  Descargando página {pagina} de {paginas}...")
        entradas.extend(pedir_pagina(pagina)["data"])

    # Comprobación: deben haber llegado todas las entradas, sin repetidas
    ids = [e["id"] for e in entradas]
    if len(ids) != total or len(set(ids)) != total:
        print(f"ADVERTENCIA: se esperaban {total} entradas únicas y llegaron "
              f"{len(ids)} ({len(set(ids))} únicas).")
    return entradas


# ---------------------------------------------------------------------------
# 2. Aplanar el JSON a una tabla en formato largo
# ---------------------------------------------------------------------------
def obtener(diccionario, *claves):
    """Lee un valor anidado sin fallar si falta algún nivel (devuelve None)."""
    for clave in claves:
        if not isinstance(diccionario, dict):
            return None
        diccionario = diccionario.get(clave)
    return diccionario


def unir(valores):
    """Une una lista de textos con '|', como en la exportación de la página."""
    valores = [str(v) for v in valores if v not in (None, "")]
    return "|".join(valores) if valores else None


def compuestos_por_rol(especies, rol):
    return unir(obtener(e, "compound", "name") for e in especies if e.get("role") == rol)


def aplanar(entradas):
    """Convierte cada entrada en una fila por parámetro (formato largo)."""
    filas = []
    for e in entradas:
        especies = obtener(e, "reaction", "species") or []
        proteinas = obtener(e, "enzyme_description", "proteins") or []

        # Información común a todos los parámetros de la entrada
        comun = {
            "EntryID": e.get("id"),
            "Organism": obtener(e, "general", "organism", "name"),
            "NCBITaxonomyID": obtener(e, "general", "organism", "ncbi_taxonomy_id"),
            "Tissue": obtener(e, "general", "tissue", "name"),
            "Strain": obtener(e, "general", "strain", "name"),
            "Comment": obtener(e, "general", "comment"),
            "ExperimentType": obtener(e, "general", "experiment_type"),
            "InsertDate": obtener(e, "general", "date_inserted"),
            "Reaction": obtener(e, "reaction", "equation"),
            "IsTransport": obtener(e, "reaction", "is_transport"),
            "Substrates": compuestos_por_rol(especies, "Substrate"),
            "Products": compuestos_por_rol(especies, "Product"),
            "Inhibitors": compuestos_por_rol(especies, "Inhibitor"),
            "Activators": compuestos_por_rol(especies, "Activator"),
            "Cofactors": compuestos_por_rol(especies, "Cofactor"),
            "OtherModifiers": compuestos_por_rol(especies, "OtherModifier"),
            "Location": unir(sorted({obtener(s, "location", "name") for s in especies} - {None})),
            "Pathway": unir(obtener(p, "name") for p in (obtener(e, "reaction", "pathways") or [])),
            "KinlawType": obtener(e, "kineticlaw", "kinlaw_type", "name"),
            "KinlawFormula": obtener(e, "kineticlaw", "formula"),
            "ECNumber": obtener(e, "enzyme_description", "ec_number"),
            "EnzymeName": obtener(e, "enzyme_description", "enzyme_name"),
            "EnzymeType": obtener(e, "enzyme_description", "wildtype"),
            "Specification": obtener(e, "enzyme_description", "mutant_spec"),
            "IsRecombinant": obtener(e, "enzyme_description", "is_recombinant"),
            "ExpressedIn": obtener(e, "enzyme_description", "expressed_in"),
            "UniprotIDs": unir(p.get("uniprot_id") for p in proteinas),
            "pH": obtener(e, "experimental_conditions", "envvar_ph", "start_value"),
            "Temperature": obtener(e, "experimental_conditions", "envvar_temperature", "start_value"),
            "TemperatureUnit": obtener(e, "experimental_conditions", "envvar_temperature", "unit"),
            "Buffer": obtener(e, "experimental_conditions", "buffer"),
            "PubMedID": obtener(e, "publication", "pubmed_id"),
            "PublicationTitle": obtener(e, "publication", "title"),
            "Author": obtener(e, "publication", "author"),
            "Journal": obtener(e, "publication", "journal"),
            "Year": obtener(e, "publication", "year"),
        }

        parametros = obtener(e, "kineticlaw", "parameter") or []
        if not parametros:
            # Entrada sin parámetros: se conserva igual, con los campos del parámetro vacíos
            filas.append(comun)
            continue

        for p in parametros:
            fila = dict(comun)
            fila.update({
                "ParameterName": p.get("name"),
                "ParameterType": obtener(p, "parameter_type", "name"),
                "ParameterSpecies": obtener(p, "species", "species_key"),
                "ParameterStartValue": p.get("start_value"),
                "ParameterEndValue": p.get("end_value"),
                "ParameterStandDev": p.get("standard_deviation"),
                "ParameterUnit": obtener(p, "unit", "name"),
            })
            filas.append(fila)

    return pd.DataFrame(filas)


# ---------------------------------------------------------------------------
# Programa principal
# ---------------------------------------------------------------------------
def main():
    CARPETA_CRUDOS.mkdir(parents=True, exist_ok=True)

    entradas = descargar_todo()

    with open(ARCHIVO_JSON, "w", encoding="utf-8") as f:
        json.dump(entradas, f, ensure_ascii=False)
    print(f"JSON crudo guardado en: {ARCHIVO_JSON}")

    tabla = aplanar(entradas)
    tabla.to_csv(ARCHIVO_CSV, index=False, encoding="utf-8")
    print(f"CSV guardado en: {ARCHIVO_CSV}")
    print(f"  {tabla['EntryID'].nunique()} entradas -> {len(tabla)} filas "
          f"(una por parámetro), {tabla.shape[1]} columnas.")


if __name__ == "__main__":
    main()
