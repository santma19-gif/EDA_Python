# EDA_Python: cinética de β-lactamasas

Análisis exploratorio de datos (EDA) con Python sobre parámetros cinéticos de **β-lactamasas**, las enzimas con las que las bacterias destruyen la penicilina y sus derivados. El trabajo va desde la descarga y limpieza de los datos hasta la comunicación de los hallazgos con gráficos.

> Proyecto en curso. Este README se irá completando a medida que avance el análisis.

## De qué trata

Las β-lactamasas son una de las principales causas de resistencia a antibióticos. Cada medición de cinética enzimática dice qué tan rápido (kcat) y con qué afinidad (Km) una de estas enzimas destruye un antibiótico, o qué tan fuerte la bloquea un inhibidor (Ki). Reuniendo miles de esas mediciones, publicadas en artículos científicos a lo largo de varias décadas, se puede explorar la "carrera" entre los antibióticos y la resistencia.

## Preguntas

1. **¿Qué tan eficiente es cada clase de β-lactamasa contra cada tipo de antibiótico?** Clase molecular de la enzima (A, B, C, D) × tipo de antibiótico (penicilinas, cefalosporinas, carbapenémicos, monobactámicos), medido con kcat/Km.
2. **¿Los inhibidores funcionan igual contra todas las clases?** Comparación del Ki de los inhibidores (ácido clavulánico, tazobactam, avibactam...) entre clases de enzimas.
3. **¿Cómo ha cambiado en el tiempo lo que se estudia?** Evolución, según el año de publicación, de los antibióticos medidos y del estudio de enzimas mutantes.

## Datos

| | |
|---|---|
| **Fuente** | [SABIO-RK](https://sabio.h-its.org), base de datos de cinética de reacciones bioquímicas curada a partir de la literatura (HITS gGmbH) |
| **Acceso** | Export API de SABIO-RK (`https://sabio.h-its.org/export-api/sabio/`) |
| **Consulta** | `ECNumber:3.5.2.6` (β-lactamasa) |
| **Fecha de descarga** | 2026-11-19 |
| **Tamaño crudo** | 2,341 entradas → _[N]_ filas (una por parámetro), 43 columnas |
| **Licencia** | Uso no comercial ([términos de SABIO-RK](https://sabio.h-its.org)) |

Por la licencia de SABIO-RK, **los datos no se incluyen en este repositorio**. El script `scripts/descargar_datos.py` los descarga directamente de la fuente (ver *Cómo reproducir*).

Cada entrada de SABIO-RK es un experimento: una enzima (silvestre o mutante) con un sustrato o inhibidor, en unas condiciones de pH y temperatura, con los parámetros cinéticos reportados en un artículo.

## Estructura del repositorio

```
EDA_Python/
├── README.md
├── datos/
│   ├── crudos/       ← datos tal como se descargan (no se suben a GitHub)
│   └── procesados/   ← datos después de la limpieza
├── scripts/
│   └── descargar_datos.py
└── figuras/          ← gráficos generados por el análisis
```

## Cómo reproducir

Requisitos: Python 3 con `pandas` y `requests`.

```bash
git clone https://github.com/santma19-gif/EDA_Python.git
cd EDA_Python
pip install pandas requests
python scripts/descargar_datos.py
```

El script consulta la API de SABIO-RK, guarda la respuesta original en `datos/crudos/` (JSON) y la convierte a un CSV en formato largo, sin ninguna limpieza. Como la base se actualiza, una descarga en otra fecha puede traer más entradas que las reportadas aquí.

## Herramientas

- **pandas** para el manejo y la limpieza de los datos.
- **requests** para la descarga desde la API.
- Biblioteca de gráficos: _por definir_.

## Uso de inteligencia artificial

En este trabajo se usó Claude (Anthropic) como apoyo para elegir el conjunto de datos, plantear las preguntas, escribir el script de descarga y revisar la redacción. El detalle se documenta en la bitácora del trabajo.

## Licencia

El código de este repositorio está bajo la licencia MIT. Los datos pertenecen a SABIO-RK y se rigen por sus propios términos de uso.