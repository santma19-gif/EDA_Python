# %% Importar librerias
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# %% Cargar datos
df = pd.read_csv('datos/crudos/sabio_betalactamasas_2026-11-19.csv', low_memory=False, )

print(df.shape)
df.head()