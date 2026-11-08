#!/usr/bin/env python
"""El suelo: UNA funcion de prediccion entre las que acepta el arnes
(suelo_arnes.build_fold_adatas). Para cada etiqueta pedida (perturbaciones
del panel comun + el control "NT"), devuelve 400 celulas muestreadas CON
REEMPLAZO del perfil de control de la linea, sin modificar ningun conteo.
Nunca recibe el nombre ni la identidad de la linea -- solo su matriz de
control (control_pool_X).
"""
import numpy as np

N_PER_LABEL = 400


def predict_floor(control_pool_X, labels, rng):
    """control_pool_X: csr_matrix (n_control, n_genes). labels: lista de
    etiquetas a predecir. rng: numpy Generator con semilla ya fijada por el
    arnes. Devuelve dict[label] -> csr_matrix (400, n_genes), copiado tal
    cual de las filas de control elegidas (con reemplazo)."""
    n_control = control_pool_X.shape[0]
    out = {}
    for label in labels:
        chosen = rng.integers(0, n_control, size=N_PER_LABEL)
        out[label] = control_pool_X[chosen]
    return out
