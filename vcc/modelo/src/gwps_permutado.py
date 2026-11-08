#!/usr/bin/env python
"""EL CONTROL PERMUTADO de la via de estructura externa: un DESARREGLO, no una
emision nueva.

🔴 LA PREGUNTA. `fid` pregunta "de los genes que de verdad se mueven, acerte el
signo?". NO pregunta si acerte el signo DE ESE BLANCO. Si gwps emite una
direccion COMPARTIDA entre los 33 (un eje comun), `fid` la premia igual. Este
brazo mide cuanto del `fid` sobrevive cuando el perfil se asigna al blanco
EQUIVOCADO.

LA CONSTRUCCION:
  * para la perturbacion `h` se emite el perfil de gwps de OTRO blanco `g != h`
  * el desarreglo es SOLO sobre los 27 blancos CON perfil
  * 🔴 los 6 replegados (DTX3L, H1F0, HERC5, HERC6, MAP3K14, ZBP1) siguen
    replegados al suelo en LOS DOS BRAZOS, asi que su aporte es IDENTICO y la
    comparacion queda limpia sobre los 27
  * 🔴 TODO LO DEMAS SIN TOCAR: la misma regla de escala (orden y signo, `r` en
    [0,2]) de `gwps_emision.construye_ratios`, el mismo `_reescala_entero`, el
    mismo arnes, la misma poblacion real (`20260001 + fold`), el mismo bundle.
    Lo unico que cambia es QUE perfil recibe cada etiqueta.

🔴 POR QUE ESTE CONTROL Y NO OTRO: el brazo permutado emite el MISMO numero de
genes con la MISMA distribucion de magnitud -- es el mismo conjunto de 27
vectores `r`, solo reasignado. Esta IGUALADO EN POTENCIA por construccion.

SEMILLAS:
  * desarreglo: `20272500 + j`, j = 0..7. **LOS MISMOS 8 desarreglos en las seis
    lineas**, para que las lineas sean comparables. Se comprueba y se reporta
    que ninguno tiene punto fijo.
  * remuestreo de la emision: familia nueva `20272000 + 100*fold + j`. Sin
    colision con 2026000x-2026600x ni 2027100x.

P1 en los dos brazos, igual que la corrida principal. No se promedia P1 con P2.
"""
import numpy as np

from gwps_emision import construye_ratios, SIN_GWPS
from cand01_covar import _reescala_entero, N_PER_LABEL, CONTROL_LABEL

SEED_DESARREGLO = 20272500


def desarreglo(etiquetas, j):
    """Permutacion SIN PUNTO FIJO de `etiquetas`, con semilla 20272500 + j.

    Rechazo: se sortea una permutacion y se descarta si deja algun punto fijo.
    Devuelve (sigma, diag) con `sigma[h] = g` el blanco cuyo perfil recibe `h`.
    """
    et = list(etiquetas)
    rng = np.random.default_rng(SEED_DESARREGLO + j)
    intentos = 0
    while True:
        intentos += 1
        perm = rng.permutation(len(et))
        if not any(perm[i] == i for i in range(len(et))):
            break
        if intentos > 10000:
            raise RuntimeError("no se encontro desarreglo")
    sigma = {et[i]: et[perm[i]] for i in range(len(et))}
    fijos = [h for h, g in sigma.items() if h == g]
    return sigma, {"j": j, "semilla": SEED_DESARREGLO + j, "intentos": intentos,
                   "n": len(et), "puntos_fijos": fijos,
                   "sin_punto_fijo": len(fijos) == 0,
                   "sigma": {h: g for h, g in sorted(sigma.items())}}


def hace_predict_fn_perm(var_names, registro, j):
    """`predict_fn` del brazo PERMUTADO. El arnes no se modifica."""
    cache = {}

    def predict_perm(control_pool_X, labels, rng):
        if "r" not in cache:
            targets = [l for l in labels if l != CONTROL_LABEL]
            ratios, diag = construye_ratios(var_names, targets)   # SIN TOCAR
            conper = [t for t in targets if t not in SIN_GWPS]
            sigma, dd = desarreglo(conper, j)
            cache["r"] = (ratios, sigma)
            if registro is not None:
                registro["diag"] = diag
                registro["desarreglo"] = dd
        ratios, sigma = cache["r"]
        n_control = control_pool_X.shape[0]
        out = {}
        for label in labels:
            chosen = rng.integers(0, n_control, size=N_PER_LABEL)
            bloque = control_pool_X[chosen]
            if label == CONTROL_LABEL:
                out[label] = bloque.tocsr()
            else:
                # 🔴 los replegados conservan su r = 1 (suelo); los 27 reciben
                # el perfil de OTRO blanco
                fuente = sigma.get(label, label)
                out[label] = _reescala_entero(bloque, ratios[fuente], rng)
        return out

    return predict_perm
