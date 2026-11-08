#!/usr/bin/env python
"""VIA DE ESTRUCTURA EXTERNA -- la emision a partir del perfil de k562_gwps.

QUE ES LA VIA. `control_o_no.md` midio que el coeficiente correcto es distinto en
cada blanco, asi que no hay un mapa compartido que 33 blancos revelen. La via de
estructura externa trae ese mapa DE FUERA: una tabla que dice, por blanco, QUE
GENES se mueven y EN QUE SENTIDO. Aqui la tabla es el perfil de respuesta de
k562_gwps.

🔴 DOS POBLACIONES DISTINTAS. gwps = K562 SIN estimular, knockdown CRISPRi.
Jiang = seis lineas CON IFNB. No son la misma cosa, y nada de esta construccion
supone que lo sean.

EL PERFIL, igual que `gwps_jiang_datos.md` §1.1:
    p_g = log2( (X_g + eps) / (media de X sobre las 514 filas core_control + eps) )
    archivo /datos/ref/replogle/k562_gwps_raw_bulk.h5ad, matriz `X`, eps = 1e-8

🔴 LA MAGNITUD NO VIAJA, y sigue declarada como hueco desde `candidato_01.md`.
De `p_g` se usan SOLO EL ORDEN Y EL SIGNO. La magnitud que se aplica es una
ELECCION DE ESCALA, no una cifra medida:

    dentro de los genes seleccionados que BAJAN en gwps (m de ellos), se ordenan
    por |p| descendente y el i-esimo recibe   s = -(m - i + 1)/m
    dentro de los que SUBEN (n), igual con signo +:   s = +(n - i + 1)/n
    r_h = 1 + s_h

El gen que mas baja en gwps recibe `s = -1` y se apaga del todo; el que menos,
`s = -1/m`. `r` vive en **[0, 2]** por construccion, es **adimensional**, y es
funcion UNICAMENTE del rango -- el valor de `log2FC` de gwps no entra en ninguna
cifra, solo su orden y su signo. Es la misma clase de eleccion que la v3 de
`candidato_01` (`r = 1 - rho`, acotada en [0,2]): acotada, adimensional y sin
ninguna constante elegida aqui. **No es una medicion de cuanto baja un gen.**

🔴 SIN PUERTA DE SIGNIFICANCIA, y la razon va declarada. La via se especifica
como "el perfil da el ORDEN y el SIGNO", asi que entran los 8124 genes mapeados y
la `s` de cada uno sale de su rango dentro de su sentido. Meter una puerta
anadiria una constante de seleccion que no esta en la especificacion.

Se MIDIO la puerta analoga a la de `candidato_01` (BH al 0.05 del preset sobre
`z = p_g / sd_nulo`, con `sd_nulo` la sd de las 514 filas `core_control`) y el
resultado esta en el reporte: selecciona de **0 a 829 genes** por blanco de los
8248, con **mediana 1**, **7 de los 27 blancos en cero** y solo 4 por encima de
90 (`CEBPG` 829, `ETS2` 630, `BRD9` 127, `TYK2` 96). **Esa variante NO se
puntuo**, y eso queda declarado como no medido -- elegir entre las dos mirando
el `s` esta prohibido.

🔴 EL REPLIEGUE AL SUELO, declarado, no elegido. Los 6 blancos del panel que NO
estan en gwps (DTX3L, H1F0, HERC5, HERC6, MAP3K14, ZBP1) reciben `r_h = 1` en
todos los genes, que es EXACTAMENTE la emision del suelo (control remuestreado
sin tocar). No hay con que emitirlos por esta via.

🔴 ID2 y RAP1GAP tienen P1 y P2 en gwps. **No se promedian.** Se usa **P1** y se
declara; `gwps_emision_p2()` da la variante con P2 de esas dos filas, que se
mide y se reporta aparte.

LOS GENES: gwps tiene 8248 columnas, 2 nombres duplicados (`HSPA14`, `TBCE`) que
se excluyen enteros. De los 8246 restantes, 8124 caen en los 34025 `var_names`
de Jiang. Los genes de Jiang que no estan en gwps reciben `r_h = 1`.

EL CONTRATO DE EMISION: conteos ENTEROS muestreados, nunca un perfil suave
(`paso1_contrato.md`, HARD_CAP_NNZ). Se reusa `_reescala_entero` de
`cand01_covar` SIN MODIFICAR: binomial para `r<=1`, Poisson aditivo para `r>1`,
esperanza `x*r` en las dos ramas y enteros en las dos.

NO TRANSDUCTIVO. `predict_fn` recibe `(control_pool_X, labels, rng)` y nada mas.
La tabla de gwps entra por clausura y es la ESTRUCTURA EXTERNA que la via trae a
proposito: no es el lado perturbado de Jiang, ni ninguna cifra del `b` o del
anchor, ni nada de `~/vcc/benchmark/`. Los ratios NO dependen del pool de
control -- son los mismos en las seis lineas -- asi que esta emision no estima
nada en contexto, y eso queda declarado como propiedad de la via.
"""
import numpy as np
import scipy.sparse as sp

from cand01_covar import _bh_mascara, _reescala_entero, N_PER_LABEL, CONTROL_LABEL

GWPS = "/datos/ref/replogle/k562_gwps_raw_bulk.h5ad"
EPS = 1e-8
FDR_Q = 0.05          # de.p_adj_threshold del preset vcc2026
SIN_GWPS = ("DTX3L", "H1F0", "HERC5", "HERC6", "MAP3K14", "ZBP1")


def _carga_gwps():
    """Devuelve (perfiles, sd_nulo, nombres) de gwps. `perfiles[(gen,constr)]`."""
    import h5py
    f = h5py.File(GWPS, "r")
    gt = np.array([x.decode() if isinstance(x, bytes) else x
                   for x in f["obs/gene_transcript"][:]])
    partes = [x.split("_") for x in gt]
    g_gen = np.array([p[1] if len(p) > 1 else "" for p in partes])
    g_con = np.array([p[2] if len(p) > 2 else "" for p in partes])
    core = f["obs/core_control"][:]
    cats = np.array([x.decode() if isinstance(x, bytes) else x
                     for x in f["var/__categories/gene_name"][:]])
    vn = cats[f["var/gene_name"][:]]
    u, c = np.unique(vn, return_counts=True)
    dup = set(u[c > 1])
    ok = np.array([g not in dup for g in vn])

    X = f["X"]
    idx_c = np.where(core)[0]
    # media de control y sd del nulo, de las PROPIAS filas de control
    bloque = np.empty((len(idx_c), X.shape[1]), dtype=np.float64)
    for s in range(0, len(idx_c), 64):
        sl = idx_c[s:s + 64]
        bloque[s:s + len(sl)] = X[sl, :].astype(np.float64)
    ctrl = bloque.mean(axis=0)
    lr = np.log2((bloque + EPS) / (ctrl + EPS))
    sd_nulo = lr.std(axis=0, ddof=1)
    del bloque, lr

    return f, X, g_gen, g_con, ctrl, sd_nulo, vn, ok


def construye_ratios(var_names, targets, usar_p2=()):
    """`ratios[g]` = vector (len(var_names),) de factores `r_h` en [0,2].

    `usar_p2` es el conjunto de blancos para los que se toma la fila P2 en vez
    de la P1. Por defecto se usa P1 en los dos blancos que tienen las dos.
    NO se promedia nunca.
    """
    f, X, g_gen, g_con, ctrl, sd_nulo, gvn, ok_var = _carga_gwps()
    n_j = len(var_names)
    pos_j = {g: i for i, g in enumerate(var_names)}

    # mapa gwps -> Jiang, solo nombres no duplicados presentes en los dos
    col_g, col_j = [], []
    for i, g in enumerate(gvn):
        if ok_var[i] and g in pos_j:
            col_g.append(i); col_j.append(pos_j[g])
    col_g = np.array(col_g); col_j = np.array(col_j)

    ratios, diag = {}, {}
    for g in targets:
        if g in SIN_GWPS:
            ratios[g] = np.ones(n_j)
            diag[g] = {"repliegue_suelo": True, "motivo": "no esta en gwps",
                       "construccion": None, "n_selec": 0,
                       "n_baja": 0, "n_sube": 0, "r_min": 1.0, "r_max": 1.0}
            continue
        cand = np.where(g_gen == g)[0]
        cons = [g_con[i] for i in cand]
        quiere = "P2" if g in usar_p2 else "P1"
        if len(cand) == 1:
            fila = int(cand[0]); usada = cons[0]
        else:
            # tiene P1 y P2: se elige una, NUNCA se promedia
            k = cons.index(quiere) if quiere in cons else 0
            fila = int(cand[k]); usada = cons[k]

        x = X[fila, :].astype(np.float64)
        p = np.log2((x + EPS) / (ctrl + EPS))

        # sin puerta: entran todos los genes con perfil finito y no nulo
        pasa = np.isfinite(p)

        # ORDEN y SIGNO: s por rango normalizado dentro de cada sentido
        s = np.zeros(len(gvn))
        for sentido, sgn in ((p < 0, -1.0), (p > 0, +1.0)):
            m = pasa & sentido
            k = int(m.sum())
            if k == 0:
                continue
            ii = np.where(m)[0]
            orden = ii[np.argsort(-np.abs(p[ii]), kind="stable")]
            s[orden] = sgn * (np.arange(k, 0, -1) / float(k))

        rat = np.ones(n_j)
        rat[col_j] = 1.0 + s[col_g]
        rat = np.clip(rat, 0.0, 2.0)
        ratios[g] = rat
        sel_j = (rat != 1.0)
        diag[g] = {
            "repliegue_suelo": False, "construccion": usada,
            "n_selec_gwps": int(pasa.sum()),
            "n_selec_en_jiang": int(sel_j.sum()),
            "n_baja": int((rat < 1.0).sum()), "n_sube": int((rat > 1.0).sum()),
            "r_min": float(rat.min()), "r_max": float(rat.max()),
            "gen_blanco_en_gwps": bool(g in set(gvn[col_g])),
        }
    diag["_gwps"] = {
        "n_genes_gwps": int(len(gvn)), "n_duplicados": int(len(gvn) - int(ok_var.sum())),
        "n_mapeados_a_jiang": int(len(col_j)), "n_var_jiang": int(n_j),
        "n_core_control": int(len(ctrl)) and int(sd_nulo.size),
        "blancos_replegados": list(SIN_GWPS),
        "p2_usado_en": sorted(usar_p2),
    }
    f.close()
    return ratios, diag


def hace_predict_fn(var_names, registro=None, usar_p2=()):
    """`predict_fn` que el arnes acepta SIN MODIFICARLO."""
    cache = {}

    def predict_gwps(control_pool_X, labels, rng):
        if "r" not in cache:
            targets = [l for l in labels if l != CONTROL_LABEL]
            cache["r"] = construye_ratios(var_names, targets, usar_p2)
        ratios, diag = cache["r"]
        if registro is not None:
            registro["diag"] = diag
        n_control = control_pool_X.shape[0]
        out = {}
        for label in labels:
            chosen = rng.integers(0, n_control, size=N_PER_LABEL)
            bloque = control_pool_X[chosen]
            if label == CONTROL_LABEL:
                out[label] = bloque.tocsr()       # el control predicho = el suelo
            else:
                out[label] = _reescala_entero(bloque, ratios[label], rng)
        return out

    return predict_gwps
