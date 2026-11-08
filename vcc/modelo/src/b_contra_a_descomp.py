#!/usr/bin/env python
"""Descomposicion de la varianza del delta REAL: compartida entre blancos contra
especifica del blanco. Sin suponer ningun modelo.

POR QUE ESTA MEDICION Y NO LA HIPOTESIS. `candidato_01.md` uso `B > A` para
afirmar que el `b` "captura la componente dominante de la varianza". Si `B > A`
cae a profundidad igualada, esa afirmacion se queda sin su apoyo y hay que
medirla directamente. Esto la mide.

COMO, SIN MODELO. Lo unico que se supone es que las dos mitades de un blanco
tienen ruido INDEPENDIENTE (son celulas disjuntas) y de media cero. Con
`h1_i`, `h2_i` el delta de las dos mitades del blanco `i`, y
`m1 = media_i h1_i`, `m2 = media_i h2_i`:

    V_tot = (1/n) Σ_i <h1_i, h2_i> / G        -> la senal reproducible total
    V_sh  = <m1, m2> / G                      -> la parte COMPARTIDA
    V_sp  = (1/n) Σ_i <h1_i-m1, h2_i-m2> / G  -> la parte ESPECIFICA

`V_tot = V_sh + V_sp` es una identidad algebraica exacta, no un ajuste: los dos
terminos cruzados se anulan porque Σ_i (h1_i - m1) = 0. Y como el ruido de las
dos mitades es independiente, su esperanza no lleva ruido: es el producto de las
senales verdaderas. Por eso se mide con el producto CRUZADO de mitades y no con
la varianza de una sola mitad -- esa si llevaria el ruido dentro.

EL CONTROL TIENE QUE IR PARTIDO, y no es un detalle. Todo delta es
`log2(pb) - log2(cr)`. Si las dos mitades se dividen por el MISMO `cr`, el
termino `-log2(cr)` es identico en las dos y cae entero dentro de `m1` y `m2`
con la MISMA realizacion de ruido: `<m1, m2>` entonces cuenta el ruido de
muestreo del control como senal compartida. Con el control partido en dos
mitades independientes (`cr1` para `h1`, `cr2` para `h2`) ese ruido es
independiente entre `m1` y `m2` y no cuenta. Se reportan las dos.

Lectura y medicion. No construye nada, no puntua nada.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from suelo_arnes import (IFNB_PATH, build_fold_adatas, load_global_obs,
                         load_var_names)
from cand01_covar import hace_predict_fn

LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
PANEL_PATH = "/home/ubuntu/vcc/modelo/out/suelo_panel.json"
OUT = "/home/ubuntu/vcc/modelo/out/b_contra_a_descomp.json"
SEED_BASE_REAL = 20260001
SEED_BASE_EMISION = 20266000
SEED_MITADES = 777
SEED_CTRL_SPLIT = 20267999
UMBRAL_EXPR = 1.0


def descompone(H1, H2):
    """H1, H2: (G, n) deltas de las dos mitades. Centra cada columna sobre los
    genes (quita el offset escalar de la columna, que no es perfil) y devuelve
    las tres cantidades y las fracciones."""
    G, n = H1.shape
    A = H1 - H1.mean(axis=0, keepdims=True)
    B = H2 - H2.mean(axis=0, keepdims=True)
    m1 = A.mean(axis=1); m2 = B.mean(axis=1)
    v_tot = float(np.einsum("gi,gi->", A, B) / (n * G))
    v_sh = float(m1 @ m2 / G)
    Ac = A - m1[:, None]; Bc = B - m2[:, None]
    v_sp = float(np.einsum("gi,gi->", Ac, Bc) / (n * G))
    return {
        "V_tot": v_tot, "V_compartida": v_sh, "V_especifica": v_sp,
        "identidad_resid": float(v_tot - (v_sh + v_sp)),
        "frac_compartida": (v_sh / v_tot) if v_tot != 0 else None,
        "frac_especifica": (v_sp / v_tot) if v_tot != 0 else None,
        "n_blancos": int(n), "n_genes": int(G),
        "por_blanco_especifica": [
            float(Ac[:, i] @ Bc[:, i] / G) for i in range(n)],
        "por_blanco_total": [float(A[:, i] @ B[:, i] / G) for i in range(n)],
    }


def main():
    panel = json.load(open(PANEL_PATH))["panel_comun"]
    cell_type, gene = load_global_obs(IFNB_PATH)
    var_names = load_var_names(IFNB_PATH)
    n_genes = len(var_names)
    salida = {"meta": {"panel": panel, "semilla_mitades": SEED_MITADES,
                       "semilla_control_partido": SEED_CTRL_SPLIT}, "lineas": {}}

    for line in LINES:
        fold_idx = LINES.index(line)
        reg = {}
        pf = hace_predict_fn(var_names, reg)
        real, pred, meta = build_fold_adatas(
            line, panel, pf, var_names, SEED_BASE_REAL + fold_idx,
            SEED_BASE_EMISION + 100 * fold_idx, cell_type, gene, n_genes,
            IFNB_PATH)
        Xr = real.X.tocsr(); lr = np.array(real.obs["target_gene"])
        Xp = pred.X.tocsr(); lp = np.array(pred.obs["target_gene"])

        def pb(M, idx):
            s = np.asarray(M[idx].sum(axis=0)).ravel()
            tot = s.sum()
            return s / tot * 5e4 if tot > 0 else s

        ctrl_idx = np.where(lr == "NT")[0]
        cr = pb(Xr, ctrl_idx)
        cp = pb(Xp, np.where(lp == "NT")[0])
        exp = (cr >= UMBRAL_EXPR) & (cp >= UMBRAL_EXPR)

        rng_c = np.random.default_rng(SEED_CTRL_SPLIT + fold_idx)
        permc = rng_c.permutation(ctrl_idx)
        mc = len(permc) // 2
        cr1 = pb(Xr, permc[:mc]); cr2 = pb(Xr, permc[mc:])

        rng = np.random.default_rng(SEED_MITADES)
        p1, p2 = [], []
        for g in panel:
            ix = np.where(lr == g)[0]
            perm = rng.permutation(ix)
            m = len(perm) // 2
            p1.append(perm[:m]); p2.append(perm[m:])

        def D(cols, denom):
            with np.errstate(divide="ignore", invalid="ignore"):
                M = np.column_stack([np.log2(pb(Xr, c)[exp] / denom[exp])
                                     for c in cols])
            M[~np.isfinite(M)] = 0.0
            return M

        comun = descompone(D(p1, cr), D(p2, cr))
        indep = descompone(D(p1, cr1), D(p2, cr2))
        salida["lineas"][line] = {
            "n_genes_expresados": int(exp.sum()),
            "n_control_cells": meta["n_control_cells"],
            "n_control_por_mitad": [int(mc), int(len(permc) - mc)],
            "control_comun": comun, "control_independiente": indep,
        }
        sys.stderr.write(
            "%-6s  comun: comp=%.3f esp=%.3f (V=%.4f)   indep: comp=%.3f esp=%.3f (V=%.4f)\n"
            % (line, comun["frac_compartida"], comun["frac_especifica"],
               comun["V_tot"], indep["frac_compartida"],
               indep["frac_especifica"], indep["V_tot"]))
        sys.stderr.flush()

    with open(OUT, "w") as fh:
        json.dump(salida, fh, indent=1)
    sys.stderr.write("escrito %s\n" % OUT)


if __name__ == "__main__":
    main()
