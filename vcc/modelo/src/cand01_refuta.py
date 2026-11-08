#!/usr/bin/env python
"""La refutacion del candidato 01, medida en las SEIS lineas con filas por
perturbacion.

LA PREGUNTA: el orden de coexpresion del pool de control, ¿predice que genes se
mueven cuando se perturba el blanco?

LOS CUATRO NUMEROS, todos Spearman del delta de pseudobulk contra el control,
sobre los genes expresados (>= 1 en CP50k en los DOS controles):
  A) mitad-contra-mitad del REAL -- el TECHO alcanzable a 400 celulas. Si esto
     fuera ~0 el test no tendria potencia y no se podria concluir nada.
  B) real de un blanco contra el real de OTRO -- la respuesta GENERICA
     compartida. Es la cantidad que el `b` del marcador captura.
  C) la prediccion del candidato contra el real del MISMO blanco -- lo que se
     quiere medir.
  D) la prediccion del candidato contra el real de OTRO blanco -- el control
     negativo. `C - D` es lo ESPECIFICO DEL BLANCO que el candidato captura.

Y una medida de signo que no depende de la escala: en los genes que el lado real
mueve de verdad (`|dlog_real| > 0.25`), la fraccion en que el signo predicho
coincide con el real. 0.5 es azar exacto.

Lectura y medicion. No ajusta nada, no puntua nada, no toca ningun bundle.
"""
import itertools
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
OUT = "/home/ubuntu/vcc/modelo/out/cand01_refuta.json"
SEED_BASE_REAL = 20260001
SEED_BASE_EMISION = 20266000
SEED_MITADES = 777          # solo para partir el real en mitades, no toca nada mas
UMBRAL_EXPR = 1.0           # CP50k, el mismo orden que el bulk_target_sum del preset
UMBRAL_SIGNO = 0.25         # |dlog_real| por encima del cual el signo es medible


def rankea(M):
    """Rangos por columna para hacer Spearman como Pearson sobre rangos: con
    3168 pares es lo unico que hace la medicion barata."""
    orden = np.argsort(M, axis=0, kind="stable")
    r = np.empty_like(M, dtype=np.float64)
    n = M.shape[0]
    ar = np.arange(n, dtype=np.float64)
    for k in range(M.shape[1]):
        r[orden[:, k], k] = ar
    return r


def normaliza(R):
    R = R - R.mean(axis=0, keepdims=True)
    s = np.sqrt((R * R).sum(axis=0))
    s[s == 0] = 1.0
    return R / s


def main():
    panel = json.load(open(PANEL_PATH))["panel_comun"]
    cell_type, gene = load_global_obs(IFNB_PATH)
    var_names = load_var_names(IFNB_PATH)
    n_genes = len(var_names)
    salida = {}

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

        cr = pb(Xr, np.where(lr == "NT")[0])
        cp = pb(Xp, np.where(lp == "NT")[0])
        exp = (cr >= UMBRAL_EXPR) & (cp >= UMBRAL_EXPR)
        ne = int(exp.sum())

        rng = np.random.default_rng(SEED_MITADES)
        D_real, D_pred, D_h1, D_h2 = [], [], [], []
        for g in panel:
            ix = np.where(lr == g)[0]
            perm = rng.permutation(ix)
            mitad = len(perm) // 2
            D_real.append(np.log2(pb(Xr, ix)[exp] / cr[exp]))
            D_h1.append(np.log2(pb(Xr, perm[:mitad])[exp] / cr[exp]))
            D_h2.append(np.log2(pb(Xr, perm[mitad:])[exp] / cr[exp]))
            D_pred.append(np.log2(pb(Xp, np.where(lp == g)[0])[exp] / cp[exp]))
        R = np.column_stack(D_real); P = np.column_stack(D_pred)
        H1 = np.column_stack(D_h1); H2 = np.column_stack(D_h2)
        for A in (R, P, H1, H2):
            A[~np.isfinite(A)] = 0.0

        nR = normaliza(rankea(R)); nP = normaliza(rankea(P))
        nH1 = normaliza(rankea(H1)); nH2 = normaliza(rankea(H2))

        A = [float(nH1[:, i] @ nH2[:, i]) for i in range(len(panel))]
        C = [float(nP[:, i] @ nR[:, i]) for i in range(len(panel))]
        MC = nP.T @ nR          # fila = prediccion de a, columna = real de b
        MR = nR.T @ nR
        B, D = [], []
        for i, k in itertools.combinations(range(len(panel)), 2):
            B.append(float(MR[i, k]))
        for i in range(len(panel)):
            for k in range(len(panel)):
                if i != k:
                    D.append(float(MC[i, k]))

        signo = []
        for i in range(len(panel)):
            grande = np.abs(R[:, i]) > UMBRAL_SIGNO
            signo.append(float(np.mean(np.sign(P[grande, i])
                                       == np.sign(R[grande, i])))
                         if grande.sum() else None)

        diag = reg.get("diag", {})
        pool = diag.pop("_pool", {})
        salida[line] = {
            "fold_idx": fold_idx, "n_genes_expresados": ne,
            "n_control_cells": meta["n_control_cells"],
            "pool": pool,
            "A_mitades_real": A, "B_cruzado_real": B,
            "C_pred_vs_mismo": C, "D_pred_vs_otro": D,
            "concordancia_signo": signo,
            "por_pert": {
                g: {
                    "A_mitades": A[i], "C_pred_vs_mismo": C[i],
                    "concordancia_signo": signo[i],
                    "n_real_movidos": int((np.abs(R[:, i]) > UMBRAL_SIGNO).sum()),
                    "n_pred_movidos": int((np.abs(P[:, i]) > UMBRAL_SIGNO).sum()),
                    "n_responde": diag.get(g, {}).get("n_responde"),
                    "corr_abs_max_pool": diag.get(g, {}).get("corr_abs_max"),
                    "r_min": diag.get(g, {}).get("r_min"),
                    "r_max": diag.get(g, {}).get("r_max"),
                    "r_mediano": diag.get(g, {})
                        .get("r_mediano_de_los_que_responden"),
                }
                for i, g in enumerate(panel)
            },
        }
        sys.stderr.write(
            "%-6s expresados=%5d  A=%+.4f  B=%+.4f  C=%+.4f  D=%+.4f  signo=%.4f\n"
            % (line, ne, float(np.median(A)), float(np.median(B)),
               float(np.median(C)), float(np.median(D)),
               float(np.nanmean([s for s in signo if s is not None]))))

    with open(OUT, "w") as fh:
        json.dump(salida, fh, indent=1)
    sys.stderr.write(f"\nescrito {OUT}\n")


if __name__ == "__main__":
    main()
