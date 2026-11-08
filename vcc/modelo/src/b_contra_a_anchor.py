#!/usr/bin/env python
"""`B - A` a profundidad igualada Y con CONTROL POR MITAD, con barra pareada.

EL SEGUNDO CONFUNDIDOR, que no estaba en el encargo y aparecio al declarar la
construccion. En `cand01_refuta.py` las dos columnas de CADA correlacion se
dividen por el MISMO pseudobulk de control `cr`:

    A  = corr( log2(mitad_1/cr), log2(mitad_2/cr) )
    B  = corr( log2(blanco_p/cr), log2(blanco_q/cr) )

`cell-eval2/src/cell_eval2/anchor.py:9-18` declara que eso sesga la fiabilidad
mitad-contra-mitad HACIA ARRIBA -- "shares the control's sampling noise between
the two quantities whose agreement is being measured and biases the anchor
upward" -- lo midio en 0.5-2.3 % sobre `de_lfc_nmae`, y cita a `ceiling.py`
convirtiendo una fiabilidad de 0.54 en 0.74 por la misma causa. Y
`anchor.py:768-779` REHUSA un ancla que no venga con control por mitad: "a
replicate anchor is only valid under per-half controls ('pred')". No es una
preferencia: es el estimando que el marcador declara no negociable.

Asi que la version conforme al ancla es:

    A_anc = corr( log2(mitad_1(p)/c1), log2(mitad_2(p)/c2) )
    B_anc = corr( log2(mitad_1(p)/c1), log2(mitad_2(q)/c2) )

con `c1`, `c2` las dos mitades disjuntas del pool de control. Las dos
cantidades: 200 celulas por lado, controles disjuntos entre los dos lados, y la
MISMA estructura mitad_1-contra-mitad_2 en las dos. Profundidad igualada y ruido
de control no compartido en ninguna de las dos.

La barra es la misma del encargo: bootstrap PAREADO sobre los 33 BLANCOS,
semillas 20267000 + 100*fold + replica, particion nueva por replica, y
`B - A` por replica.

Lectura y medicion. No construye nada, no puntua nada.
"""
import itertools
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from suelo_arnes import (IFNB_PATH, build_fold_adatas, load_global_obs,
                         load_var_names)
from cand01_covar import hace_predict_fn

LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
PANEL_PATH = "/home/ubuntu/vcc/modelo/out/suelo_panel.json"
OUT = "/home/ubuntu/vcc/modelo/out/b_contra_a_anchor.json"
SEED_BASE_REAL = 20260001
SEED_BASE_EMISION = 20266000
SEED_MITADES = 777
SEED_FAMILIA_BOOT = 20267000
SEED_CTRL_SPLIT = 20267999
N_REPLICAS = 200
UMBRAL_EXPR = 1.0


def rankea(M):
    orden = np.argsort(M, axis=0, kind="stable")
    r = np.empty_like(M, dtype=np.float64)
    ar = np.arange(M.shape[0], dtype=np.float64)
    for k in range(M.shape[1]):
        r[orden[:, k], k] = ar
    return r


def normaliza(R):
    R = R - R.mean(axis=0, keepdims=True)
    s = np.sqrt((R * R).sum(axis=0))
    s[s == 0] = 1.0
    return R / s


def nrank(M):
    return normaliza(rankea(M))


def main():
    n_rep = int(sys.argv[1]) if len(sys.argv) > 1 else N_REPLICAS
    panel = json.load(open(PANEL_PATH))["panel_comun"]
    nb = len(panel)
    cell_type, gene = load_global_obs(IFNB_PATH)
    var_names = load_var_names(IFNB_PATH)
    n_genes = len(var_names)
    salida = {"meta": {"n_blancos": nb, "n_replicas": n_rep,
                       "semilla_familia": SEED_FAMILIA_BOOT,
                       "semilla_control_partido": SEED_CTRL_SPLIT,
                       "panel": panel}, "lineas": {}}

    for line in LINES:
        t0 = time.time()
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
        ne = int(exp.sum())
        idx_por_blanco = [np.where(lr == g)[0] for g in panel]

        def D(cols, denom):
            with np.errstate(divide="ignore", invalid="ignore"):
                M = np.column_stack([np.log2(pb(Xr, c)[exp] / denom[exp])
                                     for c in cols])
            M[~np.isfinite(M)] = 0.0
            return M

        # ---- punto, con la 777 y el control partido con su propia semilla ----
        rng_c = np.random.default_rng(SEED_CTRL_SPLIT + fold_idx)
        permc = rng_c.permutation(ctrl_idx)
        mc = len(permc) // 2
        cr1 = pb(Xr, permc[:mc]); cr2 = pb(Xr, permc[mc:])
        rng777 = np.random.default_rng(SEED_MITADES)
        p1, p2 = [], []
        for ix in idx_por_blanco:
            perm = rng777.permutation(ix)
            m = len(perm) // 2
            p1.append(perm[:m]); p2.append(perm[m:])
        n1 = nrank(D(p1, cr1)); n2 = nrank(D(p2, cr2))
        M = n1.T @ n2
        A_pt = [float(M[i, i]) for i in range(nb)]
        B_pt = [float(M[i, k]) for i, k in itertools.combinations(range(nb), 2)] \
             + [float(M[k, i]) for i, k in itertools.combinations(range(nb), 2)]

        boot = {k: [] for k in ("A_med", "B_med", "A_mean", "B_mean",
                                "d_med", "d_mean", "n_pares", "n_distintos")}
        for rep in range(n_rep):
            rng = np.random.default_rng(SEED_FAMILIA_BOOT + 100 * fold_idx + rep)
            boot_idx = rng.integers(0, nb, size=nb)
            q1, q2 = [], []
            for ix in idx_por_blanco:
                perm = rng.permutation(ix)
                m = len(perm) // 2
                q1.append(perm[:m]); q2.append(perm[m:])
            # el control se reparte de nuevo en cada replica, con la misma rng
            pc = rng.permutation(ctrl_idx)
            k = len(pc) // 2
            d1 = pb(Xr, pc[:k]); d2 = pb(Xr, pc[k:])
            g1 = nrank(D(q1, d1)); g2 = nrank(D(q2, d2))
            Mg = g1.T @ g2
            a = Mg[boot_idx, boot_idx]
            jj, kk = np.triu_indices(nb, k=1)
            ok = boot_idx[jj] != boot_idx[kk]
            bi, bk = boot_idx[jj][ok], boot_idx[kk][ok]
            b = np.concatenate([Mg[bi, bk], Mg[bk, bi]])
            boot["A_med"].append(float(np.median(a)))
            boot["A_mean"].append(float(np.mean(a)))
            boot["B_med"].append(float(np.median(b)))
            boot["B_mean"].append(float(np.mean(b)))
            boot["d_med"].append(float(np.median(b) - np.median(a)))
            boot["d_mean"].append(float(np.mean(b) - np.mean(a)))
            boot["n_distintos"].append(int(len(np.unique(boot_idx))))
            boot["n_pares"].append(int(len(b)))

        salida["lineas"][line] = {
            "fold_idx": fold_idx, "n_genes_expresados": ne,
            "n_control_cells": meta["n_control_cells"],
            "n_control_por_mitad": [int(mc), int(len(permc) - mc)],
            "A_anclado": A_pt, "B_anclado": B_pt, "bootstrap": boot,
            "segundos": round(time.time() - t0, 1),
        }
        d = np.array(boot["d_med"])
        lo, hi = np.percentile(d, [2.5, 97.5])
        sys.stderr.write(
            "%-6s ne=%5d A=%+.4f B=%+.4f d=%+.4f IC95=[%+.4f,%+.4f] %s  %.0fs\n"
            % (line, ne, float(np.median(A_pt)), float(np.median(B_pt)),
               float(np.median(B_pt)) - float(np.median(A_pt)), lo, hi,
               "EXCLUYE_0" if (lo > 0 or hi < 0) else "incluye_0",
               time.time() - t0))
        sys.stderr.flush()
        with open(OUT, "w") as fh:
            json.dump(salida, fh, indent=1)
    sys.stderr.write("escrito %s\n" % OUT)


if __name__ == "__main__":
    main()
