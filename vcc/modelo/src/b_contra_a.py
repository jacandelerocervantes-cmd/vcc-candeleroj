#!/usr/bin/env python
"""B contra A a PROFUNDIDAD IGUALADA, con barra pareada por bootstrap sobre BLANCOS.

LA HIPOTESIS QUE SE CIERRA. `candidato_01.md` §1 midio `B > A` en cinco de las
seis lineas y lo dejo como hipotesis por dos razones: 5 de 6 bajo una moneda
tiene p = 7/64 = 0.109, y ninguna de las dos cantidades llevaba dispersion.

EL CONFUNDIDOR QUE VA PRIMERO. En `cand01_refuta.py`:
  A  se calcula entre MITADES del mismo blanco -- 200 celulas por lado.
  B  se calcula entre poblaciones COMPLETAS de dos blancos -- 400 por lado.
`A` esta atenuada por profundidad respecto de `B` POR CONSTRUCCION. Aqui NO se
corrige `A` hacia arriba (esa es la correccion Spearman-Brown que el anchor
rechaza explicitamente): se IGUALA la profundidad.

  A_mitades  : mitad_1 del blanco p  vs  mitad_2 del blanco p     (200 vs 200)
  B_completa : poblacion de p        vs  poblacion de q           (400 vs 400)
  B_mitades  : mitad_1 del blanco p  vs  mitad_1 del blanco q     (200 vs 200)

LA BARRA. Bootstrap PAREADO sobre los 33 BLANCOS, no sobre los 528 pares: dos
pares que comparten un blanco no son independientes y remuestrear pares daria
una barra falsamente estrecha. En cada replica se remuestrean los 33 blancos con
reemplazo, se vuelve a partir en mitades con semilla propia, y se recalculan `A`
y `B` SOBRE LOS MISMOS blancos remuestreados. La cantidad de interes es
`B - A` por replica, PAREADA.

Pares de un blanco consigo mismo (un blanco que sale repetido en la replica) se
EXCLUYEN de `B`: a profundidad igualada `mitad_1(p)` contra `mitad_1(p)` es el
mismo vector y daria 1.0 exacto.

SEMILLAS: familia nueva 20267000 + 100*fold + replica. La particion en mitades
del punto conserva la 777 de `cand01_refuta.py` para que el punto REPRODUZCA
candidato_01 exactamente antes de anadir nada.

Lectura y medicion. No construye ningun candidato, no toca ningun arnes, no
puntua nada, no escribe ningun bundle.
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
OUT = "/home/ubuntu/vcc/modelo/out/b_contra_a.json"
REF = "/home/ubuntu/vcc/modelo/out/cand01_refuta.json"

SEED_BASE_REAL = 20260001        # identica a cand01_refuta.py -> mismo lado real
SEED_BASE_EMISION = 20266000     # identica -> misma mascara `exp`
SEED_MITADES = 777               # identica -> el punto reproduce candidato_01
SEED_FAMILIA_BOOT = 20267000     # familia NUEVA de este encargo
SEED_CTRL_SPLIT = 20267999       # control partido, solo para la sensibilidad
N_REPLICAS = 200
UMBRAL_EXPR = 1.0


def rankea(M):
    """Identica a cand01_refuta.py:48. Rangos por columna via argsort estable:
    los empates NO reciben rango promedio, reciben rangos consecutivos en orden
    de indice de gen. Se replica tal cual para reproducir el punto; la
    consecuencia se mide aparte (n_cero_forzado)."""
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


def nrank(M):
    return normaliza(rankea(M))


def pares_distintos(boot_idx):
    """Pares (j,k), j<k, de ranuras de la replica cuyos blancos subyacentes son
    DISTINTOS. Excluye el par de un blanco consigo mismo."""
    n = len(boot_idx)
    jj, kk = np.triu_indices(n, k=1)
    ok = boot_idx[jj] != boot_idx[kk]
    return boot_idx[jj][ok], boot_idx[kk][ok]


def main():
    lineas = sys.argv[1].split(",") if len(sys.argv) > 1 else LINES
    n_rep = int(sys.argv[2]) if len(sys.argv) > 2 else N_REPLICAS
    panel = json.load(open(PANEL_PATH))["panel_comun"]
    nb = len(panel)
    ref = json.load(open(REF)) if os.path.exists(REF) else {}
    cell_type, gene = load_global_obs(IFNB_PATH)
    var_names = load_var_names(IFNB_PATH)
    n_genes = len(var_names)
    salida = {"meta": {
        "n_blancos": nb, "n_pares_completo": nb * (nb - 1) // 2,
        "n_replicas": n_rep, "semilla_familia": SEED_FAMILIA_BOOT,
        "semilla_mitades_punto": SEED_MITADES,
        "semilla_real": SEED_BASE_REAL, "semilla_emision": SEED_BASE_EMISION,
        "panel": panel, "n_celulas_por_blanco": 400,
    }, "lineas": {}}

    for line in lineas:
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
        # mascara IDENTICA a cand01_refuta.py (depende de `cp`, que es el control
        # predicho = remuestreo de 400 celulas del control real, sin modificar)
        exp = (cr >= UMBRAL_EXPR) & (cp >= UMBRAL_EXPR)
        # mascara de solo-verdad, para la verificacion de que la conclusion no
        # cuelga de la emision del candidato
        exp_v = cr >= UMBRAL_EXPR
        ne = int(exp.sum()); ne_v = int(exp_v.sum())

        idx_por_blanco = [np.where(lr == g)[0] for g in panel]

        # control partido en dos mitades independientes, solo para sensibilidad
        rng_c = np.random.default_rng(SEED_CTRL_SPLIT + fold_idx)
        permc = rng_c.permutation(ctrl_idx)
        mc = len(permc) // 2
        cr1 = pb(Xr, permc[:mc]); cr2 = pb(Xr, permc[mc:])

        def deltas(mask, denom, cols):
            with np.errstate(divide="ignore", invalid="ignore"):
                M = np.column_stack([np.log2(pb(Xr, c)[mask] / denom[mask])
                                     for c in cols])
            nofin = int((~np.isfinite(M)).sum())
            M[~np.isfinite(M)] = 0.0
            return M, nofin

        # ---- punto: profundidad completa (no depende de la particion) ----
        R, nofin_R = deltas(exp, cr, idx_por_blanco)
        nR = nrank(R)
        MR = nR.T @ nR
        B_comp = [float(MR[i, k]) for i, k in itertools.combinations(range(nb), 2)]

        # ---- punto: mitades con la semilla 777 de cand01_refuta ----
        rng777 = np.random.default_rng(SEED_MITADES)
        p1, p2 = [], []
        for ix in idx_por_blanco:
            perm = rng777.permutation(ix)
            m = len(perm) // 2
            p1.append(perm[:m]); p2.append(perm[m:])
        H1, nofin_H1 = deltas(exp, cr, p1)
        H2, nofin_H2 = deltas(exp, cr, p2)
        nH1 = nrank(H1); nH2 = nrank(H2)
        A_pt = [float(nH1[:, i] @ nH2[:, i]) for i in range(nb)]
        M1 = nH1.T @ nH1
        B_mit = [float(M1[i, k]) for i, k in itertools.combinations(range(nb), 2)]

        # verificacion contra candidato_01
        verif = {}
        if line in ref:
            verif = {
                "A_max_abs_dif": float(np.max(np.abs(
                    np.array(A_pt) - np.array(ref[line]["A_mitades_real"])))),
                "B_completa_max_abs_dif": float(np.max(np.abs(
                    np.array(B_comp) - np.array(ref[line]["B_cruzado_real"])))),
                "n_genes_expresados_igual":
                    ne == ref[line]["n_genes_expresados"],
                "n_control_cells_igual":
                    meta["n_control_cells"] == ref[line]["n_control_cells"],
            }

        # ---- sensibilidad: mascara de solo-verdad ----
        Rv, _ = deltas(exp_v, cr, idx_por_blanco); nRv = nrank(Rv)
        MRv = nRv.T @ nRv
        H1v, _ = deltas(exp_v, cr, p1); H2v, _ = deltas(exp_v, cr, p2)
        nH1v = nrank(H1v); nH2v = nrank(H2v); M1v = nH1v.T @ nH1v
        sens_v = {
            "A": [float(nH1v[:, i] @ nH2v[:, i]) for i in range(nb)],
            "B_completa": [float(MRv[i, k])
                           for i, k in itertools.combinations(range(nb), 2)],
            "B_mitades": [float(M1v[i, k])
                          for i, k in itertools.combinations(range(nb), 2)],
            "n_genes": ne_v,
        }

        # ---- sensibilidad: control INDEPENDIENTE en los dos lados ----
        # quita la inflacion comun por compartir el mismo -log2(cr) en las dos
        # columnas de cada correlacion. Mismo 200 vs 200 en las dos cantidades.
        H1c, _ = deltas(exp, cr1, p1); H2c, _ = deltas(exp, cr2, p2)
        nH1c = nrank(H1c); nH2c = nrank(H2c)
        A_ind = [float(nH1c[:, i] @ nH2c[:, i]) for i in range(nb)]
        Mc = nH1c.T @ nH2c      # mitad_1(p)/c1 contra mitad_2(q)/c2
        B_ind = [float(Mc[i, k]) for i in range(nb) for k in range(nb) if i != k]
        sens_c = {"A": A_ind, "B_mitades": B_ind}

        # ---- bootstrap pareado sobre BLANCOS ----
        boot = {k: [] for k in ("A_med", "Bc_med", "Bm_med", "A_mean",
                                "Bc_mean", "Bm_mean", "d_med", "d_mean",
                                "dc_med", "dc_mean", "n_pares", "n_distintos")}
        for rep in range(n_rep):
            rng = np.random.default_rng(SEED_FAMILIA_BOOT + 100 * fold_idx + rep)
            boot_idx = rng.integers(0, nb, size=nb)
            q1, q2 = [], []
            for ix in idx_por_blanco:
                perm = rng.permutation(ix)
                m = len(perm) // 2
                q1.append(perm[:m]); q2.append(perm[m:])
            G1, _ = deltas(exp, cr, q1); G2, _ = deltas(exp, cr, q2)
            nG1 = nrank(G1); nG2 = nrank(G2)
            a = np.einsum("gi,gi->i", nG1[:, boot_idx], nG2[:, boot_idx])
            jj, kk = pares_distintos(boot_idx)
            Mg = nG1.T @ nG1
            bm = Mg[jj, kk]
            bc = MR[jj, kk]
            boot["A_med"].append(float(np.median(a)))
            boot["A_mean"].append(float(np.mean(a)))
            boot["Bm_med"].append(float(np.median(bm)))
            boot["Bm_mean"].append(float(np.mean(bm)))
            boot["Bc_med"].append(float(np.median(bc)))
            boot["Bc_mean"].append(float(np.mean(bc)))
            boot["d_med"].append(float(np.median(bm) - np.median(a)))
            boot["d_mean"].append(float(np.mean(bm) - np.mean(a)))
            boot["dc_med"].append(float(np.median(bc) - np.median(a)))
            boot["dc_mean"].append(float(np.mean(bc) - np.mean(a)))
            boot["n_distintos"].append(int(len(np.unique(boot_idx))))
            boot["n_pares"].append(int(len(jj)))

        salida["lineas"][line] = {
            "fold_idx": fold_idx,
            "n_genes_expresados": ne,
            "n_genes_solo_verdad": ne_v,
            "n_control_cells": meta["n_control_cells"],
            "n_celulas_por_blanco": [int(len(i)) for i in idx_por_blanco],
            "n_cero_forzado": {"completa": nofin_R, "mitad_1": nofin_H1,
                               "mitad_2": nofin_H2},
            "A_mitades": A_pt, "B_completa": B_comp, "B_mitades": B_mit,
            "verificacion_candidato_01": verif,
            "sens_mascara_solo_verdad": sens_v,
            "sens_control_independiente": sens_c,
            "bootstrap": boot,
            "segundos": round(time.time() - t0, 1),
        }
        bm_med = float(np.median(B_mit)); a_med = float(np.median(A_pt))
        d = np.array(boot["d_med"])
        lo, hi = (np.percentile(d, [2.5, 97.5]) if len(d) else (np.nan, np.nan))
        sys.stderr.write(
            "%-6s ne=%5d A=%+.4f Bc=%+.4f Bm=%+.4f  d=%+.4f IC95=[%+.4f,%+.4f] %s  %.0fs\n"
            % (line, ne, a_med, float(np.median(B_comp)), bm_med,
               bm_med - a_med, lo, hi,
               "EXCLUYE_0" if (lo > 0 or hi < 0) else "incluye_0",
               time.time() - t0))
        sys.stderr.flush()
        with open(OUT if len(lineas) == len(LINES)
                  else OUT.replace(".json", "_%s.json" % "_".join(lineas)),
                  "w") as fh:
            json.dump(salida, fh, indent=1)
    sys.stderr.write("escrito\n")


if __name__ == "__main__":
    main()
