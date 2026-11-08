#!/usr/bin/env python
"""Rasgos de SOLO CONTROL para el par (blanco g, gen h), y el objetivo: el delta
ESPECIFICO del blanco en el lado real.

🔴 DIAGNOSTICO TRANSDUCTIVO. Lo que se ajusta aqui se ajusta SOBRE LA VERDAD de
Jiang, que es exactamente lo que un envio no tiene. Ninguna cifra, ningun
coeficiente y ningun modelo de aqui entra nunca a una submission. Lo unico que
sale es la RESPUESTA: si el techo de "solo control" es alto o es cero.

LA PREGUNTA. `b_contra_a.md` §5 midio que el delta especifico del blanco es el
61-84 % de la senal reproducible -- la componente MAYORITARIA en 6 de 6 -- y que
es reproducible. El candidato 01 no capturo nada de eso con covarianza del
control. Se discrimina entre:
  (a) la covarianza era la senal equivocada y OTRA senal del control si lo capta
  (b) NADA del control lo predice

EL OBJETIVO, igual que en `b_contra_a_descomp.py` para que sea la misma cantidad:
  H1 = log2(pb(mitad_1 de g) / c1)      c1, c2 = las dos mitades DISJUNTAS
  H2 = log2(pb(mitad_2 de g) / c2)           del pool de control
  centrar cada columna sobre los genes
  m1 = media sobre los 33 blancos de H1 ;  m2 = idem H2
  E1 = H1 - m1 ;  E2 = H2 - m2          <- el delta ESPECIFICO, dos veces
El control va PARTIDO por §4 de `b_contra_a.md`: con control comun el ruido de
muestreo del control se cuenta como senal compartida y `m` se lo come.

LOS RASGOS, todos del POOL DE CONTROL y del SIMBOLO del gen. Ninguno toca el lado
perturbado. Declarados aqui ANTES de medirlos:

  del par (g,h) -- los unicos que pueden dar correlacion media positiva (ver abajo)
   1 cov              covarianza de log1p normalizado en el pool de control
   2 r_pearson        correlacion de Pearson
   3 r_spearman       correlacion de Spearman (rangos sobre celulas)
   4 r_parcial_10     Pearson tras proyectar fuera los 10 primeros PC del control
   5 r_parcial_50     idem con 50 PC
   6 abs_r_pearson    |Pearson| -- la magnitud, separada del signo
   7 abs_r_spearman   |Spearman|
   8 delta_r          r_pearson - r_parcial_50: cuanto de la correlacion marginal
                      es programa global y no relacion directa
   9 rango_cor        rango de h entre los genes expresados por |r_pearson| con g,
                      normalizado a [0,1]
  10 mismo_modulo     1 si h cae en el mismo modulo de coexpresion que g
  11 dist_emb         distancia euclidea g-h en el embedding de genes del control
  12 cos_emb          similitud coseno g-h en ese embedding
  13 misma_familia    1 si los simbolos comparten raiz alfabetica de >= 3 letras
                      (STAT1/STAT4 -> STAT). Del SIMBOLO, no de ninguna base
                      externa: importar una seria una cifra importada
  14 es_el_blanco     1 si h ES g

  solo de h -- constantes dentro de una columna g
  15 h_media   16 h_var   17 h_disp (var/media)   18 h_frac_cero   19 h_rango_expr

  solo de g -- constantes dentro de una columna g
  20 g_media   21 g_var   22 g_frac_cero

🔴 POR QUE LOS RASGOS DE SOLO-h NO PUEDEN INFLAR EL VEREDICTO, y es una propiedad
del objetivo, no una esperanza. `E` esta centrado sobre los 33 blancos por
construccion, asi que para cada gen h la suma sobre los blancos es cero:
`Σ_g E[h,g] = 0`. Una prediccion que solo dependa de h es el mismo vector `p`
para todos los blancos, y entonces `Σ_g <p, E[:,g]> = <p, Σ_g E[:,g]> = 0`: la
correlacion MEDIA sobre blancos de cualquier predictor de solo-h es exactamente
cero. Puede acertar en un blanco y fallar en otro, nunca de media. Van incluidos
porque pueden modular a los del par (un gen de varianza alta responde mas), no
porque puedan puntuar solos. Los de solo-g son constantes dentro de la columna y
no cambian una correlacion por columna en absoluto.

LO QUE SE DESCARTO POR NO SER CALCULABLE SIN EL LADO PERTURBADO, dicho con nombre:
  * pertenencia a una via conocida (IFN-alfa/beta, JAK-STAT) desde una base de
    datos externa: no es del pool de control, seria una cifra importada (PLAN
    §0.2). No se usa.
  * respuesta conocida del gen a interferon: es el lado perturbado. No se usa.
  * efecto medido del mismo blanco en OTRA linea: es lado perturbado. No se usa.
  * distancia en el genoma / cromosoma: no esta en el h5ad (solo `var/_index`),
    y sacarla de una anotacion externa seria importar. No se usa.

SEMILLAS: familia nueva 20268000 + 100*fold (PCA/KMeans). La particion en mitades
del blanco conserva la 777 y el control partido la 20267999 + fold, las mismas de
`b_contra_a_descomp.py`, para que `E1`/`E2` sean exactamente los de alli.

Lectura y medicion. No construye ningun candidato, no toca ningun arnes.
"""
import json
import os
import re
import sys
import time

import numpy as np
import scipy.sparse as sp
from scipy.stats import rankdata
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from suelo_arnes import (IFNB_PATH, build_fold_adatas, load_global_obs,
                         load_var_names)

LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
PANEL_PATH = "/home/ubuntu/vcc/modelo/out/suelo_panel.json"
OUTDIR = "/home/ubuntu/vcc/modelo/out/control_o_no"
SEED_BASE_REAL = 20260001
SEED_MITADES = 777
SEED_CTRL_SPLIT = 20267999
SEED_FAMILIA = 20268000
UMBRAL_EXPR = 1.0
N_PC = 50
K_PARCIAL = (10, 50)
N_MODULOS = 30

RASGOS = ["cov", "r_pearson", "r_spearman", "r_parcial_10", "r_parcial_50",
          "abs_r_pearson", "abs_r_spearman", "delta_r", "rango_cor",
          "mismo_modulo", "dist_emb", "cos_emb", "misma_familia",
          "es_el_blanco", "h_media", "h_var", "h_disp", "h_frac_cero",
          "h_rango_expr", "g_media", "g_var", "g_frac_cero"]
RASGOS_PAR = RASGOS[:14]
RASGOS_SOLO_H = RASGOS[14:19]
RASGOS_SOLO_G = RASGOS[19:]


def raiz(sym):
    """Raiz alfabetica inicial del simbolo, >= 3 letras o None."""
    m = re.match(r"^[A-Za-z]+", sym or "")
    r = m.group(0).upper() if m else ""
    return r if len(r) >= 3 else None


def trivial_predict_fn(n_per_label=400):
    """`build_fold_adatas` exige un predict_fn. Este devuelve el propio control
    sin tocarlo: el lado PREDICHO no se usa para nada en este encargo (la
    mascara de genes es de solo-verdad, `cr >= 1`), asi que no se gasta la
    emision del candidato. El arnes no se modifica -- es su parametro."""
    def pf(control_pool_X, labels, rng):
        n = control_pool_X.shape[0]
        return {l: control_pool_X[rng.integers(0, n, size=n_per_label)].tocsr()
                for l in labels}
    return pf


def main():
    lineas = sys.argv[1].split(",") if len(sys.argv) > 1 else LINES
    os.makedirs(OUTDIR, exist_ok=True)
    panel = json.load(open(PANEL_PATH))["panel_comun"]
    nb = len(panel)
    cell_type, gene = load_global_obs(IFNB_PATH)
    var_names = load_var_names(IFNB_PATH)
    n_genes = len(var_names)
    idx_de = {g: i for i, g in enumerate(var_names)}

    for line in lineas:
        t0 = time.time()
        fold = LINES.index(line)
        real, _, meta = build_fold_adatas(
            line, panel, trivial_predict_fn(), var_names,
            SEED_BASE_REAL + fold, 0, cell_type, gene, n_genes, IFNB_PATH)
        Xr = real.X.tocsr(); lr = np.array(real.obs["target_gene"])
        ctrl_idx = np.where(lr == "NT")[0]
        Xc = Xr[ctrl_idx]

        def pb(M, idx):
            s = np.asarray(M[idx].sum(axis=0)).ravel(); t = s.sum()
            return s / t * 5e4 if t > 0 else s

        cr = pb(Xr, ctrl_idx)
        exp_mask = cr >= UMBRAL_EXPR            # mascara de SOLO VERDAD
        idx_exp = np.where(exp_mask)[0]
        tg_cols = np.array([idx_de[g] for g in panel])
        cols = np.union1d(idx_exp, tg_cols)
        pos = {int(c): j for j, c in enumerate(cols)}
        fila_de_exp = np.array([pos[int(c)] for c in idx_exp])
        col_de_tg = np.array([pos[int(c)] for c in tg_cols])
        G = len(idx_exp)

        # ---------- pool de control normalizado, igual que cand01_covar ----------
        Xcs = Xc.tocsc()[:, cols].tocsr().astype(np.float64)
        libs = np.asarray(Xc.sum(axis=1)).ravel()
        total = float(np.median(libs[libs > 0])) if (libs > 0).any() else 1.0
        esc = np.zeros(len(libs)); nzl = libs > 0
        esc[nzl] = total / libs[nzl]
        Y = (sp.diags(esc) @ Xcs).toarray()
        np.log1p(Y, out=Y)
        n_cel = Y.shape[0]
        frac_cero = 1.0 - (np.asarray((Xcs > 0).sum(axis=0)).ravel() / n_cel)

        ybar = Y.mean(axis=0)
        Yc = Y - ybar
        var = (Yc * Yc).sum(axis=0) / max(n_cel - 1, 1)
        sd = np.sqrt(var); sd_safe = np.where(sd > 0, sd, 1.0)
        Yn = Yc / sd_safe

        # ---------- PCA del control, embedding y modulos ----------
        pca = PCA(n_components=min(N_PC, n_cel - 1, Y.shape[1] - 1),
                  svd_solver="randomized", random_state=SEED_FAMILIA + 100 * fold)
        scores = pca.fit_transform(Yc)                      # (n_cel, k)
        comp = pca.components_                              # (k, n_cols)
        emb = (comp * np.sqrt(pca.explained_variance_)[:, None]).T   # (n_cols, k)
        nrm = np.linalg.norm(emb, axis=1); nrm_s = np.where(nrm > 0, nrm, 1.0)
        embn = emb / nrm_s[:, None]
        km = KMeans(n_clusters=N_MODULOS, n_init=10,
                    random_state=SEED_FAMILIA + 100 * fold + 1).fit(embn)
        modulo = km.labels_

        # ---------- correlaciones ----------
        Rp = Yn.T @ Yn[:, col_de_tg] / max(n_cel - 1, 1)     # (n_cols, 33) Pearson
        COV = Yc.T @ Yc[:, col_de_tg] / max(n_cel - 1, 1)    # covarianza
        Yr = rankdata(Y, axis=0).astype(np.float64)
        Yrc = Yr - Yr.mean(axis=0)
        sdr = np.sqrt((Yrc * Yrc).sum(axis=0)); sdr = np.where(sdr > 0, sdr, 1.0)
        Yrn = Yrc / sdr
        Rs = Yrn.T @ Yrn[:, col_de_tg]                       # Spearman
        del Yr, Yrc, Yrn

        Rpar = {}
        for k in K_PARCIAL:
            kk = min(k, scores.shape[1])
            S = scores[:, :kk]
            Sn = S / np.where(np.linalg.norm(S, axis=0) > 0,
                              np.linalg.norm(S, axis=0), 1.0)
            Yres = Yc - Sn @ (Sn.T @ Yc)
            s2 = np.sqrt((Yres * Yres).sum(axis=0))
            s2 = np.where(s2 > 0, s2, 1.0)
            Yresn = Yres / s2
            Rpar[k] = Yresn.T @ Yresn[:, col_de_tg]
            del Yres, Yresn
        del Yn, Yc, Y

        # ---------- el objetivo: delta ESPECIFICO con control partido ----------
        rng_c = np.random.default_rng(SEED_CTRL_SPLIT + fold)
        permc = rng_c.permutation(ctrl_idx); mc = len(permc) // 2
        c1 = pb(Xr, permc[:mc]); c2 = pb(Xr, permc[mc:])
        rng_m = np.random.default_rng(SEED_MITADES)
        p1, p2 = [], []
        for g in panel:
            ix = np.where(lr == g)[0]
            perm = rng_m.permutation(ix); m = len(perm) // 2
            p1.append(perm[:m]); p2.append(perm[m:])

        def deltas(cols_cel, denom):
            with np.errstate(divide="ignore", invalid="ignore"):
                M = np.column_stack([np.log2(pb(Xr, c)[exp_mask] / denom[exp_mask])
                                     for c in cols_cel])
            nf = int((~np.isfinite(M)).sum())
            M[~np.isfinite(M)] = 0.0
            return M, nf

        H1, nf1 = deltas(p1, c1); H2, nf2 = deltas(p2, c2)
        H1 = H1 - H1.mean(axis=0, keepdims=True)
        H2 = H2 - H2.mean(axis=0, keepdims=True)
        m1 = H1.mean(axis=1); m2 = H2.mean(axis=1)
        E1 = H1 - m1[:, None]; E2 = H2 - m2[:, None]
        E1 = E1 - E1.mean(axis=0, keepdims=True)
        E2 = E2 - E2.mean(axis=0, keepdims=True)

        # ---------- ensamblar los rasgos ----------
        raices = [raiz(var_names[int(c)]) for c in idx_exp]
        raices_tg = [raiz(g) for g in panel]
        F = np.zeros((G, nb, len(RASGOS)), dtype=np.float32)
        ri = {r: i for i, r in enumerate(RASGOS)}
        h_media = ybar[fila_de_exp]; h_var = var[fila_de_exp]
        h_fc = frac_cero[fila_de_exp]
        h_disp = h_var / np.where(h_media > 0, h_media, 1.0)
        h_rango = rankdata(h_media) / G
        mod_h = modulo[fila_de_exp]; emb_h = emb[fila_de_exp]
        for j in range(nb):
            cg = col_de_tg[j]
            rp = Rp[fila_de_exp, j]; rs = Rs[fila_de_exp, j]
            F[:, j, ri["cov"]] = COV[fila_de_exp, j]
            F[:, j, ri["r_pearson"]] = rp
            F[:, j, ri["r_spearman"]] = rs
            F[:, j, ri["r_parcial_10"]] = Rpar[10][fila_de_exp, j]
            F[:, j, ri["r_parcial_50"]] = Rpar[50][fila_de_exp, j]
            F[:, j, ri["abs_r_pearson"]] = np.abs(rp)
            F[:, j, ri["abs_r_spearman"]] = np.abs(rs)
            F[:, j, ri["delta_r"]] = rp - Rpar[50][fila_de_exp, j]
            F[:, j, ri["rango_cor"]] = rankdata(np.abs(rp)) / G
            F[:, j, ri["mismo_modulo"]] = (mod_h == modulo[cg]).astype(np.float32)
            d = emb_h - emb[cg]
            F[:, j, ri["dist_emb"]] = np.linalg.norm(d, axis=1)
            F[:, j, ri["cos_emb"]] = embn[fila_de_exp] @ embn[cg]
            rtg = raices_tg[j]
            F[:, j, ri["misma_familia"]] = np.array(
                [1.0 if (rtg is not None and r == rtg) else 0.0 for r in raices],
                dtype=np.float32)
            F[:, j, ri["es_el_blanco"]] = (idx_exp == tg_cols[j]).astype(np.float32)
            F[:, j, ri["h_media"]] = h_media
            F[:, j, ri["h_var"]] = h_var
            F[:, j, ri["h_disp"]] = h_disp
            F[:, j, ri["h_frac_cero"]] = h_fc
            F[:, j, ri["h_rango_expr"]] = h_rango
            F[:, j, ri["g_media"]] = ybar[cg]
            F[:, j, ri["g_var"]] = var[cg]
            F[:, j, ri["g_frac_cero"]] = frac_cero[cg]

        # fiabilidad del delta ESPECIFICO: el techo de esta pregunta
        num = np.einsum("gi,gi->i", E1, E2)
        A_esp = num / (np.linalg.norm(E1, axis=0) * np.linalg.norm(E2, axis=0))
        np.savez_compressed(
            os.path.join(OUTDIR, "rasgos_%s.npz" % line),
            F=F, E1=E1.astype(np.float32), E2=E2.astype(np.float32),
            A_esp=A_esp, cruz=num, rasgos=np.array(RASGOS),
            panel=np.array(panel), idx_exp=idx_exp,
            genes=np.array([var_names[int(c)] for c in idx_exp]),
            blanco_expresado=np.isin(tg_cols, idx_exp),
            n_no_finitos=np.array([nf1, nf2]),
            n_control=np.array([mc, len(permc) - mc]),
            n_cel_control=n_cel, G=G,
            var_explicada=pca.explained_variance_ratio_,
            tam_modulos=np.bincount(modulo, minlength=N_MODULOS))
        sys.stderr.write(
            "%-6s G=%5d  A_esp med=%+.4f  cruz>0 en %d/%d  blanco expresado %d/%d  "
            "nofin=%d/%d  PC1-10=%.1f%%  %.0fs\n"
            % (line, G, float(np.median(A_esp)), int((num > 0).sum()), nb,
               int(np.isin(tg_cols, idx_exp).sum()), nb, nf1, nf2,
               100 * pca.explained_variance_ratio_[:10].sum(), time.time() - t0))
        sys.stderr.flush()


if __name__ == "__main__":
    main()
