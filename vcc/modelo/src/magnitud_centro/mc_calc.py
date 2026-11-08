#!/usr/bin/env python
"""PASO 2 / parte 2: las cantidades, de los deltas de `mc_deltas.py`.

LAS CANTIDADES, declaradas antes de medir
-----------------------------------------
Todo delta entra ya CENTRADO sobre los genes (lo hizo `mc_deltas.py`), asi que
`<x,y>/G` es una covarianza y `corr` es Pearson.

TRAMPA 1 -- LA ATENUACION
  fiab[L,g]    = corr(delta_g,L,A , delta_g,L,B)          fiabilidad DENTRO de linea
  cruda[L1,L2,g] = corr(delta_g,L1,A , delta_g,L2,B)      MITADES CRUZADAS
  desat        = cruda / sqrt(fiab[L1,g] * fiab[L2,g])    DESATENUADA
  🔴 La cruzada es simetrizada: media de corr(A1,B2) y corr(B1,A2). Las dos usan
     mitades de celulas disjuntas y mitades de control disjuntas, asi que el
     ruido no se comparte en ninguna de las dos.
  🔴 Si la desatenuada pasa de 1 se REPORTA tal cual, no se recorta.
  🔴 Si fiab <= 0 la raiz no existe: la celda sale `n/d` y se cuentan.

TRAMPA 2 -- EL EJE COMPARTIDO (G_mismo / G_otro)
  G_mismo[L1,L2,g] = MISMO blanco, linea distinta   -> el centro que se busca
  G_otro[L1,L2,g]  = OTRO blanco, linea distinta    -> la respuesta comun al IFNB
                     (media sobre los 32 `q != g`; se dan sd, min, max)
  diferencia       = G_mismo - G_otro               -> la parte especifica del
                     blanco que cruza contextos
  Las tres columnas, las tres desatenuadas con su propio par de fiabilidades.

LA FRACCION CENTRO (identidad de `b_contra_a.md` §5, con el indice sobre LINEAS)
  Para el blanco g y un conjunto de lineas S, con h1_L = delta_g,L,A y
  h2_L = delta_g,L,B:
      V_tot(g) = (1/|S|) Σ_L <h1_L, h2_L> / G
      V_sh(g)  = <m1, m2> / G            con m1 = mean_L h1_L, m2 = mean_L h2_L
      V_sp(g)  = (1/|S|) Σ_L <h1_L-m1, h2_L-m2> / G
  `V_tot = V_sh + V_sp` es una identidad algebraica exacta (se reporta el
  residuo). El producto CRUZADO de mitades es lo que quita el ruido sin
  modelarlo.
  🔴 DOS VARIANTES, porque la TRAMPA 2 tambien afecta aqui:
     (a) sobre el delta CRUDO -> V_sh mezcla el centro con la respuesta comun al
         IFNB, que es compartida entre lineas Y entre blancos.
     (b) sobre el delta ESPECIFICO del blanco `E_g,L = delta_g,L - mean_q delta_q,L`,
         restado DENTRO de cada linea y DENTRO de cada mitad -> la respuesta
         comun al IFNB sale del vector antes de descomponer.
  IC por bootstrap sobre PERTURBACIONES (2000 replicas, percentil 95 %).

SIGNO Y MAGNITUD, POR SEPARADO (nunca una resumiendo a la otra)
  signo:    concordancia de signo entre contextos (azar = 0.5) y Spearman
  magnitud: pendiente de la regresion entre contextos y su R^2
            pendiente = <dA, dB> / <dA, dA>;  desatenuada = pendiente / fiab[L1,g]
            (la atenuacion de una pendiente es la fiabilidad del PREDICTOR solo,
             no la raiz del producto)
"""
import json, os, sys
import numpy as np

sys.path.insert(0, "/home/ubuntu/vcc/modelo/src")
import gwps_via as GV

LINES = GV.LINES
WORK = "/datos/modelo/out/magnitud_centro"
OUTDIR = "/home/ubuntu/vcc/modelo/out/magnitud_centro"
SEED_BOOT = 20278800          # familia nueva, fuera de las 48 de particion
N_BOOT = 2000
EPS = 1e-300


def unit(M):
    """Filas normalizadas a norma 1 (ya vienen centradas)."""
    n = np.linalg.norm(M, axis=-1, keepdims=True)
    return M / np.maximum(n, EPS), n[..., 0]


def rankea(v):
    """Rangos con empates promediados."""
    o = np.argsort(v, kind="stable")
    r = np.empty(v.size, dtype=np.float64)
    r[o] = np.arange(1, v.size + 1, dtype=np.float64)
    # promediar empates
    vs = v[o]
    i = 0
    while i < vs.size:
        k = i
        while k + 1 < vs.size and vs[k + 1] == vs[i]:
            k += 1
        if k > i:
            r[o[i:k + 1]] = r[o[i:k + 1]].mean()
        i = k + 1
    return r


def main():
    meta = json.load(open(os.path.join(OUTDIR, "deltas_meta.json")))
    panel = meta["panel"]
    P = len(panel)
    D, G = {}, None
    for L in LINES:
        z = np.load(os.path.join(WORK, f"deltas_{L}.npz"), allow_pickle=True)
        D[L] = z["D"].astype(np.float64)       # (J, P, 2, G)
        assert list(z["panel"].astype(str)) == panel
        G = D[L].shape[-1]
    J = D[LINES[0]].shape[0]
    sys.stderr.write(f"cargado: {len(LINES)} lineas, J={J}, P={P}, G={G}\n")

    # ------------------------------------------------- normalizados y rangos
    NU, NORM, RK = {}, {}, {}
    for L in LINES:
        NU[L], NORM[L] = unit(D[L])
        R = np.empty_like(D[L])
        for j in range(J):
            for g in range(P):
                for h in range(2):
                    r = rankea(D[L][j, g, h])
                    R[j, g, h] = r - r.mean()
        RK[L], _ = unit(R)
        sys.stderr.write(f"  rangos {L} ok\n")

    pares = [(a, b) for i, a in enumerate(LINES) for b in LINES[i + 1:]]
    res = {"meta_espacio": meta["espacio"], "panel": panel,
           "n_particiones": J, "n_lineas": len(LINES), "n_pares": len(pares),
           "semilla_bootstrap_base": SEED_BOOT, "n_bootstrap": N_BOOT}

    # ------------------------------------------------- 4) fiabilidad dentro de linea
    fiab = np.zeros((len(LINES), J, P))
    for li, L in enumerate(LINES):
        fiab[li] = np.einsum("jpg,jpg->jp", NU[L][:, :, 0], NU[L][:, :, 1])
    res["fiabilidad"] = {
        "por_linea_y_pert": {L: {panel[g]: {
            "media": float(fiab[li, :, g].mean()),
            "sd": float(fiab[li, :, g].std(ddof=1)),
            "min": float(fiab[li, :, g].min()),
            "max": float(fiab[li, :, g].max()),
            "los_8": [float(x) for x in fiab[li, :, g]],
        } for g in range(P)} for li, L in enumerate(LINES)},
        "por_linea": {L: {
            "media": float(fiab[li].mean()), "sd": float(fiab[li].std(ddof=1)),
            "min": float(fiab[li].min()), "max": float(fiab[li].max()),
            "mediana": float(np.median(fiab[li])),
            "n_no_positivas": int((fiab[li] <= 0).sum()),
            "n_celdas": int(fiab[li].size),
        } for li, L in enumerate(LINES)},
        "global": {"media": float(fiab.mean()), "min": float(fiab.min()),
                   "max": float(fiab.max()),
                   "n_no_positivas": int((fiab <= 0).sum()),
                   "n_celdas": int(fiab.size)},
    }

    # ------------------------------------------------- 1) por par: mismo / otro
    filas_par, filas_pert = [], []
    nd_total = 0
    for (L1, L2) in pares:
        i1, i2 = LINES.index(L1), LINES.index(L2)
        mis_c = np.zeros((J, P)); otr_c = np.zeros((J, P))
        otr_sd = np.zeros((J, P)); otr_min = np.zeros((J, P)); otr_max = np.zeros((J, P))
        mis_d = np.full((J, P), np.nan); otr_d = np.full((J, P), np.nan)
        spm = np.zeros((J, P)); sgn = np.zeros((J, P))
        pend = np.zeros((J, P)); pend_d = np.full((J, P), np.nan); r2 = np.zeros((J, P))
        for j in range(J):
            M1 = NU[L1][j, :, 0] @ NU[L2][j, :, 1].T      # A1 x B2
            M2 = NU[L1][j, :, 1] @ NU[L2][j, :, 0].T      # B1 x A2
            M = 0.5 * (M1 + M2)
            mis_c[j] = np.diag(M)
            off = M.copy(); np.fill_diagonal(off, np.nan)
            otr_c[j] = np.nanmean(off, axis=1)
            otr_sd[j] = np.nanstd(off, axis=1, ddof=1)
            otr_min[j] = np.nanmin(off, axis=1)
            otr_max[j] = np.nanmax(off, axis=1)
            # desatenuadas
            f1, f2 = fiab[i1, j], fiab[i2, j]
            den = np.sqrt(f1 * f2)
            ok = (f1 > 0) & (f2 > 0)
            mis_d[j, ok] = mis_c[j, ok] / den[ok]
            nd_total += int((~ok).sum())
            # G_otro desatenuada: cada q con su propio par de fiabilidades
            den_m = np.sqrt(np.outer(np.maximum(f1, 0), np.maximum(f2, 0)))
            with np.errstate(divide="ignore", invalid="ignore"):
                Od = off / den_m
            Od[~np.isfinite(Od)] = np.nan
            otr_d[j] = np.nanmean(Od, axis=1)
            # Spearman (rangos, mitades cruzadas, simetrizado)
            S1 = RK[L1][j, :, 0] @ RK[L2][j, :, 1].T
            S2 = RK[L1][j, :, 1] @ RK[L2][j, :, 0].T
            spm[j] = np.diag(0.5 * (S1 + S2))
            # concordancia de signo (simetrizada)
            for g in range(P):
                a1, b2 = D[L1][j, g, 0], D[L2][j, g, 1]
                b1, a2 = D[L1][j, g, 1], D[L2][j, g, 0]
                s = 0.5 * (np.mean(np.sign(a1) == np.sign(b2))
                           + np.mean(np.sign(b1) == np.sign(a2)))
                sgn[j, g] = s
                # pendiente y R^2 (simetrizados sobre las dos direcciones)
                p1 = float(a1 @ b2) / float(a1 @ a1)
                p2 = float(b1 @ a2) / float(b1 @ b1)
                pend[j, g] = 0.5 * (p1 + p2)
                r2[j, g] = mis_c[j, g] ** 2
                if f1[g] > 0:
                    pend_d[j, g] = pend[j, g] / f1[g]
        dif_c = mis_c - otr_c
        dif_d = mis_d - otr_d

        def agg(A):
            v = A[np.isfinite(A)]
            return {"media": float(v.mean()) if v.size else None,
                    "sd": float(v.std(ddof=1)) if v.size > 1 else None,
                    "min": float(v.min()) if v.size else None,
                    "max": float(v.max()) if v.size else None,
                    "n_finitas": int(v.size), "n_celdas": int(A.size)}

        filas_par.append({
            "par": f"{L1}-{L2}", "L1": L1, "L2": L2,
            "mismo_cruda": agg(mis_c), "otro_cruda": agg(otr_c),
            "dif_cruda": agg(dif_c),
            "mismo_desat": agg(mis_d), "otro_desat": agg(otr_d),
            "dif_desat": agg(dif_d),
            "fiab_L1": agg(fiab[i1]), "fiab_L2": agg(fiab[i2]),
            "otro_sd_media": float(otr_sd.mean()),
            "otro_min_global": float(otr_min.min()),
            "otro_max_global": float(otr_max.max()),
            "spearman": agg(spm), "concordancia_signo": agg(sgn),
            "pendiente": agg(pend), "pendiente_desat": agg(pend_d), "R2": agg(r2),
            "n_desat_nd": int((~np.isfinite(mis_d)).sum()),
            "n_desat_mayor_1": int(np.sum(mis_d[np.isfinite(mis_d)] > 1.0)),
        })
        for g in range(P):
            filas_pert.append({
                "par": f"{L1}-{L2}", "pert": panel[g],
                "mismo_cruda": float(mis_c[:, g].mean()),
                "mismo_cruda_sd": float(mis_c[:, g].std(ddof=1)),
                "fiab_L1": float(fiab[i1, :, g].mean()),
                "fiab_L2": float(fiab[i2, :, g].mean()),
                "mismo_desat": (float(np.nanmean(mis_d[:, g]))
                                if np.isfinite(mis_d[:, g]).any() else None),
                "otro_cruda": float(otr_c[:, g].mean()),
                "otro_desat": (float(np.nanmean(otr_d[:, g]))
                               if np.isfinite(otr_d[:, g]).any() else None),
                "dif_cruda": float(dif_c[:, g].mean()),
                "dif_desat": (float(np.nanmean(dif_d[:, g]))
                              if np.isfinite(dif_d[:, g]).any() else None),
                "spearman": float(spm[:, g].mean()),
                "concordancia_signo": float(sgn[:, g].mean()),
                "pendiente": float(pend[:, g].mean()),
                "pendiente_desat": (float(np.nanmean(pend_d[:, g]))
                                    if np.isfinite(pend_d[:, g]).any() else None),
                "R2": float(r2[:, g].mean()),
            })
        sys.stderr.write(
            f"  par {L1}-{L2}: mismo={mis_c.mean():+.4f} otro={otr_c.mean():+.4f} "
            f"dif={dif_c.mean():+.4f} desat={np.nanmean(mis_d):+.4f} "
            f"nd={int((~np.isfinite(mis_d)).sum())}\n")
    res["por_par"] = filas_par
    res["por_par_y_pert"] = filas_pert
    res["n_desat_nd_total"] = nd_total

    # ------------------------------------------------- 2) fraccion CENTRO
    def descompone(Hs, idx_lineas):
        """Hs: dict L -> (J,P,2,G). Devuelve V_tot,V_sh,V_sp por (J,P)."""
        nS = len(idx_lineas)
        Vt = np.zeros((J, P)); Vs = np.zeros((J, P)); Vp = np.zeros((J, P))
        for j in range(J):
            h1 = np.stack([Hs[L][j, :, 0] for L in idx_lineas])   # (nS,P,G)
            h2 = np.stack([Hs[L][j, :, 1] for L in idx_lineas])
            m1, m2 = h1.mean(axis=0), h2.mean(axis=0)
            Vt[j] = np.einsum("spg,spg->p", h1, h2) / nS / G
            Vs[j] = np.einsum("pg,pg->p", m1, m2) / G
            Vp[j] = np.einsum("spg,spg->p", h1 - m1, h2 - m2) / nS / G
        return Vt, Vs, Vp

    # variante (b): delta ESPECIFICO del blanco, restado dentro de linea y mitad
    E = {}
    for L in LINES:
        E[L] = D[L] - D[L].mean(axis=1, keepdims=True)

    rng = np.random.default_rng(SEED_BOOT)
    def boot(Vt, Vs):
        """IC de mean_g(Vs)/mean_g(Vt), remuestreando PERTURBACIONES."""
        vt = Vt.mean(axis=0); vs = Vs.mean(axis=0)      # media sobre particiones
        fr = []
        for _ in range(N_BOOT):
            k = rng.integers(0, P, size=P)
            t = vt[k].mean()
            fr.append(vs[k].mean() / t if t != 0 else np.nan)
        fr = np.array(fr); fr = fr[np.isfinite(fr)]
        return {"punto": float(vs.mean() / vt.mean()) if vt.mean() else None,
                "ic95": [float(np.percentile(fr, 2.5)),
                         float(np.percentile(fr, 97.5))] if fr.size else None,
                "boot_media": float(fr.mean()) if fr.size else None,
                "boot_sd": float(fr.std(ddof=1)) if fr.size > 1 else None,
                "n_boot_validas": int(fr.size)}

    centro = {}
    for nombre, Hs in (("delta_crudo", D), ("delta_especifico_del_blanco", E)):
        agregada = {}
        for etiq, S in [("las_6_lineas", LINES)] + [(f"{a}-{b}", [a, b]) for a, b in pares]:
            Vt, Vs, Vp = descompone(Hs, S)
            resid = np.abs(Vt - (Vs + Vp))
            frac_j = np.where(Vt != 0, Vs / np.where(Vt == 0, 1, Vt), np.nan)
            agregada[etiq] = {
                "n_lineas": len(S),
                "V_tot_media": float(Vt.mean()), "V_sh_media": float(Vs.mean()),
                "V_sp_media": float(Vp.mean()),
                "frac_centro_agregada": float(Vs.mean() / Vt.mean()) if Vt.mean() else None,
                "frac_especifica_agregada": float(Vp.mean() / Vt.mean()) if Vt.mean() else None,
                "residuo_identidad_max": float(resid.max()),
                "residuo_identidad_rel_max": float((resid / np.maximum(np.abs(Vt), EPS)).max()),
                "frac_por_pert_media": float(np.nanmean(frac_j)),
                "frac_por_pert_min": float(np.nanmin(frac_j)),
                "frac_por_pert_max": float(np.nanmax(frac_j)),
                "n_V_tot_no_positivas": int((Vt <= 0).sum()),
                "bootstrap_sobre_perturbaciones": boot(Vt, Vs),
                "por_pert": [{"pert": panel[g],
                              "V_tot": float(Vt[:, g].mean()),
                              "V_sh": float(Vs[:, g].mean()),
                              "V_sp": float(Vp[:, g].mean()),
                              "frac_centro": (float(Vs[:, g].mean() / Vt[:, g].mean())
                                              if Vt[:, g].mean() else None),
                              "frac_centro_sd_sobre_particiones": float(np.nanstd(frac_j[:, g], ddof=1)),
                              } for g in range(P)],
            }
            sys.stderr.write(
                f"  centro[{nombre}][{etiq}]: V_tot={Vt.mean():.6f} "
                f"V_sh={Vs.mean():.6f} frac={agregada[etiq]['frac_centro_agregada']:.4f} "
                f"resid={resid.max():.2e}\n")
        centro[nombre] = agregada
    res["fraccion_centro"] = centro

    json.dump(res, open(os.path.join(OUTDIR, "calc.json"), "w"),
              indent=1, sort_keys=True, default=str)
    sys.stderr.write(f"\nlisto: {OUTDIR}/calc.json\n")


if __name__ == "__main__":
    main()
