#!/usr/bin/env python
"""PASO 1 — APRENDER LA FUNCION sobre VCC 2025. Declarado ANTES de ajustar.

EL OBJETIVO
  y(g,h) = log2( (cpm_pert(g,h) + eps) / (cpm_ctrl(h) + eps) ),  eps = 1.0 CPM
  cpm = media sobre celulas de x/L*1e6. Es la razon que el emisor aplica, no un
  estadistico de DE: `_reescala_entero` multiplica cuentas por gen.

EL MODELO, DECLARADO ANTES DE AJUSTAR
  REGRESION LINEAL CON PENALIZACION L2 (ridge) sobre los 10 rasgos de
  `datos.RASGOS`, estandarizados (media 0, sd 1 sobre el ajuste), con
  intercepto. Sin interacciones, sin terminos no lineales, sin seleccion de
  rasgos. `alpha` de {0.1, 1, 10, 100, 1000} elegida por la validacion cruzada
  PARTIDA POR BLANCO. Es el modelo mas simple que admite esos rasgos.
  🔴 Ningun rasgo usa el gen como llave: ni el nombre del blanco, ni el del gen,
  ni la verdad de 2025 por gen. Los 28 blancos compartidos con el panel 2026
  pasan por la MISMA funcion que los otros 272.
  🔴 Los rasgos de CONTEXTO (profundidad mediana, n de control) NO entran en el
  ajuste: 2025 es UN solo contexto y su coeficiente no es identificable. Se
  declara en vez de fingirlo.

EL AJUSTE
  Filas = (blanco, gen) con el gen pasando el filtro del preset
  (CPM medio del control de 2025 > 5). Los otros genes se predicen igual, pero
  no entran al ajuste: no los mira ninguna de las seis metricas.

LA VALIDACION CRUZADA
  5 pliegues PARTIDOS POR BLANCO (ningun blanco aparece en ajuste y prueba).
  Permutacion de los 300 blancos con la semilla 20284950. Se reporta, por
  pliegue: correlacion de Pearson y Spearman entre y predicha y observada sobre
  los blancos retenidos, R^2, y la DISTRIBUCION de la correlacion POR BLANCO.
  Dos controles declarados:
    (a) predecir 0 en todo gen (el suelo: ninguna respuesta)
    (b) solo el perfil de gwps, con una escala ajustada (una sola pendiente)
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, "/home/ubuntu/vcc/modelo/src/noche")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import comun as C  # noqa: E402
import datos as D  # noqa: E402
import sem2 as S  # noqa: E402

OUT = "/home/ubuntu/vcc/modelo/out/entrega_ajuste.json"
EPS_CPM = 1.0
ALPHAS = [0.1, 1.0, 10.0, 100.0, 1000.0]
N_PLIEGUES = 5


def ridge(X, y, alpha):
    """Normal equations sobre rasgos estandarizados, con intercepto aparte."""
    mu = X.mean(0)
    sd = X.std(0)
    sd[sd <= 0] = 1.0
    Z = (X - mu) / sd
    n, p = Z.shape
    A = Z.T @ Z + alpha * np.eye(p)
    b = Z.T @ (y - y.mean())
    w = np.linalg.solve(A, b)
    return {"w": w, "b0": float(y.mean()), "mu": mu, "sd": sd}


def predice(m, X):
    return ((X - m["mu"]) / m["sd"]) @ m["w"] + m["b0"]


def spearman(a, b):
    def rk(v):
        o = np.argsort(v, kind="stable")
        r = np.empty(v.size)
        r[o] = np.arange(v.size)
        return r
    ra, rb = rk(a), rk(b)
    ra -= ra.mean(); rb -= rb.mean()
    d = np.sqrt((ra @ ra) * (rb @ rb))
    return float(ra @ rb / d) if d > 0 else None


def pearson(a, b):
    a = a - a.mean(); b = b - b.mean()
    d = np.sqrt((a @ a) * (b @ b))
    return float(a @ b / d) if d > 0 else None


def main():
    t0 = time.time()
    sem = S.comprueba()
    assert sem["sin_choques"], sem
    res = {"semillas": sem, "eps_cpm": EPS_CPM, "alphas": ALPHAS,
           "rasgos": D.RASGOS, "modelo": "ridge lineal, 10 rasgos estandarizados, "
           "intercepto, sin interacciones; alpha por CV partida por blanco",
           "declarado_antes_de_ajustar": True}
    perts, genes25, cpm_pert, ncel, cpm_ctrl, det, fano, meta = D.pseudobulk_2025()
    res["pseudobulk_2025"] = meta
    res["n_perts_2025"] = int(len(perts))
    C.log(f"2025: {len(perts)} blancos, {len(genes25)} genes, "
          f"{meta['n_control_2025']:,} celulas de control")

    y_full = np.log2((cpm_pert + EPS_CPM) / (cpm_ctrl[None, :] + EPS_CPM)).astype(np.float32)
    st25 = {"cpm": cpm_ctrl, "det": det, "fano": np.where(fano < 0, 0.0, fano),
            "n": meta["n_control_2025"]}
    filtr = np.where(cpm_ctrl > D.UMBRAL_CPM)[0]
    res["n_genes_2025"] = int(len(genes25))
    res["n_genes_filtro_2025"] = int(len(filtr))
    tf = D.tf_set()
    gw_vn, gw_perf, gw_meta = D.gwps(sorted(set(perts.tolist()) | set(D.perts_2026().tolist())))
    res["gwps"] = gw_meta
    res["blancos_2025_con_gwps"] = int(sum(1 for t in perts if t in gw_perf))
    res["genes_2025_con_gwps"] = int(len(set(genes25.tolist()) & set(gw_vn.tolist())))

    # ── matriz de rasgos de 2025, SOLO sobre los genes del filtro
    g_fit = genes25[filtr]
    st_fit = {k: (v[filtr] if isinstance(v, np.ndarray) else v) for k, v in st25.items()}
    X = D.matriz_rasgos(perts, g_fit, gw_vn, gw_perf, tf, st_fit)
    Y = y_full[:, filtr].reshape(-1)
    G = len(g_fit)
    C.log(f"diseno: {X.shape[0]:,} filas x {X.shape[1]} rasgos")
    res["forma_diseno"] = list(X.shape)

    # ── CV partida por blanco
    rng = np.random.default_rng(S.S_BARAJA)
    orden = rng.permutation(len(perts))
    pliegues = np.array_split(orden, N_PLIEGUES)
    res["cv"] = {"semilla_baraja": S.S_BARAJA, "n_pliegues": N_PLIEGUES,
                 "semillas_reservadas_no_consumidas": S.SEMILLAS["cv_pliegues"],
                 "blancos_por_pliegue": [len(p) for p in pliegues],
                 "por_alpha": {}, "blancos_de_cada_pliegue":
                 [sorted(perts[p].tolist()) for p in pliegues]}
    filas_t = {i: np.arange(i * G, (i + 1) * G) for i in range(len(perts))}

    def idx_de(ts):
        return np.concatenate([filas_t[i] for i in ts])

    mejor = None
    for a in ALPHAS:
        pl = []
        for k, pr in enumerate(pliegues):
            tr = np.concatenate([p for j, p in enumerate(pliegues) if j != k])
            itr, ite = idx_de(tr), idx_de(pr)
            m = ridge(X[itr].astype(np.float64), Y[itr].astype(np.float64), a)
            yp = predice(m, X[ite].astype(np.float64))
            yo = Y[ite].astype(np.float64)
            porb = []
            for n_, t in enumerate(pr):
                s = slice(n_ * G, (n_ + 1) * G)
                porb.append(pearson(yp[s], yo[s]))
            pl.append({"pliegue": k, "n_blancos": len(pr),
                       "pearson": pearson(yp, yo), "spearman": spearman(yp, yo),
                       "R2": float(1 - ((yo - yp) ** 2).sum() / ((yo - yo.mean()) ** 2).sum()),
                       "rmse": float(np.sqrt(((yo - yp) ** 2).mean())),
                       "pearson_por_blanco_mediana": float(np.median([x for x in porb if x is not None])),
                       "pearson_por_blanco_rango": [float(np.min([x for x in porb if x is not None])),
                                                    float(np.max([x for x in porb if x is not None]))],
                       "pearson_por_blanco": [None if x is None else round(x, 4) for x in porb]})
        res["cv"]["por_alpha"][str(a)] = {
            "pliegues": pl,
            "pearson_medio": float(np.mean([p["pearson"] for p in pl])),
            "R2_medio": float(np.mean([p["R2"] for p in pl])),
            "pearson_por_blanco_mediana": float(np.mean([p["pearson_por_blanco_mediana"] for p in pl]))}
        C.log(f"alpha={a}: pearson {res['cv']['por_alpha'][str(a)]['pearson_medio']:+.4f} "
              f"R2 {res['cv']['por_alpha'][str(a)]['R2_medio']:+.4f} "
              f"por blanco {res['cv']['por_alpha'][str(a)]['pearson_por_blanco_mediana']:+.4f}")
        C.escribe_json(OUT, res)
        v = res["cv"]["por_alpha"][str(a)]["pearson_por_blanco_mediana"]
        if mejor is None or v > mejor[1]:
            mejor = (a, v)
    res["alpha_elegida"] = mejor[0]
    res["criterio_alpha"] = "la mediana de la correlacion de Pearson POR BLANCO retenido"

    # ── controles declarados
    ctr = {}
    yo_all, yp0 = Y.astype(np.float64), np.zeros(Y.size)
    ctr["a_predecir_cero"] = {"R2": float(1 - ((yo_all - yp0) ** 2).sum()
                                          / ((yo_all - yo_all.mean()) ** 2).sum()),
                              "rmse": float(np.sqrt((yo_all ** 2).mean())),
                              "pearson": None}
    sol = []
    for k, pr in enumerate(pliegues):
        tr = np.concatenate([p for j, p in enumerate(pliegues) if j != k])
        itr, ite = idx_de(tr), idx_de(pr)
        x1 = X[itr][:, [0]].astype(np.float64)
        m = ridge(x1, Y[itr].astype(np.float64), 1.0)
        yp = predice(m, X[ite][:, [0]].astype(np.float64))
        yo = Y[ite].astype(np.float64)
        sol.append({"pliegue": k, "pearson": pearson(yp, yo),
                    "R2": float(1 - ((yo - yp) ** 2).sum() / ((yo - yo.mean()) ** 2).sum())})
    ctr["b_solo_perfil_gwps"] = {"pliegues": sol,
                                 "pearson_medio": float(np.mean([s["pearson"] for s in sol])),
                                 "R2_medio": float(np.mean([s["R2"] for s in sol]))}
    res["controles"] = ctr

    # ── ajuste final con TODOS los blancos
    m = ridge(X.astype(np.float64), Y.astype(np.float64), mejor[0])
    res["coeficientes"] = {"alpha": mejor[0], "b0": m["b0"],
                           "w": {r: float(v) for r, v in zip(D.RASGOS, m["w"])},
                           "mu": {r: float(v) for r, v in zip(D.RASGOS, m["mu"])},
                           "sd": {r: float(v) for r, v in zip(D.RASGOS, m["sd"])}}
    yp = predice(m, X.astype(np.float64))
    res["ajuste_final_en_muestra"] = {
        "pearson": pearson(yp, Y.astype(np.float64)),
        "R2": float(1 - ((Y - yp) ** 2).sum() / ((Y - Y.mean()) ** 2).sum()),
        "y_observada": {"sd": float(Y.std()), "p1": float(np.percentile(Y, 1)),
                        "p50": float(np.percentile(Y, 50)), "p99": float(np.percentile(Y, 99))},
        "y_predicha": {"sd": float(yp.std()), "p1": float(np.percentile(yp, 1)),
                       "p50": float(np.percentile(yp, 50)), "p99": float(np.percentile(yp, 99))}}
    res["dt_total_s"] = time.time() - t0
    C.escribe_json(OUT, res)
    C.log(f"ajuste listo: alpha={mejor[0]} pearson_por_blanco={mejor[1]:+.4f} "
          f"en {res['dt_total_s']/60:.1f} min")


if __name__ == "__main__":
    main()
