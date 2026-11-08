#!/usr/bin/env python
"""¿El bloque de ceros estructurales contamina tambien el PEARSON (G_mismo), o
solo la concordancia de signo?

Se recalcula `G_mismo` y `G_otro` (mitades cruzadas, simetrizado) restringiendo
el espacio a los genes con delta NO EMPATADO en el bloque del cero. El criterio
es el propio delta del lado predictor: se descartan los genes cuyo |delta| cae
en el bloque modal (el valor mas repetido de |delta|, que es el del cero).
Se da ademas la version por decil superior de |delta|.
"""
import json, os, sys
import numpy as np

sys.path.insert(0, "/home/ubuntu/vcc/modelo/src")
import gwps_via as GV

LINES = GV.LINES
WORK = "/datos/modelo/out/magnitud_centro"
OUTDIR = "/home/ubuntu/vcc/modelo/out/magnitud_centro"
NUCLEO = ["IFNAR1", "IFNAR2", "JAK1", "STAT1", "TYK2"]

meta = json.load(open(os.path.join(OUTDIR, "deltas_meta.json")))
panel = meta["panel"]
D = {L: np.load(os.path.join(WORK, f"deltas_{L}.npz"),
                allow_pickle=True)["D"].astype(np.float64) for L in LINES}
J, P, _, G = D[LINES[0]].shape
pares = [(a, b) for i, a in enumerate(LINES) for b in LINES[i + 1:]]


def corr(x, y, m=None):
    if m is not None:
        x, y = x[m], y[m]
    x = x - x.mean(); y = y - y.mean()
    nx, ny = np.linalg.norm(x), np.linalg.norm(y)
    return float(x @ y / (nx * ny)) if nx > 0 and ny > 0 else np.nan


filas = []
for (L1, L2) in pares:
    for j in range(J):
        for g in range(P):
            a1, b2 = D[L1][j, g, 0], D[L2][j, g, 1]
            b1, a2 = D[L1][j, g, 1], D[L2][j, g, 0]
            q = int((g + 7) % P)                 # un "otro blanco" fijo y declarado
            ao, bo = D[L2][j, q, 1], D[L2][j, q, 0]
            # mascara: fuera del bloque modal de |delta| del predictor
            av = np.abs(a1)
            vals, cnt = np.unique(np.round(av, 12), return_counts=True)
            modal = vals[np.argmax(cnt)]
            m_nz = np.round(av, 12) != modal
            # decil superior de |delta|
            thr = np.quantile(av, 0.9)
            m_top = av >= thr
            filas.append({
                "par": f"{L1}-{L2}", "pert": panel[g], "j": j,
                "nucleo": panel[g] in NUCLEO,
                "n_total": G, "n_fuera_del_bloque_modal": int(m_nz.sum()),
                "n_decil_superior": int(m_top.sum()),
                "mismo_todos": 0.5 * (corr(a1, b2) + corr(b1, a2)),
                "mismo_sin_bloque": 0.5 * (corr(a1, b2, m_nz) + corr(b1, a2, m_nz)),
                "mismo_decil_sup": 0.5 * (corr(a1, b2, m_top) + corr(b1, a2, m_top)),
                "otro_todos": 0.5 * (corr(a1, ao) + corr(b1, bo)),
                "otro_sin_bloque": 0.5 * (corr(a1, ao, m_nz) + corr(b1, bo, m_nz)),
                "otro_decil_sup": 0.5 * (corr(a1, ao, m_top) + corr(b1, bo, m_top)),
            })

res = {"criterio_otro_blanco": "q = (g+7) mod 33, fijo y declarado (no promedia los 32)",
       "n_filas": len(filas), "subconjuntos": {}}
for nom, sel in (("todas", lambda r: True),
                 ("nucleo (5)", lambda r: r["nucleo"]),
                 ("resto (28)", lambda r: not r["nucleo"])):
    rs = [r for r in filas if sel(r)]
    d = {"n_filas": len(rs),
         "n_genes_medio_fuera_del_bloque": float(np.mean([r["n_fuera_del_bloque_modal"] for r in rs])),
         "frac_genes_en_el_bloque_modal": float(1 - np.mean([r["n_fuera_del_bloque_modal"] for r in rs]) / G),
         "n_genes_decil_superior": float(np.mean([r["n_decil_superior"] for r in rs]))}
    for k in ("mismo_todos", "mismo_sin_bloque", "mismo_decil_sup",
              "otro_todos", "otro_sin_bloque", "otro_decil_sup"):
        v = np.array([r[k] for r in rs], dtype=float); v = v[np.isfinite(v)]
        d[k] = {"media": float(v.mean()), "sd": float(v.std(ddof=1)),
                "n": int(v.size)}
    d["dif_todos"] = d["mismo_todos"]["media"] - d["otro_todos"]["media"]
    d["dif_sin_bloque"] = d["mismo_sin_bloque"]["media"] - d["otro_sin_bloque"]["media"]
    d["dif_decil_sup"] = d["mismo_decil_sup"]["media"] - d["otro_decil_sup"]["media"]
    res["subconjuntos"][nom] = d

json.dump(res, open(os.path.join(OUTDIR, "zeros.json"), "w"), indent=1, sort_keys=True)
h = "%-12s %8s %22s %10s %10s %10s %10s %10s %10s"
print(h % ("subconj", "filas", "frac en bloque modal",
           "mis_todos", "mis_s/bloq", "mis_dec10", "otr_todos", "otr_s/bloq", "otr_dec10"))
for nom, d in res["subconjuntos"].items():
    print(h % (nom, d["n_filas"], "%.4f" % d["frac_genes_en_el_bloque_modal"],
               "%+.4f" % d["mismo_todos"]["media"], "%+.4f" % d["mismo_sin_bloque"]["media"],
               "%+.4f" % d["mismo_decil_sup"]["media"], "%+.4f" % d["otro_todos"]["media"],
               "%+.4f" % d["otro_sin_bloque"]["media"], "%+.4f" % d["otro_decil_sup"]["media"]))
print()
for nom, d in res["subconjuntos"].items():
    print(f"{nom}: dif todos={d['dif_todos']:+.4f}  sin bloque={d['dif_sin_bloque']:+.4f}  "
          f"decil10={d['dif_decil_sup']:+.4f}")
print("escrito", os.path.join(OUTDIR, "zeros.json"))
