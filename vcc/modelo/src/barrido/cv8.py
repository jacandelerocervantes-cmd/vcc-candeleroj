#!/usr/bin/env python
"""Recalcula y GUARDA las listas por blanco (Pearson y R2) de los dos modelos:
era el defecto declarado de la seccion 2.2 del informe."""
import json, sys, numpy as np
for p in ("/home/ubuntu/vcc/modelo/src", "/home/ubuntu/vcc/modelo/src/entrega"):
    sys.path.insert(0, p)
import datos as D, ajuste as A

J = "/home/ubuntu/vcc/modelo/out/barrido_k.json"
r = json.load(open(J))
alpha = float(r["alpha"]); pliegues = r["cv"]["blancos_de_cada_pliegue"] if "cv" in r else None
if pliegues is None:
    pliegues = json.load(open("/home/ubuntu/vcc/modelo/out/entrega_ajuste.json"))["cv"]["blancos_de_cada_pliegue"]
perts, genes25, cpm_pert, ncel, cpm_ctrl, det, fano, meta = D.pseudobulk_2025()
y_full = np.log2((cpm_pert + A.EPS_CPM) / (cpm_ctrl[None, :] + A.EPS_CPM)).astype(np.float32)
st = {"cpm": cpm_ctrl, "det": det, "fano": np.where(fano < 0, 0.0, fano), "n": meta["n_control_2025"]}
filtr = np.where(cpm_ctrl > D.UMBRAL_CPM)[0]
tf = D.tf_set(); gw_vn, gw_perf, _ = D.gwps(sorted(set(perts.tolist())))
g_fit = genes25[filtr]
st_fit = {k: (v[filtr] if isinstance(v, np.ndarray) else v) for k, v in st.items()}
X = D.matriz_rasgos(perts, g_fit, gw_vn, gw_perf, tf, st_fit)
Y = y_full[:, filtr].reshape(-1); G = len(g_fit)
i_fuera = [D.RASGOS.index(x) for x in r["ablacion_8"]["fuera"]]
i_8 = [i for i in range(len(D.RASGOS)) if i not in i_fuera]
fil = lambda ts: np.concatenate([np.arange(i * G, (i + 1) * G) for i in ts])
for nom, cols in (("10_rasgos", list(range(len(D.RASGOS)))), ("8_rasgos", i_8)):
    crudo = {}
    for pr in pliegues:
        ipr = [int(np.where(perts == t)[0][0]) for t in sorted(pr)]
        itr = [i for i in range(len(perts)) if i not in set(ipr)]
        m = A.ridge(X[fil(itr)][:, cols].astype(np.float64), Y[fil(itr)].astype(np.float64), alpha)
        ite = fil(ipr)
        yp = A.predice(m, X[ite][:, cols].astype(np.float64)); yo = Y[ite].astype(np.float64)
        for n_, t in enumerate(sorted(pr)):
            s = slice(n_ * G, (n_ + 1) * G)
            d = yo[s] - yp[s]
            crudo[str(t)] = [A.pearson(yp[s], yo[s]),
                             float(1 - (d ** 2).sum() / ((yo[s] - yo[s].mean()) ** 2).sum())]
    r["cv_por_modelo"][nom]["por_blanco_crudo"] = crudo
    print(nom, "guardados", len(crudo), "blancos")
json.dump(r, open(J, "w"), indent=1, default=str)
