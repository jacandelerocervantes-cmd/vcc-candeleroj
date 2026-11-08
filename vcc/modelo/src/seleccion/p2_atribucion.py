#!/usr/bin/env python
"""PASO 5a PARTE 2, etapa de ANALISIS. Lee `p2_datos.json`; no lee nada del
lado perturbado salvo la PERTENENCIA a R(g,L), que es la respuesta (y).

CONSTRUCCION, DECLARADA ANTES DE MEDIR
  Unidad: (g, gen, L1, L2), L1 != L2, gen en R(g,L1). 33 perts x 30 pares
  ordenados. Respuesta y = 1 si gen NO esta en R(g,L2) (especifico de L1
  frente a L2).
  Predictores, TODOS del CONTROL de L2, en ESTE ORDEN (fijado antes de medir):
    p1 pasa el filtro de L2 (CPM medio del control de L2 > 5.0)
    p2 CPM medio del control de L2: deciles de log10(CPM + 1e-3) (bordes de
       los cuantiles sobre las filas observadas; el nulo usa los MISMOS bordes)
    p3 tasa de deteccion en el control de L2: deciles, igual
    p4 potencia de L2: log(n de control) y log(UMI mediana/celula), lineales
  Medida A (determinista): fraccion de las filas y=1 cuyo gen no pasa el
    filtro de L2. Solo existe para p1.
  Medida B (devianza): regresion logistica (sklearn, L2 debil C=1e4), D^2 =
    1 - dev(modelo)/dev(solo intercepto). Cada predictor SOLO y ACUMULADO en
    el orden p1, p1+p2, p1+p2+p3, p1+..+p4.
  NULO: 20 sorteos, semilla 20281100+k; por sorteo, una permutacion UNIFORME
    pi de las 6 lineas y las estadisticas de control de L2 se toman de pi(L2)
    (y no cambia). Toda fraccion se da tambien ANCLADA:
    (obs - media_nulo) / (1 - media_nulo).
  RESIDUO = 1 - D^2 acumulado (crudo y anclado), y para la medida A,
    1 - fraccion determinista. Se da para R CRUDO y para R IGUALADO EN
    POTENCIA (ver p2_datos.py); en el igualado p4 es constante y se omite.
"""
import json
import os
import sys
from itertools import permutations

import numpy as np
from sklearn.linear_model import LogisticRegression

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import semillas as S  # noqa: E402

DAT = "/datos/modelo/out/seleccion/p2_datos.json"
OUT = "/datos/modelo/out/seleccion/p2_mecanismo.json"
LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
ORDEN = ["p1", "p2", "p3", "p4"]
EPS = 1e-12


def filas(var, panel, gidx):
    ed = {L: var[L]["ed"] for L in LINES}
    R = []
    for g in panel:
        for L1, L2 in permutations(LINES, 2):
            s2 = set(ed[L2].get(g, []))
            for gen in ed[L1].get(g, []):
                R.append((g, gidx[gen], LINES.index(L1), LINES.index(L2), int(gen not in s2)))
    a = np.array([(r[1], r[2], r[3], r[4]) for r in R], dtype=np.int64)
    return [r[0] for r in R], a[:, 0], a[:, 1], a[:, 2], a[:, 3]


def bordes(v):
    e = np.unique(np.quantile(v, np.linspace(0, 1, 11)[1:-1]))
    return e


def dummies(v, e):
    b = np.searchsorted(e, v, side="right")
    k = len(e) + 1
    M = np.zeros((v.size, k - 1))
    for j in range(1, k):
        M[:, j - 1] = b == j
    return M


def diseno(gi, L2, stat, E, que, usar_p4):
    cols = []
    cpm = stat["cpm"][L2, gi]
    det = stat["det"][L2, gi]
    if "p1" in que:
        cols.append((cpm > stat["umbral"]).astype(float)[:, None])
    if "p2" in que:
        cols.append(dummies(np.log10(cpm + 1e-3), E["p2"]))
    if "p3" in que:
        cols.append(dummies(det, E["p3"]))
    if "p4" in que and usar_p4:
        z = np.column_stack([stat["logn"][L2], stat["logL"][L2]])
        cols.append((z - E["p4_mu"]) / E["p4_sd"])
    return np.hstack(cols) if cols else None


def d2(X, y):
    p0 = y.mean()
    dev0 = -2 * np.sum(y * np.log(p0) + (1 - y) * np.log(1 - p0))
    if X is None or X.shape[1] == 0:
        return 0.0
    m = LogisticRegression(C=1e4, max_iter=5000, tol=1e-8)
    m.fit(X, y)
    p = np.clip(m.predict_proba(X)[:, 1], EPS, 1 - EPS)
    dev = -2 * np.sum(y * np.log(p) + (1 - y) * np.log(1 - p))
    return float(1 - dev / dev0)


def modelos(gi, L2, y, stat, E, usar_p4):
    out = {}
    for p in ORDEN:
        out["solo_" + p] = d2(diseno(gi, L2, stat, E, [p], usar_p4), y)
    for k in range(1, 5):
        out["acum_" + "+".join(ORDEN[:k])] = d2(diseno(gi, L2, stat, E, ORDEN[:k], usar_p4), y)
    m1 = y == 1
    out["A_det_p1"] = float(np.mean(stat["cpm"][L2[m1], gi[m1]] <= stat["umbral"])) if m1.any() else None
    return out


def analiza(var, panel, genes, umbral, etq):
    gidx = {g: i for i, g in enumerate(genes)}
    perts, gi, L1, L2, y = filas(var, panel, gidx)
    stat = {"cpm": np.array([var[L]["cpm"] for L in LINES]),
            "det": np.array([var[L]["det"] for L in LINES]),
            "logn": np.log(np.array([var[L]["n_control"] for L in LINES], float)),
            "logL": np.log(np.array([var[L]["L_mediana"] for L in LINES], float)),
            "umbral": umbral}
    usar_p4 = (np.ptp(stat["logn"]) > 0) or (np.ptp(stat["logL"]) > 0)
    z = np.column_stack([stat["logn"][L2], stat["logL"][L2]])
    sd = z.std(axis=0); sd[sd == 0] = 1.0
    E = {"p2": bordes(np.log10(stat["cpm"][L2, gi] + 1e-3)), "p3": bordes(stat["det"][L2, gi]),
         "p4_mu": z.mean(axis=0), "p4_sd": sd}
    # chequeo: todo gen de R(g,L1) pasa el filtro de L1
    viol = int(np.sum(stat["cpm"][L1, gi] <= umbral))
    obs = modelos(gi, L2, y, stat, E, usar_p4)
    nulo = []
    for k, sd_k in enumerate(S.SEMILLAS["p2_nulo"]):
        pi = np.random.default_rng(sd_k).permutation(6)
        r = modelos(gi, pi[L2], y, stat, E, usar_p4)
        r["perm"] = pi.tolist(); r["fijos"] = int(np.sum(pi == np.arange(6))); r["semilla"] = sd_k
        nulo.append(r)
        sys.stderr.write(f"  [{etq}] nulo k={k} fijos={r['fijos']} acum4={r['acum_p1+p2+p3+p4']:.4f}\n")
    # 2.1 por par ordenado
    por_par = []
    for a, b in permutations(range(6), 2):
        m = (L1 == a) & (L2 == b)
        m1 = m & (y == 1)
        nf = int(np.sum(stat["cpm"][b, gi[m1]] <= umbral))
        por_par.append({"L1": LINES[a], "L2": LINES[b], "n_filas": int(m.sum()), "n_especificos": int(m1.sum()),
                        "n_no_pasan_L2": nf})
    # intermedios: 20 filas y=1 equiespaciadas por CPM en L2
    i1 = np.where(y == 1)[0]
    i1 = i1[np.argsort(stat["cpm"][L2[i1], gi[i1]], kind="stable")]
    sel = i1[np.linspace(0, len(i1) - 1, 20).round().astype(int)] if len(i1) else []
    inter = [{"pert": perts[j], "gen": genes[gi[j]], "L1": LINES[L1[j]], "L2": LINES[L2[j]],
              "cpm_L1": float(stat["cpm"][L1[j], gi[j]]), "cpm_L2": float(stat["cpm"][L2[j], gi[j]]),
              "det_L2": float(stat["det"][L2[j], gi[j]]), "n_L2": int(var[LINES[L2[j]]]["n_control"]),
              "Lmed_L2": float(var[LINES[L2[j]]]["L_mediana"]), "y": int(y[j])} for j in sel]
    return {"n_filas": int(y.size), "n_y1": int(y.sum()), "frac_y1": float(y.mean()),
            "viol_filtro_L1": viol, "usar_p4": bool(usar_p4),
            "bordes": {"p2": E["p2"].tolist(), "p3": E["p3"].tolist()},
            "obs": obs, "nulo": nulo, "por_par": por_par, "intermedios": inter}


def universo(d):
    J = d["jiang"]; panel = J["panel"]; U = d["umbral"]
    out = {}
    for v in ("crudo", "igualado"):
        var = J[v]
        filt_any = set()
        for L in LINES:
            filt_any |= {J["genes"][i] for i, c in enumerate(var[L]["cpm"]) if c > U}
        filas_ = []
        for g in panel:
            sets = {L: set(var[L]["ed"].get(g, [])) for L in LINES}
            u = set().union(*sets.values())
            filas_.append({"pert": g, "U": len(u), "por_linea": {L: len(sets[L]) for L in LINES},
                           "en_6": len(set.intersection(*sets.values())),
                           "n_lineas_que_aportan": sum(1 for L in LINES if sets[L])})
        out[v] = {"filtrables_alguna_linea": len(filt_any), "por_pert": filas_}
    v = d["v2025"]
    filt = {v["genes"][i] for i, c in enumerate(v["cpm"]) if c > U}
    fp = []
    for g in v["perts"]:
        r = set(v["ed"].get(g, []))
        fp.append({"pert": g, "R": len(r), "R_fuera_filtrables": len(r - filt),
                   "frac_de_filtrables": len(r & filt) / len(filt),
                   "n_celulas": next(x["tomadas"] for x in v["por_pert"] if x["pert"] == g)})
    union = set().union(*[set(v["ed"].get(g, [])) for g in v["perts"]])
    out["v2025"] = {"n_filtrables": len(filt), "por_pert": fp, "union_R": len(union),
                    "union_R_frac": len(union & filt) / len(filt)}
    return out


def main():
    d = json.load(open(DAT))
    sem = S.comprueba(); assert sem["sin_choques"]
    J = d["jiang"]
    res = {"semillas": sem, "orden_declarado": ORDEN}
    for v in ("crudo", "igualado"):
        res[v] = analiza(J[v], J["panel"], J["genes"], d["umbral"], v)
    res["universo"] = universo(d)
    json.dump(res, open(OUT + ".tmp", "w"))
    os.replace(OUT + ".tmp", OUT)


if __name__ == "__main__":
    main()
