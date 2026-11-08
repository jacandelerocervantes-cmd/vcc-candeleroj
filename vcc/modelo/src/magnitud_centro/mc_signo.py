#!/usr/bin/env python
"""Diagnostico: ¿de donde sale la concordancia de SIGNO?

Sospecha a comprobar: un gen SIN expresion en la mitad perturbada da
`delta = log1p(0) - log1p(ctrl) = -log1p(ctrl)`, que es negativo y casi el
mismo en las dos mitades y en las dos lineas. Esos genes concuerdan en signo
por ESTRUCTURA del espacio, no por la perturbacion, e inflan la concordancia
justo donde el delta es pequeno.

Se mide la concordancia de signo por DECIL de |delta| (del lado predictor) y
por nivel de expresion del gen en el control, en un par de lineas.
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

res = {"pares": {}, "global": {}}
pares = [(a, b) for i, a in enumerate(LINES) for b in LINES[i + 1:]]

# concordancia por decil de |delta| del predictor, agregada sobre pares/perts/j
acc = {k: np.zeros(10) for k in ("n", "ok")}
accN = {k: np.zeros(10) for k in ("n", "ok")}
accR = {k: np.zeros(10) for k in ("n", "ok")}
frac_mismo_signo_ambos_neg = []
for (L1, L2) in pares:
    for j in range(J):
        for g in range(P):
            a1, b2 = D[L1][j, g, 0], D[L2][j, g, 1]
            q = np.quantile(np.abs(a1), np.linspace(0, 1, 11))
            q[0] = -np.inf; q[-1] = np.inf
            ok = np.sign(a1) == np.sign(b2)
            dig = np.clip(np.digitize(np.abs(a1), q[1:-1]), 0, 9)
            for k in range(10):
                m = dig == k
                acc["n"][k] += m.sum(); acc["ok"][k] += ok[m].sum()
                tgt = accN if panel[g] in NUCLEO else accR
                tgt["n"][k] += m.sum(); tgt["ok"][k] += ok[m].sum()
            # cuanta de la concordancia es "los dos negativos"
            amb_neg = ((a1 < 0) & (b2 < 0)).sum()
            frac_mismo_signo_ambos_neg.append(amb_neg / max(ok.sum(), 1))

res["por_decil_de_abs_delta"] = [
    {"decil": k + 1, "n": int(acc["n"][k]),
     "concordancia": float(acc["ok"][k] / acc["n"][k]),
     "concordancia_nucleo": float(accN["ok"][k] / accN["n"][k]),
     "concordancia_resto": float(accR["ok"][k] / accR["n"][k])}
    for k in range(10)]
res["frac_de_la_concordancia_que_es_ambos_negativos"] = {
    "media": float(np.mean(frac_mismo_signo_ambos_neg)),
    "min": float(np.min(frac_mismo_signo_ambos_neg)),
    "max": float(np.max(frac_mismo_signo_ambos_neg))}
# cuantos genes tienen delta negativo (el signo que impone el cero)
neg = []
for L in LINES:
    neg.append(float((D[L][:, :, :, :] < 0).mean()))
res["frac_de_genes_con_delta_negativo_por_linea"] = dict(zip(LINES, neg))

json.dump(res, open(os.path.join(OUTDIR, "signo.json"), "w"), indent=1, sort_keys=True)
print("decil |delta| | n          | concord. | nucleo  | resto")
for r in res["por_decil_de_abs_delta"]:
    print(f"   {r['decil']:2d}        | {r['n']:10,} | {r['concordancia']:.4f}   | "
          f"{r['concordancia_nucleo']:.4f}  | {r['concordancia_resto']:.4f}")
print()
f = res["frac_de_la_concordancia_que_es_ambos_negativos"]
print(f"fraccion de los genes concordantes que son AMBOS NEGATIVOS: "
      f"media {f['media']:.4f} [{f['min']:.4f}, {f['max']:.4f}]")
print("fraccion de genes con delta < 0, por linea:",
      {k: round(v, 4) for k, v in res["frac_de_genes_con_delta_negativo_por_linea"].items()})
print("escrito", os.path.join(OUTDIR, "signo.json"))
