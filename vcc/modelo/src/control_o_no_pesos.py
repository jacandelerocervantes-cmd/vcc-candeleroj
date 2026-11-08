#!/usr/bin/env python
"""Que rasgos pesan, y si su signo es CONSISTENTE entre blancos.

DONDE SE MIRA Y POR QUE. El LOTO estricto sale ~0, asi que las importancias de un
modelo que no predice nada no dicen nada. La senal esta en el ajuste DENTRO del
blanco (el discriminador), asi que es ahi donde hay que leer los pesos: un ridge
por blanco, con el gen blanco fuera de los rasgos y de las filas, estandarizado,
para que los coeficientes sean comparables entre si.

Y LA PREGUNTA QUE DECIDE SI HAY VIA: no cuales pesan, sino si **pesan en el mismo
sentido en todos los blancos**. Si el coeficiente de `r_pearson` es negativo en la
mitad de los blancos y positivo en la otra, no hay un mapa compartido que
aprender, y eso explica por si solo que el LOTO salga cero aunque dentro de cada
blanco la relacion sea fuerte. Se reporta por rasgo:

  |coef| medio      cuanto pesa
  signo mediano     en que sentido
  consistencia      fraccion de blancos con el signo del mediano (0.5 = moneda)

🔴 TRANSDUCTIVO. Se ajusta sobre la verdad de Jiang. Nada de aqui sale a un envio.
"""
import json
import os
import sys

import numpy as np

LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
DIR = "/home/ubuntu/vcc/modelo/out/control_o_no"
OUT = os.path.join(DIR, "pesos.json")
ALPHA = 10.0


def main():
    sal = {}
    for line in LINES:
        d = np.load(os.path.join(DIR, "rasgos_%s.npz" % line), allow_pickle=True)
        F = d["F"].astype(np.float64)
        E1 = d["E1"].astype(np.float64); E2 = d["E2"].astype(np.float64)
        rasgos = [str(x) for x in d["rasgos"]]
        G, nb, P = F.shape
        E = (E1 + E2) / 2.0
        E = E - E.mean(axis=0, keepdims=True)
        auto = F[:, :, rasgos.index("es_el_blanco")] > 0.5
        cols = [i for i, r in enumerate(rasgos) if r != "es_el_blanco"]
        nom = [rasgos[i] for i in cols]
        W = np.zeros((nb, len(cols)))
        for j in range(nb):
            m = ~auto[:, j]
            X = F[m][:, j, :][:, cols]
            y = E[m, j]; y = y - y.mean()
            mu = X.mean(0); sg = X.std(0); sg[sg == 0] = 1.0
            Z = (X - mu) / sg
            A = Z.T @ Z + np.eye(len(cols)) * ALPHA
            W[j] = np.linalg.solve(A, Z.T @ y)
        med = np.median(W, axis=0)
        # 🔴 consistencia = fraccion de blancos con el signo MAYORITARIO, contando
        # los ceros en ninguno de los dos lados. Comparar contra `sign(mediana)`
        # estaba roto: para un rasgo con mediana exactamente 0, `sign(0)==sign(0)`
        # contaba los ceros como aciertos y daba valores por DEBAJO del 50 %, que
        # es imposible para una consistencia de signo bien definida.
        pos = (W > 0).mean(axis=0); neg = (W < 0).mean(axis=0)
        cons = np.maximum(pos, neg)
        cero = (W == 0).mean(axis=0)
        sal[line] = {"rasgos": nom, "coef_por_blanco": W.tolist(),
                     "abs_medio": np.abs(W).mean(0).tolist(),
                     "mediano": med.tolist(), "consistencia": cons.tolist(),
                     "frac_positivos": pos.tolist(), "frac_cero": cero.tolist()}
        orden = np.argsort(-np.abs(W).mean(0))[:6]
        sys.stderr.write("%-6s  " % line + "  ".join(
            "%s:%+.3f(%.0f%%)" % (nom[q], med[q], 100 * cons[q]) for q in orden)
            + "\n")
        sys.stderr.flush()
    with open(OUT, "w") as fh:
        json.dump(sal, fh)
    sys.stderr.write("escrito %s\n" % OUT)


if __name__ == "__main__":
    main()
