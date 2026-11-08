#!/usr/bin/env python
"""EL DISCRIMINADOR, en sus dos versiones: con el gen blanco y SIN el.

QUE MIDE. Se ajusta en la MITAD de los genes DEL PROPIO blanco y se predice la
otra mitad. NO es generalizacion a un blanco nuevo y no se puede leer como tal.
Responde a otra pregunta, que es la que separa (a) de (b):

    ¿llevan los rasgos de solo-control ALGUNA informacion sobre el delta
    especifico de este blanco?

  r_mismo ~ 0  ->  (b) nada del control lo predice: los rasgos no saben nada
  r_mismo > 0  ->  los rasgos SI llevan informacion; si ademas `r_loto` no la
                   alcanza, lo que falta es transferencia entre blancos, no senal

POR QUE HACEN FALTA LAS DOS VERSIONES, y es lo que invalido la primera medicion.
`r_pearson(g,g) = 1.0000` exacto mientras el maximo en los demas ~7000 genes es
+0.3031; `dist_emb(g,g) = 0` contra un p99.99 de +0.818. **El gen blanco se
identifica con certeza desde los rasgos de correlacion solos**, asi que quitar el
indicador `es_el_blanco` no quita nada: el modelo sigue sabiendo cual es. Y su
`|E|` es 10.2 veces la mediana. Entonces "el gen perturbado baja" -- que es
cierto, es de solo control y tiene valor, pero es UN gen -- puede estar
sosteniendo el numero entero.

  version CON    : todos los genes, todos los rasgos
  version SIN    : el gen blanco fuera de los rasgos, del entrenamiento Y de la
                   evaluacion

EL NULO: permutar el objetivo sobre los GENES dentro del blanco, manteniendo los
rasgos. Rompe la correspondencia gen-objetivo sin tocar nada mas.

🔴 TRANSDUCTIVO. Nada de aqui sale a un envio.
"""
import json
import os
import sys
import time

import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor

LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
DIR = "/home/ubuntu/vcc/modelo/out/control_o_no"
OUT = os.path.join(DIR, "disc.json")
SEED_FAMILIA = 20268000
N_PERM = 20
GB_KW = dict(max_iter=200, learning_rate=0.1)
ALPHA = 10.0


def resuelve(XtX, Xty, alpha):
    reg = np.eye(XtX.shape[0]) * alpha
    reg[0, 0] = 0.0
    try:
        return np.linalg.solve(XtX + reg, Xty)
    except np.linalg.LinAlgError:
        return np.linalg.lstsq(XtX + reg, Xty, rcond=None)[0]


def metrica(p, e1, e2):
    pc = p - p.mean()
    a = e1 - e1.mean(); b = e2 - e2.mean()
    c = float(a @ b); n = np.linalg.norm(pc)
    if c <= 0 or n == 0:
        return np.nan, 0.0, 0.0
    num = float(pc @ ((a + b) / 2.0))
    return num / (n * np.sqrt(c)), num, n * np.sqrt(c)


def una(X, e1, e2, mitad, seed, modelo):
    """Dos pliegues de genes; devuelve (r por pliegue, num, den) acumulados."""
    rr, N, D = [], 0.0, 0.0
    for sel in (mitad, ~mitad):
        tr, te = ~sel, sel
        y = (e1[tr] + e2[tr]) / 2.0
        y = y - y.mean()
        if modelo == "gb":
            m = HistGradientBoostingRegressor(random_state=seed, **GB_KW)
            m.fit(X[tr], y)
            p = m.predict(X[te])
        else:
            mu = X[tr].mean(0); sg = X[tr].std(0); sg[sg == 0] = 1.0
            Zt = np.column_stack([np.ones(int(tr.sum())), (X[tr] - mu) / sg])
            Ze = np.column_stack([np.ones(int(te.sum())), (X[te] - mu) / sg])
            p = Ze @ resuelve(Zt.T @ Zt, Zt.T @ y, ALPHA)
        r, n_, d_ = metrica(p, e1[te], e2[te])
        rr.append(r)
        if np.isfinite(r):
            N += n_; D += d_
    return rr, N, D


def main():
    lineas = sys.argv[1].split(",") if len(sys.argv) > 1 else LINES
    n_perm = int(sys.argv[2]) if len(sys.argv) > 2 else N_PERM
    sal = json.load(open(OUT)) if os.path.exists(OUT) else {}
    for line in lineas:
        t0 = time.time()
        fold = LINES.index(line)
        seed = SEED_FAMILIA + 100 * fold + 60
        d = np.load(os.path.join(DIR, "rasgos_%s.npz" % line), allow_pickle=True)
        F = d["F"].astype(np.float64)
        E1 = d["E1"].astype(np.float64); E2 = d["E2"].astype(np.float64)
        rasgos = [str(x) for x in d["rasgos"]]
        G, nb, P = F.shape
        auto = F[:, :, rasgos.index("es_el_blanco")] > 0.5
        cols_sin = [i for i, r in enumerate(rasgos) if r != "es_el_blanco"]
        ent = {"panel": [str(x) for x in d["panel"]], "versiones": {}}

        for ver in ("con", "sin"):
            for modelo in ("gb", "ridge"):
                por_b, N, D = [], 0.0, 0.0
                for j in range(nb):
                    m = np.ones(G, bool) if ver == "con" else ~auto[:, j]
                    cols = list(range(P)) if ver == "con" else cols_sin
                    X = F[m][:, j, :][:, cols]
                    e1 = E1[m, j]; e2 = E2[m, j]
                    rng = np.random.default_rng(seed + j)
                    mitad = rng.permutation(int(m.sum())) < int(m.sum()) // 2
                    rr, n_, d_ = una(X, e1, e2, mitad, seed + j, modelo)
                    por_b.append(float(np.nanmean(rr))); N += n_; D += d_
                obs = {"por_blanco": por_b,
                       "R_agr": (N / D) if D > 0 else np.nan,
                       "mediana": float(np.nanmedian(por_b)),
                       "media": float(np.nanmean(por_b)),
                       "n_positivos": int(np.sum(np.array(por_b) > 0))}
                nul = []
                for b in range(n_perm if modelo == "gb" else n_perm):
                    rp = np.random.default_rng(seed + 5000 + b)
                    Nn, Dn = 0.0, 0.0
                    for j in range(nb):
                        m = np.ones(G, bool) if ver == "con" else ~auto[:, j]
                        cols = list(range(P)) if ver == "con" else cols_sin
                        X = F[m][:, j, :][:, cols]
                        o = rp.permutation(int(m.sum()))
                        e1 = E1[m, j][o]; e2 = E2[m, j][o]
                        rng = np.random.default_rng(seed + j)
                        mitad = rng.permutation(int(m.sum())) < int(m.sum()) // 2
                        _, n_, d_ = una(X, e1, e2, mitad, seed + j, modelo)
                        Nn += n_; Dn += d_
                    nul.append((Nn / Dn) if Dn > 0 else np.nan)
                obs["nulo_R_agr"] = nul
                ent["versiones"]["%s_%s" % (ver, modelo)] = obs
                sd = float(np.nanstd(nul, ddof=1)) if len(nul) > 1 else np.nan
                sys.stderr.write(
                    "%-6s %-3s %-5s R_agr=%+.4f med=%+.4f pos=%d/%d | "
                    "nulo=%+.4f sd=%.4f z=%+.2f\n"
                    % (line, ver, modelo, obs["R_agr"], obs["mediana"],
                       obs["n_positivos"], nb, float(np.nanmean(nul)), sd,
                       (obs["R_agr"] - np.nanmean(nul)) / sd if sd else np.nan))
                sys.stderr.flush()
        ent["segundos"] = round(time.time() - t0, 1)
        sal[line] = ent
        with open(OUT, "w") as fh:
            json.dump(sal, fh)
    sys.stderr.write("escrito %s\n" % OUT)


if __name__ == "__main__":
    main()
