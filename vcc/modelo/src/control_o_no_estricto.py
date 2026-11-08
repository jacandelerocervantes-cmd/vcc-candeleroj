#!/usr/bin/env python
"""La version ESTRICTA: el gen blanco fuera de los rasgos, fuera de las filas de
entrenamiento y fuera de la EVALUACION. Con el centrado por pliegue corregido de
`control_o_no_loto.py`.

POR QUE HACE FALTA. El rasgo `es_el_blanco` sube `R_agr` de forma consistente
(~+0.03 a +0.05 en las seis lineas), y hay que separar dos cosas distintas:

  * "cuando se perturba `g`, `g` baja" -- cierto, predecible de solo control (se
    sabe QUE gen se perturba sin mirar nada perturbado) y de un valor practico
    real, pero es UN gen de ~7000 cuyo delta especifico es enorme, y una
    correlacion la puede mover un solo punto extremo.
  * "la estructura de coexpresion del control predice que OTROS genes se mueven"
    -- la pregunta del encargo.

Se quita de los TRES sitios, no solo de los rasgos: si se queda en las filas de
entrenamiento es un punto extremo que el modelo persigue, y si se queda en la
evaluacion puede llevarse la correlacion entera.

Y se da la version RECORTADA (fuera el 1 % de |E| mas extremo) para ver si lo que
queda lo sostiene un punado de genes o esta repartido.

TODO EN ESPACIO DE BLOQUES, tambien con mascara: `Gram[i]` y `Cross[i,k]` se
calculan sobre las filas del blanco `i` que sobreviven a SU mascara, y la metrica
sale de ahi sin volver a tocar los genes. Por eso caben 500 permutaciones.

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
OUT = os.path.join(DIR, "estricto.json")
SEED_FAMILIA = 20268000
ALPHAS = (1.0, 10.0, 100.0, 1000.0, 10000.0)
N_PERM_RIDGE = 500
N_PERM_GB = 6
GB_KW = dict(max_iter=200, learning_rate=0.1)
RECORTE = 0.01


def resuelve(XtX, Xty, alpha):
    reg = np.eye(XtX.shape[0]) * alpha
    reg[0, 0] = 0.0
    try:
        return np.linalg.solve(XtX + reg, Xty)
    except np.linalg.LinAlgError:
        return np.linalg.lstsq(XtX + reg, Xty, rcond=None)[0]


def xty(idx, j, pi, Cross, nb):
    """Centrado por pliegue: y de entrenamiento = E_{pi(i)} + E_{pi(j)}/(nb-1)."""
    return (Cross[idx, pi[idx]].sum(0)
            + Cross[idx, pi[j]].sum(0) / float(nb - 1))


def piezas(w, j, i, Gram, Cross, Gj):
    """(numerador, ||p centrada||) con la media tomada SOBRE LOS GENES EVALUADOS."""
    s = Gram[j][0]
    num = float(w @ Cross[j, i]) - float(w @ s) / Gj * float(Cross[j, i][0])
    v = float(w @ Gram[j] @ w) - float(w @ s) ** 2 / Gj
    return num, (np.sqrt(v) if v > 0 else 0.0)


def resume(nums, nrms, cz, pi):
    ok = (nrms > 0) & (cz > 0)
    den = nrms[ok] * np.sqrt(cz[ok])
    r = nums[ok] / den
    return {"R_agr": float(nums[ok].sum() / den.sum()) if den.sum() > 0 else np.nan,
            "mediana": float(np.median(r)) if r.size else np.nan,
            "media": float(np.mean(r)) if r.size else np.nan,
            "n_positivos": int((r > 0).sum()), "n_validos": int(ok.sum()),
            "por_blanco": r.tolist(), "indices": np.where(ok)[0].tolist()}


def ridge(Gram, Cross, cruz_m, Gs, nb, perm=None, alpha_fijo=None):
    pi = np.arange(nb) if perm is None else np.asarray(perm)
    todos = np.arange(nb)
    nums, nrms, alf, ws = [], [], [], []
    for j in range(nb):
        tr = todos[todos != j]
        if alpha_fijo is not None:
            a_us = alpha_fijo
        else:
            mejor, a_us = -np.inf, ALPHAS[0]
            for a in ALPHAS:
                sn, sd = 0.0, 0.0
                for gr in np.array_split(tr, 4):
                    sub = np.setdiff1d(tr, gr)
                    w = resuelve(Gram[sub].sum(0), xty(sub, j, pi, Cross, nb), a)
                    for k in gr:
                        n_, p_ = piezas(w, k, pi[k], Gram, Cross, Gs[k])
                        c = cruz_m[k, pi[k]]
                        if c > 0 and p_ > 0:
                            sn += n_; sd += p_ * np.sqrt(c)
                m = sn / sd if sd > 0 else -np.inf
                if m > mejor:
                    mejor, a_us = m, a
        w = resuelve(Gram[tr].sum(0), xty(tr, j, pi, Cross, nb), a_us)
        n_, p_ = piezas(w, j, pi[j], Gram, Cross, Gs[j])
        nums.append(n_); nrms.append(p_); alf.append(a_us); ws.append(w)
    cz = np.array([cruz_m[j, pi[j]] for j in range(nb)])
    return np.array(nums), np.array(nrms), cz, alf, ws


def metrica_recortada(p, e1, e2, recorte=RECORTE):
    y = (e1 + e2) / 2.0
    k = max(1, int(round(recorte * len(y))))
    keep = np.argsort(-np.abs(y - y.mean()))[k:]
    pc = p[keep] - p[keep].mean()
    a = e1[keep] - e1[keep].mean(); b = e2[keep] - e2[keep].mean()
    c = float(a @ b); n = np.linalg.norm(pc)
    if c <= 0 or n == 0:
        return np.nan
    return float(pc @ ((a + b) / 2.0) / (n * np.sqrt(c)))


def main():
    lineas = sys.argv[1].split(",") if len(sys.argv) > 1 else LINES
    sal = json.load(open(OUT)) if os.path.exists(OUT) else {}
    for line in lineas:
        t0 = time.time()
        fold = LINES.index(line)
        seed = SEED_FAMILIA + 100 * fold
        d = np.load(os.path.join(DIR, "rasgos_%s.npz" % line), allow_pickle=True)
        F = d["F"].astype(np.float64)
        E1 = d["E1"].astype(np.float64); E2 = d["E2"].astype(np.float64)
        rasgos = [str(x) for x in d["rasgos"]]
        G, nb, P = F.shape
        E = (E1 + E2) / 2.0
        E = E - E.mean(axis=0, keepdims=True)
        auto = F[:, :, rasgos.index("es_el_blanco")] > 0.5
        masks = [~auto[:, j] for j in range(nb)]
        cols = [i for i, r in enumerate(rasgos) if r != "es_el_blanco"]
        Gs = np.array([int(m.sum()) for m in masks])

        flat = F[:, :, cols].reshape(-1, len(cols))
        mu = flat.mean(0); sg = flat.std(0); sg[sg == 0] = 1.0
        Zl, Gram, Cross = [], [], np.zeros((nb, nb, len(cols) + 1))
        for i in range(nb):
            Xi = (F[masks[i]][:, i, :][:, cols] - mu) / sg
            Zi = np.column_stack([np.ones(len(Xi)), Xi])
            Zl.append(Zi)
            Gram.append(Zi.T @ Zi)
            Cross[i] = (Zi.T @ E[masks[i]]).T
        Gram = np.array(Gram)
        # cruz del blanco k restringido a la mascara del pliegue j
        cruz_m = np.zeros((nb, nb))
        for j in range(nb):
            a = E1[masks[j]]; b = E2[masks[j]]
            a = a - a.mean(0); b = b - b.mean(0)
            cruz_m[j] = np.einsum("gk,gk->k", a, b)

        nums, nrms, cz, alf, ws = ridge(Gram, Cross, cruz_m, Gs, nb)
        obs = resume(nums, nrms, cz, np.arange(nb))
        a_modal = max(set(alf), key=alf.count)
        rec = [metrica_recortada(Zl[j] @ ws[j], E1[masks[j], j], E2[masks[j], j])
               for j in range(nb)]
        obs["recortado"] = rec
        obs["alpha_modal"] = a_modal
        nul = {"R_agr": [], "mediana": [], "n_positivos": []}
        for b in range(N_PERM_RIDGE):
            rp = np.random.default_rng(seed + 3500 + b)
            pi = rp.permutation(nb)
            while np.any(pi == np.arange(nb)):
                pi = rp.permutation(nb)
            n_, p_, c_, _, _ = ridge(Gram, Cross, cruz_m, Gs, nb, perm=pi,
                                     alpha_fijo=a_modal)
            rr = resume(n_, p_, c_, pi)
            for k in nul:
                nul[k].append(rr[k])
        obs["nulo"] = nul

        # --- GB estricto ---
        def gb(perm=None):
            pi = np.arange(nb) if perm is None else np.asarray(perm)
            todos = np.arange(nb)
            nn, pp = [], []
            for j in range(nb):
                tr = todos[todos != j]
                Xtr = np.concatenate([F[masks[i]][:, i, :][:, cols] for i in tr])
                ytr = np.concatenate([E[masks[i], pi[i]]
                                      + E[masks[i], pi[j]] / float(nb - 1)
                                      for i in tr])
                m = HistGradientBoostingRegressor(random_state=seed + j, **GB_KW)
                m.fit(Xtr, ytr)
                p = m.predict(F[masks[j]][:, j, :][:, cols])
                pc = p - p.mean()
                ey = E[masks[j], pi[j]]; ey = ey - ey.mean()
                nn.append(float(pc @ ey)); pp.append(float(np.linalg.norm(pc)))
            cz2 = np.array([cruz_m[j, pi[j]] for j in range(nb)])
            return np.array(nn), np.array(pp), cz2, pi

        gn, gp, gc, gpi = gb()
        go = resume(gn, gp, gc, gpi)
        gnul = []
        for b in range(N_PERM_GB):
            rp = np.random.default_rng(seed + 3700 + b)
            pi = rp.permutation(nb)
            while np.any(pi == np.arange(nb)):
                pi = rp.permutation(nb)
            n_, p_, c_, pp_ = gb(perm=pi)
            gnul.append(resume(n_, p_, c_, pp_)["R_agr"])
        go["nulo_R_agr"] = gnul

        sal[line] = {"G_evaluados": Gs.tolist(), "ridge": obs, "gb": go,
                     "panel": [str(x) for x in d["panel"]],
                     "segundos": round(time.time() - t0, 1)}
        sd = float(np.std(nul["R_agr"], ddof=1))
        sys.stderr.write(
            "%-6s ESTRICTO ridge R_agr=%+.4f med=%+.4f pos=%d/%d | nulo=%+.4f "
            "sd=%.4f z=%+.2f MDE80=%.4f p=%.4f | recortado=%+.4f | GB=%+.4f\n"
            % (line, obs["R_agr"], obs["mediana"], obs["n_positivos"],
               obs["n_validos"], float(np.mean(nul["R_agr"])), sd,
               (obs["R_agr"] - np.mean(nul["R_agr"])) / sd if sd else np.nan,
               2.802 * sd,
               (1 + sum(1 for v in nul["R_agr"] if v >= obs["R_agr"]))
               / (1 + len(nul["R_agr"])),
               float(np.nanmedian(rec)), go["R_agr"]))
        sys.stderr.flush()
        with open(OUT, "w") as fh:
            json.dump(sal, fh)
    sys.stderr.write("escrito %s\n" % OUT)


if __name__ == "__main__":
    main()
