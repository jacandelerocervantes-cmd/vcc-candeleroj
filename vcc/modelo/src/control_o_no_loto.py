#!/usr/bin/env python
"""LOTO por BLANCO con el centrado corregido por pliegue. Reemplaza a
`control_o_no_ajuste.py`, cuyo titular estaba contaminado.

🔴 EL DEFECTO QUE ESTO CORRIGE, medido y no sospechado.
El objetivo es el delta ESPECIFICO `E_i = h_i - m`, con `m` el promedio sobre los
33 blancos. Eso obliga a `Σ_i E_i = 0` por construccion. Y entonces, para el
pliegue que deja fuera el blanco `j`:

    media de los 32 de entrenamiento = (Σ_{i≠j} E_i) / 32 = -E_j / 32

Medido en A549: `corr(media de entrenamiento, blanco dejado fuera) = -1.000000`
exacto, y `max|media_tr - (-E_j/32)| = 4.4e-09`. **La media del entrenamiento ES
el negativo del blanco que se dejo fuera.** Cualquier modelo que capture esa
media se lleva correlacion NEGATIVA con el blanco nuevo sin saber nada. Lo delato
el conjunto de control `solo_h`, que no puede llevar informacion especifica del
blanco y puntuaba de -0.13 a -0.22 con desviacion EXACTAMENTE cero.

EL ARREGLO, derivado y no tanteado. El estimando correcto por pliegue usa como
promedio de referencia solo los 32 de entrenamiento: `M_j = media_{k≠j} h_k`.
Con `y_i^(j) = h_i - M_j` sale, sustituyendo `h_i = E_i + S/33` y
`Σ_{k≠j} E_k = -E_j`:

    y_i^(j) = E_i + E_j / 32          <- la y de ENTRENAMIENTO del pliegue j
    Σ_{i≠j} y_i^(j) = -E_j + 32·E_j/32 = 0     <- la media de entrenamiento es
                                                   CERO: no queda fuga

y el objetivo de EVALUACION del blanco dejado fuera es
`h_j - M_j = (33/32)·E_j`, un multiplo POSITIVO del de antes: la correlacion de
evaluacion no cambia, y el `cruz` tampoco al cancelarse la escala. Es decir, el
arreglo toca SOLO las `y` de entrenamiento, que es exactamente donde estaba la
fuga. No hace falta recalcular ningun rasgo.

EL RESUMEN. `R_agr = Σ_g <p_g,y_g> / Σ_g (||p_g||·sqrt(cruz_g))`, que es una media
de los `r_true` ponderada y vale CERO EXACTO para un predictor de solo-`h`
(porque su numerador total es `<p, Σ_g E_g> = 0`). Se dan tambien la mediana, la
media sin ponderar (la que estaba sesgada) y el conteo de signos.

EL NULO por permutacion de etiquetas de blanco, con el MISMO centrado por
pliegue, sobre el MISMO estadistico agrupado.

🔴 TRANSDUCTIVO. Se ajusta sobre la verdad de Jiang. Nada de aqui sale a un envio.
"""
import json
import os
import sys
import time

import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor

LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
DIR = "/home/ubuntu/vcc/modelo/out/control_o_no"
OUT = os.path.join(DIR, "loto.json")
SEED_FAMILIA = 20268000
ALPHAS = (1.0, 10.0, 100.0, 1000.0, 10000.0)
N_PERM_RIDGE = 500
N_PERM_GB = 8
GB_KW = dict(max_iter=200, learning_rate=0.1)
POS_COV = ["cov", "r_pearson", "r_spearman", "abs_r_pearson", "abs_r_spearman"]
POS_H = ["h_media", "h_var", "h_disp", "h_frac_cero", "h_rango_expr"]


def resuelve(XtX, Xty, alpha):
    reg = np.eye(XtX.shape[0]) * alpha
    reg[0, 0] = 0.0
    try:
        return np.linalg.solve(XtX + reg, Xty)
    except np.linalg.LinAlgError:
        return np.linalg.lstsq(XtX + reg, Xty, rcond=None)[0]


def xty_pliegue(idx, j, pi, Cross, nb):
    """Σ_{i en idx} Z_i^T y_{pi(i)}^(j) con y^(j) = E + E_{pi(j)}/32, o sea el
    centrado por pliegue. `Cross[i,k] = Z_i^T E_k`."""
    base = Cross[idx, pi[idx]].sum(0)
    corr = Cross[idx, pi[j]].sum(0) / float(nb - 1)
    return base + corr


def piezas(w, j, i, Gram, Cross, G):
    num = float(w @ Cross[j, i])
    v = float(w @ Gram[j] @ w) - float(w @ Gram[j][0]) ** 2 / G
    return num, (np.sqrt(v) if v > 0 else 0.0)


def resume(nums, nrms, cruz, pi):
    ok = (nrms > 0) & (cruz[pi] > 0)
    den = nrms[ok] * np.sqrt(cruz[pi][ok])
    r = nums[ok] / den
    return {"R_agr": float(nums[ok].sum() / den.sum()) if den.sum() > 0 else np.nan,
            "mediana": float(np.median(r)) if r.size else np.nan,
            "media": float(np.mean(r)) if r.size else np.nan,
            "n_positivos": int((r > 0).sum()), "n_validos": int(ok.sum()),
            "por_blanco": r.tolist(), "indices": np.where(ok)[0].tolist()}


def ridge_loto(Gram, Cross, cruz, G, nb, perm=None, alpha_fijo=None):
    pi = np.arange(nb) if perm is None else np.asarray(perm)
    todos = np.arange(nb)
    nums, nrms, alf = [], [], []
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
                    # CV interna: el centrado del pliegue INTERNO excluye el grupo
                    w = resuelve(Gram[sub].sum(0),
                                 xty_pliegue(sub, j, pi, Cross, nb), a)
                    for k in gr:
                        n_, p_ = piezas(w, k, pi[k], Gram, Cross, G)
                        if cruz[pi[k]] > 0 and p_ > 0:
                            sn += n_; sd += p_ * np.sqrt(cruz[pi[k]])
                m = sn / sd if sd > 0 else -np.inf
                if m > mejor:
                    mejor, a_us = m, a
        w = resuelve(Gram[tr].sum(0), xty_pliegue(tr, j, pi, Cross, nb), a_us)
        n_, p_ = piezas(w, j, pi[j], Gram, Cross, G)
        nums.append(n_); nrms.append(p_); alf.append(a_us)
    return np.array(nums), np.array(nrms), alf


def gb_loto(F, E, cruz, seed, perm=None, importancia=False):
    """GB con el mismo centrado por pliegue: y de entrenamiento = E + E_j/32."""
    G, nb, P = F.shape
    pi = np.arange(nb) if perm is None else np.asarray(perm)
    todos = np.arange(nb)
    nums, nrms = [], []
    imp = np.zeros((nb, P)) if importancia else None
    rng = np.random.default_rng(seed)
    for j in range(nb):
        tr = todos[todos != j]
        Xtr = F[:, tr, :].reshape(-1, P)
        ytr = np.concatenate([E[:, pi[i]] + E[:, pi[j]] / float(nb - 1)
                              for i in tr])
        m = HistGradientBoostingRegressor(random_state=seed + j, **GB_KW)
        m.fit(Xtr, ytr)
        p = m.predict(F[:, j, :])
        pc = p - p.mean()
        nums.append(float(pc @ E[:, pi[j]])); nrms.append(float(np.linalg.norm(pc)))
        if importancia and nrms[-1] > 0 and cruz[pi[j]] > 0:
            base = nums[-1] / (nrms[-1] * np.sqrt(cruz[pi[j]]))
            Xte = F[:, j, :]
            for q in range(P):
                Xp = Xte.copy()
                Xp[:, q] = Xp[rng.permutation(G), q]
                pp = m.predict(Xp); ppc = pp - pp.mean()
                nn = np.linalg.norm(ppc)
                imp[j, q] = base - (float(ppc @ E[:, pi[j]]) /
                                    (nn * np.sqrt(cruz[pi[j]])) if nn > 0 else 0.0)
    return np.array(nums), np.array(nrms), imp


def main():
    lineas = sys.argv[1].split(",") if len(sys.argv) > 1 else LINES
    n_gb = int(sys.argv[2]) if len(sys.argv) > 2 else N_PERM_GB
    sal = json.load(open(OUT)) if os.path.exists(OUT) else {}
    for line in lineas:
        t0 = time.time()
        fold = LINES.index(line)
        seed = SEED_FAMILIA + 100 * fold
        d = np.load(os.path.join(DIR, "rasgos_%s.npz" % line), allow_pickle=True)
        F = d["F"].astype(np.float64)
        E = ((d["E1"].astype(np.float64) + d["E2"].astype(np.float64)) / 2.0)
        E = E - E.mean(axis=0, keepdims=True)
        rasgos = [str(x) for x in d["rasgos"]]
        cruz = d["cruz"].astype(np.float64)
        G, nb, P = F.shape
        conj = {"todos": list(range(P)),
                "sin_autoblanco": [i for i, r in enumerate(rasgos)
                                   if r != "es_el_blanco"],
                "solo_covarianza": [i for i, r in enumerate(rasgos)
                                    if r in POS_COV],
                "solo_h": [i for i, r in enumerate(rasgos) if r in POS_H]}
        ent = {"G": int(G), "n_blancos": int(nb),
               "A_esp": d["A_esp"].tolist(), "cruz": cruz.tolist(),
               "panel": [str(x) for x in d["panel"]],
               "blanco_expresado": d["blanco_expresado"].tolist(),
               "rasgos": rasgos, "conjuntos": {}}
        for nom, cols in conj.items():
            Fs = F[:, :, cols]
            flat = Fs.reshape(-1, len(cols))
            mu = flat.mean(0); sg = flat.std(0); sg[sg == 0] = 1.0
            Z = np.concatenate([np.ones((G, nb, 1)), (Fs - mu) / sg], axis=2)
            Gram = np.einsum("gjp,gjq->jpq", Z, Z)
            Cross = np.einsum("gjp,gi->jip", Z, E)
            nums, nrms, alf = ridge_loto(Gram, Cross, cruz, G, nb)
            obs = resume(nums, nrms, cruz, np.arange(nb))
            a_modal = max(set(alf), key=alf.count)
            nul = {"R_agr": [], "mediana": [], "n_positivos": []}
            for b in range(N_PERM_RIDGE):
                rp = np.random.default_rng(seed + 2500 + b)
                pi = rp.permutation(nb)
                while np.any(pi == np.arange(nb)):
                    pi = rp.permutation(nb)
                n_, p_, _ = ridge_loto(Gram, Cross, cruz, G, nb, perm=pi,
                                       alpha_fijo=a_modal)
                rr = resume(n_, p_, cruz, pi)
                for k in nul:
                    nul[k].append(rr[k])
            obs["alpha_modal"] = a_modal
            obs["nulo"] = nul
            ent["conjuntos"][nom] = obs
            sd = float(np.std(nul["R_agr"], ddof=1))
            sys.stderr.write(
                "%-6s %-16s R_agr=%+.4f med=%+.4f media=%+.4f pos=%d/%d | "
                "nulo=%+.4f sd=%.4f z=%+.2f MDE80=%.4f p=%.4f\n"
                % (line, nom, obs["R_agr"], obs["mediana"], obs["media"],
                   obs["n_positivos"], obs["n_validos"],
                   float(np.mean(nul["R_agr"])), sd,
                   (obs["R_agr"] - np.mean(nul["R_agr"])) / sd if sd else np.nan,
                   2.802 * sd,
                   (1 + sum(1 for v in nul["R_agr"] if v >= obs["R_agr"]))
                   / (1 + len(nul["R_agr"]))))
            sys.stderr.flush()
        for nom in ("todos", "sin_autoblanco"):
            cols = conj[nom]
            gn, gr_, imp = gb_loto(F[:, :, cols], E, cruz, seed,
                                   importancia=(nom == "todos"))
            go = resume(gn, gr_, cruz, np.arange(nb))
            gnul = []
            for b in range(n_gb if nom == "todos" else 0):
                rp = np.random.default_rng(seed + 2700 + b)
                pi = rp.permutation(nb)
                while np.any(pi == np.arange(nb)):
                    pi = rp.permutation(nb)
                n_, p_, _ = gb_loto(F[:, :, cols], E, cruz, seed, perm=pi)
                gnul.append(resume(n_, p_, cruz, pi)["R_agr"])
            go["nulo_R_agr"] = gnul
            if imp is not None:
                go["importancia"] = np.nanmean(imp, axis=0).tolist()
                go["importancia_rasgos"] = [rasgos[i] for i in cols]
            ent["conjuntos"][nom]["gb"] = go
            sys.stderr.write(
                "%-6s %-16s GB R_agr=%+.4f med=%+.4f pos=%d/%d | nulo=%s\n"
                % (line, nom, go["R_agr"], go["mediana"], go["n_positivos"],
                   go["n_validos"],
                   ("%+.4f sd=%.4f" % (np.mean(gnul), np.std(gnul, ddof=1)))
                   if gnul else "n/d"))
            sys.stderr.flush()
        ent["segundos"] = round(time.time() - t0, 1)
        sal[line] = ent
        with open(OUT, "w") as fh:
            json.dump(sal, fh)
    sys.stderr.write("escrito %s\n" % OUT)


if __name__ == "__main__":
    main()
