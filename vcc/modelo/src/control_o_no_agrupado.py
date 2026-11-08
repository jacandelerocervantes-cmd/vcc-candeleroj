#!/usr/bin/env python
"""El resumen AGRUPADO, y por que la media sin ponderar no servia.

🔴 EL DEFECTO QUE ESTO CORRIGE, declarado porque es mio y porque lo delato una
comprobacion que puse a proposito. Se declaro que un predictor que solo dependa
del gen `h` tiene correlacion MEDIA exactamente cero sobre los blancos, porque el
objetivo esta centrado sobre los 33 blancos y por tanto `Σ_g E[h,g] = 0`. Eso es
cierto de los NUMERADORES `<p, E_g>`, no del cociente: cada blanco divide por su
propio `sqrt(cruz_g)`, y una media de cocientes con denominadores distintos no
hereda la propiedad. Medido, el conjunto `solo_h` daba de -0.13 a -0.22 con
desviacion EXACTAMENTE cero -- un predictor incapaz de llevar nada especifico del
blanco puntuando fuerte y negativo. El resumen estaba roto, no el experimento.

EL RESUMEN QUE SI LA HEREDA:

    R_agr = Σ_g <p_g, y_g>  /  Σ_g ( ||p_g|| · sqrt(cruz_g) )

Es una media de los `r_true_g` ponderada por `||p_g||·sqrt(cruz_g)`, y como el
numerador total es `<p, Σ_g E_g> = 0` para un predictor de solo-`h`, vale CERO
EXACTO para ellos. Se comprueba sobre `solo_h` y tiene que dar ~0: si no da ~0,
el resumen sigue roto y no se interpreta nada.

Se reportan los tres, y el agrupado es el titular:
  R_agr      ponderado, hereda la propiedad de cero  <- el que se lee
  mediana    robusto a un blanco con `cruz` chico
  media      el roto, para que se vea cuanto movia

Y el conteo de signos sobre los 33 blancos, que no depende de ninguna
ponderacion.

Todo en espacio de BLOQUES: la metrica no vuelve a tocar los G genes, asi que las
permutaciones salen gratis y el nulo se calcula sobre el MISMO estadistico
agrupado y no sobre otro.

🔴 TRANSDUCTIVO. Nada de aqui sale a un envio.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
DIR = "/home/ubuntu/vcc/modelo/out/control_o_no"
OUT = os.path.join(DIR, "agrupado.json")
SEED_FAMILIA = 20268000
ALPHAS = (1.0, 10.0, 100.0, 1000.0, 10000.0)
N_PERM = 500
POS_COV = ["cov", "r_pearson", "r_spearman", "abs_r_pearson", "abs_r_spearman"]
POS_H = ["h_media", "h_var", "h_disp", "h_frac_cero", "h_rango_expr"]


def resuelve(XtX, Xty, alpha):
    reg = np.eye(XtX.shape[0]) * alpha
    reg[0, 0] = 0.0
    try:
        return np.linalg.solve(XtX + reg, Xty)
    except np.linalg.LinAlgError:
        return np.linalg.lstsq(XtX + reg, Xty, rcond=None)[0]


def piezas(w, j, i, Gram, Cross, cruz, G):
    """(numerador, ||p||) para la prediccion del blanco j contra el objetivo i."""
    num = float(w @ Cross[j, i])
    v = float(w @ Gram[j] @ w) - float(w @ Gram[j][0]) ** 2 / G
    return num, (np.sqrt(v) if v > 0 else 0.0)


def loto(Gram, Cross, cruz, G, nb, perm=None, alpha_fijo=None):
    """Devuelve numeradores, normas y el alpha usado, por blanco."""
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
                    w = resuelve(Gram[sub].sum(0), Cross[sub, pi[sub]].sum(0), a)
                    for k in gr:
                        n_, p_ = piezas(w, k, pi[k], Gram, Cross, cruz, G)
                        if cruz[pi[k]] > 0 and p_ > 0:
                            sn += n_; sd += p_ * np.sqrt(cruz[pi[k]])
                m = sn / sd if sd > 0 else -np.inf
                if m > mejor:
                    mejor, a_us = m, a
        w = resuelve(Gram[tr].sum(0), Cross[tr, pi[tr]].sum(0), a_us)
        n_, p_ = piezas(w, j, pi[j], Gram, Cross, cruz, G)
        nums.append(n_); nrms.append(p_); alf.append(a_us)
    return np.array(nums), np.array(nrms), alf


def resume(nums, nrms, cruz, pi):
    ok = (nrms > 0) & (cruz[pi] > 0)
    den = nrms[ok] * np.sqrt(cruz[pi][ok])
    r = nums[ok] / den
    return {"R_agr": float(nums[ok].sum() / den.sum()) if den.sum() > 0 else np.nan,
            "mediana": float(np.median(r)) if r.size else np.nan,
            "media": float(np.mean(r)) if r.size else np.nan,
            "n_positivos": int((r > 0).sum()), "n_validos": int(ok.sum()),
            "por_blanco": r.tolist(),
            "indices": np.where(ok)[0].tolist()}


def main():
    lineas = sys.argv[1].split(",") if len(sys.argv) > 1 else LINES
    sal = json.load(open(OUT)) if os.path.exists(OUT) else {}
    for line in lineas:
        fold = LINES.index(line)
        seed = SEED_FAMILIA + 100 * fold
        d = np.load(os.path.join(DIR, "rasgos_%s.npz" % line), allow_pickle=True)
        F = d["F"].astype(np.float64); E1 = d["E1"]; E2 = d["E2"]
        rasgos = [str(x) for x in d["rasgos"]]
        cruz = d["cruz"].astype(np.float64)
        G, nb, P = F.shape
        Y = ((E1 + E2) / 2.0).astype(np.float64)
        Y = Y - Y.mean(axis=0, keepdims=True)
        conj = {"todos": list(range(P)),
                "sin_autoblanco": [i for i, r in enumerate(rasgos)
                                   if r != "es_el_blanco"],
                "solo_covarianza": [i for i, r in enumerate(rasgos)
                                    if r in POS_COV],
                "solo_h": [i for i, r in enumerate(rasgos) if r in POS_H]}
        ent = {"G": int(G), "A_esp": d["A_esp"].tolist(), "cruz": cruz.tolist(),
               "panel": [str(x) for x in d["panel"]],
               "blanco_expresado": d["blanco_expresado"].tolist(),
               "conjuntos": {}}
        for nom, cols in conj.items():
            Fs = F[:, :, cols]
            flat = Fs.reshape(-1, len(cols))
            mu = flat.mean(0); sg = flat.std(0); sg[sg == 0] = 1.0
            Z = np.concatenate([np.ones((G, nb, 1)), (Fs - mu) / sg], axis=2)
            Gram = np.einsum("gjp,gjq->jpq", Z, Z)
            Cross = np.einsum("gjp,gi->jip", Z, Y)
            nums, nrms, alf = loto(Gram, Cross, cruz, G, nb)
            obs = resume(nums, nrms, cruz, np.arange(nb))
            a_modal = max(set(alf), key=alf.count)
            nul_agr, nul_med, nul_pos = [], [], []
            for b in range(N_PERM):
                rp = np.random.default_rng(seed + 1500 + b)
                pi = rp.permutation(nb)
                while np.any(pi == np.arange(nb)):
                    pi = rp.permutation(nb)
                n_, p_, _ = loto(Gram, Cross, cruz, G, nb, perm=pi,
                                 alpha_fijo=a_modal)
                rr = resume(n_, p_, cruz, pi)
                nul_agr.append(rr["R_agr"]); nul_med.append(rr["mediana"])
                nul_pos.append(rr["n_positivos"])
            obs["alpha_modal"] = a_modal
            obs["nulo_R_agr"] = nul_agr
            obs["nulo_mediana"] = nul_med
            obs["nulo_n_positivos"] = nul_pos
            ent["conjuntos"][nom] = obs
            sd = float(np.std(nul_agr, ddof=1))
            sys.stderr.write(
                "%-6s %-16s R_agr=%+.4f  mediana=%+.4f  media=%+.4f  pos=%d/%d  "
                "| nulo R_agr=%+.4f sd=%.4f  z=%+.2f  MDE80=%.4f\n"
                % (line, nom, obs["R_agr"], obs["mediana"], obs["media"],
                   obs["n_positivos"], obs["n_validos"],
                   float(np.mean(nul_agr)), sd,
                   (obs["R_agr"] - np.mean(nul_agr)) / sd if sd > 0 else np.nan,
                   2.802 * sd))
            sys.stderr.flush()
        sal[line] = ent
        with open(OUT, "w") as fh:
            json.dump(sal, fh)
    sys.stderr.write("escrito %s\n" % OUT)


if __name__ == "__main__":
    main()
