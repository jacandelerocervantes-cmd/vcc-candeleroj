#!/usr/bin/env python
"""Ajuste y validacion cruzada POR BLANCO. Discrimina (a) de (b).

🔴 TRANSDUCTIVO Y PROHIBIDO COMO VIA. Todo lo que se ajusta aqui se ajusta sobre
la verdad de Jiang. Ningun coeficiente sale de este archivo hacia nada.

TRES NUMEROS POR BLANCO, y responden a TRES preguntas distintas que no hay que
mezclar:

  r_loto   -- se ajusta en los otros 32 blancos y se predice el que se dejo
              fuera. ES EL TITULAR: "¿generaliza a un blanco NUEVO?"
  r_mismo  -- se ajusta en la MITAD de los genes DEL PROPIO blanco y se predice
              la otra mitad. NO es generalizacion a un blanco nuevo y no se
              puede leer como tal; es el DISCRIMINADOR: "¿llevan los rasgos
              alguna informacion sobre el delta especifico de este blanco?"
  A_esp    -- fiabilidad mitad-contra-mitad del propio delta especifico. El
              TECHO: nada puede correlacionar con la senal mas de lo que la
              senal correlaciona consigo misma.

  r_loto ~ 0  y  r_mismo ~ 0   -> (b) NADA del control lo predice
  r_loto ~ 0  y  r_mismo > 0   -> la senal esta en los rasgos pero no transfiere
                                   entre blancos con n = 33: no se puede decidir
  r_loto > 0                   -> (a) otra senal del control SI lo capta

LA METRICA, con el ruido fuera por producto cruzado de mitades. Para la
prediccion `p` del blanco j, con `E1`, `E2` las dos mitades del delta especifico
y `t` la senal verdadera (desconocida):
    E[<E1_j, E2_j>] = ||t_j||^2          (el ruido de las dos mitades es
    E[<p, (E1_j+E2_j)/2>] = <p, t_j>      independiente, asi que no entra)
    r_true = <p, (E1+E2)/2> / (||p|| * sqrt(<E1_j, E2_j>))
Se reporta tambien `r_raw`, la Pearson contra `(E1+E2)/2` sin corregir, que es
la que esta atenuada. `r_true` no esta definida si `<E1_j,E2_j> <= 0`.

EL NULO, por PERMUTACION y no por una formula: se permutan las ETIQUETAS de
blanco del objetivo (los rasgos del blanco j se emparejan con el `E` del blanco
pi(j)) y se repite la misma validacion cruzada entera. Eso da la distribucion de
`media(r)` bajo "nada predice" conservando toda la estructura de los rasgos y la
dependencia entre genes.

Lectura y medicion. No construye ningun candidato.
"""
import json
import os
import sys
import time

import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor

LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
DIR = "/home/ubuntu/vcc/modelo/out/control_o_no"
OUT = os.path.join(DIR, "ajuste.json")
SEED_FAMILIA = 20268000
ALPHAS = (1.0, 10.0, 100.0, 1000.0, 10000.0)
N_PERM_RIDGE = 200
N_PERM_GB = 8
N_PERM_MISMO = 30
GB_KW = dict(max_iter=200, learning_rate=0.1)

CONJUNTOS = {
    "todos": None,
    "sin_autoblanco": ["es_el_blanco"],
    "solo_covarianza": None,     # se fija abajo con la lista positiva
    "solo_h": None,
}
POS_COV = ["cov", "r_pearson", "r_spearman", "abs_r_pearson", "abs_r_spearman"]
POS_H = ["h_media", "h_var", "h_disp", "h_frac_cero", "h_rango_expr"]


def r_true(p, y, cruz):
    """Correlacion con la senal VERDADERA. `y` ya viene centrada por columna."""
    pc = p - p.mean()
    np_ = np.linalg.norm(pc)
    if np_ == 0 or cruz <= 0:
        return np.nan
    return float(pc @ y / (np_ * np.sqrt(cruz)))


def r_raw(p, y):
    pc = p - p.mean(); yc = y - y.mean()
    d = np.linalg.norm(pc) * np.linalg.norm(yc)
    return float(pc @ yc / d) if d > 0 else np.nan


def bloques(Z, Y):
    """Z: (G, nb, P+1) con intercepto en la columna 0. Y: (G, nb) CENTRADA por
    columna. Devuelve Gram[j] = Z_j^T Z_j, Cross[j,i] = Z_j^T Y_i para TODO par,
    y yy[i] = ||Y_i||^2.

    🔴 Con estos bloques la metrica NO necesita volver a tocar los G genes, y por
    eso las permutaciones salen gratis. Para `w` y el blanco j contra el objetivo
    i, con `p = Z_j w`:
        <p - media(p), Y_i> = w . Cross[j,i]        porque Y_i esta centrada
        ||p - media(p)||^2  = w^T Gram[j] w - (w . s_j)^2 / G
    donde `s_j = Z_j^T 1 = Gram[j][0, :]`, que ya esta dentro de Gram porque la
    columna 0 de Z es el intercepto. Exacto, no aproximado."""
    G, nb, P = Z.shape
    Gram = np.einsum("gjp,gjq->jpq", Z, Z)
    Cross = np.einsum("gjp,gi->jip", Z, Y)
    yy = np.einsum("gi,gi->i", Y, Y)
    return Gram, Cross, yy


def resuelve(XtX, Xty, alpha):
    P = XtX.shape[0]
    reg = np.eye(P) * alpha
    reg[0, 0] = 0.0                      # el intercepto no se penaliza
    try:
        return np.linalg.solve(XtX + reg, Xty)
    except np.linalg.LinAlgError:
        return np.linalg.lstsq(XtX + reg, Xty, rcond=None)[0]


def metrica_bloque(w, j, i, Gram, Cross, yy, cruz, G):
    """(r_true, r_raw) para la prediccion del blanco j contra el objetivo i,
    calculada solo con los bloques. Equivalente exacto a hacerlo sobre los G
    genes."""
    num = float(w @ Cross[j, i])
    v = float(w @ Gram[j] @ w) - float(w @ Gram[j][0]) ** 2 / G
    if v <= 0:
        return np.nan, np.nan
    nrm = np.sqrt(v)
    rt = num / (nrm * np.sqrt(cruz[i])) if cruz[i] > 0 else np.nan
    rr = num / (nrm * np.sqrt(yy[i])) if yy[i] > 0 else np.nan
    return rt, rr


def ridge_loto(Gram, Cross, yy, cruz, G, nb, perm=None, alphas=ALPHAS,
               alpha_fijo=None):
    """LOTO en espacio de bloques. `perm` empareja los rasgos del blanco j con el
    objetivo de perm[j]. Con `alpha_fijo` se salta la CV interna de alpha."""
    pi = np.arange(nb) if perm is None else np.asarray(perm)
    rs, rr, alf = [], [], []
    todos = np.arange(nb)
    for j in range(nb):
        tr = todos[todos != j]
        if alpha_fijo is not None:
            mejor_a = alpha_fijo
        else:
            grupos = np.array_split(tr, 4)
            mejor, mejor_a = -np.inf, alphas[0]
            for a in alphas:
                sc = []
                for gr in grupos:
                    sub = np.setdiff1d(tr, gr, assume_unique=False)
                    w = resuelve(Gram[sub].sum(0), Cross[sub, pi[sub]].sum(0), a)
                    for k in gr:
                        t, _ = metrica_bloque(w, k, pi[k], Gram, Cross, yy,
                                              cruz, G)
                        if np.isfinite(t):
                            sc.append(t)
                m = np.mean(sc) if sc else -np.inf
                if m > mejor:
                    mejor, mejor_a = m, a
        w = resuelve(Gram[tr].sum(0), Cross[tr, pi[tr]].sum(0), mejor_a)
        t, v = metrica_bloque(w, j, pi[j], Gram, Cross, yy, cruz, G)
        rs.append(t); rr.append(v); alf.append(mejor_a)
    return np.array(rs), np.array(rr), np.array(alf)


def gb_loto(F, Y, cruz, seed, perm=None, importancia=False, rasgos=None):
    G, nb, P = F.shape
    pi = np.arange(nb) if perm is None else perm
    rs, rr = [], []
    imp = np.zeros((nb, P)) if importancia else None
    rng = np.random.default_rng(seed)
    for j in range(nb):
        tr = np.array([i for i in range(nb) if i != j])
        Xtr = F[:, tr, :].reshape(-1, P)
        ytr = Y[:, pi[tr]].reshape(-1, order="F") if False else \
              np.concatenate([Y[:, pi[i]] for i in tr])
        m = HistGradientBoostingRegressor(random_state=seed + j, **GB_KW)
        m.fit(Xtr, ytr)
        Xte = F[:, j, :]
        p = m.predict(Xte)
        base = r_true(p, Y[:, pi[j]], cruz[pi[j]])
        rs.append(base); rr.append(r_raw(p, Y[:, pi[j]]))
        if importancia and np.isfinite(base):
            for q in range(P):
                Xp = Xte.copy()
                Xp[:, q] = Xp[rng.permutation(G), q]
                imp[j, q] = base - r_true(m.predict(Xp), Y[:, pi[j]], cruz[pi[j]])
    return np.array(rs), np.array(rr), imp


def mismo_blanco(F, E1, E2, seed, n_perm=N_PERM_MISMO):
    """Ajusta en la MITAD de los genes del propio blanco y predice la otra. NO es
    generalizacion a un blanco nuevo. Con su nulo por permutacion de genes."""
    G, nb, P = F.shape
    rng = np.random.default_rng(seed)
    mitad = rng.permutation(G) < (G // 2)
    out = {"gb": [], "ridge": [], "gb_nulo": [], "ridge_nulo": []}
    for j in range(nb):
        X = F[:, j, :]
        e1 = E1[:, j].astype(np.float64); e2 = E2[:, j].astype(np.float64)
        y = (e1 + e2) / 2.0
        obs_g, obs_r = [], []
        for sel in (mitad, ~mitad):
            tr, te = ~sel, sel
            yc1 = e1[te] - e1[te].mean(); yc2 = e2[te] - e2[te].mean()
            cz = float(yc1 @ yc2)
            yt = y[te] - y[te].mean()
            g = HistGradientBoostingRegressor(random_state=seed + j, **GB_KW)
            g.fit(X[tr], y[tr] - y[tr].mean())
            obs_g.append(r_true(g.predict(X[te]), yt, cz))
            Zt = np.column_stack([np.ones(tr.sum()), X[tr]])
            Ze = np.column_stack([np.ones(te.sum()), X[te]])
            mu = Zt.mean(0); mu[0] = 0.0
            sg = Zt.std(0); sg[sg == 0] = 1.0; sg[0] = 1.0
            w = resuelve(((Zt - mu) / sg).T @ ((Zt - mu) / sg),
                         ((Zt - mu) / sg).T @ (y[tr] - y[tr].mean()), 10.0)
            obs_r.append(r_true(((Ze - mu) / sg) @ w, yt, cz))
        out["gb"].append(float(np.nanmean(obs_g)))
        out["ridge"].append(float(np.nanmean(obs_r)))
    # nulo: permutar los genes del objetivo dentro del blanco
    for b in range(n_perm):
        rp = np.random.default_rng(seed + 1000 + b)
        ng, nr = [], []
        for j in range(nb):
            X = F[:, j, :]
            o = rp.permutation(G)
            e1 = E1[o, j].astype(np.float64); e2 = E2[o, j].astype(np.float64)
            y = (e1 + e2) / 2.0
            sel = mitad; tr, te = ~sel, sel
            yc1 = e1[te] - e1[te].mean(); yc2 = e2[te] - e2[te].mean()
            cz = float(yc1 @ yc2); yt = y[te] - y[te].mean()
            g = HistGradientBoostingRegressor(random_state=seed + j, **GB_KW)
            g.fit(X[tr], y[tr] - y[tr].mean())
            ng.append(r_true(g.predict(X[te]), yt, cz))
        out["gb_nulo"].append(float(np.nanmean(ng)))
    return out


def main():
    lineas = sys.argv[1].split(",") if len(sys.argv) > 1 else LINES
    n_perm_gb = int(sys.argv[2]) if len(sys.argv) > 2 else N_PERM_GB
    salida = {}
    if os.path.exists(OUT):
        salida = json.load(open(OUT))
    for line in lineas:
        t0 = time.time()
        fold = LINES.index(line)
        seed = SEED_FAMILIA + 100 * fold
        d = np.load(os.path.join(DIR, "rasgos_%s.npz" % line), allow_pickle=True)
        F = d["F"].astype(np.float64); E1 = d["E1"]; E2 = d["E2"]
        rasgos = [str(x) for x in d["rasgos"]]
        cruz = d["cruz"].astype(np.float64)
        G, nb, P = F.shape
        Y = ((E1 + E2) / 2.0).astype(np.float64)
        Y = Y - Y.mean(axis=0, keepdims=True)

        res = {"G": int(G), "n_blancos": int(nb),
               "A_esp": d["A_esp"].tolist(), "cruz": cruz.tolist(),
               "blanco_expresado": d["blanco_expresado"].tolist(),
               "rasgos": rasgos, "conjuntos": {}}

        conj = {"todos": list(range(P)),
                "sin_autoblanco": [i for i, r in enumerate(rasgos)
                                   if r != "es_el_blanco"],
                "solo_covarianza": [i for i, r in enumerate(rasgos)
                                    if r in POS_COV],
                "solo_h": [i for i, r in enumerate(rasgos) if r in POS_H]}

        for nom, cols in conj.items():
            Fs = F[:, :, cols]
            # escalado con estadisticas de los RASGOS solamente: no entra ningun
            # valor del objetivo, asi que no hay fuga
            flat = Fs.reshape(-1, len(cols))
            mu = flat.mean(0); sg = flat.std(0); sg[sg == 0] = 1.0
            Zs = (Fs - mu) / sg
            Z = np.concatenate([np.ones((G, nb, 1)), Zs], axis=2)
            Gram, Cross, yy = bloques(Z, Y)
            rs, rr, alf = ridge_loto(Gram, Cross, yy, cruz, G, nb)
            ent = {"ridge_loto_true": rs.tolist(), "ridge_loto_raw": rr.tolist(),
                   "ridge_alpha": alf.tolist()}
            # El nulo usa el alpha MODAL elegido sobre los datos reales en vez de
            # repetir la CV interna en cada permutacion. Queda declarado: afinar
            # alpha dentro del nulo lo haria MAS optimista (mas libertad bajo el
            # nulo), asi que fijarlo es la direccion conservadora para el nulo y
            # anticonservadora para nada.
            a_modal = ALPHAS[int(np.bincount(
                [list(ALPHAS).index(a) for a in alf]).argmax())]
            ent["ridge_alpha_modal"] = a_modal
            nul = []
            for b in range(N_PERM_RIDGE):
                rp = np.random.default_rng(seed + 500 + b)
                pi = rp.permutation(nb)
                while np.any(pi == np.arange(nb)) and nb > 2:
                    pi = rp.permutation(nb)
                v, _, _ = ridge_loto(Gram, Cross, yy, cruz, G, nb, perm=pi,
                                     alpha_fijo=a_modal)
                nul.append(float(np.nanmean(v)))
            ent["ridge_nulo_media"] = nul
            res["conjuntos"][nom] = ent
            sys.stderr.write(
                "  %-16s ridge_loto media=%+.4f  nulo media=%+.4f sd=%.4f\n"
                % (nom, float(np.nanmean(rs)), float(np.mean(nul)),
                   float(np.std(nul, ddof=1))))
            sys.stderr.flush()

        # --- GB sobre el conjunto completo y sin autoblanco ---
        for nom in ("todos", "sin_autoblanco"):
            cols = conj[nom]
            Fs = F[:, :, cols]
            gs, gr, imp = gb_loto(Fs, Y, cruz, seed, importancia=(nom == "todos"))
            res["conjuntos"][nom]["gb_loto_true"] = gs.tolist()
            res["conjuntos"][nom]["gb_loto_raw"] = gr.tolist()
            if imp is not None:
                res["conjuntos"][nom]["gb_importancia"] = \
                    np.nanmean(imp, axis=0).tolist()
                res["conjuntos"][nom]["gb_importancia_rasgos"] = \
                    [rasgos[i] for i in cols]
            nul = []
            for b in range(n_perm_gb if nom == "todos" else 0):
                rp = np.random.default_rng(seed + 700 + b)
                pi = rp.permutation(nb)
                while np.any(pi == np.arange(nb)) and nb > 2:
                    pi = rp.permutation(nb)
                v, _, _ = gb_loto(Fs, Y, cruz, seed, perm=pi)
                nul.append(float(np.nanmean(v)))
            res["conjuntos"][nom]["gb_nulo_media"] = nul
            sys.stderr.write(
                "  %-16s gb_loto    media=%+.4f  nulo media=%+.4f sd=%.4f (%d perm)\n"
                % (nom, float(np.nanmean(gs)), float(np.mean(nul)),
                   float(np.std(nul, ddof=1)), len(nul)))
            sys.stderr.flush()

        res["mismo_blanco"] = mismo_blanco(F, E1, E2, seed + 50)
        sys.stderr.write(
            "  mismo_blanco     gb media=%+.4f  nulo media=%+.4f sd=%.4f\n"
            % (float(np.nanmean(res["mismo_blanco"]["gb"])),
               float(np.mean(res["mismo_blanco"]["gb_nulo"])),
               float(np.std(res["mismo_blanco"]["gb_nulo"], ddof=1))))
        res["segundos"] = round(time.time() - t0, 1)
        salida[line] = res
        with open(OUT, "w") as fh:
            json.dump(salida, fh)
        sys.stderr.write("%-6s listo en %.0fs\n\n" % (line, time.time() - t0))
        sys.stderr.flush()


if __name__ == "__main__":
    main()
