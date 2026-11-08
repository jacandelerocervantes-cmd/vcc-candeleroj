#!/usr/bin/env python
"""Reproducibilidad del COEFICIENTE dentro de un blanco (W) contra entre blancos (E).

EL CONTROL QUE FALTA. control_o_no.md §7 concluye "el mapa es distinto en cada
blanco" porque la consistencia de signo de los coeficientes entre blancos no pasa
del 64 %. Con n=33 eso no excluye una moneda, y hay dos explicaciones
indistinguibles con lo medido ahi:
  (i) los coeficientes DIFIEREN de verdad entre blancos
  (ii) son LOS MISMOS pero se estiman con tanto ruido que el signo baila
Falta la linea base de reproducibilidad DENTRO de un blanco. Es el mismo control
que `A` en b_contra_a.md, construido alli para los deltas y no para los coeficientes.

🔴 A PROFUNDIDAD IGUALADA, que es la leccion de b_contra_a.md. Los ~7000 genes de
cada blanco se parten en dos mitades disjuntas y LAS DOS cantidades se calculan
sobre ajustes de MEDIA MITAD de genes:
  W (dentro) : coef(mitad_A del blanco g)  vs  coef(mitad_B del blanco g)
  E (entre)  : coef(mitad_A del blanco g)  vs  coef(mitad_A del blanco q)   g!=q
El E a genes completos va en su propia columna y declarado, no mezclado con W.

CONDICIONES, las mismas de control_o_no.md:
  * version ESTRICTA: el gen blanco fuera de los rasgos, de las filas y de la
    evaluacion (§4.3). Se hereda el mismo `es_el_blanco` de los rasgos.
  * los mismos 22 rasgos, estandarizados dentro de cada mitad como en §7.
  * 🔴 los 3 rasgos de solo-`g` salen a CERO EXACTO por colinealidad con el
    intercepto: siguen en el ajuste (mismo ridge) y salen del conteo.
  * mismo ridge que `control_o_no_pesos.py`. El `alpha` se elige por CV interna
    sobre genes DENTRO de cada mitad: si difiere entre las dos mitades eso ES
    ruido de estimacion y se reporta. Se da tambien la version con alpha=10 fijo,
    que es el de §7, para comparabilidad.

LA MEDIDA: consistencia de signo pareada (dos vectores, fraccion de rasgos en que
coinciden en signo; moneda = 50 %) y correlacion de Pearson entre los dos vectores
de coeficientes (usa la magnitud). Las dos, por si discrepan.

LAS BARRAS: bootstrap pareado sobre blancos, 200 replicas, como b_contra_a.md §2.
Por replica: nueva particion de genes, remuestreo de los 33 blancos con reemplazo,
y W y E recalculados SOBRE LOS MISMOS. La cantidad guardada es W - E, PAREADA.
Semillas: familia nueva 20269000 + 100*fold + replica.

🔴 DIAGNOSTICO TRANSDUCTIVO. Todo se ajusta sobre la verdad de Jiang. Ninguna cifra
de aqui entra a una submission, ningun coeficiente es una via.
"""
import json
import os
import sys
import time

import numpy as np

LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
DIR = "/home/ubuntu/vcc/modelo/out/control_o_no"
OUT = os.path.join(DIR, "coef_dentro_entre.json")
SEED_FAMILIA = 20269000
ALPHAS = (1.0, 10.0, 100.0, 1000.0, 10000.0)
ALPHA_FIJO = 10.0          # el de §7
N_CV = 4
N_BOOT = 200
SOLO_G = ("g_media", "g_var", "g_frac_cero")


def ajusta(Z, y, folds, alpha_fijo=None):
    """Ridge estandarizado. Devuelve (w, alpha_usado). alpha por CV sobre genes."""
    P = Z.shape[1]
    G_tot = Z.T @ Z
    b_tot = Z.T @ y
    if alpha_fijo is not None:
        a_us = alpha_fijo
    else:
        Gf = [Z[f].T @ Z[f] for f in folds]
        bf = [Z[f].T @ y[f] for f in folds]
        mejor, a_us = -np.inf, ALPHAS[0]
        for a in ALPHAS:
            pred = np.empty_like(y)
            for k, f in enumerate(folds):
                A = G_tot - Gf[k] + np.eye(P) * a
                try:
                    w = np.linalg.solve(A, b_tot - bf[k])
                except np.linalg.LinAlgError:
                    w = np.linalg.lstsq(A, b_tot - bf[k], rcond=None)[0]
                pred[f] = Z[f] @ w
            sp, sy = pred.std(), y.std()
            s = float(np.mean((pred - pred.mean()) * (y - y.mean())) / (sp * sy)) \
                if sp > 0 and sy > 0 else -np.inf
            if s > mejor:
                mejor, a_us = s, a
    A = G_tot + np.eye(P) * a_us
    try:
        w = np.linalg.solve(A, b_tot)
    except np.linalg.LinAlgError:
        w = np.linalg.lstsq(A, b_tot, rcond=None)[0]
    return w, a_us


def coefs(F, E, auto, cols, filas, rng, alpha_fijo=None):
    """Coeficientes por blanco sobre el subconjunto de genes `filas` (bool, G).

    El gen blanco sale de las filas (version ESTRICTA) y `es_el_blanco` no esta
    en `cols`. Estandarizacion y centrado de `y` DENTRO del subconjunto.
    """
    nb = F.shape[1]
    W = np.zeros((nb, len(cols)))
    alf = np.zeros(nb)
    n_fil = np.zeros(nb, dtype=np.int64)
    for j in range(nb):
        m = filas & ~auto[:, j]
        X = F[m][:, j, :][:, cols].astype(np.float64)
        y = E[m, j].astype(np.float64)
        y = y - y.mean()
        mu = X.mean(0); sg = X.std(0); sg[sg == 0] = 1.0
        Z = (X - mu) / sg
        idx = rng.permutation(Z.shape[0])
        folds = np.array_split(idx, N_CV) if alpha_fijo is None else None
        W[j], alf[j] = ajusta(Z, y, folds, alpha_fijo)
        n_fil[j] = Z.shape[0]
    return W, alf, n_fil


def signo_par(a, b):
    """Fraccion de rasgos en que dos vectores de coeficientes coinciden en signo."""
    sa = np.sign(a); sb = np.sign(b)
    return float((sa == sb).mean())


def corr_par(a, b):
    if a.std() == 0 or b.std() == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def dentro(WA, WB):
    """W por blanco: signo y correlacion entre mitad_A y mitad_B del MISMO blanco."""
    s = np.array([signo_par(WA[j], WB[j]) for j in range(WA.shape[0])])
    c = np.array([corr_par(WA[j], WB[j]) for j in range(WA.shape[0])])
    return s, c


def entre(WA, ids):
    """E por par: signo y correlacion entre mitad_A de dos blancos DISTINTOS.

    `ids` es el blanco ORIGINAL de cada fila; los pares con el mismo original se
    excluyen (en bootstrap un blanco puede salir repetido y se compararia consigo
    mismo, que es 100 % por construccion).
    """
    n = WA.shape[0]
    s, c = [], []
    for g in range(n):
        for q in range(g + 1, n):
            if ids[g] == ids[q]:
                continue
            s.append(signo_par(WA[g], WA[q]))
            c.append(corr_par(WA[g], WA[q]))
    return np.array(s), np.array(c)


def por_rasgo_dentro(WA, WB):
    """Por rasgo: fraccion de BLANCOS en que las dos mitades coinciden en signo."""
    return (np.sign(WA) == np.sign(WB)).mean(axis=0)


def por_rasgo_entre(WA, ids):
    n, P = WA.shape
    ac = np.zeros(P); k = 0
    for g in range(n):
        for q in range(g + 1, n):
            if ids[g] == ids[q]:
                continue
            ac += (np.sign(WA[g]) == np.sign(WA[q])).astype(float); k += 1
    return ac / k if k else np.full(P, np.nan)


def mayoria(W):
    """El estadistico de §7: fraccion de blancos con el signo MAYORITARIO."""
    pos = (W > 0).mean(axis=0); neg = (W < 0).mean(axis=0)
    return np.maximum(pos, neg)


def main():
    sal = {"n_boot": N_BOOT, "semilla_familia": SEED_FAMILIA,
           "alphas": list(ALPHAS), "alpha_fijo_s7": ALPHA_FIJO, "lineas": {}}
    for fold, line in enumerate(LINES):
        t0 = time.time()
        d = np.load(os.path.join(DIR, "rasgos_%s.npz" % line), allow_pickle=True)
        F = d["F"]
        E1 = d["E1"].astype(np.float64); E2 = d["E2"].astype(np.float64)
        rasgos = [str(x) for x in d["rasgos"]]
        panel = [str(x) for x in d["panel"]]
        A_esp = d["A_esp"]; cruz = d["cruz"]
        G, nb, _ = F.shape
        E = (E1 + E2) / 2.0
        E = E - E.mean(axis=0, keepdims=True)
        auto = F[:, :, rasgos.index("es_el_blanco")] > 0.5
        cols = [i for i, r in enumerate(rasgos) if r != "es_el_blanco"]
        nom = [rasgos[i] for i in cols]
        cuenta = np.array([r not in SOLO_G for r in nom])     # 18 de 21
        nom_c = [r for r in nom if r not in SOLO_G]

        # ---------- punto: replica 0 ----------
        rng = np.random.default_rng(SEED_FAMILIA + 100 * fold + 0)
        perm = rng.permutation(G)
        mA = np.zeros(G, bool); mA[perm[:G // 2]] = True
        mB = np.zeros(G, bool); mB[perm[G // 2:]] = True
        todo = np.ones(G, bool)
        ids = np.arange(nb)

        res = {}
        for etq, af in (("alpha_cv", None), ("alpha_10", ALPHA_FIJO)):
            WA, aA, nA = coefs(F, E, auto, cols, mA, rng, af)
            WB, aB, nB = coefs(F, E, auto, cols, mB, rng, af)
            WF, aF, nF = coefs(F, E, auto, cols, todo, rng, af)
            KA = WA[:, cuenta]; KB = WB[:, cuenta]; KF = WF[:, cuenta]
            ws, wc = dentro(KA, KB)
            es, ec = entre(KA, ids)
            fs, fc = entre(KF, ids)
            res[etq] = {
                "W_signo": float(ws.mean()), "W_corr": float(np.nanmean(wc)),
                "W_signo_por_blanco": ws.tolist(), "W_corr_por_blanco": wc.tolist(),
                "E_signo": float(es.mean()), "E_corr": float(np.nanmean(ec)),
                "E_signo_mediana": float(np.median(es)),
                "Efull_signo": float(fs.mean()), "Efull_corr": float(np.nanmean(fc)),
                "W_por_rasgo": por_rasgo_dentro(KA, KB).tolist(),
                "E_por_rasgo": por_rasgo_entre(KA, ids).tolist(),
                "Efull_por_rasgo": por_rasgo_entre(KF, ids).tolist(),
                "mayoria_A_por_rasgo": mayoria(KA).tolist(),
                "mayoria_full_por_rasgo": mayoria(KF).tolist(),
                # perfil: corr SOBRE BLANCOS del coeficiente de cada rasgo entre mitades
                "r_perfil_por_rasgo": [corr_par(KA[:, p], KB[:, p])
                                       for p in range(KA.shape[1])],
                "alpha_A": aA.tolist(), "alpha_B": aB.tolist(),
                "alpha_full": aF.tolist(),
                "alpha_discrepa": float((aA != aB).mean()),
                "n_filas_A": nA.tolist(), "n_filas_B": nB.tolist(),
                "n_filas_full": nF.tolist(),
                "abs_medio_A": np.abs(KA).mean(0).tolist(),
                "abs_medio_full": np.abs(KF).mean(0).tolist(),
            }
        sys.stderr.write("%-6s punto  W=%.3f E=%.3f (cv)  W=%.3f E=%.3f (a=10)  %.0fs\n" % (
            line, res["alpha_cv"]["W_signo"], res["alpha_cv"]["E_signo"],
            res["alpha_10"]["W_signo"], res["alpha_10"]["E_signo"], time.time() - t0))
        sys.stderr.flush()

        # ---------- bootstrap pareado sobre blancos ----------
        dif_s, dif_c, bw_s, be_s, bw_c, be_c = [], [], [], [], [], []
        for rep in range(1, N_BOOT + 1):
            rb = np.random.default_rng(SEED_FAMILIA + 100 * fold + rep)
            pr = rb.permutation(G)
            bA = np.zeros(G, bool); bA[pr[:G // 2]] = True
            bB = np.zeros(G, bool); bB[pr[G // 2:]] = True
            WA, _, _ = coefs(F, E, auto, cols, bA, rb, None)
            WB, _, _ = coefs(F, E, auto, cols, bB, rb, None)
            KA = WA[:, cuenta]; KB = WB[:, cuenta]
            sel = rb.integers(0, nb, nb)
            ws, wc = dentro(KA[sel], KB[sel])
            es, ec = entre(KA[sel], sel)
            bw_s.append(float(ws.mean())); be_s.append(float(es.mean()))
            bw_c.append(float(np.nanmean(wc))); be_c.append(float(np.nanmean(ec)))
            dif_s.append(bw_s[-1] - be_s[-1]); dif_c.append(bw_c[-1] - be_c[-1])
            if rep % 50 == 0:
                sys.stderr.write("  %s rep %d  W-E signo %.4f +- %.4f\n" % (
                    line, rep, np.mean(dif_s), np.std(dif_s, ddof=1)))
                sys.stderr.flush()
        for etq, v in (("signo", dif_s), ("corr", dif_c)):
            v = np.array(v)
            res["boot_" + etq] = {
                "media": float(v.mean()), "sd": float(v.std(ddof=1)),
                "ic95": [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))],
                "frac_positivas": float((v > 0).mean())}
        res["boot_W_signo_media"] = float(np.mean(bw_s))
        res["boot_E_signo_media"] = float(np.mean(be_s))
        res["boot_W_corr_media"] = float(np.mean(bw_c))
        res["boot_E_corr_media"] = float(np.mean(be_c))
        res["rasgos_contados"] = nom_c
        res["rasgos_ajuste"] = nom
        res["panel"] = panel
        res["A_esp"] = A_esp.tolist(); res["cruz"] = cruz.tolist()
        res["blanco_expresado"] = d["blanco_expresado"].tolist()
        res["G"] = int(G)
        sal["lineas"][line] = res
        sys.stderr.write("%-6s  W-E signo = %+.4f  IC [%+.4f, %+.4f]   %.0fs\n" % (
            line, res["boot_signo"]["media"], res["boot_signo"]["ic95"][0],
            res["boot_signo"]["ic95"][1], time.time() - t0))
        sys.stderr.flush()
        with open(OUT, "w") as fh:
            json.dump(sal, fh)

    # nulos de los dos estadisticos, exactos, para que los numeros se lean
    from math import comb
    n = 33
    ev = sum(comb(n, k) * max(k, n - k) for k in range(n + 1)) / (2.0 ** n) / n
    sal["nulo_mayoria_n33"] = float(ev)
    sal["nulo_signo_pareado"] = 0.5
    with open(OUT, "w") as fh:
        json.dump(sal, fh)
    sys.stderr.write("escrito %s | nulo mayoria n=33 = %.4f\n" % (OUT, ev))


if __name__ == "__main__":
    main()
