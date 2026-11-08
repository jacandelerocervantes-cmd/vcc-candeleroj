#!/usr/bin/env python
"""PASO 5a PARTE 1: los genes filtrables de cada contexto del reto, SOLO de sus
controles. Superficie de ENTREGA (D/E/F) y de validacion (A/B/C).

CONSTRUCCION, DECLARADA ANTES DE MEDIR
  * Filtro = el del preset `vcc2026`, leido de `EvalConfig`: un gen pasa si el
    CPM MEDIO (media aritmetica sobre celulas de x/L*1e6) en el grupo de
    REFERENCIA supera `filter_gene_min_cpm_cell` (estricto, `>`). En la
    puntuacion real la referencia es el pool de control, asi que se computa de
    los 18,400 controles de cada contexto, todos, sin submuestrear.
  * Comprobacion EN CODIGO: `compute_de` del preset se corre con referencia =
    los 18,400 controles y blanco = 400 controles (semilla 20281001+i); el
    conjunto de genes probados tiene que ser IDENTICO al de esta cuenta.
  * Cantidades por gen y contexto, todas del control: CPM medio, tasa de
    deteccion (fraccion de celulas con cuenta > 0), Fano de las CUENTAS CRUDAS
    (var muestral ddof=1 / media; NaN si media = 0).
  * 1.4 fragilidad: para los genes que pasan en X y no en Y, CPM_Y / 5.0.
    CONTROL: el mismo conteo entre DOS MITADES del MISMO contexto (9,200 vs
    9,200), (a) mitades aleatorias de celulas (semilla 20281021+i) y (b)
    mitades por guia `ntc_id` (las 23 primeras por nombre vs las 23 ultimas),
    deterministas.
"""
import json
import os
import sys
import time

import h5py
import numpy as np
import scipy.sparse as sp

sys.path.insert(0, "/home/ubuntu/vcc/modelo/src/ntc_nulo")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from puerta import corre_puerta, CFG, params  # noqa: E402
import semillas as S  # noqa: E402

RETO = "/datos/reto"
CTX = [("A", "val"), ("B", "val"), ("C", "val"), ("D", "final"), ("E", "final"), ("F", "final")]
UMBRAL = float(CFG.filter.filter_gene_min_cpm_cell)
OUT = "/datos/modelo/out/seleccion/p1_filtrables.json"
N_VERIF = 400


def carga(letra, part):
    p = os.path.join(RETO, part, f"context_{letra}.h5ad")
    with h5py.File(p, "r") as f:
        shape = tuple(int(x) for x in f["X"].attrs["shape"])
        X = sp.csr_matrix((f["X/data"][:], f["X/indices"][:], f["X/indptr"][:]), shape=shape)
        dec = lambda a: np.array([x.decode() if isinstance(x, bytes) else x for x in a])
        vn = dec(f["var/_index/values"][:])
        cats = dec(f["obs/ntc_id/categories"][:])
        cods = f["obs/ntc_id/codes"][:].astype(np.int64)
        tcat = dec(f["obs/target_gene/categories"][:])
    return X, vn, cats, cods, tcat


def stats(X):
    Xw = X.astype(np.float64).tocsr()
    n = Xw.shape[0]
    L = np.asarray(Xw.sum(axis=1)).ravel()
    Lc = L.copy(); Lc[Lc <= 0] = 1.0
    cpm = np.asarray(Xw.T @ ((1e6 / Lc) / n)).ravel()
    det = np.diff(Xw.tocsc().indptr) / n
    mu = np.asarray(Xw.mean(axis=0)).ravel()
    X2 = Xw.copy(); X2.data **= 2
    m2 = np.asarray(X2.mean(axis=0)).ravel()
    var = (m2 - mu ** 2) * n / (n - 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        fano = np.where(mu > 0, var / mu, np.nan)
    return {"cpm": cpm, "det": det, "mu": mu, "fano": fano,
            "L_mediana": float(np.median(L)), "nnz_mediana": float(np.median(np.diff(Xw.indptr))),
            "n": int(n)}


def main():
    t0 = time.time()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    sem = S.comprueba()
    assert sem["sin_choques"], sem
    res = {"preset": params(), "umbral": UMBRAL, "semillas": sem, "ctx": {}}
    vn0 = None
    for i, (L, part) in enumerate(CTX):
        t = time.time()
        X, vn, cats, cods, tcat = carga(L, part)
        assert list(tcat) == ["non-targeting"], tcat
        if vn0 is None:
            vn0 = vn
            gn = np.loadtxt(os.path.join(RETO, part, "gene_names.csv"), dtype=str,
                            skiprows=1, delimiter=",", comments=None)
            res["genes_igual_gene_names_csv"] = bool(np.array_equal(gn, vn0))
        assert np.array_equal(vn, vn0), f"orden de genes distinto en {L}"
        st = stats(X)
        pasa = st["cpm"] > UMBRAL
        # --- comprobacion en codigo contra compute_de del preset
        rng = np.random.default_rng(S.SEMILLAS["p1_verif_filtro"][i])
        it = np.sort(rng.choice(X.shape[0], N_VERIF, replace=False))
        _ns, npr, _ex, tab = corre_puerta(X, X[it], vn, devuelve_tabla=True)
        probados = set(tab["feature"].to_list())
        mios = set(vn[pasa].tolist())
        verif = {"n_probados_compute_de": int(npr), "n_mios": len(mios),
                 "identicos": probados == mios,
                 "solo_compute_de": len(probados - mios), "solo_mios": len(mios - probados),
                 "semilla": S.SEMILLAS["p1_verif_filtro"][i]}
        # --- control de fragilidad: mitades del MISMO contexto
        N = X.shape[0]
        rng2 = np.random.default_rng(S.SEMILLAS["p1_mitades"][i])
        perm = rng2.permutation(N)
        h1, h2 = np.sort(perm[: N // 2]), np.sort(perm[N // 2:])
        orden = np.argsort(cats)
        g1 = set(orden[:23].tolist())
        k1 = np.where(np.isin(cods, list(g1)))[0]
        k2 = np.where(~np.isin(cods, list(g1)))[0]
        mit = {}
        for etq, a, b in (("aleatoria", h1, h2), ("por_guia", k1, k2)):
            ca, cb = stats(X[a])["cpm"], stats(X[b])["cpm"]
            mit[etq] = {"n_a": int(a.size), "n_b": int(b.size),
                        "cpm_a": ca.tolist(), "cpm_b": cb.tolist(),
                        "pasa_a": int((ca > UMBRAL).sum()), "pasa_b": int((cb > UMBRAL).sum())}
        res["ctx"][L] = {"particion": part, "n": st["n"], "L_mediana": st["L_mediana"],
                         "nnz_mediana": st["nnz_mediana"], "n_pasa": int(pasa.sum()),
                         "cpm": st["cpm"].tolist(), "det": st["det"].tolist(),
                         "fano": [None if not np.isfinite(x) else float(x) for x in st["fano"]],
                         "genes_pasan": vn[pasa].tolist(), "verif": verif, "mitades": mit,
                         "guias_mitad_a": sorted(cats[orden[:23]].tolist()),
                         "dt_s": time.time() - t}
        sys.stderr.write(f"[{L}] pasa={int(pasa.sum())} verif={verif['identicos']} "
                         f"({verif['n_probados_compute_de']}) dt={time.time()-t:.0f}s\n")
        del X
    res["genes"] = vn0.tolist()
    res["dt_total_s"] = time.time() - t0
    tmp = OUT + ".tmp"
    json.dump(res, open(tmp, "w"))
    os.replace(tmp, OUT)


if __name__ == "__main__":
    main()
