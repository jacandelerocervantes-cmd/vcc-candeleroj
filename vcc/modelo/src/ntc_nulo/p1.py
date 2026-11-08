#!/usr/bin/env python
"""PARTE 1 (y los insumos de la PARTE 2): el NULO de ntc_id, por contexto.

LA CONSTRUCCION, declarada antes de medir
-----------------------------------------
Cada context_{A..F}.h5ad trae 18,400 celulas NO DIRIGIDAS, 46 guias `ntc_id` de
400 celulas cada una. El nulo se construye partiendo las 46 GUIAS en dos mitades
de 23, y comparando mitad contra mitad con la MISMA puerta del preset vcc2026
(`puerta.py`, que lee cada parametro del propio preset).

  * LA UNIDAD ES LA GUIA, NO LA CELULA. Las 400 celulas de una guia van enteras
    a un lado. Ninguna particion parte una guia -- se verifica en cada particion
    (`guia_partida=false`, conjuntos disjuntos, union = las 46).
  * SEMILLA: `20274000 + 100*idx_contexto + j`, con idx A=0 B=1 C=2 D=3 E=4 F=5
    y j = 0..7 (ocho particiones por contexto).
  * `rng.permutation(46)` sobre la lista de las 46 categorias ORDENADA
    alfabeticamente; los primeros 23 son la mitad de REFERENCIA, los ultimos 23
    la mitad BLANCO.

TRES BRAZOS DE `n`, declarados aparte
-------------------------------------
  1. `completo` : 9,200 vs 9,200 -- las 23 guias enteras de cada lado.
  2. `n400_400` : 400 vs 400 -- el `n` POR PERTURBACION del reto, a los dos
     lados. Es el brazo que pide el encargo.
  3. `n400_ref_completa` : 400 (blanco) vs 9,200 (referencia) -- AÑADIDO y
     declarado. El reto no compara 400 contra 400: compara las 400 celulas de
     una perturbacion contra el pool de control del contexto
     (`manifest.json`: ground_truth_cells = 138,400 = 300x400 + 18,400). Este
     brazo es el que se parece al regimen de la entrega; el 2 es el que pide el
     encargo. Los tres se reportan por separado y no se mezclan.
  El submuestreo a 400 se hace SIN REEMPLAZO sobre las 9,200 celulas de la
  mitad, con el mismo `rng` y en orden declarado (referencia primero, blanco
  despues). No parte guias: las guias siguen enteras de un lado.

QUE SE GUARDA. Por corrida: n_sig, n_probados, p_min, cuantiles de p, |log2FC|,
y la lista de genes llamados. Del brazo 2 y 3 se guardan ademas los vectores
crudos (p_value, log2FC, CPM medio de la referencia) para el piso de deteccion
de la Parte 2.

NO se toca el lado perturbado -- no existe en 2026. No se lee nada de
~/vcc/benchmark/ ni de /datos/benchmark_VEDADO_*.
"""
import json, os, sys, time
import numpy as np
import scipy.sparse as sp
import h5py
import polars as pl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from puerta import corre_puerta, params, CFG

RETO = "/datos/reto"
CTX = [("A", "val", 0), ("B", "val", 1), ("C", "val", 2),
       ("D", "final", 3), ("E", "final", 4), ("F", "final", 5)]
SEM_BASE = 20274000
N_PART = 8
N_ENTREGA = 400
OUTDIR = "/home/ubuntu/vcc/modelo/out/ntc_nulo"
SOLO = os.environ.get("NULO_CTX")          # p.ej. "D" para una sola corrida


def carga(letra, part):
    p = os.path.join(RETO, part, f"context_{letra}.h5ad")
    with h5py.File(p, "r") as f:
        shape = tuple(int(x) for x in f["X"].attrs["shape"])
        X = sp.csr_matrix((f["X/data"][:], f["X/indices"][:], f["X/indptr"][:]),
                          shape=shape)
        vn = np.array([x.decode() if isinstance(x, bytes) else x
                       for x in f["var/_index/values"][:]])
        cats = np.array([x.decode() if isinstance(x, bytes) else x
                         for x in f["obs/ntc_id/categories"][:]])
        cods = f["obs/ntc_id/codes"][:].astype(np.int64)
    return X, vn, cats, cods


def cpm_medio(X):
    """Media aritmetica, sobre estas celulas, del CPM por celula x/L*1e6.
    Es la cantidad con la que `_apply_cpm_filter` decide el conjunto de genes."""
    Xw = X.astype(np.float64)
    L = np.asarray(Xw.sum(axis=1)).ravel()
    L[L <= 0] = 1.0
    return np.asarray(Xw.T @ ((1e6 / L) / X.shape[0])).ravel()


def resumen_p(tabla):
    p = tabla["p_value"].to_numpy().astype(float)
    lf = np.abs(tabla["log2_fold_change"].to_numpy().astype(float))
    pv = p[~np.isnan(p)]
    m = pv.size
    # KS contra uniforme, calculado aqui (es un diagnostico del nulo, no la puerta)
    ps = np.sort(pv)
    i = np.arange(1, m + 1)
    ks = float(max(np.max(i / m - ps), np.max(ps - (i - 1) / m))) if m else None
    return {
        "n_probados": int(tabla.height), "n_p_validos": int(m),
        "p_min": float(ps[0]) if m else None,
        "p_min_esperado_uniforme": float(1.0 / (m + 1)) if m else None,
        "frac_p_menor_05": float(np.mean(pv < 0.05)) if m else None,
        "frac_p_menor_01": float(np.mean(pv < 0.01)) if m else None,
        "frac_p_menor_001": float(np.mean(pv < 0.001)) if m else None,
        "p_mediana": float(np.median(pv)) if m else None,
        "ks_uniforme": ks,
        "abs_lfc_mediana": float(np.nanmedian(lf)),
        "abs_lfc_p99": float(np.nanpercentile(lf, 99)),
        "abs_lfc_max": float(np.nanmax(lf)),
        "n_lfc_nan": int(np.isnan(tabla["log2_fold_change"].to_numpy()).sum()),
        "umbral_bh_rango1": float(CFG.de.p_adj_threshold / max(m, 1)),
    }


def una_corrida(X, vn, idx_ref, idx_tgt, etiqueta, guarda=None):
    t0 = time.time()
    ns, npr, ex, tabla = corre_puerta(X[idx_ref], X[idx_tgt], vn,
                                      devuelve_tabla=True)
    r = {"brazo": etiqueta, "n_ref": int(idx_ref.size), "n_tgt": int(idx_tgt.size),
         "n_sig": ns, "segundos": round(time.time() - t0, 1)}
    r.update(resumen_p(tabla))
    sig = tabla.filter(pl.col("p_adj") < CFG.de.p_adj_threshold)
    r["genes_sig"] = sig["feature"].to_list()
    r["p_adj_min"] = float(np.nanmin(tabla["p_adj"].to_numpy()))
    if guarda is not None:
        rc = cpm_medio(X[idx_ref])
        pos = {g: i for i, g in enumerate(vn)}
        gi = np.array([pos[g] for g in tabla["feature"].to_list()])
        np.savez_compressed(
            guarda,
            feature_idx=gi.astype(np.int32),
            p_value=tabla["p_value"].to_numpy().astype(np.float64),
            p_adj=tabla["p_adj"].to_numpy().astype(np.float64),
            log2fc=tabla["log2_fold_change"].to_numpy().astype(np.float64),
            ref_cpm=rc[gi].astype(np.float64),
        )
    return r


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    salida = {"puerta": params(), "construccion": {
        "semilla": "20274000 + 100*idx_contexto + j, j=0..7",
        "idx_contexto": {l: k for l, _, k in CTX},
        "n_particiones": N_PART, "guias_por_mitad": 23,
        "unidad": "la guia ntc_id entera (400 celulas), nunca la celula",
        "brazos": ["completo 9200v9200", "n400_400 400v400",
                   "n400_ref_completa 400v9200 (anadido, declarado)"],
    }, "corridas": [], "particiones": []}

    for letra, part, k in CTX:
        if SOLO and letra != SOLO:
            continue
        sys.stderr.write(f"\n===== contexto {letra} ({part}) =====\n")
        X, vn, cats, cods = carga(letra, part)
        orden = np.argsort(cats)                 # categorias ordenadas alfabeticamente
        cats_ord = cats[orden]
        for j in range(N_PART):
            sem = SEM_BASE + 100 * k + j
            rng = np.random.default_rng(sem)
            perm = rng.permutation(46)
            g1 = set(orden[perm[:23]].tolist())   # codigos de la mitad REFERENCIA
            g2 = set(orden[perm[23:]].tolist())
            i1 = np.where(np.isin(cods, list(g1)))[0]
            i2 = np.where(np.isin(cods, list(g2)))[0]
            info = {
                "contexto": letra, "j": j, "semilla": sem,
                "n_guias_ref": len(g1), "n_guias_tgt": len(g2),
                "guias_ref": sorted(cats[list(g1)].tolist()),
                "guias_tgt": sorted(cats[list(g2)].tolist()),
                "disjuntas": len(g1 & g2) == 0,
                "union_46": len(g1 | g2) == 46,
                "celulas_ref": int(i1.size), "celulas_tgt": int(i2.size),
                "guia_partida": False,
                "celulas_por_guia_ref_ok": bool(
                    all(int((cods == c).sum()) == 400 for c in g1)),
                "celulas_por_guia_tgt_ok": bool(
                    all(int((cods == c).sum()) == 400 for c in g2)),
            }
            assert info["disjuntas"] and info["union_46"]
            assert i1.size == 9200 and i2.size == 9200
            salida["particiones"].append(info)

            s400_ref = np.sort(rng.choice(i1, N_ENTREGA, replace=False))
            s400_tgt = np.sort(rng.choice(i2, N_ENTREGA, replace=False))

            for etq, ir, it, gfile in (
                ("completo", i1, i2, None),
                ("n400_400", s400_ref, s400_tgt,
                 os.path.join(OUTDIR, f"crudo_{letra}_{j}_n400_400.npz")),
                ("n400_ref_completa", i1, s400_tgt,
                 os.path.join(OUTDIR, f"crudo_{letra}_{j}_n400_ref.npz")),
            ):
                r = una_corrida(X, vn, ir, it, etq, guarda=gfile)
                r.update({"contexto": letra, "j": j, "semilla": sem,
                          "particion": "aleatoria_de_guias"})
                salida["corridas"].append(r)
                sys.stderr.write(
                    f"  [{letra} j={j} {etq:18s}] n_sig={r['n_sig']:4d} "
                    f"probados={r['n_probados']} p_min={r['p_min']:.3e} "
                    f"({r['segundos']}s)\n")
        del X
        with open(os.path.join(OUTDIR, "p1.json"), "w") as fh:
            json.dump(salida, fh, indent=1, sort_keys=True, default=str)

    with open(os.path.join(OUTDIR, "p1.json"), "w") as fh:
        json.dump(salida, fh, indent=1, sort_keys=True, default=str)
    sys.stderr.write(f"\nlisto: {OUTDIR}/p1.json\n")


if __name__ == "__main__":
    main()
