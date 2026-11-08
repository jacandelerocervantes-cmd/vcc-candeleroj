#!/usr/bin/env python
"""PASO 2 — EMITIR sobre D/E/F. Un solo h5ad, escrito por bloques y REANUDABLE.

LA EMISION, declarada
  * Por contexto se leen sus 18,400 celulas de control (la ENTRADA del modelo).
    De ahi salen los rasgos del gen (CPM, deteccion, Fano, filtro) y la
    profundidad: 🔴 ADAPTATIVO AL CONTEXTO, se recalculan con ese control.
  * Para cada una de las 300 perturbaciones se predice `y(g,h)` con los
    coeficientes de `entrega_ajuste.json`, se RECORTA a [-2, +2] en log2 (razon
    0.25x a 4x) y se convierte en razones `2^y`.
  * Las 400 celulas se toman del pool de control CON REEMPLAZO (semilla
    `s_control(ctx)`) y se reescalan con `cand01_covar._reescala_entero` SIN
    MODIFICAR (semilla `s_ruido(ctx)`), dentro de un `predict_fn(control_pool_X,
    labels, rng)` que respeta el contrato del arnes SIN MODIFICARLO.
  * 🔴 NINGUNA fila de control va a la prediccion. 300 x 400 x 3 = 360,000
    filas, en el orden D, E, F y, dentro de cada contexto, el de
    `pert_counts.csv`.
  * Los 23 blancos sin perfil de gwps: REPLIEGUE DECLARADO -- los cuatro rasgos
    de gwps entran en su valor de NO-INFORMACION (perfil 0, rango 0, signo 0,
    `tiene_gwps_par = 0`), asi que reciben la parte de la funcion que depende
    del gen y del contexto. Es la misma funcion, con el rasgo ausente marcado.

El h5ad se escribe con h5py por bloques (datasets redimensionables): nunca se
tiene la matriz entera en memoria. Si se corta, se reanuda por `progreso.json`.
"""
import json
import os
import sys
import time

import h5py
import numpy as np
import scipy.sparse as sp

sys.path.insert(0, "/home/ubuntu/vcc/modelo/src/noche")
sys.path.insert(0, "/home/ubuntu/vcc/modelo/src")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import comun as C  # noqa: E402
import datos as D  # noqa: E402
import sem2 as S  # noqa: E402
from cand01_covar import _reescala_entero  # noqa: E402  (sin modificar)
from suelo_arnes import decode_categorical  # noqa: E402

WORK = "/datos/modelo/out/entrega"
H5 = os.path.join(WORK, "prediction.h5ad")
PROG = os.path.join(WORK, "progreso.json")
AJUSTE = "/home/ubuntu/vcc/modelo/out/entrega_ajuste.json"
RETO = "/datos/reto/final"
N_POR = 400
RECORTE = 2.0
CTX = ["D", "E", "F"]


def libre():
    import shutil
    return shutil.disk_usage("/datos")[2] / (1024 ** 3)


# ───────────────────── el h5ad por bloques ─────────────────────
def crea(path, n_filas, n_genes):
    with h5py.File(path, "w") as f:
        f.attrs["encoding-type"] = "anndata"
        f.attrs["encoding-version"] = "0.1.0"
        X = f.create_group("X")
        X.attrs["encoding-type"] = "csr_matrix"
        X.attrs["encoding-version"] = "0.1.0"
        X.attrs["shape"] = np.array([n_filas, n_genes], dtype="int64")
        X.create_dataset("data", shape=(0,), maxshape=(None,), dtype="float32",
                         chunks=(1 << 21,))
        X.create_dataset("indices", shape=(0,), maxshape=(None,), dtype="int32",
                         chunks=(1 << 21,))
        X.create_dataset("indptr", shape=(1,), maxshape=(None,), dtype="int64",
                         data=np.zeros(1, dtype="int64"), chunks=(1 << 14,))


def anexa(path, M):
    """Anexa un bloque CSR al final."""
    M = M.tocsr()
    with h5py.File(path, "a") as f:
        X = f["X"]
        d, i, p = X["data"], X["indices"], X["indptr"]
        n0, nnz0 = p.shape[0] - 1, d.shape[0]
        d.resize((nnz0 + M.nnz,)); d[nnz0:] = M.data.astype("float32")
        i.resize((nnz0 + M.nnz,)); i[nnz0:] = M.indices.astype("int32")
        p.resize((n0 + M.shape[0] + 1,))
        p[n0 + 1:] = (M.indptr[1:].astype("int64") + nnz0)
        return int(p.shape[0] - 1), int(d.shape[0])


def cadena(g, nombre, vals):
    dt = h5py.string_dtype(encoding="utf-8")
    ds = g.create_dataset(nombre, data=np.asarray(vals, dtype=object), dtype=dt)
    ds.attrs["encoding-type"] = "string-array"
    ds.attrs["encoding-version"] = "0.2.0"
    return ds


def cierra_obs_var(path, filas_tg, filas_ctx, genes):
    with h5py.File(path, "a") as f:
        for k in ("obs", "var"):
            if k in f:
                del f[k]
        o = f.create_group("obs")
        o.attrs["encoding-type"] = "dataframe"
        o.attrs["encoding-version"] = "0.2.0"
        o.attrs["_index"] = "_index"
        o.attrs["column-order"] = np.asarray(["target_gene", "context"],
                                             dtype=h5py.string_dtype(encoding="utf-8"))
        cadena(o, "_index", [f"c{i}" for i in range(len(filas_tg))])
        cadena(o, "target_gene", filas_tg)
        cadena(o, "context", filas_ctx)
        v = f.create_group("var")
        v.attrs["encoding-type"] = "dataframe"
        v.attrs["encoding-version"] = "0.2.0"
        v.attrs["_index"] = "_index"
        v.attrs["column-order"] = np.asarray([], dtype=h5py.string_dtype(encoding="utf-8"))
        cadena(v, "_index", list(map(str, genes)))


# ───────────────────── el predict_fn del contrato ─────────────────────
def hace_predict_fn(razones):
    """predict_fn(control_pool_X, labels, rng) -> dict[label] = csr (400, G).
    `razones[label]` es el vector de razones por gen. El contrato del arnes no
    se modifica: la funcion solo ve la matriz de control, las etiquetas y rng."""
    def predict_entrega(control_pool_X, labels, rng):
        n = control_pool_X.shape[0]
        out = {}
        for lab in labels:
            ch = rng.integers(0, n, size=N_POR)
            bloque = control_pool_X[ch]
            out[lab] = _reescala_entero(bloque, razones[lab], rng)
        return out
    return predict_entrega


def lee_control(ctx):
    p = os.path.join(RETO, f"context_{ctx}.h5ad")
    with h5py.File(p, "r") as f:
        shape = tuple(int(x) for x in f["X"].attrs["shape"])
        X = sp.csr_matrix((f["X"]["data"][:].astype(np.int64), f["X"]["indices"][:],
                           f["X"]["indptr"][:]), shape=shape)
        tg = np.asarray(decode_categorical(f["obs"]["target_gene"]), dtype=str)
        node = f["var"]["_index"]
        vals = node["values"][:] if hasattr(node, "keys") else node[:]
        vn = np.asarray([v.decode() if isinstance(v, bytes) else v for v in vals], dtype=str)
    assert set(tg.tolist()) == {"non-targeting"}, sorted(set(tg.tolist()))[:3]
    return X, vn


def main():
    os.makedirs(WORK, exist_ok=True)
    aj = json.load(open(AJUSTE))
    co = aj["coeficientes"]
    w = np.asarray([co["w"][r] for r in D.RASGOS])
    mu = np.asarray([co["mu"][r] for r in D.RASGOS])
    sd = np.asarray([co["sd"][r] for r in D.RASGOS])
    b0 = float(co["b0"])
    genes = D.genes_2026()
    perts = D.perts_2026()
    G, P = len(genes), len(perts)
    n_filas = P * N_POR * len(CTX)
    tf = D.tf_set()
    gw_vn, gw_perf, _ = D.gwps(sorted(set(perts.tolist())))

    prog = json.load(open(PROG)) if os.path.exists(PROG) else {
        "bloques": [], "n_filas_escritas": 0, "nnz": 0, "forma_objetivo": [n_filas, G],
        "recorte_log2": RECORTE, "sin_gwps": sorted(set(perts.tolist()) - set(gw_perf)),
        "por_contexto": {}}
    if not os.path.exists(H5) or prog["n_filas_escritas"] == 0:
        crea(H5, n_filas, G)
        prog["bloques"] = []; prog["n_filas_escritas"] = 0; prog["nnz"] = 0
    hechos = {tuple(b) for b in prog["bloques"]}

    for ctx in CTX:
        pend = [p for p in perts if (ctx, str(p)) not in hechos]
        if not pend:
            C.log(f"{ctx}: ya estaba completo"); continue
        if libre() < C.MIN_FREE_GIB:
            raise SystemExit(f"PARAR: {libre():.1f} GiB libres")
        t0 = time.time()
        Xc, vn = lee_control(ctx)
        assert list(vn) == list(genes), "el orden de var del control no es el de gene_names.csv"
        st = D.stats_control(Xc)
        R = D.matriz_rasgos(perts, genes, gw_vn, gw_perf, tf, st)
        Y = (((R - mu) / sd) @ w + b0).reshape(P, G)
        n_rec = int((np.abs(Y) > RECORTE).sum())
        Y = np.clip(Y, -RECORTE, RECORTE)
        razones = {str(p): np.exp2(Y[i]).astype(np.float64) for i, p in enumerate(perts)}
        prog["por_contexto"][ctx] = {
            "n_control": st["n"], "profundidad_mediana": st["profundidad_mediana"],
            "n_genes_filtro": int((st["cpm"] > D.UMBRAL_CPM).sum()),
            "y_sd": float(Y.std()), "y_p1": float(np.percentile(Y, 1)),
            "y_p50": float(np.percentile(Y, 50)), "y_p99": float(np.percentile(Y, 99)),
            "n_recortados": n_rec, "frac_recortados": float(n_rec / Y.size),
            "razon_min": float(np.exp2(Y.min())), "razon_max": float(np.exp2(Y.max())),
            "semilla_control": S.s_control(ctx), "semilla_ruido": S.s_ruido(ctx),
            "dt_rasgos_s": time.time() - t0}
        C.escribe_json(PROG, prog)
        C.log(f"{ctx}: rasgos y prediccion listos en {time.time()-t0:.0f}s | "
              f"y sd {Y.std():.4f} p1 {np.percentile(Y,1):+.3f} p99 {np.percentile(Y,99):+.3f} "
              f"| recortados {n_rec}")
        fn = hace_predict_fn(razones)
        rng_c = np.random.default_rng(S.s_control(ctx))
        rng_r = np.random.default_rng(S.s_ruido(ctx))
        class RNG:                     # el contrato pasa UN rng; se compone
            def __init__(s, a, b): s.a, s.b = a, b
            def integers(s, *A, **K): return s.a.integers(*A, **K)
            def binomial(s, *A, **K): return s.b.binomial(*A, **K)
            def poisson(s, *A, **K): return s.b.poisson(*A, **K)
        rng = RNG(rng_c, rng_r)
        t1 = time.time()
        for j, p in enumerate(perts):
            if (ctx, str(p)) in hechos:
                continue
            M = fn(Xc, [str(p)], rng)[str(p)]
            assert M.shape == (N_POR, G), M.shape
            nf, nnz = anexa(H5, M)
            prog["bloques"].append([ctx, str(p)])
            prog["n_filas_escritas"] = nf
            prog["nnz"] = nnz
            if j % 25 == 0 or j == P - 1:
                C.escribe_json(PROG, prog)
                C.log(f"{ctx}: {j+1}/{P} perts | filas {nf:,} | nnz {nnz:,} | "
                      f"{(time.time()-t1)/max(j+1,1):.2f} s/pert | {libre():.0f} GiB")
                if libre() < C.MIN_FREE_GIB:
                    C.escribe_json(PROG, prog)
                    raise SystemExit(f"PARAR: {libre():.1f} GiB libres")
        C.escribe_json(PROG, prog)
        del Xc, R, Y, razones
        C.log(f"{ctx} COMPLETO en {(time.time()-t0)/60:.1f} min")

    # ── obs y var, al final
    filas_tg, filas_ctx = [], []
    for ctx in CTX:
        for p in perts:
            filas_tg += [str(p)] * N_POR
            filas_ctx += [ctx] * N_POR
    cierra_obs_var(H5, filas_tg, filas_ctx, genes)
    prog["obs_var_escritos"] = True
    prog["tam_bytes"] = os.path.getsize(H5)
    C.escribe_json(PROG, prog)
    C.log(f"h5ad cerrado: {prog['n_filas_escritas']:,} filas, {prog['nnz']:,} nnz, "
          f"{prog['tam_bytes']/2**30:.2f} GiB")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "prueba":
        import anndata as ad
        p = "/tmp/claude-1000/-home-ubuntu-vcc/ab167ee1-8cb2-4956-a8fb-3e6d93e8518a/scratchpad/mini.h5ad"
        crea(p, 6, 4)
        for k in range(3):
            M = sp.csr_matrix(np.array([[1, 0, 2, 0], [0, 3, 0, 4]], dtype=np.float32))
            print("anexado", anexa(p, M))
        cierra_obs_var(p, ["G1"] * 2 + ["G2"] * 2 + ["G3"] * 2, ["D"] * 6,
                       ["a", "b", "c", "d"])
        a = ad.read_h5ad(p)
        print("LEIDO OK", a)
        print(a.X.toarray())
        print(a.obs)
        print(list(a.var_names))
    else:
        main()
