#!/usr/bin/env python
"""Los datos y los RASGOS de la entrega. Nada aqui usa el gen como llave.

FUENTES, declaradas
  * VCC 2025: los TRES splits (train 150 + validation 50 + test 100 = 300
    perturbaciones CON VERDAD), 18,080 genes. Se usa para APRENDER LA FUNCION.
  * k562_gwps bulk (`/datos/ref/replogle/k562_gwps_raw_bulk.h5ad`): perfil de
    respuesta por par (blanco, gen). 11,258 filas x 8,248 genes; 514 filas
    `core_control`. El perfil es `log2((x + eps)/(control + eps))`, el mismo
    que `gwps_jiang_datos.py`. Cubre 277 de las 300 finales y 7,681 de los
    18,533 genes de 2026.
  * TFactprofiler (`/datos/ref/tfact/TFactprofiler_all.csv`): el conjunto de
    `source` son 2,017 factores de transcripcion; 80 de las 300 lo son.
  * Controles de D/E/F: los 18,400 de cada contexto. 🔴 Los rasgos del gen y la
    profundidad SE RECALCULAN con el control del contexto que toque.

🔴 LA VERDAD DEL GEN NO VIAJA: de 2025 solo sale la funcion (los coeficientes).
Ningun delta de 2025, ni de los 28 blancos compartidos, entra en la emision.
"""
import json
import os
import sys

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, "/home/ubuntu/vcc/modelo/src")
sys.path.insert(0, "/home/ubuntu/vcc/modelo/src/noche")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from suelo_arnes import decode_categorical  # noqa: E402
import comun as C  # noqa: E402

WORK = "/datos/modelo/out/entrega"
GWPS = "/datos/ref/replogle/k562_gwps_raw_bulk.h5ad"
TFACT = "/datos/ref/tfact/TFactprofiler_all.csv"
RETO = "/datos/reto/final"
V25 = "/datos/vcc2025"
SPLITS = {"train": "train/adata_Training.h5ad",
          "validation": "validation/adata_Validation.h5ad",
          "test": "test/adata_Test.h5ad"}
CONTROL = "non-targeting"
EPS = 1e-8
UMBRAL_CPM = 5.0
RASGOS = ["gwps_perfil", "gwps_rango", "gwps_signo", "tiene_gwps_par",
          "log10_cpm_h", "det_h", "log_fano_h", "pasa_filtro_h",
          "g_es_TF", "log10_cpm_g"]


def log(m):
    C.log(m)


def genes_2026():
    d = pd.read_csv(os.path.join(RETO, "gene_names.csv"))
    col = "gene_name" if "gene_name" in d.columns else d.columns[0]
    return np.asarray(d[col].astype(str).tolist(), dtype=str)


def perts_2026():
    return np.asarray(pd.read_csv(os.path.join(RETO, "pert_counts.csv"))
                      ["target_gene"].astype(str).tolist(), dtype=str)


def genes_2025():
    return np.asarray(pd.read_csv(os.path.join(V25, "gene_names.csv"), header=None)[0]
                      .astype(str).tolist(), dtype=str)


def tf_set():
    s = pd.read_csv(TFACT, usecols=["source"])["source"].astype(str)
    return set(s.unique().tolist())


# ───────────────────────── 2025: la verdad para aprender ─────────────────────
def pseudobulk_2025(bloque=40000):
    """Media de CPM por (blanco, gen) sobre los tres splits, y las sumas del
    control para sus rasgos. CPM por celula = x/L*1e6; la media es sobre
    celulas. Cacheado."""
    cache = os.path.join(WORK, "pb2025.npz")
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        return (z["perts"].astype(str), z["genes"].astype(str), z["cpm_pert"],
                z["n_cel"], z["cpm_ctrl"], z["det_ctrl"], z["fano_ctrl"],
                json.loads(str(z["meta"])))
    os.makedirs(WORK, exist_ok=True)
    genes = genes_2025()
    G = len(genes)
    acc = {}
    n_cel = {}
    c_sum = np.zeros(G, dtype=np.float64)     # suma de CPM del control
    c_raw = np.zeros(G, dtype=np.float64)     # suma de cuentas crudas
    c_raw2 = np.zeros(G, dtype=np.float64)
    c_nnz = np.zeros(G, dtype=np.float64)
    n_ctrl = 0
    meta = {"splits": {}, "n_celulas": 0}
    for split, rel in SPLITS.items():
        path = os.path.join(V25, rel)
        with h5py.File(path, "r") as f:
            shape = [int(x) for x in f["X"].attrs["shape"]]
            tg = np.asarray(decode_categorical(f["obs"]["target_gene"]), dtype=object).astype(str)
            indptr = f["X"]["indptr"][:].astype(np.int64)
            assert shape[1] == G, (shape, G)
            meta["splits"][split] = {"shape": shape, "n_perts": len(set(tg) - {CONTROL}),
                                     "n_control": int((tg == CONTROL).sum())}
            meta["n_celulas"] += shape[0]
            for i0 in range(0, shape[0], bloque):
                i1 = min(i0 + bloque, shape[0])
                a, b = int(indptr[i0]), int(indptr[i1])
                data = f["X"]["data"][a:b].astype(np.float64)
                ind = f["X"]["indices"][a:b]
                ip = indptr[i0:i1 + 1] - indptr[i0]
                Xb = sp.csr_matrix((data, ind, ip), shape=(i1 - i0, G))
                L = np.asarray(Xb.sum(axis=1)).ravel()
                L[L <= 0] = 1.0
                Xn = Xb.multiply((1e6 / L)[:, None]).tocsr()
                lab = tg[i0:i1]
                for t in np.unique(lab):
                    m = np.where(lab == t)[0]
                    s = np.asarray(Xn[m].sum(axis=0)).ravel()
                    if t == CONTROL:
                        c_sum += s
                        Xr = Xb[m]
                        c_raw += np.asarray(Xr.sum(axis=0)).ravel()
                        X2 = Xr.copy(); X2.data **= 2
                        c_raw2 += np.asarray(X2.sum(axis=0)).ravel()
                        c_nnz += np.diff(Xr.tocsc().indptr)
                        n_ctrl += len(m)
                    else:
                        acc[t] = acc.get(t, 0) + s
                        n_cel[t] = n_cel.get(t, 0) + len(m)
            log(f"2025 {split}: {shape[0]:,} celulas leidas")
    perts = np.asarray(sorted(acc), dtype=str)
    cpm_pert = np.stack([acc[t] / n_cel[t] for t in perts])
    ncel = np.asarray([n_cel[t] for t in perts], dtype=np.int64)
    cpm_ctrl = c_sum / n_ctrl
    mu = c_raw / n_ctrl
    var = c_raw2 / n_ctrl - mu ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        fano = np.where(mu > 0, var / mu, np.nan)
    det = c_nnz / n_ctrl
    meta["n_control_2025"] = int(n_ctrl)
    meta["n_perts"] = int(len(perts))
    np.savez_compressed(cache, perts=perts, genes=genes, cpm_pert=cpm_pert,
                        n_cel=ncel, cpm_ctrl=cpm_ctrl, det_ctrl=det,
                        fano_ctrl=np.nan_to_num(fano, nan=-1.0),
                        meta=json.dumps(meta))
    return perts, genes, cpm_pert, ncel, cpm_ctrl, det, np.nan_to_num(fano, nan=-1.0), meta


# ─────────────────────────────── gwps ───────────────────────────────
def gwps(targets):
    """Perfil log2 por (blanco, gen gwps) para los blancos pedidos. Si un
    blanco tiene varias construcciones, se promedian sus perfiles (declarado).
    Devuelve (var_names_gwps, {blanco: perfil}, meta). Cacheado."""
    cache = os.path.join(WORK, "gwps.npz")
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        vn = z["vn"].astype(str)
        P = z["P"]
        keys = z["keys"].astype(str)
        return vn, {k: P[i] for i, k in enumerate(keys)}, json.loads(str(z["meta"]))
    os.makedirs(WORK, exist_ok=True)
    with h5py.File(GWPS, "r") as f:
        gt = np.asarray([x.decode() if isinstance(x, bytes) else x
                         for x in f["obs/gene_transcript"][:]], dtype=str)
        g_gen = np.asarray([x.split("_")[1] if len(x.split("_")) > 1 else "" for x in gt])
        core = f["obs/core_control"][:].astype(bool)
        cats = np.asarray([x.decode() if isinstance(x, bytes) else x
                           for x in f["var/__categories/gene_name"][:]], dtype=str)
        vn = cats[f["var/gene_name"][:]]
        X = f["X"]
        idx_c = np.where(core)[0]
        ctrl = np.zeros(X.shape[1], dtype=np.float64)
        for s in range(0, len(idx_c), 64):
            ctrl += X[idx_c[s:s + 64], :].astype(np.float64).sum(0)
        ctrl /= len(idx_c)
        quiero = set(map(str, targets))
        filas = {}
        for i, g in enumerate(g_gen):
            if g in quiero and not core[i]:
                filas.setdefault(g, []).append(i)
        perf = {}
        for g, ii in filas.items():
            acc = np.zeros(X.shape[1], dtype=np.float64)
            for i in ii:
                x = X[i, :].astype(np.float64)
                acc += np.log2((x + EPS) / (ctrl + EPS))
            perf[g] = (acc / len(ii)).astype(np.float32)
    # nombres duplicados: fuera (como gwps_jiang_datos)
    u, cnt = np.unique(vn, return_counts=True)
    dup = set(u[cnt > 1].tolist())
    keep = np.asarray([g not in dup for g in vn])
    vn2 = vn[keep]
    perf = {g: p[keep] for g, p in perf.items()}
    meta = {"n_filas_gwps": int(len(gt)), "n_core_control": int(core.sum()),
            "n_genes_gwps": int(len(vn)), "n_genes_sin_duplicados": int(len(vn2)),
            "n_duplicados_fuera": int(len(dup)),
            "blancos_pedidos": int(len(set(map(str, targets)))),
            "blancos_con_perfil": int(len(perf)),
            "construcciones_por_blanco": {g: len(v) for g, v in filas.items() if len(v) > 1},
            "promedio_de_construcciones": "si (media de los perfiles log2)"}
    keys = np.asarray(sorted(perf), dtype=str)
    np.savez_compressed(cache, vn=vn2, P=np.stack([perf[k] for k in keys]),
                        keys=keys, meta=json.dumps(meta))
    return vn2, perf, meta


def rango_centrado(v):
    """Rango medio normalizado a [-0.5, 0.5] (empates promediados)."""
    n = v.size
    o = np.argsort(v, kind="stable")
    r = np.empty(n, dtype=np.float64)
    sv = v[o]
    i = 0
    while i < n:
        j = i
        while j + 1 < n and sv[j + 1] == sv[i]:
            j += 1
        r[o[i:j + 1]] = (i + j) / 2.0
        i = j + 1
    return r / max(n - 1, 1) - 0.5


# ─────────────────── rasgos del gen, DEL CONTROL DEL CONTEXTO ───────────────
def stats_control(Xc):
    """CPM medio, tasa de deteccion y Fano por gen, del pool de control dado."""
    Xw = Xc.tocsr().astype(np.float64)
    n = Xw.shape[0]
    L = np.asarray(Xw.sum(axis=1)).ravel()
    Lc = L.copy(); Lc[Lc <= 0] = 1.0
    cpm = np.asarray(Xw.T @ ((1e6 / Lc) / n)).ravel()
    det = np.diff(Xw.tocsc().indptr) / n
    mu = np.asarray(Xw.mean(axis=0)).ravel()
    X2 = Xw.copy(); X2.data **= 2
    var = np.asarray(X2.mean(axis=0)).ravel() - mu ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        fano = np.where(mu > 0, var / mu, 0.0)
    return {"cpm": cpm, "det": det, "fano": fano, "n": int(n),
            "profundidad_mediana": float(np.median(L))}


def matriz_rasgos(targets, genes, gw_vn, gw_perf, tf, st):
    """(n_targets*n_genes, n_rasgos) en el orden de RASGOS. Filas por bloques de
    blanco: fila = t*len(genes) + h."""
    G = len(genes)
    pos_gw = {g: i for i, g in enumerate(gw_vn)}
    col = np.asarray([pos_gw.get(g, -1) for g in genes])
    hay = col >= 0
    log10cpm = np.log10(st["cpm"] + 1e-3)
    det = st["det"]
    logfano = np.log10(np.maximum(st["fano"], 1e-3))
    pasa = (st["cpm"] > UMBRAL_CPM).astype(np.float64)
    pos_gen = {g: i for i, g in enumerate(genes)}
    X = np.zeros((len(targets) * G, len(RASGOS)), dtype=np.float32)
    for ti, t in enumerate(targets):
        a = ti * G
        p = np.zeros(G); r = np.zeros(G); s = np.zeros(G); tiene = np.zeros(G)
        if t in gw_perf:
            pf = gw_perf[t]
            v = pf[col[hay]]
            p[hay] = v
            rk = rango_centrado(pf.astype(np.float64))
            r[hay] = rk[col[hay]]
            s[hay] = np.sign(v)
            tiene[hay] = 1.0
        ig = pos_gen.get(t, -1)
        cpm_g = log10cpm[ig] if ig >= 0 else np.log10(1e-3)
        X[a:a + G, 0] = p
        X[a:a + G, 1] = r
        X[a:a + G, 2] = s
        X[a:a + G, 3] = tiene
        X[a:a + G, 4] = log10cpm
        X[a:a + G, 5] = det
        X[a:a + G, 6] = logfano
        X[a:a + G, 7] = pasa
        X[a:a + G, 8] = 1.0 if t in tf else 0.0
        X[a:a + G, 9] = cpm_g
    return X


if __name__ == "__main__":
    if sys.argv[1] == "pb":
        r = pseudobulk_2025()
        print("perts", r[0].shape, "genes", r[1].shape, "cpm_pert", r[2].shape)
        print(json.dumps(r[7], indent=1)[:800])
    elif sys.argv[1] == "gwps":
        t = sorted(set(perts_2026().tolist()) | set(pseudobulk_2025()[0].tolist()))
        vn, perf, meta = gwps(t)
        print("gwps genes", vn.shape, "perfiles", len(perf))
        print(json.dumps(meta, indent=1)[:900])
