#!/usr/bin/env python
"""PARTE 3: el regimen del reto contra el de Jiang, en cifras.

QUE SE MIDE, con la MISMA definicion en los dos lados:
  * profundidad por celula = suma de cuentas de la fila (UMIs). En Jiang se
    cruza contra `obs/nCount_RNA`; en el reto se suma del propio CSR.
  * nnz por celula = `diff(indptr)`. En Jiang se cruza contra `obs/nFeature_RNA`.
  * n de celulas de control: 18,400 por contexto en el reto; las celulas `NT`
    de cada linea en Jiang.
  * tamano del panel probable: nº de genes cuyo CPM MEDIO POR CELULA en el pool
    de control supera 5.0 -- la MISMA cantidad con la que `_apply_cpm_filter`
    del preset decide que genes entran en BH (`filter_gene_min_cpm_cell`).

🔴 NO se concluye si la verdad vacia de Jiang transfiere. Se dan las cantidades.
NO se toca el lado perturbado de nada: en Jiang solo se leen las celulas `NT`.
"""
import json, os, sys
import numpy as np
import scipy.sparse as sp
import h5py

IFNB = "/datos/ref/convertido/Seurat_object_IFNB_Perturb_seq.h5ad"
RETO = "/datos/reto"
CTX = [("A", "val"), ("B", "val"), ("C", "val"),
       ("D", "final"), ("E", "final"), ("F", "final")]
OUT = "/home/ubuntu/vcc/modelo/out/ntc_nulo/p3.json"
UMBRAL_CPM = 5.0


def dec(node):
    cods = node["codes"][:]
    cats = [c.decode() if isinstance(c, bytes) else c for c in node["categories"][:]]
    return np.array([cats[c] if c >= 0 else None for c in cods], dtype=object)


def runs(idx):
    """Agrupa indices ordenados en tramos contiguos para leer en bloque."""
    out, i = [], 0
    while i < idx.size:
        j = i
        while j + 1 < idx.size and idx[j + 1] == idx[j] + 1:
            j += 1
        out.append((int(idx[i]), int(idx[j]) + 1))
        i = j + 1
    return out


def lee_filas(f, idx, n_genes):
    indptr = f["X/indptr"][:].astype(np.int64)
    D, I = f["X/data"], f["X/indices"]
    datos, indices, nuevo = [], [], [0]
    for a, b in runs(np.sort(idx)):
        lo, hi = int(indptr[a]), int(indptr[b])
        d = D[lo:hi]; ii = I[lo:hi]
        cortes = indptr[a:b + 1] - lo
        datos.append(d); indices.append(ii)
        nuevo.extend((nuevo[-1] + cortes[1:]).tolist())
    return sp.csr_matrix(
        (np.concatenate(datos).astype(np.float64), np.concatenate(indices),
         np.array(nuevo, dtype=np.int64)), shape=(idx.size, n_genes))


def stats(v, nombre):
    v = np.asarray(v, dtype=np.float64)
    return {f"{nombre}_media": float(v.mean()), f"{nombre}_mediana": float(np.median(v)),
            f"{nombre}_p1": float(np.percentile(v, 1)),
            f"{nombre}_p99": float(np.percentile(v, 99)),
            f"{nombre}_min": float(v.min()), f"{nombre}_max": float(v.max()),
            f"{nombre}_sd": float(v.std(ddof=1))}


def panel_probable(X):
    L = np.asarray(X.sum(axis=1)).ravel(); L[L <= 0] = 1.0
    cpm = np.asarray(X.T @ ((1e6 / L) / X.shape[0])).ravel()
    return int((cpm > UMBRAL_CPM).sum()), cpm


res = {"umbral_cpm_gen": UMBRAL_CPM, "reto": [], "jiang": [],
       "jiang_archivo": IFNB, "jiang_n_genes": None}

# ---------------------------------------------------------------- reto
for letra, part in CTX:
    p = os.path.join(RETO, part, f"context_{letra}.h5ad")
    with h5py.File(p, "r") as f:
        shape = tuple(int(x) for x in f["X"].attrs["shape"])
        indptr = f["X/indptr"][:].astype(np.int64)
        data = f["X/data"][:]
    nnz = np.diff(indptr)
    prof = np.add.reduceat(data.astype(np.float64), indptr[:-1])
    prof[nnz == 0] = 0.0
    with h5py.File(p, "r") as f:
        X = sp.csr_matrix((f["X/data"][:], f["X/indices"][:], f["X/indptr"][:]),
                          shape=shape)
    npan, _ = panel_probable(X)
    del X
    d = {"fuente": "reto", "unidad": letra, "particion": part,
         "n_celulas_control": int(shape[0]), "n_genes_panel": int(shape[1]),
         "nnz_total": int(indptr[-1]),
         "n_genes_cpm_mayor_5": npan,
         "frac_panel_probable": round(npan / shape[1], 4)}
    d.update(stats(prof, "profundidad")); d.update(stats(nnz, "nnz"))
    d["nnz_sobre_panel"] = round(float(nnz.mean()) / shape[1], 4)
    res["reto"].append(d)
    sys.stderr.write(f"  reto {letra}: prof_med={np.median(prof):.0f} "
                     f"nnz_med={np.median(nnz):.0f} panel={npan}\n")

# ---------------------------------------------------------------- jiang
with h5py.File(IFNB, "r") as f:
    n_genes = int(f["X"].attrs["shape"][1])
    res["jiang_n_genes"] = n_genes
    ct = dec(f["obs/cell_type"]); gn = dec(f["obs/gene"])
    indptr = f["X/indptr"][:].astype(np.int64)
    ncount = f["obs/nCount_RNA"][:].astype(np.float64)
    nfeat = f["obs/nFeature_RNA"][:].astype(np.float64)
    nnz_all = np.diff(indptr)
    res["jiang_nfeature_igual_indptr"] = bool(np.array_equal(nfeat, nnz_all))
    res["jiang_n_celulas_total"] = int(ct.size)
    for linea in sorted(set(ct.tolist())):
        m = (ct == linea) & (gn == "NT")
        idx = np.where(m)[0]
        X = lee_filas(f, idx, n_genes)
        prof_real = np.asarray(X.sum(axis=1)).ravel()
        npan, _ = panel_probable(X)
        nnz = np.diff(X.indptr)
        d = {"fuente": "jiang_IFNB", "unidad": linea,
             "n_celulas_control_NT": int(idx.size),
             "n_celulas_linea_total": int((ct == linea).sum()),
             "n_genes_panel": n_genes,
             "n_genes_cpm_mayor_5": npan,
             "frac_panel_probable": round(npan / n_genes, 4),
             "nCount_RNA_igual_suma_fila": bool(
                 np.allclose(ncount[idx], prof_real)),
             "nnz_total": int(X.nnz)}
        d.update(stats(prof_real, "profundidad")); d.update(stats(nnz, "nnz"))
        d["nnz_sobre_panel"] = round(float(nnz.mean()) / n_genes, 4)
        res["jiang"].append(d)
        sys.stderr.write(f"  jiang {linea}: NT={idx.size} prof_med={np.median(prof_real):.0f} "
                         f"nnz_med={np.median(nnz):.0f} panel={npan}\n")
        del X

json.dump(res, open(OUT, "w"), indent=1, sort_keys=True, default=str)
sys.stderr.write(f"\nlisto: {OUT}\n")
