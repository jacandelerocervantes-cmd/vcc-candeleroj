#!/usr/bin/env python
"""Paso 1, Parte 1 — verifica lectura de /datos/reto/{val,final}/context_*.h5ad.

Modo backed: abre cada .h5ad con h5py, lee /obs, /var y los arrays crudos de
/X (data, indices, indptr) de un CSR ya existente en disco. No escribe nada,
no modifica nada, no toca benchmark/.

Uso:
    obs_venv/bin/python -I paso1_leer_entrada.py > /tmp/paso1_raw.json
"""
import sys
import json
import h5py
import numpy as np
import csv

FILES = {
    "A": "/datos/reto/val/context_A.h5ad",
    "B": "/datos/reto/val/context_B.h5ad",
    "C": "/datos/reto/val/context_C.h5ad",
    "D": "/datos/reto/final/context_D.h5ad",
    "E": "/datos/reto/final/context_E.h5ad",
    "F": "/datos/reto/final/context_F.h5ad",
}
GENE_NAMES = {
    "val": "/datos/reto/val/gene_names.csv",
    "final": "/datos/reto/final/gene_names.csv",
}
PERT_COUNTS = {
    "val": "/datos/reto/val/pert_counts.csv",
    "final": "/datos/reto/final/pert_counts.csv",
}
MANIFEST = {
    "val": "/datos/reto/val/manifest.json",
    "final": "/datos/reto/final/manifest.json",
}


def read_csv_col(path, col=0, skip_header=True):
    out = []
    with open(path, newline="") as fh:
        r = csv.reader(fh)
        rows = list(r)
    start = 1 if skip_header else 0
    for row in rows[start:]:
        out.append(row[col])
    return out


def decode_str_array(ds):
    if hasattr(ds, "keys") and "values" in ds.keys():
        # nullable-string-array: {values, mask}; mask=True significa NaN
        vals = ds["values"][:]
        mask = ds["mask"][:] if "mask" in ds.keys() else None
        decoded = [v.decode() if isinstance(v, bytes) else v for v in vals]
        if mask is not None:
            decoded = [None if m else v for v, m in zip(decoded, mask)]
        return decoded
    vals = ds[:]
    if vals.dtype.kind in ("S", "O"):
        return [v.decode() if isinstance(v, bytes) else v for v in vals]
    return list(vals)


def read_categorical(group, name):
    """Lee una columna categórica guardada como grupo codes/categories, o como
    array plano si no es categórica."""
    node = group[name]
    if hasattr(node, "keys"):
        codes = node["codes"][:]
        cats = decode_str_array(node["categories"])
        values = [cats[c] if c >= 0 else None for c in codes]
        return values, cats
    else:
        values = decode_str_array(node)
        return values, sorted(set(values))


def inspect_file(letter, path):
    out = {"letter": letter, "path": path}
    with h5py.File(path, "r") as f:
        xattrs = dict(f["X"].attrs)
        shape = [int(v) for v in xattrs["shape"]]
        out["X_encoding"] = xattrs.get("encoding-type")
        out["shape"] = shape
        n_obs, n_vars = shape

        data = f["X"]["data"]
        indices = f["X"]["indices"]
        indptr = f["X"]["indptr"][:]
        out["X_data_dtype"] = str(data.dtype)
        out["X_indices_dtype"] = str(indices.dtype)
        out["nnz_total"] = int(len(data))

        data_vals = data[:]  # ~109M float32 => ~440MB, fits in RAM easily
        out["X_min"] = float(data_vals.min())
        out["X_max"] = float(data_vals.max())
        # enteros: comprobar sobre TODO el array (no es demasiado grande)
        is_int = np.all(data_vals == np.round(data_vals))
        out["X_all_values_integer"] = bool(is_int)
        out["X_sample_first10"] = [float(v) for v in data_vals[:10]]

        nnz_per_cell = np.diff(indptr).astype(np.int64)
        # suma por celda (UMIs) a partir de data+indptr, sin densificar
        sums = np.add.reduceat(data_vals, indptr[:-1]) if n_obs > 0 else np.array([])
        # reduceat falla si alguna fila tiene 0 nnz (indptr repetido) -> corregir
        zero_rows = nnz_per_cell == 0
        if zero_rows.any():
            sums = sums.copy()
            sums[zero_rows] = 0.0

        out["umi_median"] = float(np.median(sums))
        out["umi_mean"] = float(np.mean(sums))
        out["nnz_median"] = float(np.median(nnz_per_cell))
        out["nnz_p95"] = float(np.percentile(nnz_per_cell, 95))
        out["nnz_max"] = int(nnz_per_cell.max())
        out["nnz_min"] = int(nnz_per_cell.min())

        # obs
        obs = f["obs"]
        ctx_vals, ctx_cats = read_categorical(obs, "context")
        ntc_vals, ntc_cats = read_categorical(obs, "ntc_id")
        tg_vals, tg_cats = read_categorical(obs, "target_gene")

        out["context_unique"] = sorted(set(ctx_vals))
        out["context_n_unique"] = len(set(ctx_vals))
        out["target_gene_unique"] = sorted(set(tg_vals))
        out["target_gene_n_unique"] = len(set(tg_vals))
        out["ntc_id_n_unique"] = len(set(ntc_vals))
        out["ntc_id_categories"] = sorted(set(ntc_vals))

        from collections import Counter
        ntc_counts = Counter(ntc_vals)
        out["ntc_id_counts_min"] = min(ntc_counts.values())
        out["ntc_id_counts_max"] = max(ntc_counts.values())
        out["ntc_id_counts_all_400"] = bool(all(v == 400 for v in ntc_counts.values()))

        # var / gene order
        var_index = decode_str_array(f["var"]["_index"])
        out["n_vars"] = len(var_index)
        out["var_index_first5"] = var_index[:5]
        out["var_index_last5"] = var_index[-5:]
        out["_var_index_full"] = var_index  # para comparar fuera, se quita luego

    return out


def main():
    results = {}
    for letter, path in FILES.items():
        sys.stderr.write(f"leyendo {letter}: {path}\n")
        results[letter] = inspect_file(letter, path)

    # gene_names.csv: comparar orden exacto, val y final
    gene_names = {}
    for part, p in GENE_NAMES.items():
        gene_names[part] = read_csv_col(p, col=0, skip_header=True)

    gene_order_check = {}
    for letter, res in results.items():
        part = "val" if letter in "ABC" else "final"
        ref = gene_names[part]
        var_idx = res.pop("_var_index_full")
        same_len = len(var_idx) == len(ref)
        same_order = same_len and all(a == b for a, b in zip(var_idx, ref))
        first_mismatch = None
        if not same_order:
            for i, (a, b) in enumerate(zip(var_idx, ref)):
                if a != b:
                    first_mismatch = {"index": i, "var": a, "gene_names_csv": b}
                    break
        gene_order_check[letter] = {
            "n_var_h5ad": len(var_idx),
            "n_gene_names_csv": len(ref),
            "same_length": same_len,
            "identical_order": same_order,
            "first_mismatch": first_mismatch,
        }

    gene_names_val_vs_final_identical = gene_names["val"] == gene_names["final"]

    # pert_counts.csv: contar perturbaciones, comparar val vs final
    pert_counts = {}
    for part, p in PERT_COUNTS.items():
        pert_counts[part] = read_csv_col(p, col=0, skip_header=True)

    pert_check = {
        "val_n_perturbaciones": len(pert_counts["val"]),
        "final_n_perturbaciones": len(pert_counts["final"]),
        "val_has_duplicates": len(pert_counts["val"]) != len(set(pert_counts["val"])),
        "final_has_duplicates": len(pert_counts["final"]) != len(set(pert_counts["final"])),
        "val_eq_final_as_sets": set(pert_counts["val"]) == set(pert_counts["final"]),
        "n_overlap": len(set(pert_counts["val"]) & set(pert_counts["final"])),
    }

    # manifest.json
    manifests = {}
    for part, p in MANIFEST.items():
        with open(p) as fh:
            manifests[part] = json.load(fh)

    # ntc_id: val vs final, mismas guias o distintas
    ntc_val_union = set()
    ntc_final_union = set()
    for letter, res in results.items():
        if letter in "ABC":
            ntc_val_union |= set(res["ntc_id_categories"])
        else:
            ntc_final_union |= set(res["ntc_id_categories"])
    ntc_check = {
        "val_union": sorted(ntc_val_union),
        "final_union": sorted(ntc_final_union),
        "val_eq_final": ntc_val_union == ntc_final_union,
        "n_overlap": len(ntc_val_union & ntc_final_union),
        "n_val_only": len(ntc_val_union - ntc_final_union),
        "n_final_only": len(ntc_final_union - ntc_val_union),
    }

    out = {
        "per_file": results,
        "gene_order_check": gene_order_check,
        "gene_names_val_vs_final_identical": gene_names_val_vs_final_identical,
        "pert_check": pert_check,
        "ntc_check": ntc_check,
        "manifests": manifests,
    }
    print(json.dumps(out, indent=2, default=str))


if __name__ == "__main__":
    main()
