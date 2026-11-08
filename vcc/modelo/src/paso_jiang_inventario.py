#!/usr/bin/env python
"""Inventario de pliegues sobre Jiang (parte 1). Lee SOLO /obs (y el atributo
`shape`/`var._index` para contar genes) de los cinco convertidos de
/datos/ref/convertido/. Nunca toca /X. No lee nada de benchmark/.
"""
import sys
import json
import h5py
import numpy as np
from collections import Counter

FILES = {
    "IFNB": "/datos/ref/convertido/Seurat_object_IFNB_Perturb_seq.h5ad",
    "IFNG": "/datos/ref/convertido/Seurat_object_IFNG_Perturb_seq.h5ad",
    "INS":  "/datos/ref/convertido/Seurat_object_INS_Perturb_seq.h5ad",
    "TGFB": "/datos/ref/convertido/Seurat_object_TGFB_Perturb_seq.h5ad",
    "TNFA": "/datos/ref/convertido/Seurat_object_TNFA_Perturb_seq.h5ad",
}


def decode_node(node):
    """Devuelve (values_array_as_list_or_ndarray, dtype_str, encoding)."""
    if hasattr(node, "keys"):
        attrs = dict(node.attrs)
        enc = attrs.get("encoding-type")
        if enc == "categorical":
            codes = node["codes"][:]
            cats_raw = node["categories"][:]
            cats = [c.decode() if isinstance(c, bytes) else c for c in cats_raw]
            values = [cats[c] if c >= 0 else None for c in codes]
            return values, f"categorical<{node['categories'].dtype}>", enc
        elif enc == "nullable-string-array":
            vals = node["values"][:]
            mask = node["mask"][:] if "mask" in node.keys() else None
            decoded = [v.decode() if isinstance(v, bytes) else v for v in vals]
            if mask is not None:
                decoded = [None if m else v for v, m in zip(decoded, mask)]
            return decoded, "nullable-string-array", enc
        else:
            # nullable-integer-array / nullable-boolean-array u otros
            vals = node["values"][:] if "values" in node.keys() else None
            return list(vals), f"nullable-array<{enc}>", enc
    else:
        vals = node[:]
        return vals, str(node.dtype), "plain"


def main():
    out = {"files": {}}

    for label, path in FILES.items():
        sys.stderr.write(f"leyendo obs de {label}: {path}\n")
        finfo = {"path": path}
        with h5py.File(path, "r") as f:
            xattrs = dict(f["X"].attrs)
            shape = [int(v) for v in xattrs["shape"]]
            finfo["n_obs"] = shape[0]
            finfo["n_genes"] = shape[1]

            var_index_node = f["var"]["_index"]
            var_vals, _, _ = decode_node(var_index_node)
            finfo["var_index"] = list(var_vals)  # se compara fuera, se recorta luego

            obs = f["obs"]
            cols = {}
            for name in obs.keys():
                if name == "_index":
                    continue
                vals, dtype_str, enc = decode_node(obs[name])
                vals_arr = np.array(vals, dtype=object)
                try:
                    n_unique = len(set(v for v in vals if v is not None))
                except TypeError:
                    n_unique = int(len(np.unique(np.asarray(vals))))
                cols[name] = {
                    "dtype": dtype_str,
                    "encoding": enc,
                    "n_unique": n_unique,
                }
                # guardamos los valores crudos SOLO de las columnas que se procesan
                # fuera de este script (b-e); pathway/sample/sample_ID no se usan ahi
                if name in ("cell_type", "gene", "guide"):
                    cols[name]["_raw"] = vals
            finfo["obs_columns"] = cols
        out["files"][label] = finfo
        sys.stderr.write(f"  n_obs={finfo['n_obs']} n_genes={finfo['n_genes']} "
                          f"n_obs_cols={len(finfo['obs_columns'])}\n")

    # ---- f) eje de genes: mismo n y mismo orden entre los 5 ----
    labels = list(FILES.keys())
    ref_label = labels[0]
    ref_genes = out["files"][ref_label]["var_index"]
    gene_axis_check = {}
    for label in labels:
        genes = out["files"][label]["var_index"]
        same_len = len(genes) == len(ref_genes)
        same_order = same_len and genes == ref_genes
        gene_axis_check[label] = {
            "n_genes": len(genes),
            "same_length_as_ref": same_len,
            "identical_order_as_ref": same_order,
        }
    out["gene_axis_check"] = gene_axis_check
    out["gene_axis_ref"] = ref_label

    # quitar var_index crudo para no inflar el json de salida
    for label in labels:
        out["files"][label].pop("var_index")

    # ---- a) columnas por archivo, con roles identificados empíricamente ----
    role_guess = {}
    for label in labels:
        cols = out["files"][label]["obs_columns"]
        # cell_type: pocas categorias, valores tipo nombre de linea celular
        # gene: categorias que incluyen simbolos de genes + un control
        # guide: muchas mas categorias que 'gene' (varias guias por gen)
        role_guess[label] = {
            "linea_celular_col": "cell_type" if "cell_type" in cols else None,
            "perturbacion_col": "gene" if "gene" in cols else None,
            "guia_col": "guide" if "guide" in cols else None,
            "cell_type_n_unique": cols.get("cell_type", {}).get("n_unique"),
            "gene_n_unique": cols.get("gene", {}).get("n_unique"),
            "guide_n_unique": cols.get("guide", {}).get("n_unique"),
        }
    out["role_guess"] = role_guess

    # separar los _raw para guardarlos aparte (no al stdout principal, muy grandes)
    raw_blob = {}
    for label in labels:
        cols = out["files"][label]["obs_columns"]
        raw_blob[label] = {}
        for name in ("cell_type", "gene", "guide"):
            if name in cols:
                raw_blob[label][name] = cols[name].pop("_raw")

    with open("/home/ubuntu/vcc/modelo/out/jiang_raw_cols.json", "w") as fh:
        json.dump(raw_blob, fh)

    print(json.dumps(out, indent=2, default=str))


if __name__ == "__main__":
    main()
