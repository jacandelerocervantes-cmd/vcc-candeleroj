#!/usr/bin/env python
"""Panel común del cribado: perturbaciones de IFNB con >=400 celulas en las
SEIS lineas (A549, BXPC3, HAP1, HT29, K562, MCF7), no solo en la linea
excluida de turno. Lee obs (cell_type, gene) directo del h5ad de IFNB via
h5py -- nunca toca /X. No lee nada de benchmark/.
"""
import sys
import json
from collections import Counter, defaultdict

import h5py
import numpy as np

IFNB_PATH = "/datos/ref/convertido/Seurat_object_IFNB_Perturb_seq.h5ad"
LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
CONTROL_LABEL = "NT"
MIN_CELLS = 400
OUT_PATH = "/home/ubuntu/vcc/modelo/out/suelo_panel.json"


def decode_categorical(node):
    codes = node["codes"][:]
    cats_raw = node["categories"][:]
    cats = [c.decode() if isinstance(c, bytes) else c for c in cats_raw]
    return [cats[c] if c >= 0 else None for c in codes]


def main():
    sys.stderr.write(f"leyendo obs de IFNB: {IFNB_PATH}\n")
    with h5py.File(IFNB_PATH, "r") as f:
        n_obs = int(f["X"].attrs["shape"][0])
        cell_type = decode_categorical(f["obs"]["cell_type"])
        gene = decode_categorical(f["obs"]["gene"])
        guide = decode_categorical(f["obs"]["guide"])
    assert len(cell_type) == len(gene) == len(guide) == n_obs
    sys.stderr.write(f"  n_obs={n_obs}\n")

    cell_type = np.array(cell_type, dtype=object)
    gene = np.array(gene, dtype=object)
    guide = np.array(guide, dtype=object)

    lines_present = sorted(set(cell_type.tolist()))
    assert lines_present == sorted(LINES), (lines_present, LINES)

    # cuenta (linea, pert) para toda pert != control, en cada linea
    per_line_counts = {}  # line -> {pert: n}
    per_line_control = {}  # line -> {"n_cells": int, "n_guides": int}
    for line in LINES:
        mask = cell_type == line
        genes_l = gene[mask]
        guides_l = guide[mask]
        pert_counter = Counter(g for g in genes_l if g != CONTROL_LABEL)
        per_line_counts[line] = dict(pert_counter)
        ctrl_mask = genes_l == CONTROL_LABEL
        n_control_cells = int(ctrl_mask.sum())
        n_control_guides = len(set(guides_l[ctrl_mask].tolist()))
        per_line_control[line] = {
            "n_control_cells": n_control_cells,
            "n_control_guides": n_control_guides,
        }
        sys.stderr.write(
            f"  {line}: total={mask.sum()} n_perts={len(pert_counter)} "
            f"n_control_cells={n_control_cells} n_control_guides={n_control_guides}\n"
        )

    # union de todas las perturbaciones vistas en IFNB
    all_perts = sorted(set().union(*[set(d.keys()) for d in per_line_counts.values()]))

    panel = []
    descartes_por_linea = defaultdict(list)  # line -> [pert con <400 en esa linea]
    per_pert_counts_table = []  # intermedio: una fila por perturbacion x linea
    for pert in all_perts:
        counts = {line: per_line_counts[line].get(pert, 0) for line in LINES}
        min_count = min(counts.values())
        in_panel = min_count >= MIN_CELLS
        if in_panel:
            panel.append(pert)
        for line in LINES:
            if counts[line] < MIN_CELLS:
                descartes_por_linea[line].append(pert)
            per_pert_counts_table.append({
                "pert": pert, "linea": line, "n_cells": counts[line],
                "cumple_400": counts[line] >= MIN_CELLS,
            })

    panel = sorted(panel)

    out = {
        "condicion": "IFNB",
        "min_cells": MIN_CELLS,
        "lineas": LINES,
        "n_perts_union_ifnb": len(all_perts),
        "n_panel_comun": len(panel),
        "panel_comun": panel,
        "descartes_por_linea": {
            line: {"n_descartadas": len(ps), "perts": sorted(ps)}
            for line, ps in descartes_por_linea.items()
        },
        "per_line_control": per_line_control,
        "per_pert_counts_table": per_pert_counts_table,
    }

    if len(panel) < 20:
        sys.stderr.write(
            f"PARAR: panel comun tiene solo {len(panel)} perturbaciones (<20), "
            "las metricas de DE no informan con un panel tan chico.\n"
        )
        with open(OUT_PATH, "w") as fh:
            json.dump(out, fh, indent=2)
        sys.exit(2)

    with open(OUT_PATH, "w") as fh:
        json.dump(out, fh, indent=2)
    sys.stderr.write(f"panel comun: {len(panel)} perturbaciones. Escrito en {OUT_PATH}\n")
    print(json.dumps({k: v for k, v in out.items() if k != "per_pert_counts_table"}, indent=2))


if __name__ == "__main__":
    main()
