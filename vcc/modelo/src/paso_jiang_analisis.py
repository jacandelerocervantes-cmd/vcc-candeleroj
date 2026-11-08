#!/usr/bin/env python
"""Parte 1 b-e sobre Jiang: controles, perturbaciones, solape entre lineas,
cruces (linea,pert) con >=400 celulas. Lee jiang_raw_cols.json (cell_type,
gene, guide crudos por archivo) y jiang_parte1_raw.json (n_obs). No toca
/X ni benchmark/.
"""
import json
import numpy as np
from collections import Counter, defaultdict

RAW_COLS_PATH = "/home/ubuntu/vcc/modelo/out/jiang_raw_cols.json"
PARTE1_PATH = "/home/ubuntu/vcc/modelo/out/jiang_parte1_raw.json"
CONTROL_LABEL = "NT"  # regla citada: gene == "NT" (confirmado via guide, NTg1..NTg14)

RETO_CONTROL_CELLS = 18400
RETO_N_GUIDES = 46
RETO_CELLS_PER_GUIDE = 400
RETO_CELLS_PER_PERT = 400


def main():
    raw = json.load(open(RAW_COLS_PATH))
    parte1 = json.load(open(PARTE1_PATH))
    labels = list(raw.keys())

    out = {"per_file_line": {}, "overlap": {}, "cells_per_pert_threshold": {}, "control_guides_detail": {}}

    for label in labels:
        cell_type = raw[label]["cell_type"]
        gene = raw[label]["gene"]
        guide = raw[label]["guide"]
        n_obs = parte1["files"][label]["n_obs"]
        assert len(cell_type) == len(gene) == len(guide) == n_obs

        lines = sorted(set(cell_type))
        out["per_file_line"][label] = {}

        # perturbacion -> set de lineas que la tienen (para d)
        pert_lines = defaultdict(set)

        for line in lines:
            idx = [i for i, c in enumerate(cell_type) if c == line]
            genes_this_line = [gene[i] for i in idx]
            guides_this_line = [guide[i] for i in idx]

            total_cells = len(idx)
            gene_counter = Counter(genes_this_line)
            n_control_cells = gene_counter.get(CONTROL_LABEL, 0)
            pert_counter = {g: c for g, c in gene_counter.items() if g != CONTROL_LABEL}
            n_perts = len(pert_counter)
            cells_per_pert_vals = np.array(list(pert_counter.values()))

            for g in pert_counter:
                pert_lines[g].add(line)

            # control: guias distintas y celulas por guia, SOLO para gene==NT
            control_guide_counter = Counter(
                gu for g, gu in zip(genes_this_line, guides_this_line) if g == CONTROL_LABEL
            )
            n_control_guides = len(control_guide_counter)
            control_guide_vals = np.array(list(control_guide_counter.values())) if control_guide_counter else np.array([0])

            out["per_file_line"][label][line] = {
                "total_cells": total_cells,
                "n_perturbaciones_distintas": n_perts,
                "cells_per_pert_median": float(np.median(cells_per_pert_vals)) if n_perts else None,
                "cells_per_pert_min": int(cells_per_pert_vals.min()) if n_perts else None,
                "cells_per_pert_p5": float(np.percentile(cells_per_pert_vals, 5)) if n_perts else None,
                "n_control_cells": n_control_cells,
                "n_control_guides": n_control_guides,
                "control_cells_per_guide_median": float(np.median(control_guide_vals)),
                "control_cells_per_guide_min": int(control_guide_vals.min()),
                "control_cells_per_guide_max": int(control_guide_vals.max()),
                # contraste explicito contra el reto
                "reto_control_cells": RETO_CONTROL_CELLS,
                "reto_n_guides": RETO_N_GUIDES,
                "reto_cells_per_guide": RETO_CELLS_PER_GUIDE,
            }

        # d) solape de perturbaciones entre lineas, dentro de este archivo
        n_lines = len(lines)
        overlap_hist = Counter(len(s) for s in pert_lines.values())
        out["overlap"][label] = {
            "n_lines": n_lines,
            "n_perts_total_union": len(pert_lines),
            "histograma_en_cuantas_lineas": {str(k): v for k, v in sorted(overlap_hist.items())},
            "perts_en_todas_las_lineas": sorted(g for g, s in pert_lines.items() if len(s) == n_lines),
        }

        # e) cruces (linea,pert) con >=400 celulas, para este archivo
        n_ge_400 = 0
        n_lt_400 = 0
        for line in lines:
            idx = [i for i, c in enumerate(cell_type) if c == line]
            genes_this_line = [gene[i] for i in idx]
            pert_counter = Counter(g for g in genes_this_line if g != CONTROL_LABEL)
            for g, c in pert_counter.items():
                if c >= RETO_CELLS_PER_PERT:
                    n_ge_400 += 1
                else:
                    n_lt_400 += 1
        out["cells_per_pert_threshold"][label] = {
            "n_cruces_linea_pert_total": n_ge_400 + n_lt_400,
            "n_cruces_con_400_o_mas": n_ge_400,
            "n_cruces_con_menos_de_400": n_lt_400,
        }

    print(json.dumps(out, indent=2, default=str))


if __name__ == "__main__":
    main()
