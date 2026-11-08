#!/usr/bin/env python
"""La componente GENUINA de `r` para `jac`, por particion.

POR QUE HACE FALTA RECALCULAR: el bundle persiste solo el agregado por
particion (`anchor_splits.parquet`: split_index, seed, metric, value,
n_perturbations). El agregado MEZCLA la componente de convencion (las uniones
vacias que cobran 1.0) con la genuina, asi que la descomposicion de `r` no se
puede hacer leyendo artefactos. Hay que recomputar las particiones.

QUE SE REPRODUCE, exactamente:
  * la particion: `cell_eval2.ceiling._disjoint_halves` (la misma funcion que
    usa `anchor.py:117`), sobre la MISMA poblacion y con las MISMAS semillas
    derivadas que el bundle estampo en `anchor_splits.parquet`
  * la puntuacion: `half_b` como prediccion contra `half_a` como verdad
    (`anchor.py:120`), con `control_source=pred` (`anchor.py:103`)
  * la etiqueta de control `non-targeting`, sin `--control`: necesaria para que
    el `groupby` de `_disjoint_halves` reparta las celulas igual que el anchor

LINEAS: A549 y HAP1. Criterio declarado ANTES de medir: son las dos celdas que
el encargo nombra por su span/SE sobre el crudo (1.19 y 1.16). No se eligieron
por su valor ni por su signo.

No modifica nada: `_disjoint_halves` y `contar_de_por_perturbacion` se usan tal
cual.
"""
import json
import os
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from suelo_arnes import (
    IFNB_PATH, CONTROL_LABEL, CELLEVAL2_PY,
    build_fold_adatas, load_global_obs, load_var_names,
)
from suelo_floor import predict_floor
from techo_arnes import contar_de_por_perturbacion

FOLD_IDX = {"A549": 0, "BXPC3": 1, "HAP1": 2, "HT29": 3, "K562": 4, "MCF7": 5}
# Por defecto las dos que el encargo de `descomposicion.md` nombraba. Las
# lineas a correr se pueden pasar por argv (hueco 1 de `huecos.md`: completar
# las cuatro que faltaban, con el MISMO metodo, sin tocar nada mas).
LINEAS = {k: FOLD_IDX[k] for k in (sys.argv[1:] or ["A549", "HAP1"])}
PANEL_PATH = "/home/ubuntu/vcc/modelo/out/suelo_panel.json"
BUNDLES = "/datos/modelo/out/nativas_comp"
OUT_JSON = "/home/ubuntu/vcc/modelo/out/descomp_rgen.json"
WORK = "/datos/modelo/out/descomp_rgen"
CONTROL_COMP = "non-targeting"
SEED_BASE_REAL = 20260001
SEED_BASE_PRED = 20261001


def semillas_del_bundle(line):
    """Las semillas derivadas que el bundle estampo -- no se recalculan, se leen,
    para que la reproduccion no dependa de reimplementar SEED_DERIVATION."""
    import polars as pl
    p = os.path.join(BUNDLES, f"BUNDLE_{line}", "anchor_splits.parquet")
    d = pl.read_parquet(p).filter(
        pl.col("metric") == "de_wilcoxon_sig_jaccard"
    ).sort("split_index")
    return list(zip(d["split_index"].to_list(), d["seed"].to_list(),
                    d["value"].to_list()))


def corre(cmd, log):
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    with open(log, "w") as fh:
        fh.write("CMD: " + " ".join(cmd) + "\n\nSTDOUT:\n" + p.stdout
                 + "\n\nSTDERR:\n" + p.stderr
                 + f"\nEXIT:{p.returncode}\nDUR:{time.time()-t0:.1f}\n")
    if p.returncode != 0:
        raise RuntimeError(f"fallo: {' '.join(cmd[:4])}... ver {log}")
    return p


def main():
    os.makedirs(WORK, exist_ok=True)
    panel = json.load(open(PANEL_PATH))["panel_comun"]
    cell_type, gene = load_global_obs(IFNB_PATH)
    var_names = load_var_names(IFNB_PATH)
    n_genes = len(var_names)

    from cell_eval2.ceiling import _disjoint_halves

    salida = json.load(open(OUT_JSON)) if os.path.exists(OUT_JSON) else {}

    for line, fold_idx in LINEAS.items():
        if line in salida:
            sys.stderr.write(f"{line}: ya hecho\n")
            continue
        sys.stderr.write(f"\n=== {line} ===\n")
        real_ad, _desc, _ = build_fold_adatas(
            line, panel, predict_floor, var_names,
            SEED_BASE_REAL + fold_idx, SEED_BASE_PRED + fold_idx,
            cell_type, gene, n_genes, IFNB_PATH,
        )
        del _desc
        real_ad.obs["target_gene"] = real_ad.obs["target_gene"].replace(
            {CONTROL_LABEL: CONTROL_COMP}
        )

        splits = []
        for idx, seed, valor_bundle in semillas_del_bundle(line):
            t0 = time.time()
            half_a, half_b = _disjoint_halves(real_ad, "target_gene", CONTROL_COMP, seed)
            ra = os.path.join(WORK, f"{line}_s{idx}_a.h5ad")
            rb = os.path.join(WORK, f"{line}_s{idx}_b.h5ad")
            half_a.write_h5ad(ra)
            half_b.write_h5ad(rb)
            na, nb = half_a.n_obs, half_b.n_obs
            del half_a, half_b

            outdir = os.path.join(WORK, f"{line}_s{idx}_run")
            cmd = [
                CELLEVAL2_PY, "run",
                "-ap", rb, "-ar", ra,                 # half_b como PREDICCION
                "--preset", "vcc2026", "--profile", "vcc2026",
                "--pert-col", "target_gene",
                "--write-degenes",
                "--set", "control_source=pred",
                "-o", outdir,
            ]
            sys.stderr.write(f"  split {idx} (seed {seed})... ")
            sys.stderr.flush()
            corre(cmd, os.path.join(WORK, f"{line}_s{idx}.log"))

            cuentas = contar_de_por_perturbacion(outdir, control=CONTROL_COMP)
            import pandas as pd
            res = pd.read_csv(os.path.join(outdir, "results.csv"))
            jac = {str(r.perturbation): float(r.value)
                   for r in res[res.metric == "de_wilcoxon_sig_jaccard"].itertuples()}
            agg = pd.read_csv(os.path.join(outdir, "agg_results.csv"))
            rep = float(agg.loc[agg["statistic"] == "mean"].iloc[0]
                        ["de_wilcoxon_sig_jaccard"])

            vacias = [p for p, c in cuentas.items() if c["union_vacia"]]
            genuinos = [jac[p] for p in jac if p not in vacias]
            E = len(vacias)
            sal = {
                "split_index": idx, "seed": seed,
                "valor_del_bundle": valor_bundle,
                "valor_recomputado": rep,
                "dif_vs_bundle": abs(rep - valor_bundle),
                "n_half_a": int(na), "n_half_b": int(nb),
                "E_union_vacia": E,
                "n_genuinos": len(genuinos),
                "media_genuinos": float(np.mean(genuinos)) if genuinos else float("nan"),
                "suma_genuinos": float(np.sum(genuinos)) if genuinos else 0.0,
                "componente_convencion": E / 33.0,
                "componente_genuina": (float(np.sum(genuinos)) / 33.0) if genuinos else 0.0,
                "perts_union_vacia": sorted(vacias),
                # Hueco 2: la cohorte del LADO REAL de esta particion. El lado
                # real del anchor es `half_a` -- media poblacion -- asi que su
                # cohorte no tiene por que coincidir con la del run completo.
                # Se guarda por perturbacion para poder compararla, no inferirla.
                "por_pert": {
                    p: {
                        "n_DE_real": c["n_DE_real"],
                        "n_DE_pred": c["n_DE_pred"],
                        "n_DE_union": c["n_DE_union"],
                        "jac": jac.get(p),
                    }
                    for p, c in sorted(cuentas.items())
                },
            }
            splits.append(sal)
            sys.stderr.write(
                f"rep={rep:.6f} bundle={valor_bundle:.6f} dif={sal['dif_vs_bundle']:.1e} "
                f"E={E} gen_med={sal['media_genuinos']:.6f} ({time.time()-t0:.0f}s)\n"
            )

            for q in (ra, rb,
                      os.path.join(outdir, "de_real.parquet"),
                      os.path.join(outdir, "de_pred.parquet")):
                if os.path.exists(q):
                    os.remove(q)

        del real_ad
        gen = [s["componente_genuina"] for s in splits]
        conv = [s["componente_convencion"] for s in splits]
        crudo = [s["valor_recomputado"] for s in splits]
        salida[line] = {
            "splits": splits,
            "r_crudo_media": float(np.mean(crudo)),
            "r_crudo_sd": float(np.std(crudo, ddof=0)),
            "r_genuina_media": float(np.mean(gen)),
            "r_genuina_sd": float(np.std(gen, ddof=0)),
            "r_convencion_media": float(np.mean(conv)),
            "r_convencion_sd": float(np.std(conv, ddof=0)),
            "max_dif_vs_bundle": max(s["dif_vs_bundle"] for s in splits),
        }
        with open(OUT_JSON, "w") as fh:
            json.dump(salida, fh, indent=2, sort_keys=True)
        sys.stderr.write(
            f"  [{line}] r crudo={salida[line]['r_crudo_media']:.6f}"
            f" sd={salida[line]['r_crudo_sd']:.6f} | "
            f"r genuina={salida[line]['r_genuina_media']:.6f}"
            f" sd={salida[line]['r_genuina_sd']:.6f} | "
            f"max dif vs bundle={salida[line]['max_dif_vs_bundle']:.2e}\n"
        )

    sys.stderr.write(f"\nlisto: {OUT_JSON}\n")


if __name__ == "__main__":
    main()
