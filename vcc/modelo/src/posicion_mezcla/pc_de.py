#!/usr/bin/env python
"""PARTE C / paso 1: recomputar las tablas DE REALES de Jiang, con `p` por gen.

🔴 POR QUE HACE FALTA. El cache `/datos/modelo/out/pds_obj/ed_real.json` guarda
SOLO las listas de genes del conjunto ED, ordenadas alfabeticamente: no trae `p`
ni el orden por `p`. La version IGUALADA EN POTENCIA de la Parte C necesita
recortar cada linea a los `m` genes de `p` MAS BAJA, asi que hay que recomputar.

🔴 TRANSDUCTIVO. Lee el lado real perturbado de Jiang. Va al estante del `b`:
ninguna cifra de aqui entra en una submission ni en un rasgo de un candidato.

LA CONSTRUCCION, declarada antes de medir
-----------------------------------------
* Poblacion real: `build_fold_adatas` con `seed_real = 20260001 + fold`, la
  MISMA que uso `pds_objetivo` -- se verifica por huella sha256 contra el cache.
* La puerta: cada parametro leido de `EvalConfig.from_preset("vcc2026")`, por el
  mismo camino que el marcador en CPU (`compute_de` -> `prep_de_side`), igual
  que en `ntc_nulo.md`. Backend `auto` -> `pdex` (sin CUDA).
* Conjunto ED = filas con `p_adj < 0.05`, **excluyendo el propio gen blanco** de
  cada perturbacion, que es lo que hace `de_sig_jaccard` (#172).
* VERIFICACION: el conjunto ED recomputado se compara gen por gen contra
  `ed_real.json`. Si no coincide, se reporta la discrepancia como dato.
"""
import json, os, sys
import numpy as np
import polars as pl

sys.path.insert(0, "/home/ubuntu/vcc/modelo/src")
from suelo_arnes import (IFNB_PATH, build_fold_adatas, load_global_obs,
                         load_var_names)
import gwps_via as GV
from cell_eval2.config import EvalConfig
from cell_eval2.de_compute import compute_de
from cell_eval2.de import prep_de_side

CFG = EvalConfig.from_preset("vcc2026")
LINES = GV.LINES
PANEL_PATH = GV.PANEL_PATH
SEED_BASE_REAL = 20260001
ED_CACHE = "/datos/modelo/out/pds_obj/ed_real.json"
WORK = "/datos/modelo/out/posicion_mezcla"
OUTDIR = "/home/ubuntu/vcc/modelo/out/posicion_mezcla"


def predict_trivial():
    def fn(control_pool_X, labels, rng):
        n = control_pool_X.shape[0]
        return {l: control_pool_X[rng.integers(0, n, size=400)].tocsr()
                for l in labels}
    return fn


def main():
    panel = sorted(map(str, json.load(open(PANEL_PATH))["panel_comun"]))
    cell_type, gene = load_global_obs(IFNB_PATH)
    var_names = np.asarray(load_var_names(IFNB_PATH)).astype(str)
    n_genes = len(var_names)
    cache = json.load(open(ED_CACHE))
    c = CFG
    meta = {"puerta": {"preset": "vcc2026", "p_adj_threshold": c.de.p_adj_threshold,
                       "backend_pedido": c.de.backend, "input_type": c.input_type,
                       "target_sum": c.target_sum, "mean_calc": c.de.mean_calc,
                       "epsilon": c.de.epsilon, "fdr_scope": c.de.fdr_scope,
                       "filter_gene_min_cpm_cell": c.filter.filter_gene_min_cpm_cell,
                       "nan_lfc_policy": c.de.nan_lfc_policy,
                       "min_abs_log2fc": c.de.min_abs_log2fc},
            "panel": panel, "lineas": {}}

    for line in (sys.argv[1:] or LINES):
        fold = LINES.index(line)
        sys.stderr.write(f"\n===== {line} (fold {fold}) =====\n")
        real_ad, _, _ = build_fold_adatas(
            line, panel, predict_trivial(), var_names.tolist(),
            SEED_BASE_REAL + fold, 999_600 + fold,
            cell_type, gene, n_genes, IFNB_PATH)
        GV.renombra(real_ad)
        huella = GV.huella(real_ad)
        hcache = cache.get(f"{line}__huella_real")
        sys.stderr.write(f"  huella {huella[:16]}... cache {str(hcache)[:16]}... "
                         f"IGUAL={huella == hcache}\n")

        df = compute_de(
            real_ad, backend=c.de.backend, groupby="target_gene",
            reference=GV.CONTROL_COMP, mean_calc=c.de.mean_calc,
            epsilon=c.de.epsilon, input_type=c.input_type,
            target_sum=c.target_sum, clip_value=c.de.clip_value,
            filter_gene_min_cpm_cell=c.filter.filter_gene_min_cpm_cell,
            fdr_scope=c.de.fdr_scope, threads=c.num_threads,
            replicate_col=c.de.replicate_col, device="cpu")
        del real_ad
        tabla, _ = prep_de_side(df, name="real", sort_by=c.de.sort_by,
                                nan_lfc_policy=c.de.nan_lfc_policy,
                                min_abs_log2fc=c.de.min_abs_log2fc)
        # conjunto ED: p_adj < T, y FUERA el propio gen blanco (#172)
        sig = (tabla.filter(pl.col("p_adj") < c.de.p_adj_threshold)
               .filter(pl.col("target") != pl.col("feature")))
        tabla.select("target", "feature", "p_value", "p_adj",
                     "log2_fold_change").write_parquet(
            os.path.join(WORK, f"de_{line}.parquet"))

        # verificacion contra el cache
        mio = {t: set(g) for t, g in
               sig.group_by("target").agg(pl.col("feature")).iter_rows()}
        suyo = {k: set(v) for k, v in cache[line].items()}
        comp = []
        for p in panel:
            a, b = mio.get(p, set()), suyo.get(p, set())
            comp.append({"pert": p, "n_mio": len(a), "n_cache": len(b),
                         "igual": a == b,
                         "solo_mio": len(a - b), "solo_cache": len(b - a)})
        n_ig = sum(1 for x in comp if x["igual"])
        meta["lineas"][line] = {
            "fold": fold, "seed_real": SEED_BASE_REAL + fold,
            "huella_real_sha256": huella,
            "huella_coincide_con_cache": bool(huella == hcache),
            "n_genes_probados": int(tabla["feature"].n_unique()),
            "n_filas_tabla": int(tabla.height),
            "n_perts_ED_identico_al_cache": n_ig,
            "n_perts": len(panel),
            "verificacion_por_pert": comp,
        }
        sys.stderr.write(
            f"  genes probados {tabla['feature'].n_unique()}, filas {tabla.height}, "
            f"ED identico al cache en {n_ig}/{len(panel)} perts\n")
        json.dump(meta, open(os.path.join(OUTDIR, "de_meta.json"), "w"),
                  indent=1, sort_keys=True, default=str)

    json.dump(meta, open(os.path.join(OUTDIR, "de_meta.json"), "w"),
              indent=1, sort_keys=True, default=str)
    sys.stderr.write(f"\nlisto: {WORK}/de_*.parquet + {OUTDIR}/de_meta.json\n")


if __name__ == "__main__":
    main()
