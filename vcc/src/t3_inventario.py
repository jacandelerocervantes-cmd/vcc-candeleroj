"""T3 -- construye corpus/inventario.parquet: una fila por
(fuente, archivo, linea_celular_cruda, contexto, target_crudo, construccion),
sobre TODO /datos/ref, /datos/reto, /datos/vcc2025. Sin colapsar simbolos,
sin normalizar nombres de linea, sin inferir tipo_valor_en_X (se lee).
"""
import json
import re
import sys
import time

import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

OUT_PARQUET = "/home/ubuntu/vcc/corpus/inventario.parquet"
ROWS = []


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_legacy_categorical(obs_grp, col):
    codes = obs_grp[col][()]
    cats = obs_grp["__categories"][col][()]
    cats = np.array([c.decode() if isinstance(c, bytes) else c for c in cats])
    return cats[codes]


def read_new_categorical(col_obj):
    cats = col_obj["categories"][()]
    cats = np.array([c.decode() if isinstance(c, bytes) else c for c in cats])
    codes = col_obj["codes"][()]
    out = np.where(codes >= 0, cats[np.clip(codes, 0, len(cats) - 1)], None)
    return out


def sample_dense_X_stats(f, max_rows=20000):
    X = f["X"]
    n = X.shape[0]
    rows = min(max_rows, n)
    sample = X[0:rows, :]
    sample = np.asarray(sample)
    finite = sample[np.isfinite(sample)]
    is_int = bool(np.all(finite == np.round(finite))) if finite.size else None
    return {
        "dtype": str(X.dtype),
        "min": float(finite.min()) if finite.size else None,
        "max": float(finite.max()) if finite.size else None,
        "todo_entero_en_muestra": is_int,
        "n_muestra": int(finite.size),
    }


def get_x_shape(f):
    X = f["X"]
    if "shape" in X.attrs:
        return tuple(int(v) for v in X.attrs["shape"])
    return X.shape


def sample_csr_X_stats(f, max_elems=5_000_000):
    data = f["X"]["data"]
    n = min(max_elems, data.shape[0])
    sample = data[0:n]
    sample = np.asarray(sample)
    finite = sample[np.isfinite(sample)]
    is_int = bool(np.all(finite == np.round(finite))) if finite.size else None
    return {
        "dtype": str(data.dtype),
        "min": float(finite.min()) if finite.size else None,
        "max": float(finite.max()) if finite.size else None,
        "todo_entero_en_muestra": is_int,
        "n_muestra": int(finite.size),
    }


def tipo_valor_desde_stats(stats):
    if stats is None:
        return None
    if stats["todo_entero_en_muestra"]:
        return "conteo_entero"
    return "normalizado_o_log_no_determinado_por_dtype"


def gem_group_stats(gem_group_arr):
    if gem_group_arr is None:
        return None, None
    vals, counts = np.unique(gem_group_arr, return_counts=True)
    return int(len(vals)), float(np.median(counts))


def add_row(**kw):
    base = dict(
        fuente=None, archivo=None, linea_celular_cruda=None, contexto=None,
        target_crudo=None, construccion=None, es_control=None, regla_control=None,
        n_celulas=None, n_gem_groups=None, mediana_celulas_por_gem_group=None,
        grano=None, tipo_valor_en_X=None, x_dtype=None, x_min_muestra=None,
        x_max_muestra=None, x_todo_entero_muestra=None,
        espacio_de_simbolo=None, n_genes_en_eje=None,
    )
    base.update(kw)
    ROWS.append(base)


# ============================================================ REPLOGLE (sc) ==
def do_replogle_sc(source_key, path, linea_cruda):
    log(f"replogle sc: {source_key}")
    with h5py.File(path, "r") as f:
        obs = f["obs"]
        gene = load_legacy_categorical(obs, "gene")
        gt = load_legacy_categorical(obs, "gene_transcript")
        gem = obs["gem_group"][()] if "gem_group" in obs else None
        n_vars = get_x_shape(f)[1]
        xstats = sample_dense_X_stats(f)

    is_ctrl = gene == "non-targeting"
    df = pd.DataFrame({"gene": gene, "gene_transcript": gt, "gem_group": gem, "es_control": is_ctrl})
    for (g_target, g_constr, ctrl), grp in df.groupby(["gene", "gene_transcript", "es_control"], dropna=False):
        n_gg, med_gg = gem_group_stats(grp["gem_group"].values if gem is not None else None)
        add_row(
            fuente=source_key, archivo=path, linea_celular_cruda=linea_cruda, contexto=None,
            target_crudo=g_target, construccion=g_constr,
            es_control=bool(ctrl), regla_control='gene == "non-targeting" (obs/gene, legacy categorical)',
            n_celulas=int(len(grp)), n_gem_groups=n_gg, mediana_celulas_por_gem_group=med_gg,
            grano="celula_individual", tipo_valor_en_X=tipo_valor_desde_stats(xstats),
            x_dtype=xstats["dtype"], x_min_muestra=xstats["min"], x_max_muestra=xstats["max"],
            x_todo_entero_muestra=xstats["todo_entero_en_muestra"],
            espacio_de_simbolo="ambos", n_genes_en_eje=n_vars,
        )


def do_replogle_bulk(source_key, path, linea_cruda):
    log(f"replogle bulk: {source_key}")
    with h5py.File(path, "r") as f:
        obs = f["obs"]
        idx_name = obs.attrs["_index"]
        idx = obs[idx_name][()]
        idx = np.array([v.decode() if isinstance(v, bytes) else v for v in idx])
        ncf = obs["num_cells_filtered"][()]
        n_vars = get_x_shape(f)[1]
        xstats = sample_dense_X_stats(f)

    is_ctrl = pd.Series(idx).str.contains("non-targeting").values
    for i in range(len(idx)):
        add_row(
            fuente=source_key, archivo=path, linea_celular_cruda=linea_cruda, contexto=None,
            target_crudo=idx[i], construccion=idx[i],
            es_control=bool(is_ctrl[i]),
            regla_control='"non-targeting" in gene_transcript (indice de obs, substring -- '
                          'core_control NO basta: tiene 71 falsos negativos en gwps, ver gwps_fiabilidad.md)',
            n_celulas=float(ncf[i]) if np.isfinite(ncf[i]) else None,
            n_gem_groups=None, mediana_celulas_por_gem_group=None,
            grano="pseudobulk", tipo_valor_en_X=tipo_valor_desde_stats(xstats),
            x_dtype=xstats["dtype"], x_min_muestra=xstats["min"], x_max_muestra=xstats["max"],
            x_todo_entero_muestra=xstats["todo_entero_en_muestra"],
            espacio_de_simbolo="ambos", n_genes_en_eje=n_vars,
        )


# ============================================================ NADIG (sc) ==
def do_nadig(source_key, path, linea_cruda):
    log(f"nadig: {source_key}")
    with h5py.File(path, "r") as f:
        obs = f["obs"]
        gene = load_legacy_categorical(obs, "gene")
        gt = load_legacy_categorical(obs, "gene_transcript")
        gem = obs["gem_group"][()] if "gem_group" in obs else None
        n_vars = get_x_shape(f)[1]
        xstats = sample_dense_X_stats(f)

    is_ctrl = gene == "non-targeting"
    df = pd.DataFrame({"gene": gene, "gene_transcript": gt, "gem_group": gem, "es_control": is_ctrl})
    for (g_target, g_constr, ctrl), grp in df.groupby(["gene", "gene_transcript", "es_control"], dropna=False):
        n_gg, med_gg = gem_group_stats(grp["gem_group"].values if gem is not None else None)
        add_row(
            fuente=source_key, archivo=path, linea_celular_cruda=linea_cruda, contexto=None,
            target_crudo=g_target, construccion=g_constr,
            es_control=bool(ctrl), regla_control='gene == "non-targeting" (obs/gene, legacy categorical)',
            n_celulas=int(len(grp)), n_gem_groups=n_gg, mediana_celulas_por_gem_group=med_gg,
            grano="celula_individual", tipo_valor_en_X=tipo_valor_desde_stats(xstats),
            x_dtype=xstats["dtype"], x_min_muestra=xstats["min"], x_max_muestra=xstats["max"],
            x_todo_entero_muestra=xstats["todo_entero_en_muestra"],
            espacio_de_simbolo="ambos", n_genes_en_eje=n_vars,
        )


# ============================================================ ZENODO convertido ==
def do_convertido(source_key, path):
    log(f"convertido: {source_key}")
    with h5py.File(path, "r") as f:
        obs = f["obs"]
        cell_type = read_new_categorical(obs["cell_type"])
        gene = read_new_categorical(obs["gene"])
        guide = read_new_categorical(obs["guide"])
        batch = read_new_categorical(obs["Batch_info"]) if "Batch_info" in obs else None
        n_vars = get_x_shape(f)[1]
        xstats = sample_csr_X_stats(f)

    is_ctrl = gene == "NT"
    df = pd.DataFrame({"cell_type": cell_type, "gene": gene, "guide": guide, "batch": batch, "es_control": is_ctrl})
    for (linea, g_target, g_constr, ctrl), grp in df.groupby(["cell_type", "gene", "guide", "es_control"], dropna=False):
        n_gg, med_gg = gem_group_stats(grp["batch"].values if batch is not None else None)
        add_row(
            fuente=source_key, archivo=path, linea_celular_cruda=linea, contexto=None,
            target_crudo=g_target, construccion=g_constr,
            es_control=bool(ctrl), regla_control='gene == "NT" (obs/gene, categorical nuevo estilo)',
            n_celulas=int(len(grp)), n_gem_groups=n_gg, mediana_celulas_por_gem_group=med_gg,
            grano="celula_individual", tipo_valor_en_X=tipo_valor_desde_stats(xstats),
            x_dtype=xstats["dtype"], x_min_muestra=xstats["min"], x_max_muestra=xstats["max"],
            x_todo_entero_muestra=xstats["todo_entero_en_muestra"],
            espacio_de_simbolo="simbolo", n_genes_en_eje=n_vars,
        )


# ============================================================ RETO contexto (solo control) ==
def do_reto_context(source_key, path, contexto_esperado):
    log(f"reto context: {source_key}")
    with h5py.File(path, "r") as f:
        obs = f["obs"]
        target_gene = read_new_categorical(obs["target_gene"])
        ntc_id = read_new_categorical(obs["ntc_id"])
        context = read_new_categorical(obs["context"])
        n_vars = get_x_shape(f)[1]
        xstats = sample_csr_X_stats(f)

    assert set(np.unique(context)) == {contexto_esperado}, f"contexto inesperado en {path}: {np.unique(context)}"
    assert set(np.unique(target_gene)) == {"non-targeting"}, f"target_gene no es solo control en {path}"

    df = pd.DataFrame({"ntc_id": ntc_id, "target_gene": target_gene})
    n_gg, med_gg = gem_group_stats(ntc_id)
    add_row(
        fuente=source_key, archivo=path, linea_celular_cruda=None, contexto=contexto_esperado,
        target_crudo="non-targeting", construccion=None,
        es_control=True, regla_control='target_gene == "non-targeting" (unica categoria presente -- '
                                        "este archivo NO contiene celulas perturbadas, solo el pool de control; "
                                        "ver manifest.json: control_label, y ground_truth_cells no esta en el .h5ad)",
        n_celulas=int(len(df)), n_gem_groups=n_gg, mediana_celulas_por_gem_group=med_gg,
        grano="celula_individual", tipo_valor_en_X=tipo_valor_desde_stats(xstats),
        x_dtype=xstats["dtype"], x_min_muestra=xstats["min"], x_max_muestra=xstats["max"],
        x_todo_entero_muestra=xstats["todo_entero_en_muestra"],
        espacio_de_simbolo="simbolo", n_genes_en_eje=n_vars,
    )


def do_reto_pert_counts(source_key, path, partition_contexts):
    log(f"reto pert_counts: {source_key}")
    df = pd.read_csv(path)
    for g in df["target_gene"]:
        add_row(
            fuente=source_key, archivo=path, linea_celular_cruda=None,
            contexto=f"aplica_a_{'/'.join(partition_contexts)}_por_manifest_no_desglosado_por_contexto",
            target_crudo=g, construccion=None,
            es_control=False, regla_control="no aplica -- este csv solo lista genes blanco, no celulas; "
                                             "NO tiene columna n_cells (a diferencia de vcc2025)",
            n_celulas=None, n_gem_groups=None, mediana_celulas_por_gem_group=None,
            grano=None, tipo_valor_en_X=None, x_dtype=None, x_min_muestra=None, x_max_muestra=None,
            x_todo_entero_muestra=None, espacio_de_simbolo="simbolo", n_genes_en_eje=None,
        )


# ============================================================ VCC2025 ==
def do_vcc2025(source_key, path):
    log(f"vcc2025: {source_key}")
    with h5py.File(path, "r") as f:
        obs = f["obs"]
        target_gene = read_new_categorical(obs["target_gene"])
        guide_id = read_new_categorical(obs["guide_id"])
        batch = read_new_categorical(obs["batch"])
        n_vars = get_x_shape(f)[1]
        xstats = sample_csr_X_stats(f)

    is_ctrl = target_gene == "non-targeting"
    df = pd.DataFrame({"target_gene": target_gene, "guide_id": guide_id, "batch": batch, "es_control": is_ctrl})
    for (g_target, g_constr, ctrl), grp in df.groupby(["target_gene", "guide_id", "es_control"], dropna=False):
        n_gg, med_gg = gem_group_stats(grp["batch"].values)
        add_row(
            fuente=source_key, archivo=path, linea_celular_cruda=None, contexto=None,
            target_crudo=g_target, construccion=g_constr,
            es_control=bool(ctrl), regla_control='target_gene == "non-targeting" (categorical nuevo estilo); '
                                                  "linea celular NO declarada en ninguna columna de obs",
            n_celulas=int(len(grp)), n_gem_groups=n_gg, mediana_celulas_por_gem_group=med_gg,
            grano="celula_individual", tipo_valor_en_X=tipo_valor_desde_stats(xstats),
            x_dtype=xstats["dtype"], x_min_muestra=xstats["min"], x_max_muestra=xstats["max"],
            x_todo_entero_muestra=xstats["todo_entero_en_muestra"],
            espacio_de_simbolo="simbolo", n_genes_en_eje=n_vars,
        )


def do_vcc2025_pert_counts(source_key, path):
    log(f"vcc2025 pert_counts: {source_key}")
    df = pd.read_csv(path)
    for _, row in df.iterrows():
        add_row(
            fuente=source_key, archivo=path, linea_celular_cruda=None, contexto=None,
            target_crudo=row["target_gene"], construccion=None,
            es_control=(row["target_gene"] == "non-targeting"),
            regla_control='target_gene == "non-targeting" (si aparece); csv de resumen, no celulas',
            n_celulas=float(row["n_cells"]), n_gem_groups=None, mediana_celulas_por_gem_group=None,
            grano=None, tipo_valor_en_X=None, x_dtype=None, x_min_muestra=None, x_max_muestra=None,
            x_todo_entero_muestra=None, espacio_de_simbolo="simbolo", n_genes_en_eje=None,
        )


# ============================================================ ORION ==
def do_orion(source_key, path, linea_cruda):
    log(f"orion: {source_key}")
    t = pq.read_table(path, columns=["target", "feature", "target_ncells"])
    df = t.to_pandas()
    n_genes_eje = df["feature"].nunique()
    per_target = df.groupby("target")["target_ncells"].agg(["nunique", "first"])
    n_inconsistent = int((per_target["nunique"] > 1).sum())
    feat_sample = df["feature"].drop_duplicates().head(2000)
    looks_ensembl = feat_sample.str.match(r"^ENSG\d+$")
    espacio = "ambos" if looks_ensembl.any() and (~looks_ensembl).any() else ("ensembl" if looks_ensembl.all() else "simbolo")
    for target, row in per_target.iterrows():
        add_row(
            fuente=source_key, archivo=path, linea_celular_cruda=linea_cruda, contexto=None,
            target_crudo=target, construccion=None,
            es_control=False, regla_control="no aplica -- cada fila ya es perturbacion vs. referencia "
                                             "contrastada (target_mean vs ref_mean); no hay celulas de "
                                             "control expuestas por separado",
            n_celulas=float(row["first"]), n_gem_groups=None, mediana_celulas_por_gem_group=None,
            grano="resumen_de_pares(target,feature)_no_es_celula_ni_pseudobulk",
            tipo_valor_en_X=None, x_dtype=None, x_min_muestra=None, x_max_muestra=None, x_todo_entero_muestra=None,
            espacio_de_simbolo=espacio, n_genes_en_eje=int(n_genes_eje),
        )
    log(f"  orion {source_key}: n_targets={len(per_target)} target_ncells inconsistentes entre features={n_inconsistent} espacio={espacio}")


# ============================================================ FENG ==
def do_feng_simple(source_key, path, linea_cruda):
    log(f"feng: {source_key}")
    df = pd.read_csv(path, sep="\t")
    n_genes_eje = df["Expressed_Gene_Symbol"].nunique()
    for target in df["Target"].unique():
        add_row(
            fuente=source_key, archivo=path, linea_celular_cruda=linea_cruda, contexto=None,
            target_crudo=target, construccion=None,
            es_control=False, regla_control="no aplica -- tabla LFC por gen, ya contrastada; sin celulas",
            n_celulas=None, n_gem_groups=None, mediana_celulas_por_gem_group=None,
            grano="resumen_LFC_por_gen_no_es_celula_ni_pseudobulk",
            tipo_valor_en_X=None, x_dtype=None, x_min_muestra=None, x_max_muestra=None, x_todo_entero_muestra=None,
            espacio_de_simbolo="ambos", n_genes_en_eje=int(n_genes_eje),
        )


def do_feng_perline(source_key, path):
    log(f"feng perline: {source_key}")
    df = pd.read_csv(path, sep="\t")
    n_genes_eje = df["Expressed_Gene_Symbol"].nunique()
    for (linea, target), grp in df.groupby(["Cell_Line", "Target"]):
        add_row(
            fuente=source_key, archivo=path, linea_celular_cruda=linea, contexto=None,
            target_crudo=target, construccion=None,
            es_control=False, regla_control="no aplica -- tabla LFC por gen y linea, ya contrastada; sin celulas",
            n_celulas=None, n_gem_groups=None, mediana_celulas_por_gem_group=None,
            grano="resumen_LFC_por_gen_no_es_celula_ni_pseudobulk",
            tipo_valor_en_X=None, x_dtype=None, x_min_muestra=None, x_max_muestra=None, x_todo_entero_muestra=None,
            espacio_de_simbolo="ambos", n_genes_en_eje=int(n_genes_eje),
        )


def main():
    t0 = time.time()

    do_replogle_sc("replogle_k562_essential_sc", "/datos/ref/replogle/k562_essential_raw_sc.h5ad", "k562_essential")
    do_replogle_bulk("replogle_k562_essential_bulk", "/datos/ref/replogle/k562_essential_raw_bulk.h5ad", "k562_essential")
    do_replogle_sc("replogle_k562_gwps_sc", "/datos/ref/replogle/k562_gwps_raw_sc.h5ad", "k562_gwps")
    do_replogle_bulk("replogle_k562_gwps_bulk", "/datos/ref/replogle/k562_gwps_raw_bulk.h5ad", "k562_gwps")
    do_replogle_sc("replogle_rpe1_sc", "/datos/ref/replogle/rpe1_raw_sc.h5ad", "rpe1")
    do_replogle_bulk("replogle_rpe1_bulk", "/datos/ref/replogle/rpe1_raw_bulk.h5ad", "rpe1")

    do_nadig("nadig_hepg2", "/datos/ref/nadig/hepg2_raw_sc.h5ad", "hepg2")
    do_nadig("nadig_jurkat", "/datos/ref/nadig/jurkat_raw_sc.h5ad", "jurkat")

    for name in ["IFNB", "IFNG", "INS", "TGFB", "TNFA"]:
        do_convertido(f"zenodo_{name}", f"/datos/ref/convertido/Seurat_object_{name}_Perturb_seq.h5ad")

    do_orion("orion_hct116", "/datos/ref/orion/hct116.parquet", "hct116")
    do_orion("orion_hek293t", "/datos/ref/orion/hek293t.parquet", "hek293t")

    do_feng_simple("feng_genomewide", "/datos/ref/feng/GenomeWideScreen_LFC_byGene.tsv.gz", None)
    do_feng_simple("feng_targeted", "/datos/ref/feng/TargetedScreen_LFC_byGene.tsv.gz", None)
    do_feng_perline("feng_targeted_perline", "/datos/ref/feng/TargetedScreen_LFC_byGene-perLine.tsv.gz")

    do_reto_context("reto_val_context_A", "/datos/reto/val/context_A.h5ad", "A")
    do_reto_context("reto_val_context_B", "/datos/reto/val/context_B.h5ad", "B")
    do_reto_context("reto_val_context_C", "/datos/reto/val/context_C.h5ad", "C")
    do_reto_context("reto_final_context_D", "/datos/reto/final/context_D.h5ad", "D")
    do_reto_context("reto_final_context_E", "/datos/reto/final/context_E.h5ad", "E")
    do_reto_context("reto_final_context_F", "/datos/reto/final/context_F.h5ad", "F")
    do_reto_pert_counts("reto_val_pert_counts", "/datos/reto/val/pert_counts.csv", ["A", "B", "C"])
    do_reto_pert_counts("reto_final_pert_counts", "/datos/reto/final/pert_counts.csv", ["D", "E", "F"])

    do_vcc2025("vcc2025_train", "/datos/vcc2025/train/adata_Training.h5ad")
    do_vcc2025("vcc2025_validation", "/datos/vcc2025/validation/adata_Validation.h5ad")
    do_vcc2025("vcc2025_test", "/datos/vcc2025/test/adata_Test.h5ad")
    do_vcc2025_pert_counts("vcc2025_train_pert_counts", "/datos/vcc2025/train/pert_counts_Training.csv")
    do_vcc2025_pert_counts("vcc2025_validation_pert_counts", "/datos/vcc2025/validation/pert_counts_Validation.csv")
    do_vcc2025_pert_counts("vcc2025_test_pert_counts", "/datos/vcc2025/test/pert_counts_Test.csv")

    df = pd.DataFrame(ROWS)
    log(f"total filas: {len(df)}")
    import os
    os.makedirs("/home/ubuntu/vcc/corpus", exist_ok=True)
    df.to_parquet(OUT_PARQUET, index=False)
    log(f"escrito {OUT_PARQUET} -- elapsed={time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
