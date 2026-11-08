#!/usr/bin/env python
"""PARTE 0 del encargo `ntc_nulo`: las tres verificaciones de validez.

a) controls-final bajado y completo; sha256, formas, conteos por contexto;
   manifest.json JUNTO a pert_counts.csv.
b) cobertura POR SIMBOLO contra gwps (k562_gwps_raw_bulk.h5ad).
c) nnz por celula de control por contexto contra el tope 4.75e9/360000, y
   medicion del Poisson aditivo de `cand01_covar._reescala_entero`.

Nada de ~/vcc/benchmark/ ni de /datos/benchmark_VEDADO_*. Sin envios.
"""
import hashlib, json, os, sys
import numpy as np
import scipy.sparse as sp
import h5py

sys.path.insert(0, "/home/ubuntu/vcc/modelo/src")

RETO = "/datos/reto"
GWPS = "/datos/ref/replogle/k562_gwps_raw_bulk.h5ad"
OUT = "/home/ubuntu/vcc/modelo/out/ntc_nulo/p0.json"

CTX = [("A", "val"), ("B", "val"), ("C", "val"),
       ("D", "final"), ("E", "final"), ("F", "final")]

# tope duro del contrato: HARD_CAP_NNZ / filas de la entrega
HARD_CAP_NNZ = 4_750_000_000
FILAS_ENTREGA = 300 * 400 * 3          # 360000
TOPE_POR_CELULA = HARD_CAP_NNZ / FILAS_ENTREGA


def sha256(p, buf=1 << 22):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        while True:
            b = fh.read(buf)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def cats(f, col):
    g = f[f"obs/{col}"]
    c = np.array([x.decode() if isinstance(x, bytes) else x
                  for x in g["categories"][:]])
    return c, g["codes"][:]


def var_names_h5(f):
    v = f["var/_index/values"][:]
    return np.array([x.decode() if isinstance(x, bytes) else x for x in v])


# ---------------------------------------------------------------- a)
def parte_a():
    res = {"archivos": [], "contextos": [], "manifest": {}, "disparadores": []}
    archivos = []
    for part in ("val", "final"):
        d = os.path.join(RETO, part)
        for n in sorted(os.listdir(d)):
            archivos.append((part, n, os.path.join(d, n)))
    for part, n, p in archivos:
        st = os.stat(p)
        res["archivos"].append({
            "particion": part, "archivo": n, "bytes": st.st_size,
            "sha256": sha256(p),
        })
        sys.stderr.write(f"  sha256 {part}/{n}\n")

    for letra, part in CTX:
        p = os.path.join(RETO, part, f"context_{letra}.h5ad")
        with h5py.File(p, "r") as f:
            shape = tuple(int(x) for x in f["X"].attrs["shape"])
            nnz_total = int(f["X/data"].shape[0])
            indptr = f["X/indptr"][:].astype(np.int64)
            tg_c, tg_k = cats(f, "target_gene")
            cx_c, cx_k = cats(f, "context")
            nt_c, nt_k = cats(f, "ntc_id")
            vn = var_names_h5(f)
        porg = np.bincount(nt_k, minlength=len(nt_c))
        res["contextos"].append({
            "contexto": letra, "particion": part, "archivo": os.path.basename(p),
            "shape": list(shape), "n_celulas": int(shape[0]), "n_genes": int(shape[1]),
            "nnz_total": nnz_total,
            "nnz_consistente_con_indptr": bool(indptr[-1] == nnz_total),
            "target_gene_categorias": tg_c.tolist(),
            "context_categorias": cx_c.tolist(),
            "ntc_id_n_categorias": int(len(nt_c)),
            "ntc_id_celulas_min": int(porg.min()), "ntc_id_celulas_max": int(porg.max()),
            "ntc_id_todas_400": bool((porg == 400).all()),
            "ntc_id_suma": int(porg.sum()),
            "primer_gen": str(vn[0]), "ultimo_gen": str(vn[-1]),
        })
        sys.stderr.write(f"  contexto {letra}: {shape} nnz={nnz_total}\n")

    # manifest JUNTO a pert_counts.csv, en las dos particiones
    for part in ("val", "final"):
        d = os.path.join(RETO, part)
        mp, pp = os.path.join(d, "manifest.json"), os.path.join(d, "pert_counts.csv")
        res["manifest"][part] = {
            "manifest_path": mp, "pert_counts_path": pp,
            "mismo_directorio": os.path.dirname(mp) == os.path.dirname(pp),
            "manifest_existe": os.path.exists(mp),
            "pert_counts_existe": os.path.exists(pp),
            "contenido": json.load(open(mp)),
        }

    # los 46 ntc_id: mismo conjunto en los seis?
    conj = {}
    for letra, part in CTX:
        with h5py.File(os.path.join(RETO, part, f"context_{letra}.h5ad"), "r") as f:
            c, _ = cats(f, "ntc_id")
        conj[letra] = sorted(c.tolist())
    base = conj["A"]
    res["ntc_id_identicos_los_6"] = all(conj[k] == base for k in conj)
    res["ntc_id_lista"] = base
    return res


# ---------------------------------------------------------------- b)
def carga_gwps_simbolos():
    with h5py.File(GWPS, "r") as f:
        gt = np.array([x.decode() if isinstance(x, bytes) else x
                       for x in f["obs/gene_transcript"][:]])
        core = f["obs/core_control"][:]
        catg = np.array([x.decode() if isinstance(x, bytes) else x
                         for x in f["var/__categories/gene_name"][:]])
        vn = catg[f["var/gene_name"][:]]
    partes = [x.split("_") for x in gt]
    g_gen = np.array([p[1] if len(p) > 1 else "" for p in partes])
    g_con = np.array([p[2] if len(p) > 2 else "" for p in partes])
    u, c = np.unique(vn, return_counts=True)
    dup = sorted(u[c > 1].tolist())
    return {"gene_transcript": gt, "g_gen": g_gen, "g_con": g_con,
            "core_control": core, "var_gene_name": vn, "dup": dup}


def parte_b(g):
    res = {}
    gn = [l.strip() for l in open(os.path.join(RETO, "final", "gene_names.csv"))
          if l.strip()]
    gn_val = [l.strip() for l in open(os.path.join(RETO, "val", "gene_names.csv"))
              if l.strip()]
    res["gene_names_final_n"] = len(gn)
    res["gene_names_val_n"] = len(gn_val)
    res["gene_names_val_igual_final"] = gn == gn_val
    res["gene_names_duplicados_internos"] = int(len(gn) - len(set(gn)))

    med_todos = set(g["var_gene_name"].tolist())
    med_unicos = med_todos - set(g["dup"])
    res["gwps_var_columnas"] = int(len(g["var_gene_name"]))
    res["gwps_var_nombres_unicos"] = int(len(med_todos))
    res["gwps_var_nombres_duplicados"] = g["dup"]

    en = [s for s in gn if s in med_todos]
    en_u = [s for s in gn if s in med_unicos]
    res["gene_names_en_gwps_medidos"] = len(en)
    res["gene_names_en_gwps_medidos_sin_dup"] = len(en_u)
    res["gene_names_fuera_de_gwps"] = len(gn) - len(en)

    # perfil = fila de perturbacion en gwps
    perf = set(g["g_gen"].tolist()) - {""}
    res["gwps_filas_totales"] = int(len(g["gene_transcript"]))
    res["gwps_filas_core_control"] = int(g["core_control"].sum())
    res["gwps_blancos_con_perfil_unicos"] = int(len(perf))

    for part in ("val", "final"):
        import csv
        pth = os.path.join(RETO, part, "pert_counts.csv")
        with open(pth) as fh:
            rd = list(csv.DictReader(fh))
        col = "target_gene" if "target_gene" in rd[0] else list(rd[0])[0]
        perts = [r[col].strip() for r in rd]
        sinp = sorted(p for p in perts if p not in perf)
        conp = [p for p in perts if p in perf]
        # cuantas tienen 1 o 2 construcciones
        nconstr = {}
        for p in conp:
            nconstr[p] = int(np.sum(g["g_gen"] == p))
        res[f"pert_{part}"] = {
            "columna_leida": col,
            "n_perts": len(perts),
            "n_duplicadas": len(perts) - len(set(perts)),
            "con_perfil_en_gwps": len(conp),
            "sin_perfil_en_gwps": len(sinp),
            "lista_sin_perfil": sinp,
            "con_2_construcciones": sorted(k for k, v in nconstr.items() if v > 1),
            "en_gene_names": int(sum(1 for p in perts if p in set(gn))),
            "fuera_de_gene_names": sorted(p for p in perts if p not in set(gn)),
        }
    return res


# ---------------------------------------------------------------- c)
def parte_c_nnz():
    res = []
    for letra, part in CTX:
        with h5py.File(os.path.join(RETO, part, f"context_{letra}.h5ad"), "r") as f:
            indptr = f["X/indptr"][:].astype(np.int64)
        nnz = np.diff(indptr)
        res.append({
            "contexto": letra, "n_celulas": int(nnz.size),
            "nnz_media": float(nnz.mean()), "nnz_mediana": float(np.median(nnz)),
            "nnz_p99": float(np.percentile(nnz, 99)), "nnz_max": int(nnz.max()),
            "nnz_min": int(nnz.min()), "nnz_sd": float(nnz.std(ddof=1)),
            "tope_por_celula": TOPE_POR_CELULA,
            "margen_sobre_media": float(TOPE_POR_CELULA / nnz.mean()),
            "margen_sobre_max": float(TOPE_POR_CELULA / nnz.max()),
            "frac_tope_media": float(nnz.mean() / TOPE_POR_CELULA),
            "frac_tope_max": float(nnz.max() / TOPE_POR_CELULA),
        })
        sys.stderr.write(f"  nnz {letra}: media={nnz.mean():.1f} max={nnz.max()}\n")
    return res


def lee_bloque_csr(path, filas):
    """Lee las filas `filas` (ordenadas) de un CSR en disco, sin densificar."""
    filas = np.asarray(sorted(filas), dtype=np.int64)
    with h5py.File(path, "r") as f:
        shape = tuple(int(x) for x in f["X"].attrs["shape"])
        indptr = f["X/indptr"][:].astype(np.int64)
        D, I = f["X/data"], f["X/indices"]
        datos, indices, nuevo = [], [], [0]
        for r in filas:
            a, b = int(indptr[r]), int(indptr[r + 1])
            datos.append(D[a:b]); indices.append(I[a:b]); nuevo.append(nuevo[-1] + (b - a))
    return sp.csr_matrix(
        (np.concatenate(datos), np.concatenate(indices), np.array(nuevo, dtype=np.int64)),
        shape=(len(filas), shape[1]))


def parte_c_poisson():
    """Aplica `_reescala_entero` a 400 celulas de control con las dos emisiones
    de densidad que se venian usando, y cuenta el nnz resultante."""
    from cand01_covar import _reescala_entero, construye_respuesta
    from gwps_emision import construye_ratios

    LETRA, PART = "D", "final"
    path = os.path.join(RETO, PART, f"context_{LETRA}.h5ad")
    with h5py.File(path, "r") as f:
        vn = var_names_h5(f)
        ncel = int(f["X"].attrs["shape"][0])
    import csv
    with open(os.path.join(RETO, PART, "pert_counts.csv")) as fh:
        rd = list(csv.DictReader(fh))
    col = "target_gene" if "target_gene" in rd[0] else list(rd[0])[0]
    perts = [r[col].strip() for r in rd]
    vnset = set(vn.tolist())

    SEM = 20274999                      # declarada, fija, para esta medicion sola
    rng = np.random.default_rng(SEM)
    filas = np.sort(rng.choice(ncel, size=400, replace=False))
    bloque = lee_bloque_csr(path, filas)
    nnz0 = int(bloque.nnz)
    sys.stderr.write(f"  muestra 400 celulas ctx {LETRA}: nnz0={nnz0}\n")

    out = {"contexto": LETRA, "semilla_muestra": SEM, "n_celulas": 400,
           "nnz_antes": nnz0,
           "nnz_antes_por_celula_media": nnz0 / 400.0,
           "tope_por_celula": TOPE_POR_CELULA, "brazos": []}

    # ---- emision 1: cand01 (ratios estimados EN CONTEXTO sobre el pool entero)
    pool = lee_bloque_csr(path, np.arange(ncel))
    tg = [p for p in perts if p in vnset][:3]
    sys.stderr.write(f"  cand01 construye_respuesta sobre pool {pool.shape} blancos={tg}\n")
    ratios_c, diag_c = construye_respuesta(pool, vn, tg)
    del pool
    for g in tg:
        rat = ratios_c[g]
        r2 = np.random.default_rng(SEM + 1)
        nb = _reescala_entero(bloque, rat, r2)
        out["brazos"].append({
            "emision": "cand01_covar.construye_respuesta (en contexto)",
            "blanco": g,
            "n_genes_r_mayor_1": int((rat > 1.0).sum()),
            "n_genes_r_menor_1": int((rat < 1.0).sum()),
            "n_genes_r_igual_1": int((rat == 1.0).sum()),
            "r_min": float(rat.min()), "r_max": float(rat.max()),
            "nnz_despues": int(nb.nnz),
            "delta_nnz": int(nb.nnz) - nnz0,
            "nnz_despues_por_celula_media": float(nb.nnz) / 400.0,
            "nnz_max_por_celula": int(np.diff(nb.indptr).max()),
            "frac_tope_max": float(np.diff(nb.indptr).max() / TOPE_POR_CELULA),
            "revienta_tope": bool(np.diff(nb.indptr).max() > TOPE_POR_CELULA),
            "suma_antes": float(bloque.sum()), "suma_despues": float(nb.sum()),
        })
        sys.stderr.write(f"    cand01 {g}: nnz {nnz0} -> {nb.nnz} (r>1: {(rat>1).sum()})\n")
        del nb

    # ---- emision 2: gwps (ratios externos, los que mas genes suben)
    tg2 = [p for p in perts if p in vnset][:3]
    sys.stderr.write(f"  gwps construye_ratios blancos={tg2}\n")
    ratios_g, diag_g = construye_ratios(vn, tg2)
    for g in tg2:
        rat = ratios_g[g]
        r2 = np.random.default_rng(SEM + 2)
        nb = _reescala_entero(bloque, rat, r2)
        out["brazos"].append({
            "emision": "gwps_emision.construye_ratios (externa)",
            "blanco": g,
            "n_genes_r_mayor_1": int((rat > 1.0).sum()),
            "n_genes_r_menor_1": int((rat < 1.0).sum()),
            "n_genes_r_igual_1": int((rat == 1.0).sum()),
            "r_min": float(rat.min()), "r_max": float(rat.max()),
            "nnz_despues": int(nb.nnz),
            "delta_nnz": int(nb.nnz) - nnz0,
            "nnz_despues_por_celula_media": float(nb.nnz) / 400.0,
            "nnz_max_por_celula": int(np.diff(nb.indptr).max()),
            "frac_tope_max": float(np.diff(nb.indptr).max() / TOPE_POR_CELULA),
            "revienta_tope": bool(np.diff(nb.indptr).max() > TOPE_POR_CELULA),
            "suma_antes": float(bloque.sum()), "suma_despues": float(nb.sum()),
        })
        sys.stderr.write(f"    gwps {g}: nnz {nnz0} -> {nb.nnz} (r>1: {(rat>1).sum()})\n")
        del nb

    # ---- cota estructural: el brazo Poisson escribe SOLO en entradas ya guardadas
    out["cota_estructural"] = {
        "afirmacion": ("_reescala_entero opera sobre b.data (entradas ya guardadas) y "
                       "llama eliminate_zeros(); el brazo Poisson no puede crear un "
                       "no-cero nuevo, asi que nnz_despues <= nnz_antes siempre"),
        "verificado_en_todos_los_brazos": all(b["delta_nnz"] <= 0 for b in out["brazos"]),
    }
    return out


def main():
    sys.stderr.write("PARTE 0a ...\n"); a = parte_a()
    sys.stderr.write("PARTE 0b ...\n"); g = carga_gwps_simbolos(); b = parte_b(g)
    sys.stderr.write("PARTE 0c nnz ...\n"); c1 = parte_c_nnz()
    sys.stderr.write("PARTE 0c poisson ...\n"); c2 = parte_c_poisson()
    todo = {"a_validez_datos": a, "b_cobertura_simbolos": b,
            "c_nnz_por_contexto": c1, "c_poisson": c2,
            "constantes": {"HARD_CAP_NNZ": HARD_CAP_NNZ,
                           "filas_entrega": FILAS_ENTREGA,
                           "tope_por_celula": TOPE_POR_CELULA}}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(todo, open(OUT, "w"), indent=2, sort_keys=True, default=str)
    sys.stderr.write(f"\nlisto: {OUT}\n")


if __name__ == "__main__":
    main()
