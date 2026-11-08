#!/usr/bin/env python
"""PARTE 0b corregida: `gene_names.csv` trae CABECERA (`gene_name`), asi que la
lectura por lineas de `p0.py` contaba 18534 simbolos. Aqui se lee con csv.reader
saltando la cabecera, y se cruza el orden contra `var/_index` del h5ad."""
import csv, json, os
import numpy as np, h5py

RETO = "/datos/reto"
GWPS = "/datos/ref/replogle/k562_gwps_raw_bulk.h5ad"
OUT = "/home/ubuntu/vcc/modelo/out/ntc_nulo/p0b.json"


def lee_col(p):
    with open(p) as fh:
        rd = csv.reader(fh)
        cab = next(rd)
        return cab, [r[0].strip() for r in rd if r and r[0].strip()]


with h5py.File(GWPS, "r") as f:
    gt = np.array([x.decode() if isinstance(x, bytes) else x
                   for x in f["obs/gene_transcript"][:]])
    core = f["obs/core_control"][:]
    catg = np.array([x.decode() if isinstance(x, bytes) else x
                     for x in f["var/__categories/gene_name"][:]])
    vn_g = catg[f["var/gene_name"][:]]
partes = [x.split("_") for x in gt]
g_gen = np.array([p[1] if len(p) > 1 else "" for p in partes])
u, c = np.unique(vn_g, return_counts=True)
dup = set(u[c > 1].tolist())
medidos = set(vn_g.tolist())
medidos_sin_dup = medidos - dup
perfiles = set(g_gen.tolist()) - {""}

res = {"gwps": {"filas": int(len(gt)), "filas_core_control": int(core.sum()),
                "var_columnas": int(len(vn_g)),
                "var_nombres_unicos": int(len(medidos)),
                "var_nombres_duplicados": sorted(dup),
                "blancos_con_perfil_unicos": int(len(perfiles))}}

for part, letras in (("val", "ABC"), ("final", "DEF")):
    cab, gn = lee_col(os.path.join(RETO, part, "gene_names.csv"))
    with h5py.File(os.path.join(RETO, part, f"context_{letras[0]}.h5ad"), "r") as f:
        v = f["var/_index/values"][:]
        vn_h = [x.decode() if isinstance(x, bytes) else x for x in v]
    en = [s for s in gn if s in medidos]
    en_sd = [s for s in gn if s in medidos_sin_dup]
    cab2, perts = lee_col(os.path.join(RETO, part, "pert_counts.csv"))
    sinp = sorted(p for p in perts if p not in perfiles)
    res[part] = {
        "gene_names_cabecera": cab,
        "gene_names_n_simbolos": len(gn),
        "gene_names_unicos": len(set(gn)),
        "orden_identico_a_var_index": gn == vn_h,
        "n_genes_h5ad": len(vn_h),
        "en_gwps_medidos": len(en),
        "en_gwps_medidos_sin_duplicados": len(en_sd),
        "fuera_de_gwps": len(gn) - len(en),
        "cobertura_pct": round(100.0 * len(en) / len(gn), 4),
        "pert_cabecera": cab2,
        "pert_n": len(perts),
        "pert_con_perfil_en_gwps": len(perts) - len(sinp),
        "pert_sin_perfil_en_gwps": len(sinp),
        "pert_lista_sin_perfil": sinp,
        "pert_todas_en_gene_names": all(p in set(gn) for p in perts),
    }

json.dump(res, open(OUT, "w"), indent=2, sort_keys=True)
for k in res:
    print(k, json.dumps({a: b for a, b in res[k].items()
                         if not isinstance(b, list) or len(b) < 5}, sort_keys=True))
print("\nsin perfil val :", res["val"]["pert_lista_sin_perfil"])
print("sin perfil final:", res["final"]["pert_lista_sin_perfil"])
print("\nescrito", OUT)
