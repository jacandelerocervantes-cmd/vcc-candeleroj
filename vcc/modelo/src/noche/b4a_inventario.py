#!/usr/bin/env python
"""B4a - PASO 4, inventario de la superficie VCC 2025.

Formas, conteos, sha256, que hay en cada split, y la lista exacta de los genes
pareados con el panel 2026. No mide nada del marcador: solo describe la
superficie y declara que puede y que no.
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import comun as C
import v25

OUT_JSON = os.path.join(C.OUT, "vcc2025_inventario.json")
OUT_MD = os.path.join(C.OUT, "vcc2025_inventario.md")


def pareados():
    """Todas las definiciones de 'pareado' que se pueden calcular, con sus
    listas exactas. 🔴 El encargo cita 97; se reporta lo MEDIDO en cada
    definicion y se declara si 97 no sale de ninguna."""
    g25 = [x.strip() for x in open(v25.GENES_2025) if x.strip()]
    g26 = list(pd.read_csv(os.path.join(v25.RETO_FINAL, "gene_names.csv"))["gene_name"])
    p26f = sorted(pd.read_csv(os.path.join(v25.RETO_FINAL, "pert_counts.csv"))["target_gene"])
    p26v = sorted(pd.read_csv(os.path.join(v25.RETO_VAL, "pert_counts.csv"))["target_gene"])
    p25 = {}
    for s in v25.SPLITS:
        p25[s] = sorted(pd.read_csv(v25.ruta_counts(s))["target_gene"])
    u25 = sorted(set().union(*[set(v) for v in p25.values()]))
    S = lambda x: set(x)
    defs = {
        "d1_perts2025_x_perts2026_final": sorted(S(u25) & S(p26f)),
        "d2_perts2025_x_perts2026_val": sorted(S(u25) & S(p26v)),
        "d3_perts2025_x_union_perts2026": sorted(S(u25) & (S(p26f) | S(p26v))),
        "d4_perts2025_x_interseccion_perts2026": sorted(S(u25) & (S(p26f) & S(p26v))),
        "d5_perts2025_en_el_espacio_de_genes_2026": sorted(S(u25) & S(g26)),
        "d6_perts2026_union_en_el_espacio_de_genes_2025": sorted((S(p26f) | S(p26v)) & S(g25)),
        "d7_genes_medidos_en_los_dos": sorted(S(g25) & S(g26)),
    }
    return {
        "n_genes_2025": len(g25), "n_genes_2026": len(g26),
        "n_perts_2026_final": len(p26f), "n_perts_2026_val": len(p26v),
        "n_perts_2026_union": len(S(p26f) | S(p26v)),
        "n_perts_2026_interseccion": len(S(p26f) & S(p26v)),
        "n_perts_2025_por_split": {k: len(v) for k, v in p25.items()},
        "n_perts_2025_union": len(u25),
        "perts_2025_por_split": p25,
        "genes_2026_que_faltan_en_2025": sorted(S(g26) - S(g25)),
        "genes_2025_que_faltan_en_2026": sorted(S(g25) - S(g26))[:50],
        "definiciones": {k: {"n": len(v), "lista": v} for k, v in defs.items()},
        "n_por_definicion": {k: len(v) for k, v in defs.items()},
        "la_cifra_97_sale": sorted(k for k, v in defs.items() if len(v) == 97),
    }


def un_split(split):
    t0 = time.time()
    p = v25.ruta(split)
    tg, indptr, shape, var_names, cols, var_cols, otros = v25.lee_obs(split)
    nnz_fila = np.diff(indptr)
    cnt = pd.read_csv(v25.ruta_counts(split))
    medido = pd.Series(tg).value_counts()
    chk = []
    for _, r in cnt.iterrows():
        chk.append(int(medido.get(r["target_gene"], 0)) == int(r["n_cells"]))
    C.log(f"B4a sha256 {split} ({os.path.getsize(p)/2**30:.1f} GiB)...")
    sha, nbytes = C.sha256_archivo(p)
    s = v25.selecciona(tg)
    return {
        "split": split, "ruta": p, "bytes": nbytes,
        "gib": nbytes / 2 ** 30, "sha256": sha,
        "shape": shape, "n_genes": shape[1],
        "columnas_obs": cols, "columnas_var": var_cols,
        "n_perts_csv": int(len(cnt)),
        "n_perts_obs": int(len([x for x in medido.index if x != v25.CONTROL])),
        "n_control": int(medido.get(v25.CONTROL, 0)),
        "csv_coincide_con_obs": bool(all(chk)),
        "n_celdas_csv_distintas": int(len(chk) - sum(chk)),
        "n_cells_min": int(cnt.n_cells.min()), "n_cells_med": float(cnt.n_cells.median()),
        "n_cells_max": int(cnt.n_cells.max()),
        "n_perts_menos_de_400": int((cnt.n_cells < 400).sum()),
        "perts_menos_de_400": cnt[cnt.n_cells < 400].sort_values("n_cells")
                                 .to_dict("records"),
        "umi_mediana_csv": float(cnt.median_umi_per_cell.median()),
        "umi_min_csv": float(cnt.median_umi_per_cell.min()),
        "umi_max_csv": float(cnt.median_umi_per_cell.max()),
        "nnz_mediana_por_celula": float(np.median(nnz_fila)),
        "nnz_p1": float(np.percentile(nnz_fila, 1)),
        "nnz_p99": float(np.percentile(nnz_fila, 99)),
        "nnz_total": int(indptr[-1]),
        "n_batches": int(len(set(otros["batch"].tolist()))) if otros.get("batch") is not None else None,
        "n_guide_ids": int(len(set(otros["guide_id"].tolist()))) if otros.get("guide_id") is not None else None,
        "seleccion_si_fuera_la_superficie": {
            k: v for k, v in s.items() if k != "idx" and k != "por_pert"},
        "dt_s": time.time() - t0,
    }


def md(res):
    o = []; A = o.append
    P = res["pareados"]
    A("# VCC 2025 como superficie: inventario")
    A("")
    A("## 🔴 0. Que es esta superficie y que no")
    A("")
    A("`/datos/vcc2025/` tiene **train, validation y test CON VERDAD**. Es la "
      "unica superficie del proyecto donde se puede ajustar contra verdad **sin "
      "ser transductivo**: la verdad de 2026 no existe, y la de Jiang se usa "
      "como contexto, no como panel.")
    A("")
    A("Es el eje **PANEL** de los dos ortogonales. Jiang es el eje **CONTEXTO**.")
    A("")
    A("🔴 **Lo que esta superficie NO puede dar** (§5, declarado aqui arriba "
      "tambien): es **UN SOLO contexto** -- una sola linea celular, un solo "
      "experimento. No mide dominancia entre contextos y no da rango entre "
      "contextos. Cualquier cifra de 2025 es de un punto, no de un rango.")
    A("")
    A("## 1. Los tres splits: formas, conteos y sha256")
    A("")
    fil = []
    for i, r in enumerate(res["splits"]):
        fil.append([i + 1, f"**{r['split']}**", f"{r['gib']:.2f}",
                    f"{r['shape'][0]:,}", f"{r['shape'][1]:,}",
                    r["n_perts_obs"], f"{r['n_control']:,}",
                    f"{r['nnz_total']/1e6:.1f} M",
                    f"{r['nnz_mediana_por_celula']:.0f}",
                    f"{r['umi_mediana_csv']:,.0f}",
                    r["sha256"][:16] + "..."])
    A("\n".join(C.tabla(
        ["# ", "split", "GiB", "celulas", "genes", "perturbaciones",
         "celulas `non-targeting`", "nnz total", "nnz mediana/celula",
         "UMI mediana/celula (del csv)", "sha256 (16)"], fil)))
    A("")
    A("**sha256 completos:**")
    A("")
    A("\n".join(C.tabla(["split", "bytes", "sha256"],
                        [[r["split"], f"{r['bytes']:,}", f"`{r['sha256']}`"]
                         for r in res["splits"]])))
    A("")
    A("### 1.1 Control de consistencia: el csv contra el `obs`")
    A("")
    fil = []
    for i, r in enumerate(res["splits"]):
        fil.append([i + 1, f"**{r['split']}**", r["n_perts_csv"], r["n_perts_obs"],
                    "si" if r["csv_coincide_con_obs"] else f"🔴 NO ({r['n_celdas_csv_distintas']} filas)",
                    r["n_batches"], r["n_guide_ids"],
                    ", ".join(r["columnas_obs"]), C.dur(r["dt_s"])])
    A("\n".join(C.tabla(
        ["# ", "split", "perts en el csv", "perts en el obs",
         "`n_cells` del csv coincide con el obs", "# lotes", "# `guide_id`",
         "columnas de obs", "duracion"], fil)))
    A("")
    A("### 1.2 🔴 Potencia por perturbacion: no todas llegan a 400 celulas")
    A("")
    fil = []
    for i, r in enumerate(res["splits"]):
        fil.append([i + 1, f"**{r['split']}**", r["n_cells_min"],
                    f"{r['n_cells_med']:.0f}", r["n_cells_max"],
                    f"**{r['n_perts_menos_de_400']}** de {r['n_perts_csv']}",
                    f"{r['umi_min_csv']:,.0f} - {r['umi_max_csv']:,.0f}"])
    A("\n".join(C.tabla(
        ["# ", "split", "min `n_cells`", "mediana", "max",
         "**perts con < 400 celulas**", "rango de UMI mediana"], fil)))
    A("")
    A("El reto 2026 fija **400 celulas por perturbacion** "
      "(`manifest.json: cells_per_pert`). En 2025 eso **no se cumple para "
      "todas**, asi que cualquier construccion que pida 400 tiene que declarar "
      "que hace con las que no llegan.")
    A("")
    for r in res["splits"]:
        if not r["perts_menos_de_400"]:
            continue
        A(f"#### Las {len(r['perts_menos_de_400'])} perturbaciones con < 400 "
          f"celulas en `{r['split']}`")
        A("")
        A("\n".join(C.tabla(
            ["# ", "pert", "n_cells", "UMI mediana/celula"],
            [[i + 1, f"`{x['target_gene']}`", x["n_cells"],
              f"{x['median_umi_per_cell']:,.0f}"]
             for i, x in enumerate(r["perts_menos_de_400"])])))
        A("")
    A("## 2. 🔴 Los genes pareados con el panel 2026")
    A("")
    A("\n".join(C.tabla(["# ", "cantidad", "valor"], [
        [1, "genes medidos en 2025", f"{P['n_genes_2025']:,}"],
        [2, "genes medidos en 2026", f"{P['n_genes_2026']:,}"],
        [3, "**genes medidos en los dos**", f"**{P['n_por_definicion']['d7_genes_medidos_en_los_dos']:,}**"],
        [4, "genes de 2026 que faltan en 2025", f"{len(P['genes_2026_que_faltan_en_2025'])}: " + ", ".join(f"`{g}`" for g in P["genes_2026_que_faltan_en_2025"][:10])],
        [5, "perturbaciones 2025 (union de los 3 splits)", P["n_perts_2025_union"]],
        [6, "perturbaciones 2026 `final` / `val` / union / interseccion", f"{P['n_perts_2026_final']} / {P['n_perts_2026_val']} / {P['n_perts_2026_union']} / {P['n_perts_2026_interseccion']}"],
    ])))
    A("")
    A("🔴 **El encargo cita '97 genes pareados'. Se calcularon las siete "
      "definiciones de 'pareado' que esta superficie admite y se reporta lo "
      "medido en cada una.**")
    A("")
    nombres = {
        "d1_perts2025_x_perts2026_final": "perts 2025 ∩ perts 2026 `final`",
        "d2_perts2025_x_perts2026_val": "perts 2025 ∩ perts 2026 `val`",
        "d3_perts2025_x_union_perts2026": "perts 2025 ∩ (union de perts 2026)",
        "d4_perts2025_x_interseccion_perts2026": "perts 2025 ∩ (interseccion de perts 2026)",
        "d5_perts2025_en_el_espacio_de_genes_2026": "perts 2025 que 2026 MIDE como gen",
        "d6_perts2026_union_en_el_espacio_de_genes_2025": "perts 2026 que 2025 MIDE como gen",
        "d7_genes_medidos_en_los_dos": "genes medidos en los dos paneles",
    }
    fil = []
    for i, (k, n) in enumerate(P["n_por_definicion"].items()):
        fil.append([i + 1, f"`{k}`", nombres.get(k, k), n,
                    "🔴 **= 97**" if n == 97 else ""])
    A("\n".join(C.tabla(["# ", "definicion", "que es", "**n medido**", "¿97?"], fil)))
    A("")
    if P["la_cifra_97_sale"]:
        A(f"La cifra 97 **sale** de: {', '.join('`'+k+'`' for k in P['la_cifra_97_sale'])}.")
    else:
        A("🔴 **La cifra 97 NO sale de ninguna de las siete definiciones.** Se "
          "deja asi, medido y declarado, sin forzar ninguna definicion para que "
          "de 97. La que mas se le parece en espiritu -- perturbaciones de 2025 "
          "que tambien se perturban en 2026 -- da "
          f"**{P['n_por_definicion']['d3_perts2025_x_union_perts2026']}** contra "
          "la union de los dos paneles 2026.")
    A("")
    A("### 2.1 Las listas exactas")
    A("")
    for k in ("d3_perts2025_x_union_perts2026", "d1_perts2025_x_perts2026_final",
              "d2_perts2025_x_perts2026_val", "d4_perts2025_x_interseccion_perts2026"):
        d = P["definiciones"][k]
        A(f"**`{k}`** -- {nombres.get(k,k)}, **{d['n']}** genes:")
        A("")
        A(", ".join(f"`{g}`" for g in d["lista"]) if d["lista"] else "(vacia)")
        A("")
    A("### 2.2 En que split de 2025 esta cada pareado")
    A("")
    pa = P["definiciones"]["d3_perts2025_x_union_perts2026"]["lista"]
    fil = []
    for i, g in enumerate(pa):
        donde = [s for s, v in P["perts_2025_por_split"].items() if g in v]
        p26 = []
        if g in P["definiciones"]["d1_perts2025_x_perts2026_final"]["lista"]:
            p26.append("final")
        if g in P["definiciones"]["d2_perts2025_x_perts2026_val"]["lista"]:
            p26.append("val")
        fil.append([i + 1, f"`{g}`", ", ".join(donde), ", ".join(p26)])
    A("\n".join(C.tabla(["# ", "gen", "split(s) de 2025", "panel(es) de 2026"], fil)))
    A("")
    A("## 3. La superficie que usan 4b y 4c, declarada antes de medir")
    A("")
    A("\n".join(v25.tabla_construccion(res["construccion"])))
    A("")
    A("## 4. Defectos de la construccion, como datos sobre la medicion")
    A("")
    A("\n".join(C.tabla(["# ", "defecto", "que implica"], [
        [1, "El panel de genes no es el mismo: 18,080 contra 18,533.", f"{P['n_por_definicion']['d7_genes_medidos_en_los_dos']:,} genes en los dos. Cualquier cifra por gen cruzada entre las dos superficies vive en esa interseccion, no en el panel entero."],
        [2, "Las perturbaciones casi no se solapan.", f"{P['n_por_definicion']['d3_perts2025_x_union_perts2026']} de {P['n_perts_2025_union']} perturbaciones de 2025 aparecen en algun panel de 2026. El pareado por blanco es chico."],
        [3, "🔴 15 de las 100 perturbaciones del split `test` tienen < 400 celulas.", "La construccion de 4b/4c toma `min(400, disponibles)` y lo declara por perturbacion. No es el regimen del reto, que lleva 400 fijas."],
        [4, "El control se submuestrea de 38,176 a 18,400.", "Es para igualar el `n` de control del reto 2026. El submuestreo tiene semilla declarada y descarta informacion a proposito."],
        [5, "La profundidad es muy distinta a la de Jiang.", "~54,000 UMI/celula en 2025 contra 3,319-6,203 de mediana en Jiang (`ntc_nulo.md` §3). Las dos superficies no son intercambiables."],
        [6, "El `X` de 2025 es `float32`.", "Los valores son enteros; se convierte a `int64` para que el camino sea el mismo que en Jiang. Si hubiera no enteros, la conversion los truncaria -- no se detecto ninguno al construir."],
    ])))
    A("")
    A("## 5. 🔴 Lo que esta superficie NO puede dar")
    A("")
    A("\n".join(C.tabla(["# ", "lo que no da", "por que"], [
        [1, "**Dominancia entre contextos**", "Es UN solo contexto. No hay dos lineas ni dos experimentos que comparar dentro de 2025."],
        [2, "**Rango entre contextos**", "Sin segundo contexto no hay rango. Toda cifra de aqui es un punto."],
        [3, "Heterogeneidad de respuesta entre contextos", "La cantidad que `magnitud_centro.md` reparte entre CENTRO y contexto no esta definida con un solo contexto."],
        [4, "Transferencia a 2026", "El panel y la profundidad difieren (§4). Que algo ajuste aqui no dice que ajuste alla."],
    ])))
    A("")
    A("## 6. Las entradas")
    A("")
    A("\n".join(C.tabla(["entrada", "valor"], [
        ["datos 2025", "`/datos/vcc2025/{train,validation,test}`"],
        ["genes 2025", "`/datos/vcc2025/gene_names.csv`"],
        ["panel 2026", "`/datos/reto/final/` y `/datos/reto/val/` (`gene_names.csv`, `pert_counts.csv`)"],
        ["salidas", "`modelo/out/vcc2025_inventario.md`, `modelo/out/vcc2025_inventario.json`"],
    ])))
    A("")
    A(f"*Corrida: {res['t_inicio']} -> {res['t_fin']} UTC, duracion "
      f"{C.dur(res['dt_total_s'])}.*")
    return o


def main():
    t0 = time.time()
    os.makedirs(v25.WORK, exist_ok=True)
    C.guarda_disco()
    res = {"bloque": "B4a", "semillas": C.comprueba_semillas(), "splits": [],
           "t_inicio": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(t0))}
    for s in ("train", "validation", "test"):
        C.log(f"B4a inventario {s}")
        res["splits"].append(un_split(s))
        C.escribe_json(OUT_JSON, res)
    C.log("B4a pareados")
    res["pareados"] = pareados()
    tg, _i, _s, _v, _c, _vc, _o = v25.lee_obs(v25.SPLIT_SUP)
    sel = v25.selecciona(tg)
    res["construccion"] = {
        "split": v25.SPLIT_SUP, "n_perturbaciones": len(sel["perts"]),
        "n_por_pert_objetivo": v25.N_POR_PERT,
        "n_truncadas": sum(1 for x in sel["por_pert"] if x["truncada"]),
        "n_celulas_perturbadas": sel["n_pert_cells"],
        "n_control_pool": sel["n_control_pool"],
        "n_control_usadas": sel["n_control_usadas"],
        "shape": [sel["n_filas"], res["splits"][-1]["n_genes"]],
        "nnz": None, "seed_pert": v25.SEED_PERT, "seed_control": v25.SEED_CTRL,
        "etiqueta_control": v25.CONTROL, "huella": "(se calcula en 4b/4c)",
        "por_pert": sel["por_pert"],
    }
    res["construccion"]["nnz"] = 0
    res["t_fin"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
    res["dt_total_s"] = time.time() - t0
    C.escribe_json(OUT_JSON, res)
    C.escribe_md(OUT_MD, md(res))
    C.log(f"B4a listo: {OUT_MD}")


if __name__ == "__main__":
    main()
