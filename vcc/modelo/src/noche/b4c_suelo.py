#!/usr/bin/env python
"""B4c - PASO 4: las seis metricas sobre un SUELO de la superficie VCC 2025, con
8 semillas, igual que `dispersion_u`.

EL SUELO: `suelo_floor.predict_floor` SIN MODIFICAR -- 400 celulas remuestreadas
CON REEMPLAZO del pool de control, por etiqueta, sin tocar ningun conteo. Es la
misma emision que fijo el suelo de Jiang.

Semillas: `20280300 + j`, `j = 0..7` (familia NUEVA, comprobada en codigo).

Se puntua con el camino nativo: `cell-eval2 run --write-degenes` para el `u`, y
`cell-eval2 score --real-bundle` contra el bundle de 4b para el `s`
(`from_replicate`). Si 4b no dejo bundle, el `u` se reporta igual y el `s` queda
declarado como no calculado.
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import comun as C
import v25

sys.path.insert(0, C.SRC)

OUT_JSON = os.path.join(C.OUT, "vcc2025_suelo.json")
OUT_MD = os.path.join(C.OUT, "vcc2025_suelo.md")
BUNDLE = os.path.join(v25.WORK, "BUNDLE_2025_test")
N_SEMILLAS = 8
SEEDS = C.SEMILLAS_NUEVAS["b4_suelo"]


def lee_control_del_real(real_path):
    """El pool de control del h5ad real ya construido, sin cargar la matriz
    entera."""
    import h5py
    from suelo_arnes import read_rows_csr, decode_categorical
    with h5py.File(real_path, "r") as f:
        n_genes = int(f["X"].attrs["shape"][1])
        node = f["obs"]["target_gene"]
        lab = (np.array(decode_categorical(node), dtype=object)
               if hasattr(node, "keys")
               else np.array([x.decode() if isinstance(x, bytes) else x
                              for x in node[:]], dtype=object))
        vnode = f["var"]["_index"]
        var_names = [x.decode() if isinstance(x, bytes) else x for x in vnode[:]]
    idx = np.where(lab == v25.CONTROL)[0]
    X = read_rows_csr(real_path, np.sort(idx), n_genes)
    perts = sorted(set(lab.tolist()) - {v25.CONTROL})
    return X, perts, var_names


def una_semilla(j, real_path, control_pool, perts, var_names, hay_bundle):
    import anndata as ad
    from suelo_floor import predict_floor
    from suelo_arnes import read_agg_metrics
    from dispersion_u import leer_scored
    from techo_arnes import contar_de_por_perturbacion

    t0 = time.time()
    C.guarda_disco()
    seed = SEEDS[j]
    rng = np.random.default_rng(seed)
    etiquetas = perts + [v25.CONTROL]
    pred = predict_floor(control_pool, etiquetas, rng)
    bloques, labs = [], []
    for l in etiquetas:
        b = pred[l]
        assert b.shape[0] == 400, (l, b.shape)
        bloques.append(b)
        labs.extend([l] * b.shape[0])
    X = sp.vstack(bloques, format="csr")
    del pred, bloques
    adata = ad.AnnData(X=X, obs=pd.DataFrame({"target_gene": labs}),
                       var=pd.DataFrame(index=var_names))
    adata.obs_names_make_unique()
    hue = C.huella_adata(adata)
    pp = os.path.join(v25.WORK, f"pred_suelo_j{j}.h5ad")
    rd = os.path.join(v25.WORK, f"run_suelo_j{j}")
    shape = list(adata.shape)
    nnz = int(X.nnz)
    adata.write_h5ad(pp)
    del adata, X
    C.guarda_disco()

    cmd = [C.CELLEVAL2_PY, "run", "-ap", pp, "-ar", real_path,
           "--preset", "vcc2026", "--profile", "vcc2026",
           "--pert-col", "target_gene", "--write-degenes", "-o", rd]
    C.log(f"B4c j={j} cell-eval2 run...")
    p, dt_run = C.corre(cmd, os.path.join(v25.WORK, f"b4c_run_j{j}.log"))
    rec = {"j": j, "seed_emision": int(seed), "huella_pred": hue,
           "pred_shape": shape, "pred_nnz": nnz,
           "cmd_run": cmd, "run_ok": p.returncode == 0,
           "run_exitcode": p.returncode, "run_duration_s": dt_run,
           "run_error": (p.stderr[-3000:] if p.returncode else None)}
    if p.returncode == 0:
        try:
            rec["u"] = read_agg_metrics(os.path.join(rd, "agg_results.csv"))
        except Exception as e:
            rec["u_error"] = repr(e)
        try:
            from dispersion_u import hash_del_run
            rec["hash_run"] = hash_del_run(rd)
        except Exception as e:
            rec["hash_run_error"] = repr(e)
        try:
            cuentas = contar_de_por_perturbacion(rd, control=v25.CONTROL)
            vac = [q for q, c in cuentas.items() if c["union_vacia"]]
            conp = [q for q, c in cuentas.items() if c["n_DE_pred"] > 0]
            rec["cohortes"] = {
                "n_perts_total": len(cuentas),
                "E_union_vacia": len(vac),
                "n_perts_con_pred_sig": len(conp),
                "suma_n_DE_pred": int(sum(c["n_DE_pred"] for c in cuentas.values())),
                "suma_n_DE_real": int(sum(c["n_DE_real"] for c in cuentas.values())),
                "n_DE_real_vacio": int(sum(1 for c in cuentas.values()
                                           if c["n_DE_real"] == 0)),
            }
            rec["por_pert"] = {q: {k: int(v) if isinstance(v, (int, np.integer))
                                   else v for k, v in c.items()}
                               for q, c in cuentas.items()}
        except Exception as e:
            rec["cohortes_error"] = repr(e)
        if hay_bundle:
            out_csv = os.path.join(v25.WORK, f"scored_j{j}.csv")
            cmd2 = [C.CELLEVAL2_PY, "score",
                    "--user-agg", os.path.join(rd, "agg_results.csv"),
                    "--user-meta", os.path.join(rd, "run_meta.json"),
                    "--real-bundle", BUNDLE, "-o", out_csv]
            p2, dt2 = C.corre(cmd2, os.path.join(v25.WORK, f"b4c_score_j{j}.log"))
            rec["cmd_score"] = cmd2
            rec["score_ok"] = p2.returncode == 0
            rec["score_exitcode"] = p2.returncode
            rec["score_duration_s"] = dt2
            rec["score_stderr_literal"] = (p2.stderr[-2000:] if p2.returncode else None)
            if p2.returncode == 0:
                rec["scored"] = leer_scored(out_csv)
        else:
            rec["score_ok"] = None
            rec["score_nota"] = ("no se intento: 4b no dejo bundle en "
                                 f"{BUNDLE}")
    C.borra(pp, rd)
    rec["total_duration_s"] = time.time() - t0
    return rec


def val_s(rec, metrica):
    sc = rec.get("scored")
    if not sc:
        return None
    f = sc["filas"].get(metrica)
    if not f:
        return None
    for c in ("from_replicate", "from_baseline", "value"):
        if c in f and f[c] is not None:
            return f[c]
    return None


def md(res):
    o = []; A = o.append
    R = res["corridas"]
    ok = [r for r in R if r.get("run_ok")]
    A("# VCC 2025 como superficie: el SUELO con 8 semillas y las seis metricas")
    A("")
    A("## 🔴 0. Que es esto y que no")
    A("")
    A("El **suelo** de la superficie 2025: `suelo_floor.predict_floor` sin "
      "modificar -- 400 celulas remuestreadas **con reemplazo** del pool de "
      "control, por etiqueta, sin tocar ningun conteo. Es la **misma emision** "
      "que fijo el suelo de Jiang, sobre la otra superficie.")
    A("")
    A("Ocho semillas, `20280300 + j`, para que el suelo tenga dispersion "
      "medida y no un solo punto -- igual que `dispersion_u`.")
    A("")
    A("## 1. La construccion, declarada antes de medir")
    A("")
    A("\n".join(v25.tabla_construccion(res["real_meta"])))
    A("")
    A("\n".join(C.tabla(["# ", "cantidad", "valor"], [
        [1, "emision", "`suelo_floor.predict_floor`, SIN MODIFICAR"],
        [2, "celulas emitidas por etiqueta", "400 (con reemplazo del pool de control)"],
        [3, "etiquetas emitidas", f"{res['n_etiquetas']} ({res['n_perts']} perturbaciones + el control)"],
        [4, "semillas", f"`{SEEDS[0]} + j`, `j = 0..{N_SEMILLAS-1}` -- familia NUEVA, {len(SEEDS)} semillas, 0 choques"],
        [5, "puntuacion `u`", "`cell-eval2 run --preset vcc2026 --profile vcc2026 --write-degenes`, fila `mean` de `agg_results.csv`"],
        [6, "puntuacion `s`", ("`cell-eval2 score --real-bundle` contra el bundle de 4b"
                               if res["hay_bundle"] else
                               "🔴 **no calculada**: 4b no dejo bundle")],
        [7, "bundle", f"`{BUNDLE}`" + ("" if res["hay_bundle"] else " (**ausente**)")],
    ])))
    A("")
    A("## 2. Las seis metricas sobre el suelo 2025")
    A("")
    fil = []
    for m in C.VCC2026_METRICS:
        u = [r["u"][m] for r in ok if r.get("u") and m in r["u"]]
        s = [val_s(r, m) for r in ok]
        s = [x for x in s if x is not None]
        fil.append([f"`{C.METRICA_CORTA[m]}`", m, len(u),
                    C.f4(float(np.mean(u))) if u else "n/d",
                    C.f4(float(np.std(u, ddof=1))) if len(u) > 1 else "n/d",
                    C.rango(u),
                    C.f4(float(np.mean(s))) if s else "n/d",
                    C.rango(s) if s else "n/d"])
    A("\n".join(C.tabla(
        ["corta", "metrica", "n semillas", "**`u` media**", "`u` sd",
         "**rango de `u`**", "**`s` media**", "rango de `s`"], fil)))
    A("")
    A("🔴 **Rango, no solo media.** Con 8 semillas el rango es la cantidad que "
      "dice cuanta resolucion tiene la superficie.")
    A("")
    A("## 3. 🔴 Los intermedios crudos: las 8 corridas una por una")
    A("")
    fil = []
    for r in R:
        fil.append([r["j"], r["seed_emision"], r["huella_pred"][:12] + "...",
                    f"{r['pred_shape'][0]:,}x{r['pred_shape'][1]:,}",
                    f"{r['pred_nnz']/1e6:.1f} M",
                    "si" if r.get("run_ok") else f"🔴 {r.get('run_exitcode')}",
                    C.dur(r.get("run_duration_s") or 0),
                    ("si" if r.get("score_ok") else
                     ("n/d" if r.get("score_ok") is None else "🔴 NO")),
                    (r.get("hash_run") or {}).get("coincide"),
                    C.dur(r.get("total_duration_s") or 0)])
    A("\n".join(C.tabla(
        ["j", "semilla", "huella pred (12)", "forma pred", "nnz pred",
         "run ok", "duracion del run", "score ok", "`config_hash` coincide",
         "duracion total"], fil)))
    A("")
    A("### 3.1 El `u` crudo, corrida por corrida")
    A("")
    fil = []
    for r in R:
        row = [r["j"], r["seed_emision"]]
        for m in C.VCC2026_METRICS:
            row.append(C.f4((r.get("u") or {}).get(m)))
        fil.append(row)
    A("\n".join(C.tabla(["j", "semilla"] +
                        [f"`{C.METRICA_CORTA[m]}`" for m in C.VCC2026_METRICS], fil)))
    A("")
    if res["hay_bundle"]:
        A("### 3.2 El `s` crudo, corrida por corrida")
        A("")
        fil = []
        for r in R:
            fil.append([r["j"], r["seed_emision"]] +
                       [C.f4(val_s(r, m)) for m in C.VCC2026_METRICS] +
                       [C.f4(val_s(r, "avg_score"))])
        A("\n".join(C.tabla(["j", "semilla"] +
                            [f"`{C.METRICA_CORTA[m]}`" for m in C.VCC2026_METRICS] +
                            ["**avg_score**"], fil)))
        A("")
    A("### 3.3 Las cohortes moviles, corrida por corrida")
    A("")
    fil = []
    for r in R:
        c = r.get("cohortes") or {}
        fil.append([r["j"], c.get("n_perts_total"), c.get("E_union_vacia"),
                    c.get("n_perts_con_pred_sig"), c.get("suma_n_DE_pred"),
                    c.get("suma_n_DE_real"), c.get("n_DE_real_vacio")])
    A("\n".join(C.tabla(
        ["j", "perts puntuadas", "`E_union_vacia`", "perts con pred significativo",
         "suma `n_DE_pred`", "suma `n_DE_real`",
         "**perts con DE REAL vacio**"], fil)))
    A("")
    A("🔴 La ultima columna es del lado real, asi que **no varia con la semilla "
      "del suelo**: es una propiedad de la superficie 2025, igual que el "
      "`control_permutado.md` §7 lo es de Jiang.")
    A("")
    A("## 4. 🔴 Lo que esta superficie NO puede dar")
    A("")
    A("\n".join(C.tabla(["# ", "lo que no da", "por que"], [
        [1, "**Dominancia entre contextos**", "Es UN solo contexto. Las 8 semillas dan dispersion de la EMISION, no de contexto."],
        [2, "**Rango entre contextos**", "El rango de esta tabla es sobre semillas del suelo, no sobre contextos. En Jiang hay 6 lineas x 8 semillas y los dos rangos se distinguen; aqui solo existe el segundo."],
        [3, "Transferencia del suelo a 2026", "Panel y profundidad difieren; el valor del suelo es de esta superficie."],
    ])))
    A("")
    A("## 5. Defectos de la construccion")
    A("")
    A("\n".join(C.tabla(["# ", "defecto", "que implica"], [
        [1, "El lado real es UNA realizacion (una semilla de celulas por perturbacion).", "Todas las 8 corridas comparten el mismo lado real, a proposito: la dispersion medida es solo de la emision."],
        [2, "15 de 100 perturbaciones llevan < 400 celulas reales.", "Su potencia es menor y arrastra las metricas de DE. La tabla por perturbacion esta en el json."],
        [3, "El suelo emite 400 celulas por etiqueta aunque el real tenga menos.", "Los dos lados no tienen el mismo `n` en esas 15 perturbaciones."],
        [4, "El `s` depende del bundle de 4b.", "Si 4b fallo, el `s` no existe y queda declarado; el `u` no depende del bundle."],
        [5, "El control emitido tambien es control remuestreado.", "Igual que en Jiang: el suelo predice el control con el propio control."],
    ])))
    A("")
    A(f"*Corrida: {res['t_inicio']} -> {res['t_fin']} UTC, duracion "
      f"{C.dur(res['dt_total_s'])}.*")
    return o


def main():
    t0 = time.time()
    os.makedirs(v25.WORK, exist_ok=True)
    C.guarda_disco()
    real_path, real_meta = v25.arma_real()
    control_pool, perts, var_names = lee_control_del_real(real_path)
    hay_bundle = os.path.exists(os.path.join(BUNDLE, "anchor_agg.parquet"))
    C.log(f"B4c control {control_pool.shape} perts={len(perts)} "
          f"bundle={'si' if hay_bundle else 'NO'}")
    res = {"bloque": "B4c", "semillas": C.comprueba_semillas(),
           "real_meta": real_meta, "hay_bundle": hay_bundle,
           "n_perts": len(perts), "n_etiquetas": len(perts) + 1,
           "corridas": [],
           "t_inicio": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(t0))}
    for j in range(N_SEMILLAS):
        C.log(f"B4c semilla j={j} ({j+1}/{N_SEMILLAS})")
        rec = una_semilla(j, real_path, control_pool, perts, var_names, hay_bundle)
        res["corridas"].append(rec)
        res["t_fin"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
        res["dt_total_s"] = time.time() - t0
        C.escribe_json(OUT_JSON, res)
        C.escribe_md(OUT_MD, md(res))      # salida incremental: tras CADA semilla
        C.log(f"B4c j={j} run_ok={rec.get('run_ok')} "
              f"dt={rec.get('total_duration_s',0):.0f}s")
    C.borra(v25.REAL_H5)
    C.log(f"B4c listo: {OUT_MD}")


if __name__ == "__main__":
    main()
