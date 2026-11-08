#!/usr/bin/env python
"""B3 - 3c PARTE B: situar las seis emisiones ya medidas.

LAS SEIS, con las semillas que YA gastaron (j = 0 en todas), reconstruidas con
sus propias fabricas sin modificarlas:
  1. SUELO                 `suelo_floor.predict_floor`      20264000 + 100*fold
  2. gwps REAL             `gwps_emision.hace_predict_fn`   20271000 + 100*fold
  3. gwps PERMUTADO        `gwps_permutado.hace_predict_fn_perm` 20272000 + 100*fold
                           (desarreglo 20272500 + 0)
  4. ORACULO COMPLETO (1b) `pds_objetivo.predict_razon`     20275000 + 100*fold
  5. ORACULO DE COLA       idem, razones de cola            20276000 + 100*fold
  6. ORACULO DE BULTO      idem, razones de bulto           20277000 + 100*fold

🔴 EL ORACULO COMPLETO ES, POR DEFINICION, LA POSICION CORRECTA: lee las medias
del lado real perturbado de la misma poblacion que se puntua.

🔴 NO SE TRAE NINGUN MODELO ANTERIOR. La medicion se sostiene sin ellos.

LO QUE SE MIDE, por emision, linea y perturbacion:
  * `d_norm` = ||centroide emitido - centroide real|| / dispersion del CONJUNTO
    REAL. En unidades de la dispersion del propio conjunto real.
  * `frac_disp` = dispersion de la emision / dispersion del conjunto real. Que
    fraccion de la dispersion real cubre la emision.
"""
import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import comun as C

sys.path.insert(0, C.SRC)

OUT_JSON = os.path.join(C.OUT, "posicion_emisiones.json")
OUT_MD = os.path.join(C.OUT, "posicion_emisiones.md")
ED_CACHE = "/datos/modelo/out/pds_obj/ed_real.json"
TARGET_SUM = 50000.0                 # el mismo `bulk_target_sum` de magnitud_centro §1.2

EMISIONES = [
    ("suelo",     "SUELO (control remuestreado)",       20264000),
    ("gwps_real", "gwps REAL (estructura externa)",     20271000),
    ("gwps_perm", "gwps PERMUTADO (el desarreglo)",     20272000),
    ("completo",  "ORACULO COMPLETO por razon (1b)",    20275000),
    ("cola",      "ORACULO DE COLA",                    20276000),
    ("bulto",     "ORACULO DE BULTO",                   20277000),
]


def normaliza(X):
    """log1p(TS * x / suma_de_la_celula), por celula. El analogo celular del
    `bulk_lognorm` con `bulk_target_sum = 50000` de `magnitud_centro.md` §1.2.
    log1p(0) = 0, asi que la matriz sigue dispersa."""
    X = X.tocsr().astype(np.float64).copy()
    s = np.asarray(X.sum(axis=1)).ravel()
    inv = TARGET_SUM / np.maximum(s, 1.0)
    X = sp.diags(inv) @ X
    X = X.tocsr()
    X.data = np.log1p(X.data)
    return X, s


def estad(X, keep):
    """(centroide enmascarado, dispersion RMS) de un bloque ya normalizado.
    `dispersion^2 = E||x||^2 - ||c||^2` con la mascara de genes aplicada."""
    n = X.shape[0]
    c = np.asarray(X.sum(axis=0)).ravel() / n
    c = c * keep
    msq = float((X.data * X.data * keep[X.indices]).sum()) / n
    d2 = max(msq - float((c * c).sum()), 0.0)
    return c, float(np.sqrt(d2))


def mide_bloques(bloques, panel, keep):
    out = {}
    for h in panel:
        c, d = estad(bloques[h], keep)
        out[h] = (c, d)
    return out


def una_linea(line, fold_idx, panel, var_names, cell_type, gene, n_genes, cache):
    import pds_objetivo as PO
    from gwps_emision import hace_predict_fn
    from gwps_permutado import hace_predict_fn_perm

    t0 = time.time()
    real_ad, pred_suelo, control_pool, meta = C.arma_linea(
        line, fold_idx, panel, var_names, cell_type, gene, n_genes)
    hue = C.huella_adata(real_ad)
    var = np.asarray(var_names, dtype=str)
    keep_all = np.ones(len(var), dtype=np.float64)
    en_panel = np.isin(var, np.asarray(panel, dtype=str))
    keep_sin = (~en_panel).astype(np.float64)

    # ---- el conjunto REAL
    lab = np.asarray(real_ad.obs["target_gene"].tolist(), dtype=object)
    Xr = real_ad.X.tocsr()
    real_n, prof = {}, {}
    for h in panel:
        Xn, s = normaliza(Xr[lab == h])
        real_n[h] = Xn
        prof[h] = float(np.median(s))
    ctrl_n, s_ctrl = normaliza(control_pool)

    reales = {"todos": mide_bloques(real_n, panel, keep_all),
              "sin_panel": mide_bloques(real_n, panel, keep_sin)}
    c_ctrl = {k: estad(ctrl_n, kk) for k, kk in
              (("todos", keep_all), ("sin_panel", keep_sin))}
    del real_n, ctrl_n

    # ---- las razones de los tres oraculos, ANTES de soltar el lado real
    # (una sola lectura del bloque de la linea para todo el bloque B3)
    ed = cache.get(line, {})
    razones = {}
    for modo in ("completo", "cola", "bulto"):
        razones[modo], _d = PO.razones_oraculo(real_ad, panel, ed, modo)
    del real_ad

    filas = []
    for clave, etiqueta, base in EMISIONES:
        t1 = time.time()
        seed = base + 100 * fold_idx
        rng = np.random.default_rng(seed)
        if clave == "suelo":
            lab_p = np.asarray(pred_suelo.obs["target_gene"].tolist(), dtype=object)
            Xp = pred_suelo.X.tocsr()
            bloques = {h: Xp[lab_p == h] for h in panel}
            extra = {"nota": "la emision de suelo j=0 que devolvio el arnes con "
                             "su propia semilla"}
        else:
            if clave == "gwps_real":
                reg = {}
                fn = hace_predict_fn(var_names, reg)
            elif clave == "gwps_perm":
                reg = {}
                fn = hace_predict_fn_perm(var_names, reg, 0)
            else:
                reg = {}
                fn = PO.predict_razon(razones[clave])
            pred = fn(control_pool, panel + [C.CONTROL_LABEL], rng)
            bloques = {h: pred[h].tocsr() for h in panel}
            extra = {"desarreglo": reg.get("desarreglo")} if clave == "gwps_perm" else {}
            del pred
        norm = {h: normaliza(bloques[h])[0] for h in panel}
        del bloques
        for esp, keep in (("todos", keep_all), ("sin_panel", keep_sin)):
            med = mide_bloques(norm, panel, keep)
            for h in panel:
                cR, dR = reales[esp][h]
                cP, dP = med[h]
                dd = float(np.linalg.norm((cP - cR) * keep))
                filas.append({
                    "emision": clave, "espacio": esp, "linea": line, "pert": h,
                    "seed_emision": int(seed),
                    "d_centroide": dd,
                    "disp_real": dR, "disp_emision": dP,
                    "d_norm": (dd / dR if dR > 0 else None),
                    "frac_disp": (dP / dR if dR > 0 else None),
                    "disp_control": c_ctrl[esp][1],
                    "profundidad_mediana_real": prof[h],
                })
        del norm
        C.log(f"B3 {line} {clave}: {time.time()-t1:.0f}s")
    del control_pool, pred_suelo, razones
    return {"linea": line, "fold_idx": fold_idx, "huella_real": hue,
            "n_control_pool": int(meta["n_control_cells"]),
            "n_genes": int(len(var)),
            "n_genes_sin_panel": int(len(var) - int(en_panel.sum())),
            "disp_control": {k: c_ctrl[k][1] for k in c_ctrl},
            "filas": filas, "dt_s": time.time() - t0}


def resumen(filas, clave, esp, campo):
    v = [f[campo] for f in filas
         if f["emision"] == clave and f["espacio"] == esp and f[campo] is not None]
    return v


def md(res):
    o = []; A = o.append
    todas = [f for r in res["lineas"] for f in r["filas"]]
    A("# La posicion de las seis emisiones ya medidas: centroide y dispersion")
    A("")
    A("## 🔴 0. Que mide esto y que no")
    A("")
    A("Las seis emisiones ya tienen sus `u` y sus `s` en disco. Esto no vuelve a "
      "puntuarlas: **las situa**. Por emision, linea y perturbacion se dan dos "
      "cantidades geometricas en el espacio de la celula:")
    A("")
    A("* **`d_norm`** = `||centroide emitido - centroide real|| / dispersion del "
      "conjunto real`. Cuantas dispersiones reales de distancia hay entre los "
      "dos centros.")
    A("* **`frac_disp`** = `dispersion de la emision / dispersion del conjunto "
      "real`. Que fraccion de la dispersion real cubre la emision.")
    A("")
    A("🔴 **El ORACULO COMPLETO es, por definicion, la posicion correcta**: lee "
      "las medias del lado real perturbado de la misma poblacion. Es "
      "transductivo y no puede enviarse; esta aqui como la referencia de lo que "
      "'estar en el lugar' vale numericamente.")
    A("")
    A("🔴 **No se trae ningun modelo anterior.** Las seis emisiones de esta "
      "tabla son las que ya estaban medidas.")
    A("")
    A("## 1. La construccion, declarada antes de medir")
    A("")
    A("\n".join(C.tabla(["# ", "cantidad", "valor"], [
        [1, "espacio", f"`log1p({TARGET_SUM:g} * x / suma de la celula)` por celula -- el analogo celular del `bulk_lognorm` con `bulk_target_sum = {TARGET_SUM:g}` de `magnitud_centro.md` §1.2"],
        [2, "centroide", "media sobre las 400 celulas del bloque, en ese espacio"],
        [3, "dispersion", "`sqrt(E||x - c||^2)` sobre las celulas del bloque: el radio RMS"],
        [4, "dos espacios de genes", f"**todos** ({res['lineas'][0]['n_genes']}) y **sin los 33 blancos del panel** ({res['lineas'][0]['n_genes_sin_panel']}), que es el `exclusion_scope=\"panel\"` de `magnitud_centro.md` §1.2"],
        [5, "poblacion real", f"`seed_real = {C.SEED_BASE_REAL} + fold`, la MISMA de siempre, verificada por sha256"],
        [6, "realizaciones", "**una por emision (`j = 0`)**, con la semilla que ya gasto"],
        [7, "celulas", f"{C.N_PER_LABEL} reales y {C.N_PER_LABEL} emitidas por perturbacion"],
    ])))
    A("")
    A("🔴 **Si el espacio cambia, los numeros cambian.** Esta declarado por eso.")
    A("")
    A("### 1.1 Las seis emisiones y sus semillas ya gastadas")
    A("")
    fil = []
    for i, (clave, etq, base) in enumerate(EMISIONES):
        fil.append([i + 1, f"`{clave}`", etq, f"`{base} + 100*fold`",
                    "si" if clave in ("completo", "cola", "bulto") else "no",
                    "**la posicion correcta por definicion**" if clave == "completo" else ""])
    A("\n".join(C.tabla(["# ", "clave", "emision", "semilla (j = 0)",
                         "transductiva", "nota"], fil)))
    A("")
    A("### 1.2 Verificacion por linea")
    A("")
    fil = []
    for i, r in enumerate(res["lineas"]):
        fil.append([i + 1, f"**{r['linea']}**", r["n_control_pool"],
                    r["huella_real"][:16] + "...",
                    C.f4a(r["disp_control"]["todos"]),
                    C.f4a(r["disp_control"]["sin_panel"]),
                    C.dur(r["dt_s"])])
    A("\n".join(C.tabla(["# ", "linea", "pool de control", "huella real (16)",
                         "dispersion del control (todos)",
                         "dispersion del control (sin panel)", "duracion"], fil)))
    A("")

    # ---- 2
    for esp, nombre in (("sin_panel", "SIN los 33 blancos del panel (el espacio de `magnitud_centro`)"),
                        ("todos", "TODOS los genes")):
        A(f"## 2.{1 if esp == 'sin_panel' else 2} El agregado, espacio: {nombre}")
        A("")
        fil = []
        for clave, etq, base in EMISIONES:
            dn = resumen(todas, clave, esp, "d_norm")
            fd = resumen(todas, clave, esp, "frac_disp")
            fil.append([f"`{clave}`", len(dn), C.rango(dn),
                        C.f4a(float(np.median(dn))) if dn else "n/d",
                        C.rango(fd),
                        C.f4a(float(np.median(fd))) if fd else "n/d"])
        A("\n".join(C.tabla(
            ["emision", "celdas (linea x pert)",
             "**rango de `d_norm`**", "mediana `d_norm`",
             "**rango de `frac_disp`**", "mediana `frac_disp`"], fil)))
        A("")
    A("🔴 **Rango, no media**, por la misma razon que en los demas reportes: las "
      "lineas no son intercambiables.")
    A("")

    # ---- 3 por linea
    A("## 3. Por emision y linea")
    A("")
    for esp in ("sin_panel",):
        fil = []
        for clave, etq, base in EMISIONES:
            for r in res["lineas"]:
                dn = [f["d_norm"] for f in r["filas"]
                      if f["emision"] == clave and f["espacio"] == esp
                      and f["d_norm"] is not None]
                fd = [f["frac_disp"] for f in r["filas"]
                      if f["emision"] == clave and f["espacio"] == esp
                      and f["frac_disp"] is not None]
                fil.append([f"`{clave}`", f"**{r['linea']}**", len(dn),
                            C.f4a(float(np.median(dn))) if dn else "n/d",
                            C.rango(dn),
                            C.f4a(float(np.median(fd))) if fd else "n/d",
                            C.rango(fd)])
        A("\n".join(C.tabla(
            ["emision", "linea", "n perts", "mediana `d_norm`",
             "rango `d_norm`", "mediana `frac_disp`", "rango `frac_disp`"], fil)))
        A("")

    # ---- 4 intermedios
    A("## 4. 🔴 Los intermedios crudos: las celdas una por una")
    A("")
    A("Espacio **sin los 33 blancos del panel**. Una fila por (emision, linea, "
      f"perturbacion): {len(EMISIONES)} x {len(res['lineas'])} x "
      f"{len(res['lineas'][0]['filas'])//(2*len(EMISIONES))} celdas.")
    A("")
    for r in res["lineas"]:
        A(f"### 4.{res['lineas'].index(r)+1} {r['linea']}")
        A("")
        panel = sorted({f["pert"] for f in r["filas"]})
        fil = []
        for i, h in enumerate(panel):
            row = [i + 1, f"`{h}`"]
            for clave, _e, _b in EMISIONES:
                f = next((x for x in r["filas"] if x["emision"] == clave
                          and x["espacio"] == "sin_panel" and x["pert"] == h), None)
                row.append(C.f4a(f["d_norm"]) if f else "n/d")
            for clave, _e, _b in EMISIONES:
                f = next((x for x in r["filas"] if x["emision"] == clave
                          and x["espacio"] == "sin_panel" and x["pert"] == h), None)
                row.append(C.f4a(f["frac_disp"]) if f else "n/d")
            f0 = next((x for x in r["filas"] if x["espacio"] == "sin_panel"
                       and x["pert"] == h), None)
            row.append(C.f4a(f0["disp_real"]) if f0 else "n/d")
            fil.append(row)
        A("\n".join(C.tabla(
            ["# ", "pert"] + [f"d_norm {k}" for k, _, _ in EMISIONES] +
            [f"frac_disp {k}" for k, _, _ in EMISIONES] + ["disp. real"], fil)))
        A("")
    A("## 5. Defectos de la construccion, como datos sobre la medicion")
    A("")
    A("\n".join(C.tabla(["# ", "defecto", "que implica"], [
        [1, "El espacio es una eleccion.", f"`log1p` con `TS = {TARGET_SUM:g}` por celula. Otro espacio da otros numeros; se declara por eso, no por tramite."],
        [2, "`d_norm` divide por la dispersion del conjunto real.", "Donde el conjunto real es muy disperso, una misma distancia absoluta sale mas chica. La columna `disp. real` esta en las tablas de intermedios para poder deshacerlo."],
        [3, "`frac_disp` compara radios RMS, no formas.", "Dos nubes con el mismo radio y orientaciones distintas dan `frac_disp = 1`. No es una medida de solape."],
        [4, "Una sola realizacion por emision (`j = 0`).", "No hay barra de error sobre la semilla de emision. El suelo y gwps tienen 8 semillas en disco; aqui solo se uso la primera."],
        [5, "El centroide del bloque emitido se compara con el del bloque real de la MISMA perturbacion.", "No se mide si la emision confunde perturbaciones entre si. Eso es otra cantidad."],
        [6, "El ORACULO COMPLETO es transductivo.", "Su posicion no es alcanzable sin leer el lado real. Esta como referencia de escala, no como objetivo alcanzable."],
    ])))
    A("")
    A("## 6. Las entradas")
    A("")
    A("\n".join(C.tabla(["entrada", "valor"], [
        ["IFNB", "`/datos/ref/convertido/Seurat_object_IFNB_Perturb_seq.h5ad`"],
        ["gwps", "`/datos/ref/replogle/k562_gwps_raw_bulk.h5ad` via `gwps_emision.construye_ratios`, SIN MODIFICAR"],
        ["emisiones", "`suelo_floor`, `gwps_emision`, `gwps_permutado`, `pds_objetivo`, todas SIN MODIFICAR"],
        ["salidas", "`modelo/out/posicion_emisiones.md`, `modelo/out/posicion_emisiones.json`"],
    ])))
    A("")
    A(f"*Corrida: {res['t_inicio']} -> {res['t_fin']} UTC, duracion "
      f"{C.dur(res['dt_total_s'])}.*")
    return o


def main():
    t_ini = time.time()
    C.guarda_disco()
    from suelo_arnes import load_global_obs, load_var_names, IFNB_PATH
    panel = C.carga_panel()
    cell_type, gene = load_global_obs(IFNB_PATH)
    var_names = load_var_names(IFNB_PATH)
    n_genes = len(var_names)
    cache = json.load(open(ED_CACHE)) if os.path.exists(ED_CACHE) else {}
    res = {"bloque": "B3", "semillas": C.comprueba_semillas(),
           "target_sum": TARGET_SUM,
           "emisiones": [{"clave": k, "etiqueta": e, "base_semilla": b}
                         for k, e, b in EMISIONES],
           "lineas": [],
           "t_inicio": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(t_ini))}
    for i, line in enumerate(C.LINES):
        C.log(f"B3 {line} ({i+1}/6)")
        C.guarda_disco()
        r = una_linea(line, i, panel, var_names, cell_type, gene, n_genes, cache)
        res["lineas"].append(r)
        res["t_fin"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
        res["dt_total_s"] = time.time() - t_ini
        C.escribe_json(OUT_JSON, res)
        C.log(f"B3 {line}: listo dt={r['dt_s']:.0f}s")
    C.escribe_md(OUT_MD, md(res))
    C.log(f"B3 listo: {OUT_MD}")


if __name__ == "__main__":
    main()
