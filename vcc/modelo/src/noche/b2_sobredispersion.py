#!/usr/bin/env python
"""B2 - 3c PARTE A: hay heterogeneidad celula a celula POR ENCIMA del ruido de
conteo?

LA PREGUNTA, declarada antes de medir. Las emisiones por razon aplican UN vector
de razones a las 400 celulas y realizan el conteo con `_reescala_entero`
(Binomial para `r<=1`, Poisson aditivo para `r>1`): ruido de conteo alrededor de
UNA moda. Si las 400 celulas reales fueran mas dispersas que eso, habria un
hueco que una emision de UNA moda no puede llenar.

LA LINEA BASE. El cero correcto NO es Poisson teorico: son las celulas de
CONTROL de la misma linea, que ya son una mezcla (ciclo, profundidad, lote) y
por tanto ya estan sobredispersas. Se toman 400 controles de la misma linea y se
comparan A MEDIA IGUALADA, por banda de expresion.

🔴 EL EXCESO PUEDE SALIR CERO O NEGATIVO. Si sale, se dice.

LA CUARTA COLUMNA. La misma cantidad sobre UNA EMISION NUESTRA: el ORACULO
COMPLETO POR RAZON de `pds_objetivo.md` §3.1 (1b), j=0, con su semilla ya
gastada `20275000 + 100*fold`. Si nuestra emision tiene la misma
sobredispersion que la real, NO HAY HUECO QUE LLENAR y eso es el resultado.
"""
import json
import os
import sys
import time

import warnings

import numpy as np
import scipy.sparse as sp

warnings.filterwarnings("ignore", message="Mean of empty slice")
warnings.filterwarnings("ignore", message="All-NaN slice encountered")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import comun as C

sys.path.insert(0, C.SRC)

OUT_JSON = os.path.join(C.OUT, "sobredispersion.json")
OUT_MD = os.path.join(C.OUT, "sobredispersion.md")
ED_CACHE = "/datos/modelo/out/pds_obj/ed_real.json"
SEED_ORACULO_1B = 20275000          # ya gastada: 20275000 + 100*fold + j, j=0

# Bandas de expresion (media de cuentas crudas por celula en el CONTROL).
# Declaradas antes de medir, no ajustadas a los datos.
BORDES = [0.0, 1e-3, 1e-2, 0.1, 0.5, 1.0, 2.0, 5.0, 20.0, np.inf]
ETQ_BANDA = (["= 0"] +
             [f"({BORDES[i]:g}, {BORDES[i+1]:g}]" for i in range(len(BORDES) - 1)])
N_BOOT = 1000
N_BINS_CURVA = 60                   # bins finos para la curva de referencia


def media_var(Xb):
    """Media y varianza (ddof=0) por gen de un bloque de celulas, sin
    densificar. `var = E[x^2] - E[x]^2`."""
    n = Xb.shape[0]
    m = np.asarray(Xb.mean(axis=0)).ravel()
    e2 = np.asarray(Xb.multiply(Xb).mean(axis=0)).ravel()
    v = np.maximum(e2 - m * m, 0.0)
    return m.astype(np.float64), v.astype(np.float64), n


def banda_de(m):
    """Indice de banda: 0 para media exactamente 0, luego los bordes."""
    b = np.digitize(m, BORDES[1:], right=True) + 1
    b[m == 0.0] = 0
    return b


def curva_referencia(m_ctrl, v_ctrl):
    """Curva `var_control(media)`: bins finos en log10(media) sobre los genes de
    control con media > 0, con la MEDIA de la varianza en cada bin. Es la
    referencia que se evalua a la media de cada gen perturbado -- es decir, la
    igualacion de media."""
    vivo = m_ctrl > 0
    lm = np.log10(m_ctrl[vivo])
    lo, hi = lm.min(), lm.max()
    bordes = np.linspace(lo, hi, N_BINS_CURVA + 1)
    idx = np.clip(np.digitize(lm, bordes[1:-1]), 0, N_BINS_CURVA - 1)
    centros = 0.5 * (bordes[:-1] + bordes[1:])
    vv = v_ctrl[vivo]
    ref = np.full(N_BINS_CURVA, np.nan)
    nn = np.zeros(N_BINS_CURVA, dtype=int)
    for k in range(N_BINS_CURVA):
        s = idx == k
        nn[k] = int(s.sum())
        if nn[k]:
            ref[k] = float(vv[s].mean())
    ok = np.isfinite(ref)
    v_cero = float(v_ctrl[~vivo].mean()) if (~vivo).any() else 0.0
    return {"centros_log10": centros[ok], "ref": ref[ok], "n_por_bin": nn[ok],
            "var_en_media_cero": v_cero, "lo": float(lo), "hi": float(hi)}


def evalua_referencia(cur, m):
    """Interpola la curva de referencia en la media dada. Fuera de rango, se
    sujeta al extremo (declarado). Media 0 -> la varianza de los genes de
    control con media 0."""
    out = np.empty_like(m)
    cero = m <= 0
    out[cero] = cur["var_en_media_cero"]
    lm = np.log10(np.maximum(m[~cero], 1e-12))
    out[~cero] = np.interp(lm, cur["centros_log10"], cur["ref"])
    return out


def boot_ci(x, seed, n=N_BOOT):
    """IC 95 % por bootstrap sobre GENES (la unidad de resampleo son los genes,
    no los pares pert-gen: cada gen entra o sale completo)."""
    x = np.asarray(x, dtype=np.float64)
    x = x[np.isfinite(x)]
    if x.size < 3:
        return None, None, None
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, x.size, size=(n, x.size))
    med = x[idx].mean(axis=1)
    return (float(x.mean()), float(np.percentile(med, 2.5)),
            float(np.percentile(med, 97.5)))


def por_banda(m_ref_banda, m_p, v_p, ref_v, seed, m_ctrl, v_ctrl):
    """Resume por banda. `m_ref_banda` es la media por la que se asigna banda
    (la del gen perturbado). Agrega primero POR GEN sobre las perturbaciones y
    luego hace bootstrap sobre genes."""
    b = banda_de(m_ref_banda)                     # (n_pert, n_gen) o (n_gen,)
    exceso = v_p - ref_v
    fano_p = np.divide(v_p, m_p, out=np.full_like(v_p, np.nan), where=m_p > 0)
    fano_c = np.divide(v_ctrl, m_ctrl, out=np.full_like(v_ctrl, np.nan),
                       where=m_ctrl > 0)
    bc = banda_de(m_ctrl)
    filas = []
    for k in range(len(BORDES)):
        s = b == k
        n_celdas = int(s.sum())
        if n_celdas == 0:
            filas.append({"banda": k, "etiqueta": ETQ_BANDA[k], "n_celdas": 0})
            continue
        # agregacion por gen: media sobre las perturbaciones de ese gen
        col = s.any(axis=0) if s.ndim == 2 else s
        if s.ndim == 2:
            num = np.where(s, exceso, 0.0).sum(axis=0)
            den = s.sum(axis=0)
            ex_gen = np.divide(num, den, out=np.full(den.shape, np.nan),
                               where=den > 0)[col]
            fp = np.nanmean(np.where(s, fano_p, np.nan), axis=0)[col]
            vp_gen = np.divide(np.where(s, v_p, 0.0).sum(axis=0), den,
                               out=np.full(den.shape, np.nan), where=den > 0)[col]
            mp_gen = np.divide(np.where(s, m_p, 0.0).sum(axis=0), den,
                               out=np.full(den.shape, np.nan), where=den > 0)[col]
        else:
            ex_gen, fp, vp_gen, mp_gen = exceso[s], fano_p[s], v_p[s], m_p[s]
        media, lo, hi = boot_ci(ex_gen, seed + k)
        sc = bc == k
        filas.append({
            "banda": k, "etiqueta": ETQ_BANDA[k],
            "n_celdas": n_celdas, "n_genes": int(np.isfinite(ex_gen).sum()),
            "media_perturbada": float(np.nanmean(mp_gen)),
            "var_perturbada_media": float(np.nanmean(vp_gen)),
            "var_referencia_media": float(np.nanmean(
                np.where(s, ref_v, np.nan)) if s.ndim == 2 else np.nanmean(ref_v[s])),
            "fano_perturbada_mediana": float(np.nanmedian(fp)),
            "fano_control_mediana": (float(np.nanmedian(fano_c[sc]))
                                     if sc.any() else None),
            "n_genes_control_en_banda": int(sc.sum()),
            "exceso_medio": media, "ic95_lo": lo, "ic95_hi": hi,
            "frac_exceso_positivo": float(np.nanmean(ex_gen > 0)),
            "exceso_mediana": float(np.nanmedian(ex_gen)),
        })
    return filas


def una_linea(line, fold_idx, panel, var_names, cell_type, gene, n_genes, cache):
    import pds_objetivo as PO
    t0 = time.time()
    real_ad, _suelo, control_pool, meta = C.arma_linea(
        line, fold_idx, panel, var_names, cell_type, gene, n_genes)
    hue = C.huella_adata(real_ad)        # ANTES de renombrar: no se usa el marcador aqui
    lab = np.asarray(real_ad.obs["target_gene"].tolist(), dtype=object)
    Xr = real_ad.X.tocsr()
    G = Xr.shape[1]
    nP = len(panel)

    # 1. las 400 perturbadas reales, por perturbacion y gen
    m_real = np.zeros((nP, G)); v_real = np.zeros((nP, G))
    for i, h in enumerate(panel):
        mk = lab == h
        assert int(mk.sum()) == C.N_PER_LABEL, (h, int(mk.sum()))
        m_real[i], v_real[i], _ = media_var(Xr[mk])

    # 2. la linea base: 400 CONTROLES de la misma linea, sin reemplazo
    seed_ctrl = C.SEMILLAS_NUEVAS["b2_control"][fold_idx]
    rng = np.random.default_rng(seed_ctrl)
    n_ctrl_pool = control_pool.shape[0]
    sel = np.sort(rng.choice(n_ctrl_pool, size=C.N_PER_LABEL, replace=False))
    m_ctrl, v_ctrl, _ = media_var(control_pool[sel])
    cur = curva_referencia(m_ctrl, v_ctrl)

    # 3. la emision nuestra: ORACULO COMPLETO POR RAZON (1b), j=0
    ed = cache.get(line, {})
    ratios, diag = PO.razones_oraculo(real_ad, panel, ed, "completo")
    del real_ad
    fn = PO.predict_razon(ratios)
    seed_em = SEED_ORACULO_1B + 100 * fold_idx
    rng_em = np.random.default_rng(seed_em)
    pred = fn(control_pool, panel + [C.CONTROL_LABEL], rng_em)
    m_em = np.zeros((nP, G)); v_em = np.zeros((nP, G))
    for i, h in enumerate(panel):
        m_em[i], v_em[i], _ = media_var(pred[h].tocsr())
    del pred, ratios

    ref_real = evalua_referencia(cur, m_real)
    ref_em = evalua_referencia(cur, m_em)
    seed_boot = C.SEMILLAS_NUEVAS["b2_bootstrap"][fold_idx]
    bandas_real = por_banda(m_real, m_real, v_real, ref_real, seed_boot,
                            m_ctrl, v_ctrl)
    bandas_em = por_banda(m_em, m_em, v_em, ref_em, seed_boot + 1000,
                          m_ctrl, v_ctrl)

    # por perturbacion: los intermedios crudos
    fano_r = np.divide(v_real, m_real, out=np.full_like(v_real, np.nan),
                       where=m_real > 0)
    fano_e = np.divide(v_em, m_em, out=np.full_like(v_em, np.nan), where=m_em > 0)
    ex_r = v_real - ref_real
    ex_e = v_em - ref_em
    por_pert = []
    for i, h in enumerate(panel):
        viv = m_real[i] > 0
        por_pert.append({
            "pert": h,
            "n_genes_vivos": int(viv.sum()),
            "fano_real_mediana": float(np.nanmedian(fano_r[i][viv])),
            "fano_emision_mediana": float(np.nanmedian(fano_e[i][m_em[i] > 0])),
            "exceso_real_medio": float(np.nanmean(ex_r[i][viv])),
            "exceso_emision_medio": float(np.nanmean(ex_e[i][m_em[i] > 0])),
            "frac_exceso_real_positivo": float(np.nanmean(ex_r[i][viv] > 0)),
            "var_real_media": float(np.nanmean(v_real[i][viv])),
            "var_ref_media": float(np.nanmean(ref_real[i][viv])),
            "n_ed_real": int(len(ed.get(h, ()))),
        })
    fc = np.divide(v_ctrl, m_ctrl, out=np.full_like(v_ctrl, np.nan),
                   where=m_ctrl > 0)
    return {
        "linea": line, "fold_idx": fold_idx, "huella_real_sin_renombrar": hue,
        "n_control_pool": int(n_ctrl_pool), "seed_control_400": int(seed_ctrl),
        "seed_bootstrap": int(seed_boot), "seed_emision_1b": int(seed_em),
        "n_genes": int(G), "n_genes_control_media_cero": int((m_ctrl == 0).sum()),
        "fano_control_mediana": float(np.nanmedian(fc)),
        "fano_control_frac_mayor_1": float(np.nanmean(fc > 1)),
        "diag_oraculo_global": diag.get("_global"),
        "curva_referencia_n_bins": int(len(cur["ref"])),
        "bandas_real": bandas_real, "bandas_emision": bandas_em,
        "por_pert": por_pert, "dt_s": time.time() - t0,
    }


def md(res):
    o = []; A = o.append
    A("# Sobredispersion celula a celula: hay heterogeneidad por encima del ruido de conteo?")
    A("")
    A("## 🔴 0. Que mide esto y que no")
    A("")
    A("Las emisiones por razon de este proyecto aplican **un** vector de razones "
      "a las 400 celulas y realizan el conteo con `_reescala_entero` (Binomial "
      "para `r<=1`, Poisson aditivo para `r>1`): ruido de conteo alrededor de "
      "**una moda**. Esta medicion pregunta si las 400 celulas **reales** son "
      "mas dispersas que eso.")
    A("")
    A("🔴 **El cero correcto no es Poisson teorico.** Las celulas de control ya "
      "son una mezcla (ciclo, profundidad, lote), asi que ya estan "
      "sobredispersas. La linea base son 400 controles de la MISMA linea, "
      "comparados **a media igualada**.")
    A("")
    A("🔴 **El exceso puede salir cero o negativo.** Si sale, se dice. No hay "
      "nada en esta construccion que lo fuerce a ser positivo.")
    A("")
    A("## 1. La construccion, declarada antes de medir")
    A("")
    A("\n".join(C.tabla(["# ", "cantidad", "valor"], [
        [1, "espacio", "**cuentas crudas**, sin normalizar: es donde el ruido Binomial/Poisson tiene su forma conocida (`var = media` para Poisson)"],
        [2, "varianza", "`E[x^2] - E[x]^2` por gen sobre las celulas, `ddof = 0`"],
        [3, "lado real", f"400 celulas por perturbacion, `seed_real = {C.SEED_BASE_REAL} + fold` (la MISMA poblacion de todos los reportes)"],
        [4, "linea base", "400 celulas de CONTROL de la misma linea, **sin reemplazo**, semilla `20280401 + fold` (familia nueva)"],
        [5, "igualacion de media", f"curva `var_control(media)` en {N_BINS_CURVA} bins de `log10(media)` sobre los genes de control, evaluada por interpolacion en la media de cada gen perturbado"],
        [6, "bandas", f"{len(BORDES)} bandas de media de control, **declaradas antes de medir**: {', '.join(ETQ_BANDA)}"],
        [7, "IC", f"bootstrap de {N_BOOT} resampleos **sobre genes** (cada gen entra o sale completo), semilla `20280451 + fold`"],
        [8, "emision nuestra", "ORACULO COMPLETO POR RAZON (`pds_objetivo.md` §3.1, brazo 1b), `j = 0`, semilla ya gastada `20275000 + 100*fold`"],
        [9, "Fano", "`var / media`. Para Poisson vale **1**; para Binomial de un control remuestreado, `<= 1` mas la mezcla del propio control"],
    ])))
    A("")
    A("### 1.1 Por linea: los `n` y el control")
    A("")
    fil = []
    for i, r in enumerate(res["lineas"]):
        fil.append([i + 1, f"**{r['linea']}**", r["n_control_pool"],
                    C.N_PER_LABEL, r["n_genes"],
                    r["n_genes_control_media_cero"],
                    C.f4a(r["fano_control_mediana"]),
                    f"{100*r['fano_control_frac_mayor_1']:.1f} %",
                    r["curva_referencia_n_bins"], C.dur(r["dt_s"])])
    A("\n".join(C.tabla(
        ["# ", "linea", "pool de control", "controles usados", "genes",
         "genes con media de control 0", "**Fano del control, mediana**",
         "frac genes con Fano > 1", "bins de la curva", "duracion"], fil)))
    A("")
    A("🔴 La columna del Fano del control es la que fija la escala: **el control "
      "no es Poisson**, y por eso no se usa Poisson como cero.")
    A("")

    # ---- 2. el exceso por banda
    A("## 2. 🔴 EL EXCESO: varianza de las perturbadas menos la del control a media igualada")
    A("")
    A("Por banda de expresion. `exceso = var_perturbada - var_control(media "
      "igualada)`, agregado **primero por gen** (media sobre las 33 "
      "perturbaciones) y luego promediado sobre genes, con IC 95 % por "
      "bootstrap sobre genes.")
    A("")
    for r in res["lineas"]:
        A(f"### 2.{res['lineas'].index(r)+1} {r['linea']}")
        A("")
        fil = []
        for b in r["bandas_real"]:
            if not b.get("n_celdas"):
                continue
            e = next((x for x in r["bandas_emision"] if x["banda"] == b["banda"]
                      and x.get("n_celdas")), None)
            ic = (f"[{C.fe(b['ic95_lo'])}, {C.fe(b['ic95_hi'])}]"
                  if b.get("ic95_lo") is not None else "n/d")
            sig = ("si" if (b.get("ic95_lo") is not None and
                            b["ic95_lo"] > 0) else
                   ("🔴 NO (el IC cruza 0 o es negativo)"
                    if b.get("ic95_lo") is not None else "n/d"))
            fil.append([b["etiqueta"], b["n_celdas"], b["n_genes"],
                        C.fe(b["media_perturbada"]),
                        C.fe(b["var_perturbada_media"]),
                        C.fe(b["var_referencia_media"]),
                        C.fe(b["exceso_medio"]), ic, sig,
                        f"{100*b['frac_exceso_positivo']:.1f} %",
                        C.f4a(b["fano_perturbada_mediana"]),
                        C.f4a(b["fano_control_mediana"]),
                        C.fe(e["exceso_medio"]) if e else "n/d",
                        C.f4a(e["fano_perturbada_mediana"]) if e else "n/d"])
        A("\n".join(C.tabla(
            ["banda (media de control)", "celdas (pert x gen)", "genes",
             "media pert.", "var pert.", "var referencia",
             "**exceso medio**", "IC 95 %", "**exceso > 0?**",
             "frac genes con exceso > 0", "Fano real (med)",
             "Fano control (med)", "exceso de NUESTRA emision",
             "Fano de nuestra emision"], fil)))
        A("")

    # ---- 3. nuestra emision
    A("## 3. 🔴 La misma cantidad sobre UNA EMISION NUESTRA")
    A("")
    A("El ORACULO COMPLETO POR RAZON (1b) es la emision de **una moda** mas "
      "extrema disponible: lleva la razon real en todos los genes. Si su "
      "sobredispersion ya iguala la de las 400 reales, **no hay hueco que una "
      "emision de una moda deje sin llenar**, y eso es el resultado.")
    A("")
    fil = []
    for i, r in enumerate(res["lineas"]):
        er = [b["exceso_medio"] for b in r["bandas_real"]
              if b.get("exceso_medio") is not None]
        ee = [b["exceso_medio"] for b in r["bandas_emision"]
              if b.get("exceso_medio") is not None]
        fr = [b["fano_perturbada_mediana"] for b in r["bandas_real"]
              if b.get("n_celdas")]
        fe_ = [b["fano_perturbada_mediana"] for b in r["bandas_emision"]
               if b.get("n_celdas")]
        neg = sum(1 for b in r["bandas_real"]
                  if b.get("ic95_hi") is not None and b["ic95_hi"] < 0)
        pos = sum(1 for b in r["bandas_real"]
                  if b.get("ic95_lo") is not None and b["ic95_lo"] > 0)
        cruza = sum(1 for b in r["bandas_real"]
                    if b.get("ic95_lo") is not None
                    and b["ic95_lo"] <= 0 <= b["ic95_hi"])
        fil.append([i + 1, f"**{r['linea']}**", C.rango(er), C.rango(ee),
                    C.rango(fr), C.rango(fe_), pos, cruza, neg])
    A("\n".join(C.tabla(
        ["# ", "linea", "rango del exceso REAL sobre bandas",
         "rango del exceso de NUESTRA emision",
         "rango del Fano real", "rango del Fano de nuestra emision",
         "bandas con IC > 0", "bandas con IC que cruza 0", "bandas con IC < 0"],
        fil)))
    A("")
    A("## 4. 🔴 Los intermedios crudos: por perturbacion y linea")
    A("")
    for r in res["lineas"]:
        A(f"### 4.{res['lineas'].index(r)+1} {r['linea']} -- las 33 filas")
        A("")
        fil = []
        for i, p in enumerate(r["por_pert"]):
            fil.append([i + 1, f"`{p['pert']}`", p["n_ed_real"],
                        p["n_genes_vivos"],
                        C.f4a(p["fano_real_mediana"]),
                        C.f4a(p["fano_emision_mediana"]),
                        C.fe(p["var_real_media"]), C.fe(p["var_ref_media"]),
                        C.fe(p["exceso_real_medio"]),
                        C.fe(p["exceso_emision_medio"]),
                        f"{100*p['frac_exceso_real_positivo']:.1f} %"])
        A("\n".join(C.tabla(
            ["# ", "pert", "\\|ED real\\|", "genes con media > 0",
             "Fano real (med)", "Fano emision (med)", "var real (media)",
             "var referencia (media)", "**exceso real (media)**",
             "exceso emision (media)", "frac genes con exceso > 0"], fil)))
        A("")
    A("## 5. Defectos de la construccion, como datos sobre la medicion")
    A("")
    A("\n".join(C.tabla(["# ", "defecto", "que implica"], [
        [1, "Las 400 celulas reales de cada perturbacion incluyen **toda** la heterogeneidad tecnica de la linea (profundidad, lote), no solo la biologica.", "El exceso medido es una cota SUPERIOR de la heterogeneidad de respuesta: parte de el es tecnico."],
        [2, "Los 400 controles son **una** muestra, con una semilla.", "No hay barra de error sobre la propia linea base. El IC es sobre genes, no sobre la eleccion de los 400 controles."],
        [3, "La curva de referencia se interpola y se **sujeta** a los extremos fuera de rango.", "En los genes perturbados cuya media cae por encima de la media de control maxima, la referencia es la del extremo: el exceso ahi esta sobreestimado."],
        [4, "La banda se asigna por la media **perturbada**, no por la del control.", "Es lo que hace la igualacion de media; cambia la asignacion de banda de los genes que se mueven mucho."],
        [5, "La agregacion por gen promedia las 33 perturbaciones antes del bootstrap.", "El IC es sobre genes, asi que no habla de la variacion entre perturbaciones. Esa queda en la tabla de intermedios §4."],
        [6, "Una sola realizacion de nuestra emision (`j = 0`).", "No hay dispersion sobre semillas de emision en esta tabla."],
        [7, "`var = E[x^2] - E[x]^2` con `ddof = 0` sobre 400 celulas.", "Sesgo de `-var/400` respecto del estimador insesgado, el MISMO en los tres lados, asi que no mueve la diferencia."],
    ])))
    A("")
    A("## 6. Las entradas")
    A("")
    A("\n".join(C.tabla(["entrada", "valor"], [
        ["IFNB", "`/datos/ref/convertido/Seurat_object_IFNB_Perturb_seq.h5ad`"],
        ["panel", "`modelo/out/suelo_panel.json` (33)"],
        ["emision", "`pds_objetivo.razones_oraculo` + `pds_objetivo.predict_razon`, SIN MODIFICAR"],
        ["conteo entero", "`cand01_covar._reescala_entero`, SIN MODIFICAR"],
        ["salidas", "`modelo/out/sobredispersion.md`, `modelo/out/sobredispersion.json`"],
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

    res = {"bloque": "B2", "semillas": C.comprueba_semillas(),
           "bordes_banda": [float(x) for x in BORDES],
           "etiquetas_banda": ETQ_BANDA, "n_boot": N_BOOT,
           "lineas": [],
           "t_inicio": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(t_ini))}
    for i, line in enumerate(C.LINES):
        C.log(f"B2 {line} ({i+1}/6)")
        C.guarda_disco()
        r = una_linea(line, i, panel, var_names, cell_type, gene, n_genes, cache)
        res["lineas"].append(r)
        res["t_fin"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
        res["dt_total_s"] = time.time() - t_ini
        C.escribe_json(OUT_JSON, res)
        C.log(f"B2 {line}: fano_ctrl={r['fano_control_mediana']:.3f} dt={r['dt_s']:.0f}s")
    C.escribe_md(OUT_MD, md(res))
    C.log(f"B2 listo: {OUT_MD}")


if __name__ == "__main__":
    main()
