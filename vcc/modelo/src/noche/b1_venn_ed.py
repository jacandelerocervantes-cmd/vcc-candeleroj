#!/usr/bin/env python
"""B1 - 3c PARTE C: el Venn de los conjuntos ED entre contextos (Jiang, 6 lineas,
33 perturbaciones).

CONSTRUCCION, DECLARADA ANTES DE MEDIR
  * El conjunto ED REAL se obtiene con la MAQUINARIA DEL PROPIO MARCADOR, igual
    que `pds_objetivo.md` §3.2: una SONDA por linea corre `cell-eval2 run
    --write-degenes` sobre la poblacion real (que depende solo de
    `seed_real = 20260001 + fold`) y de ahi sale `R(L,h)` con
    `cell_eval2.de.prepare_de` + `p_adj < 0.05` + `_exclude_own_gene` -- las
    mismas tres lineas de `techo_arnes.contar_de_por_perturbacion`.
  * Se guarda TAMBIEN el `p_adj` de cada gen significativo, que el cache
    `/datos/modelo/out/pds_obj/ed_real.json` no tiene y que el recorte por
    potencia necesita.
  * La poblacion real se verifica por sha256 contra las huellas ya guardadas en
    los json de disco. El cache de conjuntos se recomprueba gen a gen.

🔴 ESTO ES UNA COTA DE JIANG. No se extrapola a D/E/F en ningun sentido.
"""
import json
import os
import sys
import time
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import comun as C

sys.path.insert(0, C.SRC)

WORK = os.path.join(C.WORK_BASE, "b1")
OUT_JSON = os.path.join(C.OUT, "venn_ed.json")
OUT_MD = os.path.join(C.OUT, "venn_ed.md")
ED_CACHE = "/datos/modelo/out/pds_obj/ed_real.json"
SEED_DESEMPATE = C.SEMILLAS_NUEVAS["b1_desempate"][0]
CORTES = [6, 5, 4, 3, 2, 1]


# ------------------------------------------------------------------ la sonda
def sonda(line, fold_idx, panel, var_names, cell_type, gene, n_genes):
    """Devuelve ({pert: [genes ED]}, {pert: {gen: p_adj}}, huella_real, dt)."""
    import polars as pl
    from cell_eval2.de import prepare_de
    from cell_eval2.metrics.de import _exclude_own_gene
    import techo_arnes as TA

    t0 = time.time()
    real_ad, pred_ad, _ctrl, meta = C.arma_linea(
        line, fold_idx, panel, var_names, cell_type, gene, n_genes)
    for ad in (real_ad, pred_ad):
        n = int((ad.obs["target_gene"] == C.CONTROL_LABEL).sum())
        ad.obs["target_gene"] = ad.obs["target_gene"].replace(
            {C.CONTROL_LABEL: C.CONTROL_COMP})
        assert int((ad.obs["target_gene"] == C.CONTROL_COMP).sum()) == n
    h_real = C.huella_adata(real_ad)

    rp = os.path.join(WORK, f"sonda_real_{line}.h5ad")
    pp = os.path.join(WORK, f"sonda_pred_{line}.h5ad")
    rd = os.path.join(WORK, f"sonda_run_{line}")
    real_ad.write_h5ad(rp)
    pred_ad.write_h5ad(pp)
    del real_ad, pred_ad
    cmd = [C.CELLEVAL2_PY, "run", "-ap", pp, "-ar", rp,
           "--preset", "vcc2026", "--profile", "vcc2026",
           "--pert-col", "target_gene", "--write-degenes", "-o", rd]
    p, dt_run = C.corre(cmd, os.path.join(WORK, f"sonda_{line}.log"))
    if p.returncode != 0:
        C.borra(rp, pp)
        raise RuntimeError(f"sonda {line} fallo (exit={p.returncode}): "
                           f"{p.stderr[-1500:]}")

    prepared = prepare_de(
        os.path.join(rd, "de_pred.parquet"), os.path.join(rd, "de_real.parquet"),
        control=C.CONTROL_COMP, sort_by=TA.DE_SORT_BY,
        p_adj_threshold=TA.DE_P_ADJ_THRESHOLD,
        nan_lfc_policy=TA.DE_NAN_LFC_POLICY, min_abs_log2fc=TA.DE_MIN_ABS_LOG2FC)
    T = prepared.p_adj_threshold
    sig = _exclude_own_gene(
        prepared.real_df.filter(pl.col("p_adj") < T)
        .select("target", "feature", "p_adj").unique(),
        prepared, who="b1_venn_ed", side="real")
    d = sig.to_dict(as_series=False)
    ed, pv = {h: [] for h in panel}, {h: {} for h in panel}
    for t, f, q in zip(d["target"], d["feature"], d["p_adj"]):
        t, f = str(t), str(f)
        if t not in ed:
            ed[t], pv[t] = [], {}
        ed[t].append(f)
        pv[t][f] = float(q)
    ed = {h: sorted(v) for h, v in ed.items()}
    C.borra(rp, pp, rd)
    return ed, pv, h_real, time.time() - t0, dt_run, float(T), meta


# ------------------------------------------------------------------ medicion
def jac(a, b):
    """Jaccard crudo. Devuelve (valor, union_vacia)."""
    if not a and not b:
        return None, True
    return len(a & b) / len(a | b), False


def cuenta_lineas(ed_por_linea, pert, lineas):
    c = Counter()
    for L in lineas:
        for g in ed_por_linea[L].get(pert, ()):
            c[g] += 1
    return c


def barrido(ed_por_linea, panel, lineas):
    """🔴 El barrido del corte de frecuencia. S_k = {genes llamados en >= k
    lineas}. Para cada k se da la media de `jac` sobre las 6 lineas con las DOS
    convenciones declaradas, y el maximo sobre k es la COTA INFERIOR del techo,
    restringida a ESTA familia de S (los cortes de frecuencia)."""
    out = {}
    for pert in panel:
        cnt = cuenta_lineas(ed_por_linea, pert, lineas)
        fila = {"n_genes_llamados_alguna": len(cnt),
                "n_por_corte": {}, "cortes": {}}
        for k in CORTES:
            S = {g for g, v in cnt.items() if v >= k}
            vals, uv = [], 0
            por_linea = {}
            for L in lineas:
                R = set(ed_por_linea[L].get(pert, ()))
                v, vacia = jac(S, R)
                por_linea[L] = v
                if vacia:
                    uv += 1
                else:
                    vals.append(v)
            # convencion del marcador: union vacia -> 1.0 (metrics/de.py)
            con_conv = vals + [1.0] * uv
            fila["n_por_corte"][k] = len(S)
            fila["cortes"][k] = {
                "n_S": len(S),
                "media_jac_conv_marcador": float(np.mean(con_conv)) if con_conv else None,
                "media_jac_sin_union_vacia": float(np.mean(vals)) if vals else None,
                "n_celdas_union_vacia": uv,
                "por_linea": por_linea,
                "min": (min(vals) if vals else None),
                "max": (max(vals) if vals else None),
            }
        for clave in ("media_jac_conv_marcador", "media_jac_sin_union_vacia"):
            cand = [(fila["cortes"][k][clave], k) for k in CORTES
                    if fila["cortes"][k][clave] is not None]
            fila[f"max_{clave}"] = (max(cand)[0] if cand else None)
            fila[f"argmax_{clave}"] = (max(cand)[1] if cand else None)
        out[pert] = fila
    return out


def recorta_por_potencia(ed_por_linea, pv_por_linea, panel, lineas):
    """🔴 IGUALACION DE POTENCIA. Las medianas de |ED| van de 1 (A549) a 197
    (MCF7), asi que el solape lo acota el conjunto mas chico. Se recortan las
    seis al tamano del conjunto ED mas pequeno NO VACIO de esa perturbacion,
    tomando los de `p_adj` mas baja. Empates: clave aleatoria por gen con
    semilla declarada SEED_DESEMPATE; los conjuntos vacios SIGUEN vacios."""
    rng = np.random.default_rng(SEED_DESEMPATE)
    ed2, diag = {L: {} for L in lineas}, {}
    for pert in panel:
        tam = [len(ed_por_linea[L].get(pert, ())) for L in lineas]
        no_vacios = [t for t in tam if t > 0]
        m = min(no_vacios) if no_vacios else 0
        diag[pert] = {"tamanos_crudos": dict(zip(lineas, tam)),
                      "m_objetivo": m,
                      "n_lineas_no_vacias": len(no_vacios),
                      "n_empates_en_el_corte": {}}
        for L in lineas:
            gs = list(ed_por_linea[L].get(pert, ()))
            if not gs or m == 0:
                ed2[L][pert] = []
                continue
            claves = rng.random(len(gs))
            orden = sorted(range(len(gs)),
                           key=lambda i: (pv_por_linea[L][pert].get(gs[i], 1.0),
                                          claves[i]))
            sel = [gs[i] for i in orden[:m]]
            if m < len(gs):
                p_corte = pv_por_linea[L][pert].get(gs[orden[m - 1]], 1.0)
                emp = sum(1 for g in gs
                          if pv_por_linea[L][pert].get(g, 1.0) == p_corte)
                diag[pert]["n_empates_en_el_corte"][L] = emp
            ed2[L][pert] = sorted(sel)
    return ed2, diag


def distribuciones(ed_por_linea, panel, lineas):
    """Dos definiciones distintas, las dos reportadas:
       A (pares pert-gen): suma de las distribuciones POR PERTURBACION.
       B (genes distintos): un gen cuenta las lineas en que aparece en ALGUNA
         perturbacion. NO son la misma cantidad."""
    A, por_pert = Counter(), {}
    for pert in panel:
        c = cuenta_lineas(ed_por_linea, pert, lineas)
        h = Counter(c.values())
        por_pert[pert] = {k: int(h.get(k, 0)) for k in CORTES}
        for k, v in h.items():
            A[k] += v
    c2 = Counter()
    for L in lineas:
        s = set()
        for pert in panel:
            s |= set(ed_por_linea[L].get(pert, ()))
        for g in s:
            c2[g] += 1
    B = Counter(c2.values())
    return ({k: int(A.get(k, 0)) for k in CORTES}, int(sum(A.values())),
            {k: int(B.get(k, 0)) for k in CORTES}, len(c2), por_pert)


# ------------------------------------------------------------------ md
def md(res):
    L = res["lineas"]
    panel = res["panel"]
    o = []
    A = o.append
    A("# El Venn de los conjuntos ED entre contextos (Jiang, 6 lineas, 33 perturbaciones)")
    A("")
    A("## 🔴 0. Encabezado obligatorio: esto es una COTA DE JIANG")
    A("")
    A("Todo lo que sigue se mide sobre las **seis lineas celulares de Jiang**, "
      "condicion IFNB, con el panel comun de 33 perturbaciones. **No se "
      "extrapola a D/E/F en ningun sentido, ni como cota superior ni como cota "
      "inferior.**")
    A("")
    A("La razon esta MEDIDA, no supuesta. La heterogeneidad que produce esta "
      "cota -- que las seis lineas tengan conjuntos ED de tamanos y contenidos "
      "muy distintos -- esta **medida como AUSENTE en el reto** "
      "(`ntc_nulo.md` §3):")
    A("")
    A(*["\n".join(C.tabla(
        ["# ", "cantidad", "reto 2026: max/min entre contextos",
         "Jiang: max/min entre lineas"],
        [[1, "celulas de control por unidad", "**x1.000**", "**x2.341**"],
         [2, "profundidad mediana (UMIs/celula)", "**x1.028**", "**x1.869**"],
         [3, "`nnz` mediana por celula", "**x1.092**", "**x1.602**"],
         [4, "panel probable (el `m` de BH)", "**x1.062**", "**x1.077**"]]))])
    A("")
    A("Y el nulo de `ntc_id` en contexto **da cero** conjunto ED falso positivo "
      "(`ntc_nulo.md` §2.1). Las dos cosas juntas dicen que el mecanismo que "
      "genera esta cota en Jiang **no esta medido como presente en el reto**. "
      "Lo que esta cota ES: una medicion de Jiang. Lo que NO es: una prediccion "
      "sobre D/E/F.")
    A("")
    A("## 1. La construccion, declarada antes de medir")
    A("")
    A("\n".join(C.tabla(["# ", "cantidad", "valor"], [
        [1, "conjunto ED real", "`cell-eval2 run --write-degenes` sobre la poblacion real + `prepare_de` + `p_adj < 0.05` + `_exclude_own_gene` -- la maquinaria del propio marcador, igual que `pds_objetivo.md` §3.2"],
        [2, "umbral", f"`p_adj < {res['p_adj_threshold']}` (el del preset `vcc2026`)"],
        [3, "poblacion real", f"`seed_real = {C.SEED_BASE_REAL} + fold`, la MISMA de siempre, verificada por sha256"],
        [4, "sondas", f"{len(L)}, una por linea (el lado real depende solo de `seed_real`)"],
        [5, "genes del espacio", res["n_genes"]],
        [6, "panel", f"{len(panel)} perturbaciones, el comun de las 6 lineas"],
        [7, "semilla de desempate del recorte", f"`{SEED_DESEMPATE}` (familia nueva 20280000+)"],
        [8, "`jac`", "interseccion / union, crudo. **Dos convenciones declaradas** para la union vacia: la del marcador (`= 1.0`) y excluyendo esas celdas"],
    ])))
    A("")
    A("### 1.1 Verificacion de la poblacion real y recomprobacion del cache")
    A("")
    fil = []
    for i, line in enumerate(L):
        v = res["verif"][line]
        fil.append([i + 1, f"**{line}**", v["huella_real"][:16] + "...",
                    "si" if v["huella_coincide_con_disco"] else "🔴 NO",
                    v["n_json_que_la_comparten"],
                    v["n_perts_identicas_al_cache"],
                    v["n_perts_distintas_al_cache"],
                    v["n_genes_solo_ahora"], v["n_genes_solo_cache"],
                    C.dur(v["dt_total_s"])])
    A("\n".join(C.tabla(
        ["# ", "linea", "huella real (sha256, 16)", "coincide con la de disco",
         "# json que la comparten", "perts identicas al cache",
         "perts distintas", "genes solo ahora", "genes solo en cache", "duracion"],
        fil)))
    A("")
    A("El cache `/datos/modelo/out/pds_obj/ed_real.json` se recomprobo **gen a "
      "gen**, no por conteo. "
      f"Perturbaciones identicas: **{res['cache_identicas_total']} de "
      f"{len(L) * len(panel)}**.")
    A("")

    # ---- 2. tamanos
    A("## 2. Tamano del conjunto ED real, por linea y por perturbacion")
    A("")
    A("🔴 **Rango, nunca media.** Con 15 de 33 perturbaciones vacias en A549, "
      "el techo ahi es **0 por construccion en casi medio panel**, y una media "
      "lo esconderia.")
    A("")
    fil = []
    for i, line in enumerate(L):
        t = [len(res["ed"][line][h]) for h in panel]
        fil.append([i + 1, f"**{line}**", sum(1 for x in t if x == 0),
                    sum(1 for x in t if x > 0), C.rango_i(t),
                    int(np.median(t)), int(sum(t))])
    A("\n".join(C.tabla(
        ["# ", "linea", "**# perts con ED vacio**", "# con ED no vacio",
         "rango [min, max] de |ED|", "mediana |ED|", "suma |ED|"], fil)))
    A("")
    A("La columna de vacios reproduce `control_permutado.md` §7 "
      f"(15, 0, 14, 6, 8, 0): **{res['vacios_vs_cp']}**.")
    A("")
    A("### 🔴 TABLA 2.1 -- los intermedios crudos: |ED| por perturbacion y linea, las 33 filas")
    A("")
    fil = []
    for i, h in enumerate(panel):
        t = [len(res["ed"][L_][h]) for L_ in L]
        nv = sum(1 for x in t if x == 0)
        fil.append([i + 1, f"`{h}`"] + t + [nv, C.rango_i(t),
                                            res["dist"]["por_pert"][h][6],
                                            res["dist"]["por_pert"][h][1]])
    A("\n".join(C.tabla(["# ", "pert"] + [f"**{x}**" for x in L] +
                        ["# lineas vacias", "rango |ED|", "genes en 6", "genes en 1"], fil)))
    A("")

    # ---- 3. distribucion
    A("## 3. Cuantos genes se llaman en 6, 5, 4, 3, 2 y 1 sola linea")
    A("")
    A("🔴 **Dos definiciones distintas, las dos reportadas, porque no son la "
      "misma cantidad.**")
    A("")
    A("* **A -- pares (perturbacion, gen)**: para cada perturbacion se cuenta "
      "en cuantas lineas se llama cada gen, y se SUMAN las 33 distribuciones. "
      "Es la definicion que cita el encargo.")
    A("* **B -- genes distintos**: un gen cuenta las lineas en que aparece en "
      "**alguna** perturbacion. Colapsa el panel, asi que da otro numero.")
    A("")
    for etq, clave in (("(a) CRUDA", "crudo"), ("(b) IGUALADA EN POTENCIA", "igualado")):
        d = res["dist"] if clave == "crudo" else res["dist_ig"]
        A(f"**{etq}**")
        A("")
        fil = []
        for k in CORTES:
            fil.append([k, d["A"][k], f"{100*d['A'][k]/max(d['A_total'],1):.2f} %",
                        d["B"][k], f"{100*d['B'][k]/max(d['B_total'],1):.2f} %"])
        A("\n".join(C.tabla(["# lineas en que se llama",
                             "**A: pares (pert, gen)**", "frac A",
                             "**B: genes distintos**", "frac B"], fil)))
        A("")
        A(f"Total A = **{d['A_total']}** pares, total B = **{d['B_total']}** genes distintos.")
        A("")
    A(f"🔴 **La recomprobacion del encargo.** El cache citado decia 19,068 de "
      f"21,004 en una sola linea y 126 en las seis. Medido ahora, definicion A: "
      f"**{res['dist']['A'][1]} de {res['dist']['A_total']}** en una sola y "
      f"**{res['dist']['A'][6]}** en las seis. "
      f"{'Coincide' if (res['dist']['A'][1], res['dist']['A_total'], res['dist']['A'][6]) == (19068, 21004, 126) else '🔴 NO coincide'}.")
    A("")

    # ---- 4. barrido
    A("## 4. 🔴 El techo de `jac` para una emision IDENTICA en las seis lineas")
    A("")
    A("**La interseccion NO es el optimo.** Con `R_1={a,b}` y `R_2={a,c}`, "
      "`S={a}` da media `0.5000` y `S={a,b}` da `0.6667`. Por eso no se toma la "
      "interseccion: se **barre** el corte de frecuencia. `S_k` = los genes "
      "llamados en **k o mas** lineas, con `k = 6, 5, 4, 3, 2, 1`.")
    A("")
    A("🔴 El maximo de la curva es una **COTA INFERIOR del techo, restringida a "
      "esta familia de `S`** (los seis cortes de frecuencia). **No es el techo "
      "exacto**: el optimo sobre TODOS los subconjuntos posibles no se calcula "
      "aqui y puede ser mayor.")
    A("")
    for etq, clave in (("(a) CRUDA", "barr"), ("(b) IGUALADA EN POTENCIA", "barr_ig")):
        b = res[clave]
        A(f"### 4.{1 if clave == 'barr' else 2} {etq} -- la curva entera, agregada sobre las 33 perturbaciones")
        A("")
        fil = []
        for k in CORTES:
            nS = [b[h]["cortes"][k]["n_S"] for h in panel]
            mc = [b[h]["cortes"][k]["media_jac_conv_marcador"] for h in panel]
            ms = [b[h]["cortes"][k]["media_jac_sin_union_vacia"] for h in panel]
            uv = sum(b[h]["cortes"][k]["n_celdas_union_vacia"] for h in panel)
            fil.append([f"`>= {k}`", C.rango_i(nS), int(np.median(nS)),
                        C.rango([x for x in mc if x is not None]),
                        C.f4a(float(np.mean([x for x in mc if x is not None]))) if any(x is not None for x in mc) else "n/d",
                        C.rango([x for x in ms if x is not None]),
                        C.f4a(float(np.mean([x for x in ms if x is not None]))) if any(x is not None for x in ms) else "n/d",
                        uv])
        A("\n".join(C.tabla(
            ["corte", "rango |S| sobre las 33 perts", "mediana |S|",
             "rango de la media de `jac` (conv. marcador)", "media de medias",
             "rango (sin union vacia)", "media de medias", "# celdas union vacia"],
            fil)))
        A("")
        mx = [(b[h]["max_media_jac_conv_marcador"], h) for h in panel
              if b[h]["max_media_jac_conv_marcador"] is not None]
        mxs = [(b[h]["max_media_jac_sin_union_vacia"], h) for h in panel
               if b[h]["max_media_jac_sin_union_vacia"] is not None]
        A(f"**El maximo de la curva, como COTA INFERIOR del techo (familia de "
          f"cortes de frecuencia):** con la convencion del marcador, el rango "
          f"sobre las 33 perturbaciones es {C.rango([x for x, _ in mx])} y el "
          f"maximo lo da `{max(mx)[1]}` con **{C.f4a(max(mx)[0])}**; sin las "
          f"celdas de union vacia, rango {C.rango([x for x, _ in mxs])}, maximo "
          f"`{max(mxs)[1]}` con **{C.f4a(max(mxs)[0])}**.")
        A("")
        A(f"#### 🔴 TABLA 4.{1 if clave == 'barr' else 2}b -- los intermedios crudos: la curva por perturbacion, las 33 filas")
        A("")
        fil = []
        for i, h in enumerate(panel):
            r = b[h]
            fil.append([i + 1, f"`{h}`", r["n_genes_llamados_alguna"]] +
                       [r["cortes"][k]["n_S"] for k in CORTES] +
                       [C.f4a(r["cortes"][k]["media_jac_conv_marcador"]) for k in CORTES] +
                       [f"`>= {r['argmax_media_jac_conv_marcador']}`"
                        if r["argmax_media_jac_conv_marcador"] else "n/d",
                        C.f4a(r["max_media_jac_conv_marcador"])])
        A("\n".join(C.tabla(
            ["# ", "pert", "genes llamados en alguna"] +
            [f"\\|S\\| >={k}" for k in CORTES] +
            [f"jac >={k}" for k in CORTES] + ["argmax", "**max (cota inf.)**"], fil)))
        A("")
    A("### 4.3 Rango POR LINEA en el corte de cada maximo, nunca media")
    A("")
    fil = []
    for i, line in enumerate(L):
        vs = []
        for h in panel:
            k = res["barr"][h]["argmax_media_jac_conv_marcador"]
            if k:
                v = res["barr"][h]["cortes"][k]["por_linea"][line]
                if v is not None:
                    vs.append(v)
        nv = sum(1 for h in panel if not res["ed"][line][h])
        fil.append([i + 1, f"**{line}**", nv, len(vs), C.rango(vs),
                    C.f4a(float(np.median(vs))) if vs else "n/d"])
    A("\n".join(C.tabla(
        ["# ", "linea", "# perts con ED real vacio", "# celdas con jac definido",
         "rango de `jac` en el corte del maximo", "mediana"], fil)))
    A("")
    A("🔴 En A549 y HAP1 el techo es **0 por construccion** en 15 y 14 de las 33 "
      "perturbaciones: no hay conjunto real con el que solapar. Esa es la razon "
      "de no dar media de panel.")
    A("")

    # ---- 5. la diferencia (a) - (b)
    A("## 5. 🔴 LA CIFRA QUE HABLA DE CONTEXTO: la diferencia entre (a) y (b)")
    A("")
    A("Las dos versiones son **construcciones distintas**, declaradas como "
      "distintas. (a) usa los conjuntos tal como salen. (b) recorta las seis al "
      "tamano del conjunto ED mas pequeno **no vacio** de esa perturbacion, "
      "tomando los de `p_adj` mas baja. Si (b) sube mucho respecto de (a), lo "
      "que limitaba el solape era el **tamano** (la potencia); si no sube, lo "
      "que limitaba era **que genes** (el contexto).")
    A("")
    fil = []
    for i, h in enumerate(panel):
        a = res["barr"][h]["max_media_jac_conv_marcador"]
        bb = res["barr_ig"][h]["max_media_jac_conv_marcador"]
        d = res["dig"][h]
        fil.append([i + 1, f"`{h}`", d["m_objetivo"], d["n_lineas_no_vacias"],
                    C.f4a(a), C.f4a(bb),
                    C.f4(bb - a) if (a is not None and bb is not None) else "n/d",
                    res["barr"][h]["n_genes_llamados_alguna"],
                    res["barr_ig"][h]["n_genes_llamados_alguna"]])
    A("\n".join(C.tabla(
        ["# ", "pert", "m (tamano del mas chico no vacio)", "# lineas no vacias",
         "max jac (a) cruda", "max jac (b) igualada", "**(b) - (a)**",
         "genes llamados (a)", "genes llamados (b)"], fil)))
    A("")
    ds = [res["barr_ig"][h]["max_media_jac_conv_marcador"] -
          res["barr"][h]["max_media_jac_conv_marcador"]
          for h in panel
          if res["barr"][h]["max_media_jac_conv_marcador"] is not None
          and res["barr_ig"][h]["max_media_jac_conv_marcador"] is not None]
    A(f"Rango de `(b) - (a)` sobre las {len(ds)} perturbaciones donde las dos "
      f"estan definidas: **{C.rango(ds)}**, mediana "
      f"**{C.f4(float(np.median(ds))) if ds else 'n/d'}**. "
      "Se presenta la diferencia; no se concluye de ella.")
    A("")
    A("### 5.1 Los intermedios del recorte: tamanos crudos, `m`, y empates en el corte")
    A("")
    fil = []
    for i, h in enumerate(panel):
        d = res["dig"][h]
        fil.append([i + 1, f"`{h}`"] + [d["tamanos_crudos"][x] for x in L] +
                   [d["m_objetivo"],
                    ", ".join(f"{k}:{v}" for k, v in d["n_empates_en_el_corte"].items()) or "-"])
    A("\n".join(C.tabla(["# ", "pert"] + [f"\\|ED\\| {x}" for x in L] +
                        ["**m**", "empates en el corte, por linea"], fil)))
    A("")

    # ---- 6. defectos
    A("## 6. Defectos de la construccion, como datos sobre la medicion")
    A("")
    A("\n".join(C.tabla(["# ", "defecto", "que implica"], [
        [1, "El maximo del barrido es una **cota inferior** del techo, restringida a la familia de los seis cortes de frecuencia.", "El optimo sobre los `2^n` subconjuntos no se calcula. El techo verdadero de esta familia de emisiones es **>=** lo reportado."],
        [2, "La convencion de union vacia mueve la media.", "Por eso se dan las DOS columnas. Con la del marcador, una perturbacion vacia en las seis lineas y `S` vacio aporta `1.0`."],
        [3, "Una sola realizacion de la poblacion real (`seed_real` fijo).", "El lado real no se remuestrea, asi que no hay barra de error sobre `R(L,h)`. Es el mismo lado real de todos los reportes previos, a proposito."],
        [4, "El recorte (b) usa `p_adj` del propio ajuste BH de cada linea.", "Los `p_adj` no son comparables entre lineas (distinto `m` de BH). El recorte ordena DENTRO de cada linea, nunca entre lineas."],
        [5, "`m` es el tamano del mas chico NO VACIO.", "Las lineas con ED vacio siguen vacias en (b): el recorte no las rellena. En A549, 15 de 33 siguen aportando `jac = 0`."],
        [6, "El espacio es el de Jiang (34,025 genes), no el del reto (18,533).", "El panel probable de BH difiere (`ntc_nulo.md` §3: x0.79 a x0.91). Otra razon para no extrapolar."],
        [7, "`_exclude_own_gene` quita el propio blanco del conjunto.", "Es la convencion del marcador, no una eleccion de este reporte."],
    ])))
    A("")
    A("## 7. Las entradas")
    A("")
    A("\n".join(C.tabla(["entrada", "valor"], [
        ["IFNB", "`/datos/ref/convertido/Seurat_object_IFNB_Perturb_seq.h5ad`"],
        ["panel", "`modelo/out/suelo_panel.json` (`panel_comun`, 33)"],
        ["cache recomprobado", "`/datos/modelo/out/pds_obj/ed_real.json`"],
        ["marcador", "`cell-eval2` 0.20.0, preset/profile `vcc2026`"],
        ["salidas", "`modelo/out/venn_ed.md`, `modelo/out/venn_ed.json`"],
    ])))
    A("")
    A(f"*Corrida: {res['t_inicio']} -> {res['t_fin']} UTC, "
      f"duracion {C.dur(res['dt_total_s'])}.*")
    return o


def main():
    t_ini = time.time()
    os.makedirs(WORK, exist_ok=True)
    C.guarda_disco()
    from suelo_arnes import load_global_obs, load_var_names, IFNB_PATH
    panel = C.carga_panel()
    cell_type, gene = load_global_obs(IFNB_PATH)
    var_names = load_var_names(IFNB_PATH)
    n_genes = len(var_names)
    cache = json.load(open(ED_CACHE)) if os.path.exists(ED_CACHE) else {}
    ref_huellas = C.huellas_guardadas()

    ed, pv, verif = {}, {}, {}
    T = None
    for i, line in enumerate(C.LINES):
        C.log(f"B1 sonda {line} ({i+1}/6)")
        C.guarda_disco()
        e, p, h, dt, dtr, T, meta = sonda(line, i, panel, var_names,
                                          cell_type, gene, n_genes)
        ed[line], pv[line] = e, p
        cen = cache.get(line, {})
        iguales = sum(1 for x in panel
                      if set(e.get(x, ())) == set(cen.get(x, ())))
        solo_ahora = sum(len(set(e.get(x, ())) - set(cen.get(x, ()))) for x in panel)
        solo_cache = sum(len(set(cen.get(x, ())) - set(e.get(x, ()))) for x in panel)
        verif[line] = {
            "huella_real": h,
            "huella_coincide_con_disco": h in ref_huellas.get(line, {}),
            "n_json_que_la_comparten": len(ref_huellas.get(line, {}).get(h, [])),
            "huella_cache_declarada": cache.get(f"{line}__huella_real"),
            "n_perts_identicas_al_cache": iguales,
            "n_perts_distintas_al_cache": len(panel) - iguales,
            "n_genes_solo_ahora": solo_ahora,
            "n_genes_solo_cache": solo_cache,
            "dt_total_s": dt, "dt_run_s": dtr,
            "n_control_cells": meta.get("n_control_cells"),
        }
        C.log(f"B1 {line}: |ED| suma={sum(len(v) for v in e.values())} "
              f"iguales_al_cache={iguales}/33 dt={dt:.0f}s")
        C.escribe_json(OUT_JSON, {"parcial": True, "verif": verif,
                                  "ed": ed, "p_adj_threshold": T})

    A, At, B, Bt, por_pert = distribuciones(ed, panel, C.LINES)
    barr = barrido(ed, panel, C.LINES)
    ed_ig, dig = recorta_por_potencia(ed, pv, panel, C.LINES)
    A2, At2, B2, Bt2, pp2 = distribuciones(ed_ig, panel, C.LINES)
    barr_ig = barrido(ed_ig, panel, C.LINES)

    vac = {L: sum(1 for h in panel if not ed[L][h]) for L in C.LINES}
    esperado = {"A549": 15, "BXPC3": 0, "HAP1": 14, "HT29": 6, "K562": 8, "MCF7": 0}
    res = {
        "bloque": "B1", "lineas": C.LINES, "panel": panel, "n_genes": n_genes,
        "p_adj_threshold": T, "semillas": C.comprueba_semillas(),
        "verif": verif,
        "cache_identicas_total": sum(v["n_perts_identicas_al_cache"]
                                     for v in verif.values()),
        "ed": ed, "ed_igualado": ed_ig, "diag_igualacion": dig, "dig": dig,
        "dist": {"A": A, "A_total": At, "B": B, "B_total": Bt, "por_pert": por_pert},
        "dist_ig": {"A": A2, "A_total": At2, "B": B2, "B_total": Bt2, "por_pert": pp2},
        "barr": barr, "barr_ig": barr_ig,
        "n_vacios_por_linea": vac,
        "vacios_vs_cp": ("coincide en las 6" if vac == esperado
                         else f"🔴 difiere: medido {vac}, cp §7 {esperado}"),
        "t_inicio": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(t_ini)),
        "t_fin": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()),
        "dt_total_s": time.time() - t_ini,
        "parcial": False,
    }
    C.escribe_json(OUT_JSON, res)
    C.escribe_md(OUT_MD, md(res))
    C.borra(WORK)
    C.log(f"B1 listo: {OUT_MD}")


if __name__ == "__main__":
    main()
