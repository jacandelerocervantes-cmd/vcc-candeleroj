#!/usr/bin/env python
"""PARTE 2a y 2b, calculadas SOBRE EL NULO de la Parte 1.

2a — TASA DE FALSOS POSITIVOS DEL CONJUNTO ED EN CONTEXTO.
   Dos cantidades distintas, que no se mezclan:
     * `n_sig` por corrida: cuantos genes pasa la puerta una comparacion
       control-contra-control. Es la cantidad que se compara con el `n_DE` de
       una perturbacion.
     * FPR POR GEN: (suma de n_sig) / (suma de genes probados). Es la tasa por
       prueba, y por construccion de BH se espera MUY por debajo de 0.05 --
       0.05 es el nivel de FDR del conjunto llamado, no la tasa por prueba.

2b — PISO DE DETECCION, calculado del propio nulo.
   La puerta llama un gen cuando `p_adj < alpha`. Para una senal AISLADA (el
   resto del panel nulo), BH en su escalon de rango 1 exige
   `p <= alpha/m`, con `m` = genes probados. Ese `p_crit` se calcula de cada
   corrida, no se importa.
   El nulo da la relacion |log2FC| <-> p MEDIDA, por banda de expresion (el
   piso depende fuerte del nivel del gen). Se reporta:
     * `|lfc|` minimo observado en el nulo para cada umbral de p alcanzable;
     * el techo de ruido: `|lfc|` maximo que el nulo produce en esa banda;
     * la inversion por tramos de |lfc| (minimo p por tramo), que es la lectura
       directa del piso cuando `p_crit` cae dentro del alcance del nulo;
     * 🔴 el ALCANCE: `m_filas * p_crit` dice cuantas filas del nulo se esperan
       por debajo de `p_crit`. Si es << 1, el piso en `p_crit` NO se puede leer
       directo y se da como RANGO con su extrapolacion DECLARADA.
"""
import glob, json, os, sys
import numpy as np

OUTDIR = "/home/ubuntu/vcc/modelo/out/ntc_nulo"
OUT = os.path.join(OUTDIR, "p2ab.json")
ALPHA = 0.05
BANDAS = [(5, 10), (10, 20), (20, 50), (50, 100), (100, 300), (300, 1000),
          (1000, float("inf"))]
UMBRALES_P = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]

p1 = json.load(open(os.path.join(OUTDIR, "p1.json")))
ALPHA = p1["puerta"]["de.p_adj_threshold"]
res = {"alpha": ALPHA, "bandas_cpm": [list(b) for b in BANDAS], "brazos": {}}

# ------------------------------------------------------------------ 2a
por_brazo = {}
for r in p1["corridas"]:
    por_brazo.setdefault(r["brazo"], []).append(r)

for brazo, rs in sorted(por_brazo.items()):
    ns = np.array([r["n_sig"] for r in rs])
    m = np.array([r["n_probados"] for r in rs])
    pmin = np.array([r["p_min"] for r in rs], dtype=float)
    porctx = {}
    for r in rs:
        porctx.setdefault(r["contexto"], []).append(r["n_sig"])
    genes = {}
    for r in rs:
        for g in (r["genes_sig"] or []):
            genes[g] = genes.get(g, 0) + 1
    res["brazos"][brazo] = {
        "n_corridas": len(rs),
        "n_ref": rs[0]["n_ref"], "n_tgt": rs[0]["n_tgt"],
        "n_sig_media": float(ns.mean()), "n_sig_sd": float(ns.std(ddof=1)),
        "n_sig_mediana": float(np.median(ns)),
        "n_sig_min": int(ns.min()), "n_sig_max": int(ns.max()),
        "n_corridas_con_alguno": int((ns > 0).sum()),
        "n_sig_total": int(ns.sum()),
        "m_min": int(m.min()), "m_max": int(m.max()), "m_media": float(m.mean()),
        "p_crit_rango": [ALPHA / int(m.max()), ALPHA / int(m.min())],
        "p_min_sobre_corridas_min": float(np.nanmin(pmin)),
        "p_min_sobre_corridas_mediana": float(np.nanmedian(pmin)),
        "fpr_por_gen": float(ns.sum() / m.sum()),
        "n_sig_por_contexto": {c: {"valores": v, "media": float(np.mean(v)),
                                   "sd": float(np.std(v, ddof=1)),
                                   "min": int(min(v)), "max": int(max(v))}
                               for c, v in sorted(porctx.items())},
        "genes_llamados_alguna_vez": len(genes),
        "genes_llamados_conteo": dict(sorted(genes.items(), key=lambda kv: -kv[1])),
    }

# ------------------------------------------------------------------ 2b
def pool(patron):
    P, L, C = [], [], []
    ms = []
    for f in sorted(glob.glob(os.path.join(OUTDIR, patron))):
        d = np.load(f)
        P.append(d["p_value"]); L.append(np.abs(d["log2fc"])); C.append(d["ref_cpm"])
        ms.append(d["p_value"].size)
    if not P:
        return None
    P = np.concatenate(P); L = np.concatenate(L); C = np.concatenate(C)
    ok = np.isfinite(P) & np.isfinite(L) & np.isfinite(C)
    return P[ok], L[ok], C[ok], np.array(ms)


def piso(P, L, C, m_media):
    p_crit = ALPHA / m_media
    out = {"p_crit": p_crit, "m_media": float(m_media),
           "n_filas_pool": int(P.size),
           "p_min_pool": float(P.min()),
           "filas_esperadas_bajo_p_crit": float(P.size * p_crit),
           "filas_observadas_bajo_p_crit": int((P <= p_crit).sum()),
           "p_crit_dentro_del_alcance": bool((P <= p_crit).sum() > 0),
           "bandas": []}
    if (P <= p_crit).sum() > 0:
        s = P <= p_crit
        out["lfc_de_las_filas_bajo_p_crit"] = sorted(
            [round(float(x), 4) for x in L[s]])
        out["cpm_de_las_filas_bajo_p_crit"] = sorted(
            [round(float(x), 2) for x in C[s]])
    for lo, hi in BANDAS:
        mb = (C > lo) & (C <= hi)
        if mb.sum() < 50:
            out["bandas"].append({"banda": [lo, hi if hi != float("inf") else None],
                                  "n_filas": int(mb.sum()), "insuficiente": True})
            continue
        pp, ll = P[mb], L[mb]
        d = {"banda": [lo, hi if hi != float("inf") else None],
             "n_filas": int(mb.sum()),
             "lfc_max_del_nulo": float(ll.max()),
             "lfc_p99_del_nulo": float(np.percentile(ll, 99)),
             "lfc_mediana_del_nulo": float(np.median(ll)),
             "p_min_de_la_banda": float(pp.min()),
             "filas_esperadas_bajo_p_crit": float(mb.sum() * p_crit),
             "filas_observadas_bajo_p_crit": int((pp <= p_crit).sum()),
             "anclas": {}}
        for t in UMBRALES_P + [p_crit]:
            s = pp <= t
            k = "p_crit" if t == p_crit else f"{t:.0e}"
            d["anclas"][k] = {
                "n": int(s.sum()),
                "lfc_min": float(ll[s].min()) if s.sum() else None,
                "lfc_mediana": float(np.median(ll[s])) if s.sum() else None,
            }
        # inversion por tramos de |lfc|: minimo p observado en cada tramo
        bordes = np.unique(np.quantile(ll, np.linspace(0, 1, 41)))
        tramos = []
        for a, b in zip(bordes[:-1], bordes[1:]):
            sm = (ll >= a) & (ll < b)
            if sm.sum() >= 20:
                tramos.append({"lfc_lo": float(a), "lfc_hi": float(b),
                               "n": int(sm.sum()), "p_min": float(pp[sm].min()),
                               "p_mediana": float(np.median(pp[sm]))})
        d["tramos_lfc"] = tramos
        cruza = [t for t in tramos if t["p_min"] <= p_crit]
        d["primer_tramo_que_cruza_p_crit"] = cruza[0] if cruza else None

        # COTAS MEDIDAS del piso F(t) = min{|lfc| : p <= t}, que es NO CRECIENTE en t.
        # p_crit < t  =>  F(p_crit) >= F(t): cada ancla alcanzada da una cota INFERIOR.
        # Una fila del nulo que SI cruzo p_crit da una cota SUPERIOR.
        alcanzadas = [(t, d["anclas"][("p_crit" if t == p_crit else f"{t:.0e}")])
                      for t in UMBRALES_P if t > p_crit
                      and d["anclas"][f"{t:.0e}"]["lfc_min"] is not None]
        d["cota_inferior_medida"] = None
        if alcanzadas:
            t_min, a_min = alcanzadas[-1]
            d["cota_inferior_medida"] = {
                "umbral_p_alcanzado": t_min, "n_filas_en_el_umbral": a_min["n"],
                "lfc_min": a_min["lfc_min"],
                "razon": f"p_crit={p_crit:.3e} < {t_min:.0e} y F es no creciente",
            }
        sb = pp <= p_crit
        d["cota_superior_medida"] = (
            {"lfc_min_de_las_filas_que_cruzaron": float(ll[sb].min()),
             "n_filas": int(sb.sum())} if sb.sum() else None)

        # Extrapolacion DECLARADA, sobre las ANCLAS (estadisticos de orden medidos),
        # no sobre los tramos: F(t) lineal en log10(t), leida en log10(p_crit).
        xs = np.array([np.log10(t) for t, _ in alcanzadas])
        ys = np.array([a["lfc_min"] for _, a in alcanzadas])
        if xs.size >= 3:
            A = np.polyfit(xs, ys, 1)
            if A[0] < 0:
                d["extrapolacion_declarada"] = {
                    "forma": "F(t)=min|lfc| con p<=t, lineal en log10(t), ajustada a "
                             "las anclas MEDIDAS de esta banda",
                    "anclas_usadas": [float(t) for t, _ in alcanzadas],
                    "lfc_de_las_anclas": [float(y) for y in ys],
                    "pendiente_por_decada": float(-A[0]),
                    "lfc_en_p_crit": float(np.polyval(A, np.log10(p_crit))),
                    "residuo_sd": float(np.std(ys - np.polyval(A, xs), ddof=1)),
                    "decadas_de_extrapolacion": float(xs.min() - np.log10(p_crit)),
                    "es_extrapolacion": True,
                }
        out["bandas"].append(d)
    return out


for brazo, patron in (("n400_400", "crudo_*_n400_400.npz"),
                      ("n400_ref_completa", "crudo_*_n400_ref.npz")):
    pl_ = pool(patron)
    if pl_ is None:
        continue
    P, L, C, ms = pl_
    mm = float(ms.mean())
    res["brazos"].setdefault(brazo, {})["piso_deteccion"] = piso(P, L, C, mm)
    res["brazos"][brazo]["n_archivos_pool"] = int(ms.size)
    sys.stderr.write(f"  {brazo}: pool {P.size} filas, p_min={P.min():.3e}, "
                     f"p_crit={ALPHA/mm:.3e}, bajo p_crit="
                     f"{int((P<=ALPHA/mm).sum())}\n")

json.dump(res, open(OUT, "w"), indent=1, sort_keys=True, default=str)
sys.stderr.write(f"\nlisto: {OUT}\n")
