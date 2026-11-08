#!/usr/bin/env python
"""PARTE C: el Venn literal de los conjuntos ED entre contextos, y el TECHO de
`jac` para una emision context-independiente.

🔴 TRANSDUCTIVO: lee el lado real perturbado de Jiang. Estante del `b`.

LAS CUATRO CORRECCIONES, aplicadas y declaradas
-----------------------------------------------
1. 🔴 **LA INTERSECCION NO ES EL TECHO.** El agregado es la MEDIA de los `jac`
   por contexto, y a un `S` le puede convenir cubrir bien un contexto y mal
   otro. Contraejemplo: `R_D={a,b}`, `R_E={a,c}` -> `S={a}` da media 0.5000 y
   `S={a,b}` da 0.6667. Asi que el techo se **MIDE** barriendo el corte de
   frecuencia: se ordenan los genes por en cuantas lineas se llaman y se corta
   en `>=6, >=5, ..., >=1`, dando la media de `jac` en cada corte. Se reporta la
   CURVA ENTERA y su maximo, declarado como **COTA INFERIOR del techo**
   (restringida a esa familia de `S`), **nunca como el techo exacto**.
2. 🔴 **EL CONTEO CRUDO ESTA CONFUNDIDO CON POTENCIA.** Las medianas de ED van
   de 1 (A549) a 197 (MCF7): el solape lo acota el conjunto mas chico. Se dan
   las dos versiones y se declaran DISTINTAS:
     (a) el conteo CRUDO
     (b) IGUALADO EN POTENCIA: las seis lineas recortadas al mismo tamano por
         perturbacion -- el del conjunto mas pequeno NO VACIO -- tomando los de
         `p` mas baja. Empates por semilla declarada.
   🔴 La cifra que habla de CONTEXTO es la DIFERENCIA entre (a) y (b).
3. El encabezado declara que esto es una cota **DE JIANG**. La heterogeneidad
   que produce el conteo esta MEDIDA COMO AUSENTE en el reto (`ntc_nulo.md` §3).
   **No se extrapola a D/E/F en ningun sentido.**
4. Se reporta **por linea y por perturbacion como RANGO**, nunca como media
   suelta.

LA METRICA, la del marcador y no otra (`cell_eval2/metrics/de.py:180`):
  `jac(S, R_L) = |S ∩ R_L| / |S ∪ R_L|` sobre pares (blanco, gen) UNICOS,
  con el **propio gen blanco fuera de los dos lados** (#172) y
  🔴 **union vacia -> 1.0** (convencion documentada). Se dan las dos lecturas:
  con la convencion, y **solo genuinos** (restringido a las lineas con union no
  vacia), como hizo `control_permutado.md` §6.
"""
import json, os, sys
from collections import Counter
import numpy as np
import polars as pl

sys.path.insert(0, "/home/ubuntu/vcc/modelo/src")
import gwps_via as GV

LINES = GV.LINES
PANEL_PATH = GV.PANEL_PATH
WORK = "/datos/modelo/out/posicion_mezcla"
OUTDIR = "/home/ubuntu/vcc/modelo/out/posicion_mezcla"
SEED_BASE = 20279000              # 20279000 + 100*fold + j, j = indice de pert
ALPHA = 0.05
FIABLES = ["IFNAR1", "IFNAR2", "JAK1", "STAT1", "TYK2"]   # de magnitud_centro.md §2.1
CORTES = [6, 5, 4, 3, 2, 1]


def jac(S, R):
    u = len(S | R)
    return (1.0, True) if u == 0 else (len(S & R) / u, False)


def main():
    panel = sorted(map(str, json.load(open(PANEL_PATH))["panel_comun"]))
    # tablas DE recomputadas: p por (blanco, gen)
    ED, PV = {}, {}
    for L in LINES:
        t = pl.read_parquet(os.path.join(WORK, f"de_{L}.parquet"))
        sig = t.filter((pl.col("p_adj") < ALPHA)
                       & (pl.col("target") != pl.col("feature")))
        ED[L] = {p: set() for p in panel}
        PV[L] = {p: {} for p in panel}
        for tgt, feat, pv in sig.select("target", "feature", "p_value").iter_rows():
            if tgt in ED[L]:
                ED[L][tgt].add(feat); PV[L][tgt][feat] = pv
    sys.stderr.write("cargado: " + ", ".join(
        f"{L}={sum(len(v) for v in ED[L].values())}" for L in LINES) + "\n")

    def recorta(ed, pv):
        """(b) igualado en potencia: cada linea al tamano del mas pequeno NO
        vacio, tomando los de `p` mas baja. Empates por semilla declarada."""
        out = {L: {} for L in LINES}
        info = {}
        for ip, p in enumerate(panel):
            tam = [len(ed[L][p]) for L in LINES]
            nz = [x for x in tam if x > 0]
            m = min(nz) if nz else 0
            info[p] = {"tamanos": dict(zip(LINES, tam)), "m": m,
                       "n_lineas_vacias": sum(1 for x in tam if x == 0),
                       "todas_vacias": not nz}
            for fold, L in enumerate(LINES):
                g = sorted(ed[L][p])
                if not g or m == 0:
                    out[L][p] = set(); continue
                sem = SEED_BASE + 100 * fold + ip
                rng = np.random.default_rng(sem)
                jitter = rng.permutation(len(g))
                orden = sorted(range(len(g)),
                               key=lambda i: (pv[L][p][g[i]], jitter[i]))
                out[L][p] = {g[i] for i in orden[:m]}
                info[p].setdefault("semillas", {})[L] = sem
        return out, info

    EDb, info_b = recorta(ED, PV)

    res = {"metrica": "de_sig_jaccard del preset: |S∩R|/|S∪R|, gen blanco fuera "
                      "de los dos lados (#172), union vacia -> 1.0",
           "alpha": ALPHA, "cortes": CORTES, "semilla_base_empates": SEED_BASE,
           "fiables": FIABLES, "panel": panel, "versiones": {}}

    for ver, ed in (("a_crudo", ED), ("b_igualado_en_potencia", EDb)):
        filas = []
        for p in panel:
            Rs = {L: ed[L][p] for L in LINES}
            cnt = Counter()
            for L in LINES:
                for g in Rs[L]:
                    cnt[g] += 1
            mult = {k: sum(1 for g, n in cnt.items() if n == k) for k in range(1, 7)}
            inter6 = {g for g, n in cnt.items() if n == 6}
            union = set(cnt)
            # el barrido del corte de frecuencia
            curva = []
            for c in CORTES:
                S = {g for g, n in cnt.items() if n >= c}
                js, vac = [], []
                for L in LINES:
                    v, e = jac(S, Rs[L]); js.append(v); vac.append(e)
                gen = [v for v, e in zip(js, vac) if not e]
                curva.append({
                    "corte": c, "n_S": len(S),
                    "jac_media_con_convencion": float(np.mean(js)),
                    "jac_por_linea": {L: float(x) for L, x in zip(LINES, js)},
                    "n_lineas_union_vacia": int(sum(vac)),
                    "jac_media_solo_genuinos": (float(np.mean(gen)) if gen else None),
                    "n_lineas_genuinas": len(gen),
                })
            best = max(curva, key=lambda d: d["jac_media_con_convencion"])
            bestg = max([d for d in curva if d["jac_media_solo_genuinos"] is not None],
                        key=lambda d: d["jac_media_solo_genuinos"], default=None)
            filas.append({
                "pert": p, "fiable": p in FIABLES,
                "tamanos_ED": {L: len(Rs[L]) for L in LINES},
                "n_lineas_vacias": sum(1 for L in LINES if not Rs[L]),
                "union": len(union), "interseccion_6": len(inter6),
                "multiplicidad": mult,
                "curva": curva,
                "corte_optimo": best["corte"],
                "techo_cota_inferior_con_convencion": best["jac_media_con_convencion"],
                "n_S_en_el_optimo": best["n_S"],
                "corte_optimo_genuinos": bestg["corte"] if bestg else None,
                "techo_cota_inferior_solo_genuinos": (
                    bestg["jac_media_solo_genuinos"] if bestg else None),
                "jac_en_corte_6": curva[0]["jac_media_con_convencion"],
                "jac_en_corte_1": curva[-1]["jac_media_con_convencion"],
            })
        # agregados
        def rango(k, sel=lambda r: True):
            v = [r[k] for r in filas if sel(r) and r[k] is not None]
            return {"min": float(min(v)), "max": float(max(v)),
                    "media": float(np.mean(v)), "mediana": float(np.median(v)),
                    "n": len(v)} if v else None
        mult_agg = {k: sum(r["multiplicidad"][k] for r in filas) for k in range(1, 7)}
        tot = sum(mult_agg.values())
        res["versiones"][ver] = {
            "multiplicidad_agregada": mult_agg,
            "total_llamadas_gen_pert": tot,
            "frac_en_una_sola_linea": (mult_agg[1] / tot) if tot else None,
            "frac_en_las_seis": (mult_agg[6] / tot) if tot else None,
            "interseccion_6_total": sum(r["interseccion_6"] for r in filas),
            "n_perts_con_interseccion_6_vacia": sum(
                1 for r in filas if r["interseccion_6"] == 0),
            "techo_rango": rango("techo_cota_inferior_con_convencion"),
            "techo_rango_genuinos": rango("techo_cota_inferior_solo_genuinos"),
            "techo_rango_fiables": rango("techo_cota_inferior_con_convencion",
                                         lambda r: r["fiable"]),
            "techo_rango_resto": rango("techo_cota_inferior_con_convencion",
                                       lambda r: not r["fiable"]),
            "techo_genuinos_fiables": rango("techo_cota_inferior_solo_genuinos",
                                            lambda r: r["fiable"]),
            "techo_genuinos_resto": rango("techo_cota_inferior_solo_genuinos",
                                          lambda r: not r["fiable"]),
            "corte_optimo_conteo": dict(Counter(r["corte_optimo"] for r in filas)),
            # la curva agregada: media sobre perturbaciones en cada corte
            "curva_agregada": [{
                "corte": c,
                "jac_media_sobre_perts_con_convencion": float(np.mean(
                    [r["curva"][i]["jac_media_con_convencion"] for r in filas])),
                "jac_rango_sobre_perts": [
                    float(min(r["curva"][i]["jac_media_con_convencion"] for r in filas)),
                    float(max(r["curva"][i]["jac_media_con_convencion"] for r in filas))],
                "n_S_medio": float(np.mean([r["curva"][i]["n_S"] for r in filas])),
                "jac_media_solo_genuinos": float(np.mean(
                    [r["curva"][i]["jac_media_solo_genuinos"] for r in filas
                     if r["curva"][i]["jac_media_solo_genuinos"] is not None])),
            } for i, c in enumerate(CORTES)],
            # por LINEA: el jac en el corte optimo de cada pert
            "por_linea_en_el_optimo": {L: {
                "min": float(min(r["curva"][CORTES.index(r["corte_optimo"])]
                                 ["jac_por_linea"][L] for r in filas)),
                "max": float(max(r["curva"][CORTES.index(r["corte_optimo"])]
                                 ["jac_por_linea"][L] for r in filas)),
                "media": float(np.mean([r["curva"][CORTES.index(r["corte_optimo"])]
                                        ["jac_por_linea"][L] for r in filas])),
            } for L in LINES},
            "tamanos_ED_por_linea": {L: {
                "min": min(r["tamanos_ED"][L] for r in filas),
                "max": max(r["tamanos_ED"][L] for r in filas),
                "mediana": float(np.median([r["tamanos_ED"][L] for r in filas])),
                "n_vacios": sum(1 for r in filas if r["tamanos_ED"][L] == 0),
            } for L in LINES},
            "por_pert": filas,
        }
        a = res["versiones"][ver]
        sys.stderr.write(
            f"  [{ver}] mult={mult_agg} frac_1_linea={a['frac_en_una_sola_linea']:.4f} "
            f"techo=[{a['techo_rango']['min']:.4f}, {a['techo_rango']['max']:.4f}] "
            f"media {a['techo_rango']['media']:.4f}\n")

    res["igualado_info"] = info_b
    # la cifra que habla de CONTEXTO: (a) - (b)
    A, B = res["versiones"]["a_crudo"], res["versiones"]["b_igualado_en_potencia"]
    res["diferencia_a_menos_b"] = {
        "frac_en_una_sola_linea": A["frac_en_una_sola_linea"] - B["frac_en_una_sola_linea"],
        "frac_en_las_seis": A["frac_en_las_seis"] - B["frac_en_las_seis"],
        "techo_media": A["techo_rango"]["media"] - B["techo_rango"]["media"],
        "techo_genuinos_media": (A["techo_rango_genuinos"]["media"]
                                 - B["techo_rango_genuinos"]["media"]),
        "interseccion_6_total": A["interseccion_6_total"] - B["interseccion_6_total"],
    }
    json.dump(res, open(os.path.join(OUTDIR, "venn.json"), "w"),
              indent=1, sort_keys=True, default=str)
    sys.stderr.write(f"\ndiferencia (a)-(b): "
                     f"{json.dumps(res['diferencia_a_menos_b'], indent=1)}\n")
    sys.stderr.write(f"listo: {OUTDIR}/venn.json\n")


if __name__ == "__main__":
    main()
