#!/usr/bin/env python
"""Pone las tres vias lado a lado y emite las tablas. Solo lectura.

SUELO        : dispersion_u.json      (ya medido)
CANDIDATO 01 : candidato_01*.json     (ya medido)
GWPS         : gwps_via.json          (esta corrida)

El margen absoluto de tres terminos y la sd_s del suelo se LEEN de la tabla §4
de `dispersion_b.md`, que es nuestra propia salida, en vez de recalcularse.
"""
import json
import os
import re
import sys

import numpy as np

OUT = "/home/ubuntu/vcc/modelo/out"
LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
M = {
    "pds_cosine": "pds_cosine",
    "expr_mse_unbiased_capped_norm": "mse_capped_norm",
    "de_wilcoxon_lfc_nmae": "lfc_nmae",
    "de_wilcoxon_direction_fidelity_yield_raw": "fid",
    "de_wilcoxon_direction_reach_raw": "reach",
    "de_wilcoxon_sig_jaccard": "jac",
}
METS = list(M.keys())


def jl(p):
    return [json.loads(l) for l in open(p) if l.strip()]


def carga(paths):
    r = []
    for p in paths:
        if os.path.exists(p):
            r += jl(p)
    return [x for x in r if x.get("score_ok")]


def s_de(rec, met):
    f = (rec.get("scored") or {}).get("filas", {}).get(met)
    if not f:
        return None
    v = f.get("from_replicate")
    return None if v is None else float(v)


def avg_de(rec):
    f = (rec.get("scored") or {}).get("filas", {}).get("avg_score")
    return None if not f else f.get("from_replicate")


def resume(recs):
    """{(linea,met): (media, sd, n)} y {(linea,'avg'): ...}"""
    out = {}
    for ln in LINES:
        rs = [x for x in recs if x["linea"] == ln]
        for met in METS + ["avg_score"]:
            v = [(avg_de(x) if met == "avg_score" else s_de(x, met)) for x in rs]
            v = [q for q in v if q is not None]
            if v:
                out[(ln, met)] = (float(np.mean(v)), float(np.std(v, ddof=1))
                                  if len(v) > 1 else 0.0, len(v))
    return out


def margen_de_dispersion_b():
    """Lee la tabla §4 de dispersion_b.md: margen COMPLETO y sd_s por celda."""
    txt = open(os.path.join(OUT, "dispersion_b.md")).read()
    bloque = txt.split("## 4. El margen completo")[1].split("## 5")[0]
    marg, sds, span = {}, {}, {}
    for linea in bloque.splitlines():
        c = [q.strip() for q in linea.strip().strip("|").split("|")]
        if len(c) < 12 or c[0] not in LINES:
            continue
        met = c[1].strip("`")
        try:
            span[(c[0], met)] = float(c[3])
            sds[(c[0], met)] = float(c[5])
            marg[(c[0], met)] = float(c[10].replace("*", ""))
        except ValueError:
            continue
    return marg, sds, span


def main():
    suelo = carga([os.path.join(OUT, "dispersion_u.json")])
    cand = carga([os.path.join(OUT, f"candidato_01{s}.json") for s in ("", "_p2", "_p3")])
    gwps = carga([os.path.join(OUT, "gwps_via.json")])
    gp2 = carga([os.path.join(OUT, "gwps_via_p2.json")])
    vias = {"SUELO": resume(suelo), "CAND01": resume(cand), "GWPS": resume(gwps)}
    marg, sds_suelo, span = margen_de_dispersion_b()
    sal = {"n_registros": {"SUELO": len(suelo), "CAND01": len(cand),
                           "GWPS": len(gwps), "GWPS_P2": len(gp2)}}

    print("## n registros con score_ok")
    for k, v in sal["n_registros"].items():
        print("  %-8s %d" % (k, v))

    print()
    print("### T1 los seis s por linea y por via (media +- sd sobre 8 semillas)")
    print("| linea | metrica | **SUELO** | sd | **CAND01** | sd | **GWPS** | sd |")
    print("|---|---|---|---|---|---|---|---|")
    for ln in LINES:
        for met in METS:
            row = [ln, "`%s`" % M[met]]
            for v in ("SUELO", "CAND01", "GWPS"):
                t = vias[v].get((ln, met))
                row += ["**%+.4f**" % t[0], "%.4f" % t[1]] if t else ["n/d", "n/d"]
            print("| " + " | ".join(row) + " |")

    print()
    print("### T2 avg_score por linea y por via")
    print("| linea | **SUELO** | sd | **CAND01** | sd | **GWPS** | sd |")
    print("|---|---|---|---|---|---|---|")
    for ln in LINES:
        row = [ln]
        for v in ("SUELO", "CAND01", "GWPS"):
            t = vias[v].get((ln, "avg_score"))
            row += ["**%+.4f**" % t[0], "%.4f" % t[1]] if t else ["n/d", "n/d"]
        print("| " + " | ".join(row) + " |")

    print()
    print("### T3 🔴 EL RANGO por metrica sobre las seis lineas: [min, max]")
    print("| metrica | **SUELO [min, max]** | **CAND01 [min, max]** | **GWPS [min, max]** |")
    print("|---|---|---|---|")
    for met in METS + ["avg_score"]:
        row = ["`%s`" % M.get(met, met)]
        for v in ("SUELO", "CAND01", "GWPS"):
            q = [vias[v][(ln, met)][0] for ln in LINES if (ln, met) in vias[v]]
            row.append("**[%+.4f, %+.4f]**" % (min(q), max(q)) if q else "n/d")
        print("| " + " | ".join(row) + " |")

    print()
    print("### T4 ¿supera cada via al b? (s > 0), con el margen absoluto de 3 terminos")
    print("| metrica | SUELO s>0 | SUELO s>margen | CAND01 s>0 | CAND01 s>margen | GWPS s>0 | GWPS s>margen |")
    print("|---|---|---|---|---|---|---|")
    for met in METS:
        row = ["`%s`" % M[met]]
        for v in ("SUELO", "CAND01", "GWPS"):
            pos = sum(1 for ln in LINES
                      if (ln, met) in vias[v] and vias[v][(ln, met)][0] > 0)
            sup = sum(1 for ln in LINES
                      if (ln, met) in vias[v]
                      and vias[v][(ln, met)][0] > marg.get((ln, M[met]), np.inf))
            row += ["%d/6" % pos, "%d/6" % sup]
        print("| " + " | ".join(row) + " |")

    print()
    print("### T5 ¿supera cada via al SUELO? resolucion_8 = sqrt((sd_a^2+sd_b^2)/8)")
    print("| metrica | CAND01 > SUELO | de ellas > resolucion_8 | GWPS > SUELO | de ellas > resolucion_8 |")
    print("|---|---|---|---|---|")
    det = []
    for met in METS:
        row = ["`%s`" % M[met]]
        for v in ("CAND01", "GWPS"):
            gt = res = 0
            for ln in LINES:
                a = vias[v].get((ln, met)); b = vias["SUELO"].get((ln, met))
                if not a or not b:
                    continue
                d = a[0] - b[0]
                r8 = np.sqrt((a[1] ** 2 + b[1] ** 2) / 8.0)
                if d > 0:
                    gt += 1
                    if d > r8:
                        res += 1
                det.append({"via": v, "linea": ln, "metrica": M[met],
                            "dif": d, "res8": float(r8)})
            row += ["%d/6" % gt, "%d/6" % res]
        print("| " + " | ".join(row) + " |")

    print()
    print("### T6 GWPS - SUELO celda por celda, con su resolucion_8")
    print("| linea | metrica | SUELO | GWPS | **GWPS - SUELO** | resolucion_8 | \\|dif\\|>res | margen abs. §4 |")
    print("|---|---|---|---|---|---|---|---|")
    for ln in LINES:
        for met in METS:
            a = vias["GWPS"].get((ln, met)); b = vias["SUELO"].get((ln, met))
            if not a or not b:
                continue
            d = a[0] - b[0]
            r8 = np.sqrt((a[1] ** 2 + b[1] ** 2) / 8.0)
            print("| %s | `%s` | %+.4f | %+.4f | **%+.4f** | %.4f | %s | %.4f |" % (
                ln, M[met], b[0], a[0], d, r8, "si" if abs(d) > r8 else "no",
                marg.get((ln, M[met]), float("nan"))))

    json.dump(sal, open("/tmp/claude-1000/-home-ubuntu-vcc/175bacfc-735d-42f8-a4ad-2ed66a1f8970/scratchpad/tv.json", "w"))


if __name__ == "__main__":
    main()
