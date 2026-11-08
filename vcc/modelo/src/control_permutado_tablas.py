#!/usr/bin/env python
"""Tablas del control permutado. Solo lectura."""
import json
import os
import numpy as np

OUT = "/home/ubuntu/vcc/modelo/out"
L = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
METS = ["pds_cosine", "expr_mse_unbiased_capped_norm", "de_wilcoxon_lfc_nmae",
        "de_wilcoxon_direction_fidelity_yield_raw",
        "de_wilcoxon_direction_reach_raw", "de_wilcoxon_sig_jaccard"]
SH = dict(zip(METS, ["pds_cosine", "mse_capped_norm", "lfc_nmae", "fid",
                     "reach", "jac"]))


def jl(p):
    return [json.loads(x) for x in open(p) if x.strip()]


def s(x, m):
    f = (x.get("scored") or {}).get("filas", {}).get(m)
    return None if not f else f.get("from_replicate")


def resume(recs):
    o = {}
    for ln in L:
        rs = [x for x in recs if x["linea"] == ln and x.get("score_ok")]
        for m in METS + ["avg_score"]:
            v = [s(x, m) for x in rs]
            v = [q for q in v if q is not None]
            if v:
                o[(ln, m)] = (float(np.mean(v)),
                              float(np.std(v, ddof=1)) if len(v) > 1 else 0.0,
                              len(v))
    return o


def margen():
    txt = open(os.path.join(OUT, "dispersion_b.md")).read()
    bl = txt.split("## 4. El margen completo")[1].split("## 5")[0]
    mg = {}
    for linea in bl.splitlines():
        c = [q.strip() for q in linea.strip().strip("|").split("|")]
        if len(c) < 12 or c[0] not in L:
            continue
        try:
            mg[(c[0], c[1].strip("`"))] = float(c[10].replace("*", ""))
        except ValueError:
            pass
    return mg


def venn_de(por_pert):
    """4 celdas + jaccard genuino de la celda ambos-llenos."""
    av = r0p1 = r1p0 = al = 0
    jac = []
    for p, c in por_pert.items():
        rv = c["n_DE_real"] == 0
        pv = c["n_DE_pred"] == 0
        if rv and pv:
            av += 1
        elif rv and not pv:
            r0p1 += 1
        elif not rv and pv:
            r1p0 += 1
        else:
            al += 1
            if c["n_DE_union"] > 0:
                jac.append(c["n_DE_interseccion"] / c["n_DE_union"])
    return av, r0p1, r1p0, al, (float(np.mean(jac)) if jac else None)


def main():
    perm = jl(os.path.join(OUT, "control_permutado.json"))
    real = jl(os.path.join(OUT, "gwps_via.json"))
    suelo = jl(os.path.join(OUT, "dispersion_u.json"))
    P, R, S = resume(perm), resume(real), resume(suelo)
    mg = margen()

    print("## verificaciones del brazo PERMUTADO")
    print("  score_ok %d/48 ; config_hash %d/48 ; poblacion_real_identica %d/48"
          % (sum(1 for x in perm if x["score_ok"]),
             sum(1 for x in perm if x["hash_run"]["coincide"]),
             sum(1 for x in perm if x["poblacion_real_identica"])))
    print("  huellas pred distintas: %d/48 ; todo_entero %d/48 ; negativos %d"
          % (len(set(x["huella_pred"] for x in perm)),
             sum(1 for x in perm if x["contrato_emision"]["todo_entero"]),
             sum(1 for x in perm if x["contrato_emision"]["hay_negativos"])))
    print("  huellas reales por linea: %s"
          % {ln: len(set(x["huella_real"] for x in perm if x["linea"] == ln))
             for ln in L})
    sinfijo = all(x["desarreglo"]["sin_punto_fijo"] for x in perm)
    print("  desarreglos sin punto fijo: %s (48/48)" % sinfijo)
    mismos = {}
    for j in range(8):
        sig = set(json.dumps(x["desarreglo"]["sigma"], sort_keys=True)
                  for x in perm if x["j"] == j)
        mismos[j] = len(sig)
    print("  sigma identico entre lineas para cada j: %s" % mismos)
    print("  huella_pred compartida con el brazo REAL: %d"
          % len(set(x["huella_pred"] for x in perm)
                & set(x["huella_pred"] for x in real)))

    print()
    print("### T1 los seis s, brazo PERMUTADO (media +- sd sobre 8 desarreglos)")
    print("| linea | " + " | ".join("`%s`" % SH[m] for m in METS) + " | `avg` |")
    print("|---|" + "---|" * 7)
    for ln in L:
        row = [ln]
        for m in METS + ["avg_score"]:
            t = P.get((ln, m))
            row.append("**%+.4f** ± %.4f" % (t[0], t[1]) if t else "n/d")
        print("| " + " | ".join(row) + " |")

    print()
    print("### T2 los tres brazos lado a lado")
    print("| linea | metrica | SUELO | sd | REAL | sd | PERMUTADO | sd |")
    print("|---|---|---|---|---|---|---|---|")
    for ln in L:
        for m in METS:
            row = [ln, "`%s`" % SH[m]]
            for D in (S, R, P):
                t = D.get((ln, m))
                row += ["%+.4f" % t[0], "%.4f" % t[1]] if t else ["n/d", "n/d"]
            print("| " + " | ".join(row) + " |")
    for ln in L:
        row = [ln, "`avg_score`"]
        for D in (S, R, P):
            t = D.get((ln, "avg_score"))
            row += ["%+.4f" % t[0], "%.4f" % t[1]] if t else ["n/d", "n/d"]
        print("| " + " | ".join(row) + " |")

    print()
    print("### T3 REAL - PERMUTADO celda por celda")
    print("| linea | metrica | REAL | PERMUTADO | **REAL - PERM** | resolucion_8 | \\|dif\\|>res |")
    print("|---|---|---|---|---|---|---|")
    sup = {m: 0 for m in METS}
    for ln in L:
        for m in METS:
            a, b = R.get((ln, m)), P.get((ln, m))
            if not a or not b:
                continue
            d = a[0] - b[0]
            r8 = np.sqrt((a[1] ** 2 + b[1] ** 2) / 8.0)
            ok = abs(d) > r8
            if ok and d > 0:
                sup[m] += 1
            print("| %s | `%s` | %+.4f | %+.4f | **%+.4f** | %.4f | %s |"
                  % (ln, SH[m], a[0], b[0], d, r8, "si" if ok else "no"))

    print()
    print("### T4 los conteos que contestan la pregunta")
    print("| metrica | REAL > PERM mas alla de resolucion_8 | REAL s>0 | PERM s>0 | REAL s>margen | **PERM s>margen** |")
    print("|---|---|---|---|---|---|")
    for m in METS:
        rp = sum(1 for ln in L if (ln, m) in R and (ln, m) in P
                 and R[(ln, m)][0] > 0)
        pp = sum(1 for ln in L if (ln, m) in P and P[(ln, m)][0] > 0)
        rm = sum(1 for ln in L if (ln, m) in R
                 and R[(ln, m)][0] > mg.get((ln, SH[m]), np.inf))
        pm = sum(1 for ln in L if (ln, m) in P
                 and P[(ln, m)][0] > mg.get((ln, SH[m]), np.inf))
        print("| `%s` | **%d/6** | %d/6 | %d/6 | %d/6 | **%d/6** |"
              % (SH[m], sup[m], rp, pp, rm, pm))

    print()
    print("### T5 Venn 2x2 por brazo y por linea")
    print("| linea | brazo | ambos vacios | real ∅ / pred lleno (FP) | real lleno / pred ∅ | ambos llenos | jaccard genuino medio |")
    print("|---|---|---|---|---|---|---|")
    Rvac = {}
    for ln in L:
        pp = [x for x in real if x["linea"] == ln][0]["cohortes"]["por_pert"]
        Rvac[ln] = sum(1 for c in pp.values() if c["n_DE_real"] == 0)
    for ln in L:
        for nom, recs in (("SUELO", suelo), ("REAL", real), ("PERMUTADO", perm)):
            rs = [x for x in recs if x["linea"] == ln]
            if nom == "SUELO":
                # derivado exacto de R, P y E_union_vacia (sin por_pert guardado)
                cel = []
                for x in rs:
                    E = x["cohortes"]["E_union_vacia"]
                    Pn = x["cohortes"]["n_perts_con_pred_sig"]
                    r = Rvac[ln]
                    cel.append((E, r - E, 33 - Pn - E, Pn - (r - E), None))
                a = np.mean([c[0] for c in cel]); b = np.mean([c[1] for c in cel])
                c_ = np.mean([c[2] for c in cel]); d_ = np.mean([c[3] for c in cel])
                jm = "n/d (sin `por_pert`)"
            else:
                cel = [venn_de(x["cohortes"]["por_pert"]) for x in rs]
                a = np.mean([c[0] for c in cel]); b = np.mean([c[1] for c in cel])
                c_ = np.mean([c[2] for c in cel]); d_ = np.mean([c[3] for c in cel])
                js = [c[4] for c in cel if c[4] is not None]
                jm = "%.4f" % np.mean(js) if js else "n/d"
            print("| %s | %s | %.2f | %.2f | %.2f | %.2f | %s |"
                  % (ln, nom, a, b, c_, d_, jm))

    print()
    print("### T6 perturbaciones con DE REAL vacio, MEDIDO, por linea")
    print("| linea | **# DE real vacio** | de 33 | identico en las 8 semillas | E_union_vacia SUELO (otra cantidad) |")
    print("|---|---|---|---|---|")
    for ln in L:
        sets = [frozenset(p for p, c in x["cohortes"]["por_pert"].items()
                          if c["n_DE_real"] == 0)
                for x in real if x["linea"] == ln]
        ev = sorted(set(x["cohortes"]["E_union_vacia"]
                        for x in suelo if x["linea"] == ln))
        print("| %s | **%d** | 33 | %s | %s |"
              % (ln, Rvac[ln], "si" if len(set(sets)) == 1 else "NO", ev))

    print()
    print("### T7 el rango [min, max] por metrica y por brazo")
    print("| metrica | SUELO | REAL | **PERMUTADO** |")
    print("|---|---|---|---|")
    for m in METS + ["avg_score"]:
        row = ["`%s`" % SH.get(m, "avg_score")]
        for D in (S, R, P):
            q = [D[(ln, m)][0] for ln in L if (ln, m) in D]
            row.append("[%+.4f, %+.4f]" % (min(q), max(q)) if q else "n/d")
        print("| " + " | ".join(row) + " |")

    print()
    print("### extra: DEpred por brazo (igualacion de potencia)")
    print("| linea | REAL DEpred | PERM DEpred | dif | REAL n_fid | PERM n_fid |")
    print("|---|---|---|---|---|---|")
    for ln in L:
        a = np.mean([x["cohortes"]["suma_n_DE_pred"] for x in real if x["linea"] == ln])
        b = np.mean([x["cohortes"]["suma_n_DE_pred"] for x in perm if x["linea"] == ln])
        fa = np.mean([x["cohortes"]["n_fid"] for x in real if x["linea"] == ln])
        fb = np.mean([x["cohortes"]["n_fid"] for x in perm if x["linea"] == ln])
        print("| %s | %.0f | %.0f | %+.0f | %.2f | %.2f |" % (ln, a, b, b - a, fa, fb))


if __name__ == "__main__":
    main()
