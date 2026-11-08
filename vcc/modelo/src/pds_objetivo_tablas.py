#!/usr/bin/env python
"""Tablas de `modelo/out/pds_objetivo.md`. Solo lee y agrega; no corre nada.

🔴 Los tres oraculos y el brazo de signos aleatorios SON TRANSDUCTIVOS o
diagnosticos: ninguna de sus cifras puede entrar en una submission ni en un
rasgo de un candidato. Esto se imprime en el encabezado del reporte.
"""
import json
import math
import os
import statistics as st
import sys

import pandas as pd

LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
OUT = "/home/ubuntu/vcc/modelo/out"
AZAR = 0.5                     # punto sin informacion, DERIVADO en §1

# `b` y `r` oficiales, de `anclas_nativas_v2.md` (baseline.py / anchor.py)
B_OF = {"A549": 0.486742, "BXPC3": 0.476326, "HAP1": 0.437500,
        "HT29": 0.501894, "K562": 0.520833, "MCF7": 0.503788}
R_OF = {"A549": 0.777841, "BXPC3": 0.955114, "HAP1": 0.789583,
        "HT29": 0.783523, "K562": 0.822727, "MCF7": 0.936174}
R_SD5 = {"A549": 0.029037, "BXPC3": 0.016829, "HAP1": 0.015251,
         "HT29": 0.017083, "K562": 0.017558, "MCF7": 0.016343}

BRAZOS = [("copia", "ORACULO COMPLETO (copia)"),
          ("completo", "ORACULO COMPLETO (razon)"),
          ("cola", "ORACULO DE COLA"),
          ("bulto", "ORACULO DE BULTO"),
          ("signos", "SIGNOS ALEATORIOS DENSOS")]
YA = [("SUELO", "dispersion_u.json"), ("REAL", "gwps_via.json"),
      ("PERMUTADO", "control_permutado.json")]


def jl(p):
    return [json.loads(l) for l in open(os.path.join(OUT, p)) if l.strip()]


def mm(xs):
    return (st.mean(xs), st.stdev(xs) if len(xs) > 1 else 0.0, len(xs))


def fmt(v, n=4, sign=True):
    if v is None:
        return "n/d"
    return f"{v:+.{n}f}" if sign else f"{v:.{n}f}"


def carga():
    P = jl("pds_objetivo.json")
    P = [r for r in P if not r.get("fallo_en")]
    d = {}
    for b, _ in BRAZOS:
        for l in LINES:
            rs = [r for r in P if r["brazo"] == b and r["linea"] == l]
            if rs:
                d[(b, l)] = rs
    ya = {}
    for nm, f in YA:
        rs = jl(f)
        for l in LINES:
            ya[(nm, l)] = [r for r in rs if r["linea"] == l]
    return d, ya


def u_pds(rs):
    return [r["u"]["pds_cosine"] for r in rs]


def s_pds(rs):
    return [r["scored"]["filas"]["pds_cosine"]["from_replicate"] for r in rs]


# ----------------------------------------------------------------- PARTE A
def parte_a():
    U = jl("dispersion_u.json")
    B = []
    for f in ("dispersion_b.json", "dispersion_b_p2.json", "dispersion_b_p3.json"):
        B += jl(f)
    print("\n### TABLA A1: la escala de `pds_cosine` contra el azar\n")
    print("| linea | azar (derivado) | `b` oficial | `b` - azar | `r` oficial | "
          "`r` - azar | span `r-b` | `sd_b` (8 real.) | `sd_u` suelo (8 sem.) | "
          "`r_sd` (5 part.) | `resolucion_8` b vs u |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    fil = {}
    for l in LINES:
        us = [r["u"]["pds_cosine"] for r in U if r["linea"] == l]
        bs = [r["b"]["pds_cosine"] for r in B if r["linea"] == l]
        sdu, sdb = st.stdev(us), st.stdev(bs)
        res = math.sqrt((sdu ** 2 + sdb ** 2) / 8)
        fil[l] = dict(sdu=sdu, sdb=sdb, res=res, u=st.mean(us), b8=st.mean(bs))
        print(f"| {l} | {AZAR:.4f} | {B_OF[l]:.6f} | **{B_OF[l]-AZAR:+.6f}** | "
              f"{R_OF[l]:.6f} | **{R_OF[l]-AZAR:+.6f}** | {R_OF[l]-B_OF[l]:.6f} | "
              f"{sdb:.6f} | {sdu:.6f} | {R_SD5[l]:.6f} | {res:.6f} |")
    n_r = sum(1 for l in LINES if R_OF[l] < AZAR)
    n_b = sum(1 for l in LINES if B_OF[l] < AZAR)
    print(f"\n**Conteos.** `r` < azar en **{n_r}/6**. `b` < azar en **{n_b}/6** "
          f"({', '.join(l for l in LINES if B_OF[l] < AZAR)}).")
    print(f"`r` - azar va de **{min(R_OF[l]-AZAR for l in LINES):+.6f}** a "
          f"**{max(R_OF[l]-AZAR for l in LINES):+.6f}**; "
          f"`b` - azar de **{min(B_OF[l]-AZAR for l in LINES):+.6f}** a "
          f"**{max(B_OF[l]-AZAR for l in LINES):+.6f}**.")
    print("\n**El `b` contra el azar con su propia resolucion** "
          "(`SE_b = sd_b/sqrt(8)`, 8 realizaciones de `b`):\n")
    print("| linea | `b` de 8 real. | `sd_b` | `SE_b` | `b` - azar | "
          "\\|`b`-azar\\| > `SE_b` | `u` suelo de 8 sem. | `u` - azar | "
          "`SE_u` | \\|`u`-azar\\| > `SE_u` |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for l in LINES:
        f = fil[l]
        seb, seu = f["sdb"] / math.sqrt(8), f["sdu"] / math.sqrt(8)
        db, du = f["b8"] - AZAR, f["u"] - AZAR
        print(f"| {l} | {f['b8']:.6f} | {f['sdb']:.6f} | {seb:.6f} | {db:+.6f} | "
              f"{'**si**' if abs(db) > seb else 'no'} | {f['u']:.6f} | {du:+.6f} | "
              f"{seu:.6f} | {'**si**' if abs(du) > seu else 'no'} |")

    print("\n### 🔴 TABLA A2: donde cae el `b` DEL BUNDLE dentro de sus 8 "
          "realizaciones\n")
    print("El `b` que define la escala `s` es **una sola realizacion**, la del "
          "bundle. `dispersion_b` midio 8 realizaciones de `b` sobre la misma "
          "poblacion real. Esta tabla sitúa la del bundle dentro de esa "
          "distribucion. 🔴 Es la cantidad que explica §6.1: donde el `b` del "
          "bundle es un valor bajo de su propia distribucion, el suelo cobra "
          "`s > 0` sin ninguna habilidad.\n")
    print("| linea | `b` del bundle | media de los 8 `b` | `sd_b` | "
          "**`z` del `b` del bundle** | `b` bundle - azar | `s` del SUELO (8 sem.) |")
    print("|---|---|---|---|---|---|---|")
    U = jl("dispersion_u.json")
    for l in LINES:
        f = fil[l]
        z = (B_OF[l] - f["b8"]) / f["sdb"]
        ss = [r["scored"]["filas"]["pds_cosine"]["from_replicate"]
              for r in U if r["linea"] == l]
        print(f"| {l} | {B_OF[l]:.6f} | {f['b8']:.6f} | {f['sdb']:.6f} | "
              f"**{z:+.2f}** | {B_OF[l]-AZAR:+.6f} | {st.mean(ss):+.4f} |")
    return fil


# ----------------------------------------------------------------- PARTE B
def parte_b(d, ya, fil):
    print("\n### TABLA B1: `pds_cosine` crudo (`u`) por brazo y por linea\n")
    cols = [b for b, _ in BRAZOS]
    print("| linea | " + " | ".join(f"**{b}**" for b in cols) +
          " | SUELO | gwps REAL | gwps PERM |")
    print("|---" * (len(cols) + 4) + "|")
    uu = {}
    for l in LINES:
        row = [l]
        for b in cols:
            rs = d.get((b, l))
            if not rs:
                row.append("n/d"); continue
            m, s, n = mm(u_pds(rs))
            uu[(b, l)] = (m, s, n)
            row.append(f"**{m:.4f}** ± {s:.4f}" if n > 1 else f"**{m:.4f}**")
        for nm, _ in YA:
            m, s, n = mm(u_pds(ya[(nm, l)]))
            uu[(nm, l)] = (m, s, n)
            row.append(f"{m:.4f} ± {s:.4f}")
        print("| " + " | ".join(row) + " |")

    print("\n### TABLA B2: `s` de `pds_cosine` (`from_replicate`) por brazo y por linea\n")
    print("| linea | " + " | ".join(f"**{b}**" for b in cols) +
          " | SUELO | gwps REAL | gwps PERM |")
    print("|---" * (len(cols) + 4) + "|")
    ss = {}
    for l in LINES:
        row = [l]
        for b in cols:
            rs = d.get((b, l))
            if not rs:
                row.append("n/d"); continue
            m, s, n = mm(s_pds(rs))
            ss[(b, l)] = (m, s, n)
            row.append(f"**{m:+.4f}** ± {s:.4f}" if n > 1 else f"**{m:+.4f}**")
        for nm, _ in YA:
            m, s, n = mm(s_pds(ya[(nm, l)]))
            ss[(nm, l)] = (m, s, n)
            row.append(f"{m:+.4f} ± {s:.4f}")
        print("| " + " | ".join(row) + " |")

    print("\n### 🔴 TABLA B3: la fraccion que alcanza el ORACULO DE COLA\n")
    print("Tres escalas, declaradas. La **anclada en el azar** es la que usa el "
          "punto sin informacion derivado en §1 (`0.5`), y es la que responde "
          "\"¿se puede mover `pds_cosine` tocando solo la cola?\".\n")
    print("| linea | `u` completo (razon) | `u` cola | `u` bulto | `u` azar | "
          "**(cola-azar)/(completo-azar)** | **(bulto-azar)/(completo-azar)** | "
          "suma cola+bulto | `s` cola / `s` completo |")
    print("|---|---|---|---|---|---|---|---|---|")
    fr = {}
    for l in LINES:
        if (("completo", l) not in uu) or (("cola", l) not in uu):
            print(f"| {l} | n/d | | | | | | | |"); continue
        uc = uu[("completo", l)][0]; ut = uu[("cola", l)][0]
        ub = uu[("bulto", l)][0] if ("bulto", l) in uu else None
        den = uc - AZAR
        ft = (ut - AZAR) / den if den else None
        fb = ((ub - AZAR) / den) if (ub is not None and den) else None
        sc = ss[("completo", l)][0]; sst = ss[("cola", l)][0]
        fs = sst / sc if sc else None
        fr[l] = (ft, fb, fs)
        print(f"| {l} | {uc:.4f} | {ut:.4f} | {fmt(ub,4,False)} | {AZAR:.4f} | "
              f"**{fmt(ft,4)}** | **{fmt(fb,4)}** | "
              f"{fmt((ft+fb) if (ft is not None and fb is not None) else None,4)} | "
              f"{fmt(fs,4)} |")
    vals = [v[0] for v in fr.values() if v[0] is not None]
    if vals:
        print(f"\n**Rango de la fraccion de cola anclada en el azar sobre las seis "
              f"lineas: [{min(vals):+.4f}, {max(vals):+.4f}].**")

    print("\n### 🔴 TABLA B3b: la misma fraccion SOLO sobre las perturbaciones "
          "con ED real NO VACIO\n")
    print("En las perturbaciones de ED real vacio el ORACULO DE COLA no emite "
          "NADA (`r = 1` en todos los genes), asi que ahi es el suelo por "
          "construccion y arrastra la media de la TABLA B3 hacia el azar. Aqui "
          "la misma fraccion contada solo donde la cola existe.\n")
    P2 = [r for r in jl("pds_objetivo.json") if not r.get("fallo_en")]
    Pp = jl("control_permutado.json")
    vac = {}
    for l in LINES:
        r = [x for x in Pp if x["linea"] == l][0]
        vac[l] = {g for g, v in r["cohortes"]["por_pert"].items()
                  if v["n_DE_real"] == 0}
    print("| linea | n con ED lleno | `pds` completo (ED lleno) | "
          "`pds` cola (ED lleno) | `pds` bulto (ED lleno) | "
          "**(cola-azar)/(completo-azar)** | `pds` cola (ED vacio) | "
          "`pds` suelo (ED vacio) |")
    print("|---|---|---|---|---|---|---|---|")
    frr = {}
    for l in LINES:
        def med(b, lleno=True):
            rs = [r for r in P2 if r["brazo"] == b and r["linea"] == l]
            if not rs:
                return None
            vs = []
            for r in rs:
                sel = [x for g, x in r["pds_por_pert"].items()
                       if (g not in vac[l]) == lleno]
                if sel:
                    vs.append(st.mean(sel))
            return st.mean(vs) if vs else None
        uc, ut, ub = med("completo"), med("cola"), med("bulto")
        utv = med("cola", lleno=False)
        if uc is None or ut is None:
            continue
        den = uc - AZAR
        f = (ut - AZAR) / den if den else None
        frr[l] = f
        # suelo en las vacias
        sv = []
        for j in range(8):
            fp = os.path.join("/datos/modelo/out/dispu", f"run_{l}_{j}",
                              "results.csv")
            if os.path.exists(fp):
                x = pd.read_csv(fp); x = x[x.metric == "pds_cosine"]
                m = x.perturbation.isin(vac[l])
                if m.any():
                    sv.append(x.loc[m, "value"].mean())
        print(f"| {l} | {33-len(vac[l])} | {uc:.4f} | {ut:.4f} | "
              f"{fmt(ub,4,False)} | **{fmt(f,4)}** | {fmt(utv,4,False)} | "
              f"{(f'{st.mean(sv):.4f}' if sv else 'n/d')} |")
    vv = [x for x in frr.values() if x is not None]
    if vv:
        print(f"\n**Rango de la fraccion de cola restringida sobre las seis "
              f"lineas: [{min(vv):+.4f}, {max(vv):+.4f}].**")

    print("\n### TABLA B4: el rango `[min, max]` de `pds_cosine` por brazo\n")
    print("| brazo | `u` [min, max] | `s` [min, max] | n lineas | n corridas |")
    print("|---|---|---|---|---|")
    for b, nm in BRAZOS + [(x, x) for x, _ in YA]:
        us = [uu[(b, l)][0] for l in LINES if (b, l) in uu]
        sv = [ss[(b, l)][0] for l in LINES if (b, l) in ss]
        nc = sum(uu[(b, l)][2] for l in LINES if (b, l) in uu)
        if not us:
            continue
        print(f"| {nm} | [{min(us):.4f}, {max(us):.4f}] | "
              f"[{min(sv):+.4f}, {max(sv):+.4f}] | {len(us)} | {nc} |")

    print("\n### TABLA B5: cada brazo contra el azar, con su resolucion\n")
    print("| linea | brazo | `u` | `sd` | n | `u` - azar | `SE = sd/sqrt(n)` | "
          "\\|`u`-azar\\| > `SE` |")
    print("|---|---|---|---|---|---|---|---|")
    for l in LINES:
        for b, nm in BRAZOS + [(x, x) for x, _ in YA]:
            if (b, l) not in uu:
                continue
            m, s, n = uu[(b, l)]
            se = s / math.sqrt(n) if n > 1 else None
            dd = m - AZAR
            marca = "n/d (n=1)" if se is None else ("**si**" if abs(dd) > se else "no")
            print(f"| {l} | {nm} | {m:.4f} | {s:.4f} | {n} | {dd:+.4f} | "
                  f"{fmt(se,4,False)} | {marca} |")
    return uu, ss, fr


# ----------------------------------------------------------------- PARTE C
def parte_c(ya):
    print("\n### TABLA C1: la copia del control contra el `b`, con las 8 semillas\n")
    U = jl("dispersion_u.json")
    UNA = {"A549": -0.0423, "BXPC3": 0.0870, "HAP1": 0.2125,
           "HT29": -0.0403, "K562": -0.0847, "MCF7": 0.0197}
    print("| linea | `s` de 1 semilla (`cable_suelo.md`) | `s` de 8 semillas | "
          "`sd` | `SE = sd/sqrt(8)` | `s` > 0 con 8 | \\|`s`\\| > `SE` | se sostiene |")
    print("|---|---|---|---|---|---|---|---|")
    n_pos = n_res = 0
    for l in LINES:
        ss = [r["scored"]["filas"]["pds_cosine"]["from_replicate"]
              for r in U if r["linea"] == l]
        m, s, n = mm(ss)
        se = s / math.sqrt(n)
        pos = m > 0; resu = abs(m) > se
        n_pos += pos; n_res += (pos and resu)
        print(f"| {l} | {UNA[l]:+.4f} | **{m:+.4f}** | {s:.4f} | {se:.4f} | "
              f"{'**si**' if pos else 'no'} | {'**si**' if resu else 'no'} | "
              f"{'**si**' if pos and resu else 'no'} |")
    print(f"\n`s > 0` en **{n_pos}/6** lineas con 8 semillas; de esas, "
          f"**{n_res}/6** con \\|`s`\\| por encima del `SE` de la media de 8.")

    print("\n### TABLA C2: con que puntuan las perturbaciones de ED real VACIO\n")
    P = jl("control_permutado.json")
    vac = {}
    for l in LINES:
        r = [x for x in P if x["linea"] == l][0]
        vac[l] = {g for g, v in r["cohortes"]["por_pert"].items()
                  if v["n_DE_real"] == 0}
    DIRS = {"SUELO": "/datos/modelo/out/dispu",
            "gwps REAL": "/datos/modelo/out/gwps_via",
            "gwps PERM": "/datos/modelo/out/gwps_perm"}
    extra = {}
    for b, _ in BRAZOS:
        extra[b] = None
    print("| linea | n ED real vacio | brazo | `pds` medio de los VACIOS | sd | "
          "`pds` medio de los LLENOS | sd | dif vacios-llenos | `u` total |")
    print("|---|---|---|---|---|---|---|---|---|")
    P2 = [r for r in jl("pds_objetivo.json") if not r.get("fallo_en")]
    for l in LINES:
        for nm, dd in DIRS.items():
            vv, ll, tt = [], [], []
            for j in range(8):
                f = os.path.join(dd, f"run_{l}_{j}", "results.csv")
                if not os.path.exists(f):
                    continue
                x = pd.read_csv(f); x = x[x.metric == "pds_cosine"]
                m = x.perturbation.isin(vac[l])
                if m.any():
                    vv.append(x.loc[m, "value"].mean())
                ll.append(x.loc[~m, "value"].mean()); tt.append(x.value.mean())
            if not tt:
                continue
            sv = f"{st.mean(vv):.4f}" if vv else "n/d"
            sdv = f"{st.stdev(vv):.4f}" if len(vv) > 1 else "-"
            dif = (st.mean(vv) - st.mean(ll)) if vv else None
            print(f"| {l} | {len(vac[l])} | {nm} | **{sv}** | {sdv} | "
                  f"{st.mean(ll):.4f} | {st.stdev(ll):.4f} | {fmt(dif,4)} | "
                  f"{st.mean(tt):.4f} |")
        for b, bn in BRAZOS:
            rs = [r for r in P2 if r["brazo"] == b and r["linea"] == l]
            if not rs:
                continue
            vv, ll, tt = [], [], []
            for r in rs:
                pp = r["pds_por_pert"]
                v = [x for g, x in pp.items() if g in vac[l]]
                w = [x for g, x in pp.items() if g not in vac[l]]
                if v:
                    vv.append(st.mean(v))
                ll.append(st.mean(w)); tt.append(st.mean(list(pp.values())))
            sv = f"{st.mean(vv):.4f}" if vv else "n/d"
            sdv = f"{st.stdev(vv):.4f}" if len(vv) > 1 else "-"
            dif = (st.mean(vv) - st.mean(ll)) if vv else None
            print(f"| {l} | {len(vac[l])} | {bn} | **{sv}** | {sdv} | "
                  f"{st.mean(ll):.4f} | "
                  f"{st.stdev(ll) if len(ll)>1 else 0.0:.4f} | {fmt(dif,4)} | "
                  f"{st.mean(tt):.4f} |")


def composicion():
    p = os.path.join(OUT, "pds_composicion.json")
    if not os.path.exists(p):
        return
    C = json.load(open(p))
    print("\n### TABLA B0: de que esta hecho el delta real que lee `pds_cosine`\n")
    print("Medido en el MISMO espacio de la metrica (`bulk_lognorm`, "
          "`bulk_target_sum = 50000.0` de `run_params.yaml`) y con los 33 genes "
          "blanco del panel FUERA del espacio de rasgos, que es lo que hace "
          "`exclusion_scope=\"panel\"`. `pds_cosine` usa distancia **cosine**, "
          "asi que la cantidad que reparte es la **norma al cuadrado** del "
          "delta.\n")
    print("| linea | genes en el espacio | del panel fuera | perts con ED vacio | "
          "**frac. de \\|\\|delta\\|\\|² en la COLA (media)** | [min, max] | "
          "frac. en el BULTO | `cos`(completo, cola) medio | "
          "`cos`(completo, bulto) medio |")
    print("|---|---|---|---|---|---|---|---|---|")
    for l in LINES:
        if l not in C:
            continue
        c = C[l]
        cc = [f["cos_completo_cola"] for f in c["por_pert"]
              if f["cos_completo_cola"] is not None]
        cb = [f["cos_completo_bulto"] for f in c["por_pert"]
              if f["cos_completo_bulto"] is not None]
        print(f"| {l} | {c['n_genes_espacio']} | {c['n_panel_fuera']} | "
              f"{c['n_perts_ed_vacio']} | **{c['frac_norma2_cola_media']:.4f}** | "
              f"[{c['frac_norma2_cola_min']:.4f}, {c['frac_norma2_cola_max']:.4f}] | "
              f"{1-c['frac_norma2_cola_media']:.4f} | {st.mean(cc):.4f} | "
              f"{st.mean(cb):.4f} |")
    fr = [C[l]["frac_norma2_cola_media"] for l in LINES if l in C]
    if fr:
        print(f"\n**La cola lleva entre {min(fr):.2%} y {max(fr):.2%} de la norma "
              f"al cuadrado del delta real; el bulto el resto.**")


def crudas(d):
    print("\n### TABLA D: las corridas crudas de los cuatro brazos, una por una\n")
    print("| brazo | linea | `j` | semilla emision | `u` `pds_cosine` | "
          "`s` `pds_cosine` | `nnz` pred | max valor | `score_ok` |")
    print("|---|---|---|---|---|---|---|---|---|")
    for b, nm in BRAZOS:
        for l in LINES:
            rs = sorted(d.get((b, l), []), key=lambda r: r["j"])
            for r in rs:
                print(f"| {nm} | {l} | {r['j']} | {r['seed_emision']} | "
                      f"**{r['u']['pds_cosine']:.6f}** | "
                      f"{r['scored']['filas']['pds_cosine']['from_replicate']:+.6f} | "
                      f"{r['contrato_emision']['nnz']} | "
                      f"{r['contrato_emision']['max_valor']:.0f} | "
                      f"{r['score_ok']} |")
            if len(rs) > 1:
                m, s, _ = mm(u_pds(rs)); m2, s2, _ = mm(s_pds(rs))
                print(f"| **{nm} media** | **{l}** | | | **{m:.6f}** (sd {s:.6f}) | "
                      f"**{m2:+.6f}** (sd {s2:.6f}) | | | |")


def verificaciones(d):
    P = [r for r in jl("pds_objetivo.json") if not r.get("fallo_en")]
    n = len(P)
    print("\n### TABLA E: verificaciones de las corridas\n")
    hp = [r["huella_pred"] for r in P]
    print("| | |")
    print("|---|---|")
    print(f"| corridas sin fallo | **{n}** |")
    print(f"| `score_ok` | **{sum(r['score_ok'] for r in P)}/{n}** |")
    print(f"| `config_hash` de competencia coincide | "
          f"**{sum(1 for r in P if r['hash_run'].get('coincide') is True)}/{n}** |")
    print(f"| un solo digest de `config_hash` | "
          f"`{sorted({r['hash_run'].get('referencia') for r in P})[0]}` |")
    print(f"| `poblacion_real_identica` | "
          f"**{sum(r['poblacion_real_identica'] for r in P)}/{n}** |")
    print(f"| huellas reales distintas por linea | "
          f"**{ {l: len({r['huella_real'] for r in P if r['linea']==l}) for l in LINES} }** |")
    print(f"| huellas `pred` distintas | **{len(set(hp))}/{n}** |")
    print(f"| `todo_entero` | "
          f"**{sum(r['contrato_emision']['todo_entero'] for r in P)}/{n}** |")
    print(f"| corridas con algun valor negativo | "
          f"**{sum(r['contrato_emision']['hay_negativos'] for r in P)}** |")
    sem = json.load(open("/datos/modelo/out/pds_obj/semillas.json"))
    print(f"| semillas distintas declaradas | **{sem['n_semillas']}** |")
    print(f"| choques internos de semilla | **{len(sem['choques_internos'])}** |")
    print(f"| choques con familias previas | "
          f"**{len(sem['choques_con_familias_previas'])}** |")
    g = [r["diag_emision"]["_global"] for r in P
         if r.get("diag_emision") and r["diag_emision"].get("_global")
         and "n_control_mean_cero" in r["diag_emision"]["_global"]]
    if g:
        ncz = sorted({x["n_control_mean_cero"] for x in g})
        nvar = sorted({x["n_var"] for x in g})
        print(f"| genes con media de control 0 (inalcanzables) | "
              f"**{min(ncz)} a {max(ncz)}** de "
              f"{nvar[0] if len(nvar) == 1 else f'{min(nvar)} a {max(nvar)}'}"
              f" (por linea abajo) |")
        fr = [x["fraccion_delta_real_en_genes_sin_control"] for x in g
              if x.get("fraccion_delta_real_en_genes_sin_control") is not None]
        if fr:
            print(f"| fraccion de \\|delta real\\| en esos genes | "
                  f"**{min(fr):.6f} a {max(fr):.6f}** |")
        # 🔴 por linea: la cifra NO es comun a las seis
        gl = {}
        for r in P:
            x = (r.get("diag_emision") or {}).get("_global")
            if x and "n_control_mean_cero" in x:
                gl.setdefault(r["linea"], x)
        print("\n**Los genes inalcanzables por `_reescala_entero`, por linea** "
              "(media de control 0, asi que ninguna razon los mueve):\n")
        print("| linea | genes con media de control 0 | de `n_var` | "
              "fraccion de \\|delta real\\| que llevan |")
        print("|---|---|---|---|")
        for l in LINES:
            x = gl.get(l)
            if not x:
                print(f"| {l} | n/d | n/d | n/d |")
                continue
            f_ = x.get("fraccion_delta_real_en_genes_sin_control")
            print(f"| {l} | **{x['n_control_mean_cero']}** | {x['n_var']} | "
                  f"{'n/d' if f_ is None else f'{f_:.6f}'} |")


def main():
    d, ya = carga()
    fil = parte_a()
    composicion()
    uu, ss, fr = parte_b(d, ya, fil)
    parte_c(ya)
    crudas(d)
    verificaciones(d)


if __name__ == "__main__":
    main()
