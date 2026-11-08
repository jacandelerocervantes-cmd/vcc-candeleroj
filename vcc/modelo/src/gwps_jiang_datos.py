#!/usr/bin/env python
"""Relacion entre el perfil de respuesta de k562_gwps y el delta especifico de Jiang.

🔴 DOS POBLACIONES DISTINTAS:
  gwps  = K562, SIN estimular, knockdown CRISPRi (Replogle)
  Jiang = seis lineas, CON IFNB
No son la misma cosa. Lo que se mide es si una estructura de la primera se
relaciona con la segunda, no si son equivalentes.

CONSTRUCCION, declarada antes de medir:

ENTRADA gwps: /datos/ref/replogle/k562_gwps_raw_bulk.h5ad
  `X` (11258 perturbaciones x 8248 genes) = pseudobulk normalizado por
  perturbacion. Es la UNICA matriz del archivo (no hay layers ni obsm).
  `obs/gene_transcript` identifica la perturbacion con el formato
  `{idx}_{GEN}_{P1|P2|P1P2}_{ENSG}`; `obs/core_control` marca los
  no-dirigidos.

PERFIL DE RESPUESTA de la perturbacion g, declarado:
  p_g = log2( (X_g + eps) / (media de X sobre core_control + eps) )
  eps = 1e-8. Es log-fold-change contra el pseudobulk de los no-dirigidos.

ALINEACION DE UNIVERSOS: `var/gene_name` (categorico) contra el array `genes`
  de `rasgos_<linea>.npz` (el universo de Jiang de control_o_no_rasgos.py).
  Los nombres de gen DUPLICADOS en gwps se excluyen enteros (no se elige uno).
  Se interseca por SIMBOLO.

E_g DE JIANG: `E1`, `E2` de `rasgos_<linea>.npz`, que son exactamente el delta
  ESPECIFICO con el control PARTIDO en dos mitades de
  `control_o_no_rasgos.py`. No se recalcula nada.

METRICA, con el producto cruzado de mitades de Jiang para que el ruido de Jiang
  no entre (como control_o_no.md §2.2), pero EN ESPACIO DE RANGOS porque lo que
  puede transferir es el orden:
     p = rango(perfil gwps), centrado
     a = rango(E1), b = rango(E2), centrados
     G = <p, (a+b)/2> / ( ||p|| * sqrt(<a,b>) )
  Los rangos de E1 y E2 son funciones de mitades de control DISJUNTAS, asi que
  su ruido sigue siendo independiente y el producto cruzado sigue valiendo.
  `<a,b> <= 0`  ->  n/d.
  Se dan en columna aparte: el Spearman INGENUO (contra (E1+E2)/2, atenuado por
  el ruido de Jiang) y el PEARSON con el mismo producto cruzado.

CUATRO CANTIDADES:
  G_mismo : perfil gwps de g            contra E_g
  G_otro  : perfil gwps de g            contra E_q, q != g   (CONTROL NEGATIVO)
  A_esp   : el techo, ya medido (control_o_no.md §8.1), leido del npz
  G_prom  : promedio de perfiles gwps   contra E_g  (para ver si lo que
            correlaciona es la respuesta generica de gwps y no la del blanco)

LAS DOS CONSTRUCCIONES: los blancos con P1 y P2 por separado se reportan POR
  SEPARADO y NO se promedian (gwps_fiabilidad.md midio que P1 domina en 67.6 %
  de 546 pares, asi que la media no es neutral).

EL GEN BLANCO: se reporta con y sin el, porque es el punto de mayor palanca en
  los dos lados (control_o_no.md §4.3).

SEMILLAS: familia nueva 20270000 + 100*fold, para el nulo por permutacion de
  etiquetas de gen del perfil de gwps.

🔴 NO se escribe ninguna conclusion. Solo numeros y sus controles.
"""
import json
import os
import sys

import h5py
import numpy as np

GWPS = "/datos/ref/replogle/k562_gwps_raw_bulk.h5ad"
DIR = "/home/ubuntu/vcc/modelo/out/control_o_no"
PANEL = "/home/ubuntu/vcc/modelo/out/suelo_panel.json"
OUT = "/home/ubuntu/vcc/modelo/out/gwps_jiang_datos.json"
LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
SEED_FAMILIA = 20270000
EPS = 1e-8
N_PERM = 200


def rangos(v):
    """Rangos medios (empates promediados), centrados."""
    n = len(v)
    o = np.argsort(v, kind="stable")
    r = np.empty(n, dtype=np.float64)
    sv = v[o]
    i = 0
    while i < n:
        j = i
        while j + 1 < n and sv[j + 1] == sv[i]:
            j += 1
        r[o[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return r - r.mean()


def cruzado(p, a, b):
    """<p,(a+b)/2> / (||p|| sqrt(<a,b>)); n/d si <a,b><=0."""
    c = float(a @ b)
    np_ = float(np.linalg.norm(p))
    if c <= 0 or np_ == 0:
        return np.nan, c
    return float(p @ ((a + b) / 2.0)) / (np_ * np.sqrt(c)), c


def main():
    panel = json.load(open(PANEL))["panel_comun"]

    # ---------- gwps: identidades, controles, perfiles ----------
    f = h5py.File(GWPS, "r")
    gt = np.array([x.decode() if isinstance(x, bytes) else x
                   for x in f["obs/gene_transcript"][:]])
    partes = [x.split("_") for x in gt]
    g_gen = np.array([p[1] if len(p) > 1 else "" for p in partes])
    g_con = np.array([p[2] if len(p) > 2 else "" for p in partes])
    core = f["obs/core_control"][:]
    cats = np.array([x.decode() if isinstance(x, bytes) else x
                     for x in f["var/__categories/gene_name"][:]])
    vn = cats[f["var/gene_name"][:]]

    # nombres duplicados fuera enteros
    u, cnt = np.unique(vn, return_counts=True)
    dup = set(u[cnt > 1])
    ok_var = np.array([g not in dup for g in vn])

    # filas de interes: los del panel + controles
    filas = {}
    for g in panel:
        for i in np.where(g_gen == g)[0]:
            filas.setdefault(g, []).append((g_con[i], int(i)))
    cubiertos = [g for g in panel if g in filas]
    no_cub = [g for g in panel if g not in filas]
    multi = {g: v for g, v in filas.items() if len(v) > 1}

    sys.stderr.write("cubiertos %d/33 ; no cubiertos %d ; multi-construccion %d\n"
                     % (len(cubiertos), len(no_cub), len(multi)))
    sys.stderr.write("  no cubiertos: %s\n" % ", ".join(no_cub))
    for g, v in multi.items():
        sys.stderr.write("  %s: %s\n" % (g, [c for c, _ in v]))

    X = f["X"]
    # media de control, en bloques
    idx_c = np.where(core)[0]
    ctrl = np.zeros(X.shape[1], dtype=np.float64)
    for s in range(0, len(idx_c), 64):
        ctrl += X[idx_c[s:s + 64], :].astype(np.float64).sum(0)
    ctrl /= len(idx_c)

    # media de TODAS las perturbaciones dirigidas (respuesta generica)
    idx_p = np.where(~core)[0]
    gen_mean = np.zeros(X.shape[1], dtype=np.float64)
    for s in range(0, len(idx_p), 512):
        gen_mean += X[idx_p[s:s + 512], :].astype(np.float64).sum(0)
    gen_mean /= len(idx_p)

    def perfil(fila_idx):
        x = X[fila_idx, :].astype(np.float64)
        return np.log2((x + EPS) / (ctrl + EPS))

    perf = {}                     # (gen, construccion) -> perfil completo
    for g in cubiertos:
        for c, i in filas[g]:
            perf[(g, c)] = perfil(i)
    perf_prom_cub = np.mean([perf[k] for k in perf], axis=0)
    perf_prom_all = np.log2((gen_mean + EPS) / (ctrl + EPS))
    sys.stderr.write("controles gwps: %d filas core_control ; perturbadas: %d\n"
                     % (len(idx_c), len(idx_p)))

    sal = {"cubiertos": cubiertos, "no_cubiertos": no_cub,
           "multi_construccion": {g: [c for c, _ in v] for g, v in multi.items()},
           "construcciones": {g: [c for c, _ in v] for g, v in filas.items()},
           "n_core_control": int(len(idx_c)), "n_perturbadas": int(len(idx_p)),
           "n_var_gwps": int(X.shape[1]), "n_var_duplicados": int(len(dup)),
           "duplicados": sorted(dup), "semilla_familia": SEED_FAMILIA,
           "n_perm": N_PERM, "lineas": {}}

    # ---------- por linea ----------
    for fold, line in enumerate(LINES):
        d = np.load(os.path.join(DIR, "rasgos_%s.npz" % line), allow_pickle=True)
        jg = np.array([str(x) for x in d["genes"]])
        E1 = d["E1"].astype(np.float64); E2 = d["E2"].astype(np.float64)
        jpanel = [str(x) for x in d["panel"]]
        A_esp = d["A_esp"]
        rasgos = [str(x) for x in d["rasgos"]]
        auto = d["F"][:, :, rasgos.index("es_el_blanco")] > 0.5

        # interseccion de universos
        pos_g = {g: i for i, g in enumerate(vn) if ok_var[i]}
        com_j, com_g = [], []
        for i, g in enumerate(jg):
            if g in pos_g:
                com_j.append(i); com_g.append(pos_g[g])
        com_j = np.array(com_j); com_g = np.array(com_g)
        nom_com = jg[com_j]

        res = {"n_comun": int(len(com_j)), "n_jiang": int(len(jg)),
               "filas": [], "por_blanco": {}}
        sys.stderr.write("%-6s genes comunes %d de %d\n" % (line, len(com_j), len(jg)))

        rng = np.random.default_rng(SEED_FAMILIA + 100 * fold)

        for g in cubiertos:
            if g not in jpanel:
                continue
            j = jpanel.index(g)
            for c, _ in filas[g]:
                p_full = perf[(g, c)][com_g]
                a_full = E1[com_j, j]; b_full = E2[com_j, j]
                for quita, etq in ((False, "con_gb"), (True, "sin_gb")):
                    if quita:
                        keep = ~auto[com_j, j] & (nom_com != g)
                    else:
                        keep = np.ones(len(com_j), bool)
                    p = rangos(p_full[keep])
                    a = rangos(a_full[keep]); b = rangos(b_full[keep])
                    G, cz = cruzado(p, a, b)
                    # Pearson con el mismo producto cruzado
                    pp = p_full[keep] - p_full[keep].mean()
                    aa = a_full[keep] - a_full[keep].mean()
                    bb = b_full[keep] - b_full[keep].mean()
                    Gp, czp = cruzado(pp, aa, bb)
                    # Spearman ingenuo contra la media de mitades
                    y = rangos(((a_full + b_full) / 2.0)[keep])
                    sn = float(p @ y / (np.linalg.norm(p) * np.linalg.norm(y)))
                    # G_otro: contra los E_q de los otros blancos
                    otros = []
                    for q in range(len(jpanel)):
                        if q == j:
                            continue
                        if quita:
                            kq = ~auto[com_j, j] & (nom_com != g)
                        else:
                            kq = keep
                        aq = rangos(E1[com_j, q][kq]); bq = rangos(E2[com_j, q][kq])
                        pq = rangos(p_full[kq])
                        Gq, _ = cruzado(pq, aq, bq)
                        if Gq == Gq:
                            otros.append(Gq)
                    # G_prom: perfiles promedio contra E_g
                    pc = rangos(perf_prom_cub[com_g][keep])
                    pa = rangos(perf_prom_all[com_g][keep])
                    Gc, _ = cruzado(pc, a, b)
                    Ga, _ = cruzado(pa, a, b)
                    # nulo por permutacion de etiquetas de gen del perfil gwps
                    nulo = []
                    for _ in range(N_PERM):
                        ps = rangos(p_full[keep][rng.permutation(keep.sum())])
                        Gn, _ = cruzado(ps, a, b)
                        nulo.append(Gn)
                    nulo = np.array(nulo)
                    res["filas"].append({
                        "blanco": g, "construccion": c, "gb": etq,
                        "n_genes": int(keep.sum()),
                        "G_mismo": G, "G_mismo_pearson": Gp,
                        "spearman_ingenuo": sn,
                        "G_otro_media": float(np.mean(otros)) if otros else np.nan,
                        "G_otro_mediana": float(np.median(otros)) if otros else np.nan,
                        "G_otro_sd": float(np.std(otros, ddof=1)) if len(otros) > 1 else np.nan,
                        "G_otro_min": float(np.min(otros)) if otros else np.nan,
                        "G_otro_max": float(np.max(otros)) if otros else np.nan,
                        "n_otros": len(otros),
                        "G_prom_cubiertos": Gc, "G_prom_todas": Ga,
                        "A_esp": float(A_esp[j]),
                        "cruz_rangos": cz,
                        "nulo_media": float(np.nanmean(nulo)),
                        "nulo_sd": float(np.nanstd(nulo, ddof=1)),
                        "z_vs_nulo": float((G - np.nanmean(nulo)) / np.nanstd(nulo, ddof=1))
                        if np.nanstd(nulo, ddof=1) > 0 else np.nan,
                    })
        sal["lineas"][line] = res
        with open(OUT, "w") as fh:
            json.dump(sal, fh)
        n = [r for r in res["filas"] if r["gb"] == "sin_gb"]
        sys.stderr.write("  %-6s sin_gb: G_mismo media %+.4f ; G_otro media %+.4f\n" % (
            line, np.nanmean([r["G_mismo"] for r in n]),
            np.nanmean([r["G_otro_media"] for r in n])))
        sys.stderr.flush()

    with open(OUT, "w") as fh:
        json.dump(sal, fh)
    sys.stderr.write("escrito %s\n" % OUT)


if __name__ == "__main__":
    main()
