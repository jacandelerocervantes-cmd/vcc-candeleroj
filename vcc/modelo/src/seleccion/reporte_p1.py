#!/usr/bin/env python
"""Tablas de la PARTE 1 (superficie de entrega y de validacion). Solo lee
`p1_filtrables.json`; no computa nada que no salga de los controles del reto."""
import json
from itertools import combinations, permutations

import numpy as np
from scipy.stats import spearmanr

P1 = "/datos/modelo/out/seleccion/p1_filtrables.json"
LISTAS = "/home/ubuntu/vcc/modelo/out/seleccion_filtrables_listas.json"
PARTS = {"D/E/F": ["D", "E", "F"], "A/B/C": ["A", "B", "C"]}


def tabla(enc, filas):
    o = ["| " + " | ".join(map(str, enc)) + " |", "|" + "|".join(["---"] * len(enc)) + "|"]
    for f in filas:
        o.append("| " + " | ".join("" if c is None else str(c) for c in f) + " |")
    return o


def q(v, ps=(0, 10, 25, 50, 75, 90, 100)):
    v = np.asarray([x for x in v if x is not None and np.isfinite(x)])
    return [f"{np.percentile(v, p):.3g}" for p in ps] if v.size else ["n/d"] * len(ps)


def lista_corta(gs, k=40):
    gs = sorted(gs)
    if len(gs) <= k:
        return ", ".join(f"`{g}`" for g in gs) if gs else "—"
    return f"{len(gs)} genes (lista en `seleccion_filtrables_listas.json`)"


def main():
    d = json.load(open(P1))
    U = d["umbral"]
    genes = np.array(d["genes"])
    G = len(genes)
    cx = d["ctx"]
    cpm = {L: np.array(cx[L]["cpm"]) for L in cx}
    det = {L: np.array(cx[L]["det"]) for L in cx}
    fano = {L: np.array([np.nan if x is None else x for x in cx[L]["fano"]]) for L in cx}
    pasa = {L: cpm[L] > U for L in cx}
    W = []
    W += ["## PARTE 1 — SUPERFICIE DE ENTREGA (D/E/F) Y DE VALIDACION (A/B/C), SOLO DE SUS CONTROLES", "",
          "Toda cifra de esta parte sale de los 18,400 controles de cada contexto. No hay verdad "
          "en ella y no la necesita: el filtro del preset lo decide la referencia, y la referencia "
          "en la puntuacion real ES este pool.", ""]
    # ---------------- 1.1
    W += ["### 1.1 Que genes pasan el filtro, por contexto", "",
          f"Filtro: CPM medio en el control `> {U}` (estricto), de `EvalConfig.from_preset('vcc2026')`. "
          "Columna *verif*: el conjunto de genes que `compute_de` del preset PRUEBA con referencia = los 18,400 "
          "controles y blanco = 400 controles (semilla en la columna), comparado gen a gen con esta cuenta.", ""]
    filas = []
    for i, L in enumerate(["A", "B", "C", "D", "E", "F"]):
        c = cx[L]; v = c["verif"]
        filas.append([i + 1, f"**{L}**", c["particion"], f"{c['n']:,}", f"{c['L_mediana']:,.0f}",
                      f"{c['nnz_mediana']:,.0f}", f"**{c['n_pasa']:,}**", f"{100*c['n_pasa']/G:.2f} %",
                      f"{v['n_probados_compute_de']:,}", "**si**" if v["identicos"] else f"NO ({v['solo_compute_de']}/{v['solo_mios']})",
                      v["semilla"]])
    W += tabla(["# ", "contexto", "particion", "n control", "UMI mediana/celula", "nnz mediana/celula",
                "**# pasan**", f"de {G:,}", "# probados por compute_de", "identicos gen a gen", "semilla verif"], filas)
    W += ["", f"Las listas completas: `{LISTAS}`. Orden de genes identico en los 6 contextos y a "
          f"`gene_names.csv`: **{'si' if d.get('genes_igual_gene_names_csv') else 'NO'}**.", ""]
    # ---------------- 1.2
    W += ["### 1.2 🔴 El Venn de los filtrables, dentro de cada particion", ""]
    venn = {}
    for nom, (a, b, c) in PARTS.items():
        A, B, Cc = pasa[a], pasa[b], pasa[c]
        reg = [(f"solo {a}", A & ~B & ~Cc), (f"solo {b}", ~A & B & ~Cc), (f"solo {c}", ~A & ~B & Cc),
               (f"{a}∩{b} sin {c}", A & B & ~Cc), (f"{a}∩{c} sin {b}", A & ~B & Cc),
               (f"{b}∩{c} sin {a}", ~A & B & Cc), (f"**{a}∩{b}∩{c}**", A & B & Cc)]
        venn[nom] = reg
        n3 = int((A & B & Cc).sum()); n0 = int((~A & ~B & ~Cc).sum())
        n1 = sum(int(m.sum()) for _, m in reg[:3]); n2 = sum(int(m.sum()) for _, m in reg[3:6])
        W += [f"**{nom}**", ""]
        W += tabla(["# ", "region", "# genes", "frac de la union", "genes (si <= 40)"],
                   [[i + 1, r, f"{int(m.sum()):,}", f"{100*m.sum()/max(1,(A|B|Cc).sum()):.2f} %", lista_corta(genes[m])]
                    for i, (r, m) in enumerate(reg)])
        W += [""]
        W += tabla(["pasan en", "# genes", f"frac de {G:,}"],
                   [["los 3", f"{n3:,}", f"{100*n3/G:.2f} %"], ["exactamente 2", f"{n2:,}", f"{100*n2/G:.2f} %"],
                    ["exactamente 1", f"{n1:,}", f"{100*n1/G:.2f} %"], ["ninguno", f"{n0:,}", f"{100*n0/G:.2f} %"],
                    ["union", f"{int((A|B|Cc).sum()):,}", f"{100*(A|B|Cc).sum()/G:.2f} %"]])
        W += [""]
    # ---------------- 1.3
    W += ["### 1.3 Cantidades por gen y contexto, todas del control", "",
          "Cuantiles sobre genes. CPM medio = media de x/L·1e6; deteccion = fraccion de celulas con cuenta > 0; "
          "Fano = var/media de las CUENTAS CRUDAS (depende de la profundidad: no es comparable entre contextos de "
          "profundidad distinta sin esa salvedad).", ""]
    for etq, sub in (("todos los genes con media > 0", None), ("solo los que pasan el filtro del contexto", "pasa")):
        filas = []
        for L in ["A", "B", "C", "D", "E", "F"]:
            m = (cpm[L] > 0) if sub is None else pasa[L]
            filas.append([f"**{L}**", f"{int(m.sum()):,}"] + q(cpm[L][m], (10, 50, 90)) +
                         q(det[L][m], (10, 50, 90)) + q(fano[L][m], (10, 50, 90)))
        W += [f"**{etq}**", ""]
        W += tabla(["ctx", "# genes", "CPM p10", "CPM p50", "CPM p90", "det p10", "det p50", "det p90",
                    "Fano p10", "Fano p50", "Fano p90"], filas)
        W += [""]
    W += ["**En que se diferencian los tres contextos de cada particion** (pares; sobre los genes que pasan en "
          "al menos uno del par)", ""]
    filas = []
    for nom, ls in PARTS.items():
        for a, b in combinations(ls, 2):
            m = pasa[a] | pasa[b]
            la, lb = np.log2(cpm[a][m] + 0.01), np.log2(cpm[b][m] + 0.01)
            rho = spearmanr(la, lb).statistic
            jac = (pasa[a] & pasa[b]).sum() / m.sum()
            dd = det[a][m] - det[b][m]
            filas.append([nom, f"{a}–{b}", f"{int(m.sum()):,}", f"{jac:.4f}", f"{int((pasa[a]&~pasa[b]).sum()):,}",
                          f"{int((pasa[b]&~pasa[a]).sum()):,}", f"{rho:.4f}", f"{np.median(np.abs(la-lb)):.3f}",
                          f"{np.percentile(np.abs(la-lb), 90):.3f}", f"{np.median(dd):+.4f}",
                          f"{cx[a]['L_mediana']/cx[b]['L_mediana']:.3f}"])
    W += tabla(["particion", "par", "# genes (union)", "jaccard filtrables", "pasan en 1º no 2º",
                "pasan en 2º no 1º", "Spearman log2 CPM", "mediana |log2 ratio CPM|", "p90 |log2 ratio|",
                "mediana dif. deteccion (1º−2º)", "ratio UMI mediana (1º/2º)"], filas)
    W += [""]
    for nom, ls in PARTS.items():
        a, b, c = ls
        m = (pasa[a] | pasa[b] | pasa[c]) & ~(pasa[a] & pasa[b] & pasa[c])
        idx = np.where(m)[0]
        idx = idx[np.argsort(cpm[a][idx])]
        sel = idx[np.linspace(0, len(idx) - 1, 20).round().astype(int)] if len(idx) else []
        W += [f"🔴 **Intermedios {nom}: 20 genes que pasan en algun contexto pero no en los 3**, "
              f"equiespaciados por CPM en {a} ({len(idx):,} candidatos)", ""]
        filas = []
        for i, j in enumerate(sel):
            fila = [i + 1, f"`{genes[j]}`"]
            for L in ls:
                fila += [f"{cpm[L][j]:.3f}", f"{det[L][j]:.4f}", "n/d" if not np.isfinite(fano[L][j]) else f"{fano[L][j]:.3f}",
                         "si" if pasa[L][j] else "no"]
            filas.append(fila)
        enc = ["# ", "gen"]
        for L in ls:
            enc += [f"CPM {L}", f"det {L}", f"Fano {L}", f"pasa {L}"]
        W += tabla(enc, filas) + [""]
    # ---------------- 1.4
    W += ["### 1.4 🔴 Por cuanto fallan: genes que pasan en X y NO en Y, su CPM en Y respecto de 5.0", "",
          "Columnas de distribucion: fraccion de esos genes con CPM_Y/5.0 en cada intervalo. Si se agolpan junto "
          "a 1, el filtro es fragil; si estan lejos, es estructural.", "",
          "🔴 **CONTROL** (filas `mitad`): el MISMO conteo entre dos mitades del MISMO contexto (9,200 vs 9,200), "
          "que no difieren en nada salvo muestreo — (a) aleatoria de celulas, (b) por guia `ntc_id` (23 vs 23). "
          "Defecto: cada mitad tiene la mitad de celulas que un contexto, asi que su ruido de muestreo es ~√2 "
          "mayor que el de un contexto entero; el control ACOTA POR ARRIBA la fragilidad por muestreo.", ""]
    BINS = [(0.95, 1.0000001), (0.9, 0.95), (0.8, 0.9), (0.5, 0.8), (0.2, 0.5), (-1, 0.2)]
    BN = ["(0.95,1]", "(0.9,0.95]", "(0.8,0.9]", "(0.5,0.8]", "(0.2,0.5]", "≤0.2"]

    def fila14(etq, cy, cxv):
        r = cy / U
        rx = cxv / U
        fr = [f"{100*np.mean((r > lo) & (r <= hi)):.1f} %" if r.size else "—" for lo, hi in BINS]
        return [etq, f"{r.size:,}"] + fr + ([f"{np.median(r):.3f}", f"{np.median(rx):.3f}"] if r.size else ["—", "—"])

    filas = []
    for nom, ls in PARTS.items():
        for x, y in permutations(ls, 2):
            m = pasa[x] & ~pasa[y]
            filas.append(fila14(f"{nom}: pasa **{x}**, no **{y}**", cpm[y][m], cpm[x][m]))
    for L in ["A", "B", "C", "D", "E", "F"]:
        for etq in ("aleatoria", "por_guia"):
            mm = cx[L]["mitades"][etq]
            ca, cb = np.array(mm["cpm_a"]), np.array(mm["cpm_b"])
            m1 = (ca > U) & ~(cb > U); m2 = (cb > U) & ~(ca > U)
            filas.append(fila14(f"mitad {L} ({etq}), ambos sentidos", np.concatenate([cb[m1], ca[m2]]),
                                np.concatenate([ca[m1], cb[m2]])))
    W += tabla(["comparacion", "# genes"] + [f"CPM_Y/5 en {b}" for b in BN] + ["mediana CPM_Y/5", "mediana CPM_X/5"], filas)
    W += [""]
    m = pasa["D"] & ~pasa["E"]
    idx = np.where(m)[0]
    idx = idx[np.argsort(cpm["E"][idx])]
    sel = idx[np.linspace(0, len(idx) - 1, 20).round().astype(int)] if len(idx) else []
    W += [f"🔴 **Intermedios: 20 de los {len(idx):,} genes que pasan en D y no en E**, equiespaciados por CPM en E", ""]
    W += tabla(["# ", "gen", "CPM D", "CPM E", "CPM_E/5", "CPM F", "det D", "det E", "pasa F"],
               [[i + 1, f"`{genes[j]}`", f"{cpm['D'][j]:.3f}", f"{cpm['E'][j]:.3f}", f"{cpm['E'][j]/U:.3f}",
                 f"{cpm['F'][j]:.3f}", f"{det['D'][j]:.4f}", f"{det['E'][j]:.4f}", "si" if pasa["F"][j] else "no"]
                for i, j in enumerate(sel)])
    W += [""]
    # ---------------- 1.5
    W += ["### 1.5 El Venn entre las dos particiones: ¿los filtrables de A/B/C son los de D/E/F?", "",
          "Tabla cruzada: en cuantos de A/B/C pasa el gen (filas) contra en cuantos de D/E/F (columnas).", ""]
    kA = pasa["A"].astype(int) + pasa["B"] + pasa["C"]
    kD = pasa["D"].astype(int) + pasa["E"] + pasa["F"]
    filas = []
    for i in range(4):
        filas.append([f"**{i}** de A/B/C"] + [f"{int(((kA==i)&(kD==j)).sum()):,}" for j in range(4)] +
                     [f"{int((kA==i).sum()):,}"])
    filas.append(["total"] + [f"{int((kD==j).sum()):,}" for j in range(4)] + [f"{G:,}"])
    W += tabla(["", "0 de D/E/F", "1 de D/E/F", "2 de D/E/F", "3 de D/E/F", "total"], filas)
    W += [""]
    nucA, nucD = kA == 3, kD == 3
    uniA, uniD = kA > 0, kD > 0
    filas = []
    for etq, a, b in (("nucleo (pasa en los 3)", nucA, nucD), ("union (pasa en >= 1)", uniA, uniD)):
        filas.append([etq, f"{int(a.sum()):,}", f"{int(b.sum()):,}", f"{int((a&b).sum()):,}",
                      f"{int((a&~b).sum()):,}", f"{int((b&~a).sum()):,}", f"{(a&b).sum()/max(1,(a|b).sum()):.4f}",
                      lista_corta(genes[a & ~b], 25), lista_corta(genes[b & ~a], 25)])
    W += tabla(["conjunto", "A/B/C", "D/E/F", "comun", "solo A/B/C", "solo D/E/F", "jaccard",
                "genes solo A/B/C (<= 25)", "genes solo D/E/F (<= 25)"], filas)
    W += ["", "**Jaccard de filtrables, contexto a contexto (A/B/C × D/E/F)**", ""]
    W += tabla(["", "D", "E", "F"],
               [[f"**{a}**"] + [f"{(pasa[a]&pasa[b]).sum()/(pasa[a]|pasa[b]).sum():.4f}" for b in "DEF"] for a in "ABC"])
    W += ["", "Referencia dentro de particion: " + ", ".join(
        f"{a}–{b} {(pasa[a]&pasa[b]).sum()/(pasa[a]|pasa[b]).sum():.4f}"
        for ls in PARTS.values() for a, b in combinations(ls, 2)), ""]
    kA_ = pasa["A"].astype(int) + pasa["B"] + pasa["C"]
    kD_ = pasa["D"].astype(int) + pasa["E"] + pasa["F"]
    json.dump({"umbral": U, "n_genes": G,
               "filtrables": {L: genes[pasa[L]].tolist() for L in sorted(cx)},
               "venn": {nom: {r.replace("*", ""): genes[mm].tolist() for r, mm in reg} for nom, reg in venn.items()},
               "entre_particiones": {"nucleo_solo_ABC": genes[(kA_ == 3) & (kD_ < 3)].tolist(),
                                     "nucleo_solo_DEF": genes[(kD_ == 3) & (kA_ < 3)].tolist(),
                                     "union_solo_ABC": genes[(kA_ > 0) & (kD_ == 0)].tolist(),
                                     "union_solo_DEF": genes[(kD_ > 0) & (kA_ == 0)].tolist()}},
              open(LISTAS, "w"), indent=0)
    return W, d


if __name__ == "__main__":
    W, _ = main()
    print("\n".join(W))
