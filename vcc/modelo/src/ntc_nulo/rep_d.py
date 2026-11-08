#!/usr/bin/env python
"""Reporte `ntc_nulo.md` -- 1.4 los intermedios crudos, corrida por corrida."""
import json, os
import numpy as np
OUTDIR = "/home/ubuntu/vcc/modelo/out/ntc_nulo"
BRAZOS = ["completo", "n400_400", "n400_ref_completa"]
CORTO = {"completo": "completo", "n400_400": "n400_400", "n400_ref_completa": "n400_ref"}


def main():
    p1 = json.load(open(os.path.join(OUTDIR, "p1.json")))
    cors = p1["corridas"]
    L = []; W = L.append

    W("### 1.4 🔴 Los intermedios crudos: las corridas una por una")
    W("")
    W("Ninguna cifra de las tablas de resumen es un agregado sin sus filas. Aqui")
    W("estan las **%d corridas** completas." % len(cors))
    W("")
    W("Columnas: `m` = genes probados tras el filtro CPM del preset; `p_crit` =")
    W("`0.05/m`, el escalon de BH de rango 1, que es lo que una senal AISLADA tiene")
    W("que batir; `frac p<0.05` = fraccion de los `m` genes con `p` cruda bajo 0.05,")
    W("que bajo un nulo uniforme vale ~0.05; `KS` = estadistico de")
    W("Kolmogorov-Smirnov de las `p` crudas contra la uniforme (diagnostico del nulo,")
    W("**no** parte de la puerta); `|lfc|` en el espacio log2 del propio marcador.")
    W("")
    W("**TABLA 1.3 — las %d corridas crudas.**" % len(cors))
    W("")
    W("| # | ctx | j | semilla | brazo | n ref | n blanco | `m` | **`n_sig`** | "
      "`p_adj` min | `p` min | `p_crit` | `p` esp. min unif. | frac `p<0.05` | "
      "frac `p<0.01` | `p` mediana | KS unif. | `|lfc|` mediana | `|lfc|` p99 | "
      "`|lfc|` max | `lfc` NaN | seg |")
    W("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    ordenado = sorted(cors, key=lambda r: (r["contexto"], r["j"],
                                           BRAZOS.index(r["brazo"])))
    for i, r in enumerate(ordenado, 1):
        W(f"| {i} | **{r['contexto']}** | {r['j']} | `{r['semilla']}` | "
          f"`{CORTO[r['brazo']]}` | {r['n_ref']:,} | {r['n_tgt']:,} | "
          f"{r['n_probados']:,} | **{r['n_sig']}** | {r['p_adj_min']:.4f} | "
          f"{r['p_min']:.3e} | {r['umbral_bh_rango1']:.3e} | "
          f"{r['p_min_esperado_uniforme']:.3e} | {r['frac_p_menor_05']:.4f} | "
          f"{r['frac_p_menor_01']:.4f} | {r['p_mediana']:.4f} | {r['ks_uniforme']:.4f} | "
          f"{r['abs_lfc_mediana']:.4f} | {r['abs_lfc_p99']:.4f} | "
          f"{r['abs_lfc_max']:.4f} | {r['n_lfc_nan']} | {r['segundos']} |")
    W("")
    # agregados de los diagnosticos por brazo
    W("**Los diagnosticos del nulo, agregados por brazo** (de las mismas filas de")
    W("arriba, para poder cruzarlos):")
    W("")
    W("| brazo | corridas | frac `p<0.05` media (sd) | frac `p<0.01` media (sd) | "
      "`p` mediana media | KS media (sd) | `p` min minimo | `|lfc|` mediana media | "
      "`|lfc|` max maximo | `lfc` NaN total |")
    W("|---|---|---|---|---|---|---|---|---|---|")
    for b in BRAZOS:
        rs = [r for r in cors if r["brazo"] == b]
        if not rs:
            continue
        g = lambda k: np.array([r[k] for r in rs], dtype=float)
        W(f"| `{CORTO[b]}` | {len(rs)} | "
          f"{g('frac_p_menor_05').mean():.4f} ({g('frac_p_menor_05').std(ddof=1):.4f}) | "
          f"{g('frac_p_menor_01').mean():.4f} ({g('frac_p_menor_01').std(ddof=1):.4f}) | "
          f"{g('p_mediana').mean():.4f} | "
          f"{g('ks_uniforme').mean():.4f} ({g('ks_uniforme').std(ddof=1):.4f}) | "
          f"{g('p_min').min():.3e} | {g('abs_lfc_mediana').mean():.4f} | "
          f"{g('abs_lfc_max').max():.4f} | {int(g('n_lfc_nan').sum())} |")
    W("")
    W("---")
    W("")
    return L


if __name__ == "__main__":
    L = main()
    open("/tmp/claude-1000/-home-ubuntu-vcc/52816e9e-e42f-4af3-b6ec-69a5b07f6ce9/scratchpad/rep_d.txt",
         "w").write("\n".join(L))
    print("ok, lineas:", len(L))
