#!/usr/bin/env python
"""Reporte `magnitud_centro.md` -- encabezado y construccion."""
import json, os
OUTDIR = "/home/ubuntu/vcc/modelo/out/magnitud_centro"
J = lambda n: json.load(open(os.path.join(OUTDIR, n)))


def main():
    sem = J("semillas.json"); met = J("deltas_meta.json")
    L = []; W = L.append
    esp = met["espacio"]

    W("# La descomposicion de la magnitud entre contextos: que fraccion es CENTRO")
    W("")
    W("## 🔴 0. Como se lee esta salida, y de que habla")
    W("")
    W("Esto PRESENTA DATOS. No lleva veredicto, no concluye, y no va dirigido a")
    W("ningun candidato siguiente.")
    W("")
    W("🔴 **De que contextos habla.** De **los seis de Jiang**, que son LINEAS")
    W("CELULARES (A549, BXPC3, HAP1, HT29, K562, MCF7) bajo IFNB. **No son los")
    W("contextos de entrega.** Y los de entrega son, ademas, contextos DISTINTOS de")
    W("los de validacion: `paso1_contrato.md` punto 37 midio que `val` (A/B/C) y")
    W("`final` (D/E/F) **no comparten el panel de 300** -- solapan 129 de 300, 171")
    W("exclusivas de cada lado. Asi que hay **dos saltos** entre lo que se mide aqui y")
    W("lo que se entrega: linea celular -> contexto del reto, y validacion -> test.")
    W("")
    W("🔴 **La salida NO es una fraccion que se transfiera.** Es un **SI O NO**")
    W("**ESTRUCTURAL** sobre si el centro domina, dado como **rango**. Una invariancia")
    W("demostrada cruza saltos; un numero no. Ninguna cifra de este reporte esta")
    W("puesta para ser reutilizada como constante.")
    W("")
    W("**Por que existe la medicion.** El Venn de los contextos: un metodo cae en el")
    W("centro D∩E∩F, en un par, o en un singlete. Todas las vias probadas emiten solo")
    W("el centro, por construccion. Esto mide si el centro vale la pena.")
    W("")
    W("🔴 **NO pasa por la puerta de expresion diferencial.** Es una medicion sobre la")
    W("MAGNITUD del delta, asi que la salvedad del motor (`pdex` contra `gpudge`) no la")
    W("toca: no hay Wilcoxon, no hay BH, no hay `p_adj` en ninguna cifra.")
    W("")
    W("**Lo que se leyo**: `coef_dentro_entre.md` (la maquinaria de desatenuacion),")
    W("`b_contra_a.md` §4-§5 (el control compartido y la identidad")
    W("`V_tot = V_sh + V_sp` por productos cruzados), `gwps_jiang_datos.md` (la forma")
    W("`G_mismo` / `G_otro`) y `pds_objetivo.md` TABLA B0 (el espacio del delta).")
    W("")
    W("**Lo que NO se uso**: nada de `~/vcc/benchmark/` ni de")
    W("`/datos/benchmark_VEDADO_*`. Sin envios, sin `vcc sample`, sin commits.")
    W("")
    W("---")
    W("")
    W("## 1. La construccion, declarada antes de medir")
    W("")
    W("### 1.1 Datos, panel y poblacion real")
    W("")
    W("Jiang IFNB, las seis lineas, el **panel COMUN de 33 perturbaciones** (100 % de")
    W("solape entre lineas), de `suelo_panel.json`. Poblacion real fija por linea con")
    W("`build_fold_adatas` y `seed_real = 20260001 + fold`: **400 celulas por**")
    W("**perturbacion** (sin reemplazo) y **TODO el pool de control** de la linea.")
    W("")
    W("**TABLA 1.1 — la poblacion real por linea, con su huella `sha256`.**")
    W("")
    W("| # | linea | fold | `seed_real` | celulas totales | celulas de control | "
      "celulas perturbadas | **huella `sha256` del objeto real** |")
    W("|---|---|---|---|---|---|---|---|")
    for i, l in enumerate(sorted(met["lineas"]), 1):
        d = met["lineas"][l]
        W(f"| {i} | **{l}** | {d['fold']} | `{d['seed_real']}` | "
          f"{int(d['n_celulas_real']):,} | {int(d['n_control_cells']):,} | "
          f"{int(d['n_celulas_real']) - int(d['n_control_cells']):,} | "
          f"`{d['huella_real_sha256']}` |")
    W("")
    W("### 1.2 El espacio en que se mide el delta")
    W("")
    W(f"**`{esp['comparator']}` con `bulk_target_sum = {esp['bulk_target_sum']}`**,")
    W("leido de `run_params.yaml` de las corridas ya hechas, via")
    W("`cell_eval2.prep.pseudobulk_bulk_lognorm`. Es el MISMO espacio de la TABLA B0 de")
    W("`pds_objetivo.md`. `bulk_lognorm_means` calcula `log1p(TS * P / sum_g P)` por")
    W("grupo, asi que depende **solo de la composicion** del pseudobulk y no de")
    W("cuantas celulas lo formen -- que es lo que hace comparable el delta de una")
    W("mitad de 200 celulas con el de las 400.")
    W("")
    W("| # | cantidad | valor |")
    W("|---|---|---|")
    W(f"| 1 | genes del archivo | {esp['n_genes_total']:,} |")
    W(f"| 2 | genes blanco del panel FUERA del espacio (`exclusion_scope=\"panel\"`) | "
      f"{esp['n_panel_fuera']} |")
    W(f"| 3 | **genes en el espacio** | **{esp['n_genes_espacio']:,}** |")
    W(f"| 4 | delta | `pb(g) - pb(control)` en ese espacio |")
    W(f"| 5 | centrado sobre los genes | {esp['centrado_sobre_genes']} "
      f"(como `b_contra_a.md` §5, asi que `<x,y>/G` es una covarianza) |")
    W("")
    W("🔴 **Si el espacio cambia, el reparto cambia.** Esta declarado por eso, no por")
    W("tramite.")
    W("")
    W("### 1.3 Las mitades, y el control PARTIDO (TRAMPA 1)")
    W("")
    W("* **Unidad: la celula.** Las 400 celulas perturbadas de `g` en `L` se permutan")
    W("  con la semilla de particion y se parten **200 / 200**. **400 es par, asi que**")
    W("  **aqui no hay desempate de impares.**")
    W("* 🔴 **El control va PARTIDO en dos mitades disjuntas, no compartido.** Es el")
    W("  confundidor que `b_contra_a.md` §4 destapo citando `anchor.py`: con un control")
    W("  comun, el ruido de muestreo del control entra identico en las dos cantidades")
    W("  cuyo acuerdo se mide y lo infla -- `ceiling.py` midio una fiabilidad de **0.54**")
    W("  convertida en **0.74** por ese error. Asi que")
    W("  `delta_g,L,A = pb(mitad_A) - pb(c1)` y `delta_g,L,B = pb(mitad_B) - pb(c2)`,")
    W("  con `c1` y `c2` disjuntas.")
    W("* **Desempate de impares del control, declarado y ALTERNADO**: si el pool de")
    W("  control tiene `n` impar, la celula sobrante va a `c1` cuando `j` es par y a")
    W("  `c2` cuando `j` es impar.")
    W("")
    W("**TABLA 1.2 — los cortes de control, por linea y particion.**")
    W("")
    W("| # | linea | control | `c1` / `c2` por `j` | desempate |")
    W("|---|---|---|---|---|")
    for i, l in enumerate(sorted(met["lineas"]), 1):
        d = met["lineas"][l]
        ps = d["particiones"]
        cs = ", ".join(f"{int(p['n_control_c1'])}/{int(p['n_control_c2'])}" for p in ps)
        des = ps[0]["desempate_impares"].split(":")[0]
        W(f"| {i} | **{l}** | {int(d['n_control_cells']):,} | {cs} | "
          f"{'par: 50/50 exacto' if 'par' == des else 'impar: alternado por paridad de `j`'} |")
    W("")
    W("🔴 **Cruce independiente de la construccion.** `b_contra_a.md` §4 reporto sus")
    W("cortes de control: A549 **689/690**, BXPC3 **1614/1614**, HAP1 **1586/1586**,")
    W("HT29 **1368/1368**, K562 **699/699**, MCF7 **1334/1335**. Los de la TABLA 1.2")
    W("coinciden, lo que confirma que se esta partiendo el mismo pool de control con el")
    W("mismo criterio.")
    W("")
    W("### 1.4 Las semillas: familia NUEVA, comprobada sin choque EN CODIGO")
    W("")
    W(f"Familia **`{sem['familia_nueva']}`** -> **{sem['n_semillas_nuevas']} semillas**,")
    W(f"rango `{sem['rango_nuevas'][0]}`..`{sem['rango_nuevas'][1]}`.")
    W("")
    W("El chequeo (`mc_semillas.py`) se hizo por **dos caminos**:")
    W("")
    W("| # | camino | resultado |")
    W("|---|---|---|")
    W(f"| 1 | expansion de los generadores de las {len(sem['bases_ya_usadas'])} bases ya "
      f"usadas (`base`, `base+fold`, `base+100*fold+j`, `base+replica`): "
      f"{sem['n_semillas_ya_usadas_expandidas']:,} semillas, "
      f"`{sem['rango_ya_usadas'][0]}`..`{sem['rango_ya_usadas'][1]}` | "
      f"choque = `{sem['choque_contra_expansion_de_generadores']}` |")
    W(f"| 2 | las {sem['n_semillas_registradas_en_los_json']} semillas REGISTRADAS en "
      f"los JSON de `modelo/out/` | choque = "
      f"`{sem['choque_contra_las_registradas']}` |")
    W("")
    W(f"**Sin choque: `{sem['sin_choque']}`.**")
    W("")
    W("🔴 **Un defecto de la PRIMERA version de este chequeo, declarado.** El camino 2")
    W("se escribio al principio con un regex `20\\d{6}` sobre el TEXTO de los JSON, y")
    W("eso dio un choque FALSO: `20278503` aparece dentro de los digitos del float")
    W("`-0.05420278503638758` de `control_o_no/pesos.json`. La version que corre")
    W("**parsea el JSON y recorre los valores**, contando solo enteros en")
    W("`[20200000, 20400000]`, y entonces no hay choque. Queda anotado porque el")
    W("chequeo de semillas es justo el sitio donde un falso positivo se cuela sin ruido.")
    W("")
    W("**Las 48 semillas de particion, escritas** (`semillas.json` las guarda, no solo")
    W("su conteo):")
    W("")
    W("| # | linea (fold) | las 8 semillas, `j = 0..7` |")
    W("|---|---|---|")
    for i, l in enumerate(sorted(met["lineas"], key=lambda x: met["lineas"][x]["fold"]), 1):
        d = met["lineas"][l]
        ss = ", ".join(f"`{int(p['semilla'])}`" for p in d["particiones"])
        W(f"| {i} | **{l}** ({d['fold']}) | {ss} |")
    W("")
    W(f"Y la del bootstrap sobre perturbaciones: **`{J('calc.json')['semilla_bootstrap_base']}`**,")
    W("de la misma familia nueva y fuera del rango de las 48 de particion.")
    W("")
    W("---")
    W("")
    return L


if __name__ == "__main__":
    open("/tmp/claude-1000/-home-ubuntu-vcc/52816e9e-e42f-4af3-b6ec-69a5b07f6ce9/scratchpad/mcrep_a.txt",
         "w").write("\n".join(main()))
    print("ok")
