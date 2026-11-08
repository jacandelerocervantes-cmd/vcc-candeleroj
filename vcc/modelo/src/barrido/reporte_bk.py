#!/usr/bin/env python
"""Escribe modelo/out/barrido_k.md desde barrido_k.json. Reejecutable: cada vez
que cae un brazo, se vuelve a correr y la tabla crece. NO concluye nada."""
import json, os, time
import numpy as np

JSON = "/home/ubuntu/vcc/modelo/out/barrido_k.json"
AJ = "/home/ubuntu/vcc/modelo/out/entrega_ajuste.json"
MD = "/home/ubuntu/vcc/modelo/out/barrido_k.md"
MET = ["pds_cosine", "expr_mse_unbiased_capped_norm", "de_wilcoxon_lfc_nmae",
       "de_wilcoxon_direction_fidelity_yield_raw", "de_wilcoxon_direction_reach_raw",
       "de_wilcoxon_sig_jaccard"]
ORDEN = ["k1.0", "k3.4", "r8_k1.0", "k2.0", "k1.5", "k2.5"]
ETQ = {"k1.0": "k=1.0", "k3.4": "k=3.4", "r8_k1.0": "8 rasgos, k=1.0",
       "k2.0": "k=2.0", "k1.5": "k=1.5", "k2.5": "k=2.5"}

# medido en sesion, bloques 0a/0b/0c
CLI = {"antes": "vcc 0.3.0", "despues": "vcc 0.3.0",
       "deps": "pandas 3.1.0→3.1.1, pydantic 2.14.1→2.14.3, zarr 3.4.1→3.4.2, "
               "scverse-misc 0.1.7→0.1.8, tzdata 2026.5→2026.6, python-dateutil 3.0.0→4.0.0 (mayor), six retirado",
       "revalidacion": "✓ dry run; 360,000 celulas; 18,533 genes; D/E/F 120,000 cada uno; "
                       "targets: verified against the official list"}
RESPALDO = {"ruta": "/datos/modelo/out/entrega/prediction.vcc", "bytes": 3186462720,
            "sha256": "66149bf6794c9f5d8cc062599eff9856451a62456e5a3f5ad710749656c5bff1",
            "mtime": "2026-11-05 19:12:05 UTC"}
PAQUETE = {"partition": "test", "panel_id": "vcc2026-test-1", "contextos": "D, E, F",
           "pert_counts.csv": "301 lineas (cabecera + 300)",
           "gene_names.csv": "18,534 lineas (cabecera + 18,533)",
           "val_hermano": '/datos/reto/val/ con "partition": "val", "panel_id": "vcc2026-val-1", contextos A/B/C'}
AVISOS = [
 "`de_lfc_nmae: omitted 1 perturbation(s) for an empty gate (e.g. PCCB). The gate is real-side only, so this set is the same for every submission.`",
 "`de_lfc_nmae: omitted 11 perturbation(s) for fewer than 10 gated genes (e.g. CD276). The gate is real-side only, so this set is the same for every submission.`",
 "`pred DE table omits 130 of the 76015 reference-significant (target, feature) pairs (0.171%), across 30 of 60 target(s); worst is 'DNMT1' with 12.` El marcador anade que el efecto NO es de un solo signo (su issue #291).",
]

def f(x, n=4):
    return "n/d" if x is None else f"{x:+.{n}f}" if isinstance(x, float) else str(x)

def tabla(cab, filas):
    s = "| " + " | ".join(cab) + " |\n|" + "|".join(["---"] * len(cab)) + "|\n"
    for r in filas:
        s += "| " + " | ".join(str(x) for x in r) + " |\n"
    return s + "\n"

def main():
    r = json.load(open(JSON))
    ok = [n for n in ORDEN if r["brazos"].get(n, {}).get("estado") == "OK"]
    otros = {n: v.get("estado") for n, v in r["brazos"].items() if v.get("estado") != "OK"}
    B = r["brazos"]
    L = []
    A = L.append
    A("# El barrido de k sobre VCC 2025, blancos retenidos\n")
    A("## 🔴 CÓMO SE LEE\n")
    A("Este documento mide **la sensibilidad del marcador al factor de des-encogimiento "
      "k** sobre una sub-superficie de VCC 2025 con verdad, en blancos que el modelo "
      "NO vio. No es una prediccion del score real de la entrega: la linea es otra, el "
      "panel es otro (18,080 genes de 2025 contra 18,533 del reto) y el control es otro. "
      "Lo unico que podria cruzar de aqui a una reemision es **k, un numero**.\n")
    A("🔴 **No se ha enviado nada** y no se ha tocado D/E/F: alli no hay verdad y no se "
      "puede puntuar nada. El `.vcc` de respaldo esta intacto (seccion 0.3).\n")
    A("🔴 **Sin veredicto.** No se elige k aqui. Las tablas se presentan con sus "
      "controles y sus defectos.\n")
    A(f"_Brazos cerrados: {len(ok)} de 6 — {', '.join(ETQ[n] for n in ok) or '(ninguno)'}._ "
      + (f"_Otros: {otros}._" if otros else "") + "\n")

    A("---\n\n## 0. Verificaciones previas\n")
    A("### 0.1 La CLI al dia\n")
    A(tabla(["", "version"], [["antes de `uv tool upgrade vcc-cli`", f"`{CLI['antes']}`"],
                              ["despues", f"`{CLI['despues']}` — la misma"]]))
    A(f"La version no cambio; el arbol de dependencias si: {CLI['deps']}. Por el salto "
      f"mayor de `python-dateutil` se revalido el `.h5ad` ya empaquetado con la CLI "
      f"nueva: {CLI['revalidacion']}\n")
    A("### 0.2 El paquete es el FINAL, no el de validacion\n")
    A(tabla(["campo", "valor"], [[f"`{k}`", f"`{v}`"] for k, v in PAQUETE.items()]))
    A("### 0.3 El respaldo, intocable\n")
    A(tabla(["", ""], [["ruta", f"`{RESPALDO['ruta']}`"],
                       ["bytes", f"**{RESPALDO['bytes']:,}**"],
                       ["sha256", f"`{RESPALDO['sha256']}`"],
                       ["mtime", RESPALDO["mtime"]]]))
    A("`vcc prep` corrio **en firme** (`-f`, sin `--dry-run`): el log del paso 5 dice "
      "`✓ prep complete` con `output: …/prediction.vcc`, frente al `✓ dry run — no .vcc "
      "written` del paso 4. El sha256 se verifico DESPUES de la revalidacion de 0.1: "
      "identico. No esta bajo `~/vcc/out` (vive en `/datos`, otro sistema de archivos, "
      "y no es enlace). `VCC_TOKEN` no esta definido en el entorno ni su nombre aparece "
      "en ningun archivo de `~/vcc` ni en los dotfiles. Sin commits: `~/vcc` no es repo git.\n")

    A("---\n\n## 1. El banco, declarado antes de correr\n")
    A(f"{r['que_es']}\n")
    lr = r.get("lado_real", {})
    A(tabla(["que", "valor"], [
        ["donde", "VCC 2025, pliegue 0 de la MISMA particion por blanco de la CV"],
        ["particion", f"leida de `entrega_ajuste.json` (semilla {r['semillas']['semillas'].get('baraja', 20284950)}), no reconstruida"],
        ["blancos retenidos", r.get("n_retenidos_pliegue")],
        ["reajuste", "sobre los otros 4 pliegues: ningun blanco puntuado estuvo en su ajuste"],
        ["alpha", r.get("alpha")],
        ["celulas por etiqueta", r["n_por_etiqueta"]],
        ["control del banco", f"{lr.get('n_control', '?')} de {lr.get('n_control_disponible', '?')} disponibles en 2025"],
        ["regimen de DE", "400 contra 9,200 — el que `ntc_nulo` declaro como el de entrega"],
        ["forma del lado real", lr.get("shape")],
        ["nnz del lado real", f"{lr.get('nnz', 0):,}"],
        ["recorte en log2", f"±{r['recorte_log2']} aplicado DESPUES de multiplicar por k"],
        ["sin modificar", "`suelo_arnes`, `cand01_covar._reescala_entero` (via `sup.aplica`), marcador oficial preset/profile `vcc2026`"],
        ["orden declarado", ", ".join(ETQ[n] for n in r["orden_declarado"])],
        ["por que ese orden", "extremos y ablacion primero, para que una tabla parcial abarque el rango"]]))
    hb = r.get("hash_bundle", {})
    if hb.get("coincide"):
        A(f"**El bundle reproduce el `config_hash` de competicion**: `{hb['obtenido']}` "
          f"= referencia. Las seis estan medidas con las reglas oficiales.\n")
    elif hb:
        A(f"🔴 hash del bundle: {hb}\n")
    A("### 1.1 Semillas: familia 20284000+, comprobadas EN CODIGO y guardadas enteras\n")
    s = r["semillas"]
    A(tabla(["uso", "semilla"], [[f"`{k}`", v] for k, v in sorted(s["semillas"].items(), key=lambda x: x[1])]))
    A(f"{s['n_nuevas']} semillas nuevas contra **{s['n_gastadas_enumeradas']:,}** "
      f"enumeradas (rangos vetados mas las de la entrega): **choques {s['choques'] or 'ninguno'}**.\n")
    A("### 1.2 Lo que el marcador avisa de ESTE lado real\n")
    for a in AVISOS:
        A(f"- {a}\n")
    A("\nEscala del lado real: **76,015 pares (blanco, gen) significativos sobre 60 "
      "blancos**, unos **1,267 por blanco** segun el aviso; la media medida en la "
      "seccion 6 es 1,950 porque alli la puerta es `p_adj` sin el recorte de "
      "`N_conf`. Son dos conteos distintos del mismo lado real, no una contradiccion.\n")

    A("---\n\n## 2. 1-bis — los rasgos de contexto y la ablacion de 8 rasgos\n")
    rc = r["rasgos_de_contexto"]
    A("🔴 **Correccion de premisa, verificada en el codigo.** `profundidad_mediana` y "
      "`n_control` **no estan entre los 10 rasgos**. La lista es:\n\n")
    A("`" + "`, `".join(rc["rasgos_del_modelo"]) + "`\n\n")
    A(f"{rc['nota']} Es decir: **el modelo empaquetado ya es el que no los usa**, y "
      "quitarlos de 10 daria los mismos 10.\n")
    A(tabla(["cantidad", "valor"], [
        ["sd de `profundidad_mediana` sobre las filas del ajuste", f"{rc['sd_profundidad_sobre_filas_del_ajuste']:.1f} (escalar unico)"],
        ["sd de `n_control` sobre las filas del ajuste", f"{rc['sd_n_control_sobre_filas_del_ajuste']:.1f} (escalar unico)"],
        ["n de control de 2025", f"{rc['n_control_2025']:,}"]]))
    ab = r["ablacion_8"]
    A(f"La ablacion que si se puede hacer, **criterio declarado antes de ajustar**: "
      f"{ab['criterio_declarado_antes_de_ajustar']}.\n\n")
    A(tabla(["", "rasgos"], [["fuera (2)", "`" + "`, `".join(ab["fuera"]) + "`"],
                             ["dentro (8)", "`" + "`, `".join(ab["dentro"]) + "`"]]))
    A("### 2.1 La CV partida por blanco, los dos modelos\n")
    cv = r.get("cv_por_modelo", {})
    if cv:
        for nom, lab in (("pearson_por_blanco", "Pearson por blanco retenido"),
                         ("R2_por_blanco", "R² por blanco retenido")):
            filas = []
            for campo in ("mediana", "media", "p10_p90", "min_max", "frac_mayor_0", "frac_mayor_0.1", "n"):
                fila = [campo]
                for m in ("10_rasgos", "8_rasgos"):
                    v = cv.get(m, {}).get(nom, {}).get(campo)
                    if isinstance(v, list):
                        fila.append(f"[{v[0]:+.4f}, {v[1]:+.4f}]")
                    elif isinstance(v, float):
                        fila.append(f"{v:+.4f}")
                    else:
                        fila.append(str(v))
                filas.append(fila)
            A(f"**{lab}** (los 300 blancos, cada uno predicho por un modelo que no lo vio):\n\n")
            A(tabla(["cantidad", "10 rasgos (empaquetado)", "8 rasgos"], filas))
        A("🔴 **El R² por blanco es negativo en los dos modelos**, con el Pearson "
          "positivo: la funcion acierta algo de la FORMA del perfil y falla la ESCALA. "
          "Es la cifra que motiva este barrido. En el informe de la entrega el R² "
          "publicado era el global (+0.0052), no el de cada blanco.\n")
    A("### 2.2 Intermedio crudo: 20 blancos retenidos, los DOS modelos lado a lado\n")
    crudo10 = cv.get("10_rasgos", {}).get("por_blanco_crudo", {})
    crudo8 = cv.get("8_rasgos", {}).get("por_blanco_crudo", {})
    if crudo10 and crudo8:
        pares = sorted(((t, v[0], v[1]) for t, v in crudo10.items() if v[0] is not None),
                       key=lambda x: x[1])
        idx = np.linspace(0, len(pares) - 1, 20).astype(int)
        filas = []
        for i, j in enumerate(idx):
            t, p10, r10 = pares[j]
            p8, r8 = crudo8.get(t, [None, None])
            filas.append([i + 1, f"`{t}`", f"{p10:+.4f}",
                          "n/d" if p8 is None else f"{p8:+.4f}",
                          f"{r10:+.4f}", "n/d" if r8 is None else f"{r8:+.4f}"])
        A(tabla(["#", "blanco", "Pearson 10r", "Pearson 8r", "R² 10r", "R² 8r"], filas))
        A(f"_De los {len(pares)} blancos retenidos, equiespaciados por su Pearson del "
          "modelo de 10 rasgos. Cada blanco predicho por un modelo que no lo vio._\n")
        d = [abs(crudo10[t][0] - crudo8[t][0]) for t in crudo10
             if crudo10[t][0] is not None and crudo8.get(t, [None])[0] is not None]
        A(f"\nDiferencia |Pearson 10r − Pearson 8r| por blanco: mediana "
          f"{np.median(d):.5f}, maxima {max(d):.5f}, sobre {len(d)} blancos. "
          "Es el tamano del efecto de quitar los dos rasgos del blanco g, blanco a blanco.\n")
    else:
        A("_(crudo no disponible)_\n")
    A("---\n\n## 3. Las seis metricas por brazo — `u` (fila `mean` de `agg_results.csv`)\n")
    if ok:
        A(tabla(["metrica"] + [ETQ[n] for n in ok],
                [[f"`{m}`"] + [f"{B[n]['u'][m]:.4f}" for n in ok] for m in MET]))
        A("## 4. `s` — normalizado contra el bundle real\n")
        A("**`from_replicate`:**\n\n")
        A(tabla(["metrica"] + [ETQ[n] for n in ok],
                [[f"`{m}`"] + [f((B[n].get('s') or {}).get(m, {}).get("from_replicate"))
                               for n in ok] for m in MET + ["avg_score"]]))
        A("**`from_baseline`:**\n\n")
        A(tabla(["metrica"] + [ETQ[n] for n in ok],
                [[f"`{m}`"] + [f((B[n].get('s') or {}).get(m, {}).get("from_baseline"))
                               for n in ok] for m in MET + ["avg_score"]]))
        A("## 5. La emision de cada brazo\n")
        filas = [["modelo"] + [B[n]["modelo"] for n in ok],
                 ["k"] + [f"{B[n]['k']:.1f}" for n in ok],
                 ["sd de la razon en log2"] + [f"{B[n]['sd_log2_razon']:.4f}" for n in ok],
                 ["sd de la razon"] + [f"{B[n]['sd_razon']:.4f}" for n in ok],
                 ["razon [min, max]"] + [f"[{B[n]['razon_min_max'][0]:.3f}, {B[n]['razon_min_max'][1]:.3f}]" for n in ok],
                 ["recortados a ±2"] + [B[n]["n_recortados"] for n in ok],
                 ["semilla"] + [B[n]["semilla"] for n in ok],
                 ["segundos"] + [B[n]["segundos"] for n in ok]]
        A(tabla([""] + [ETQ[n] for n in ok], filas))
        A("_Referencia: la sd de la `y` observada en 2025 es 0.1772; la de la entrega "
          "empaquetada, 0.0138 en muestra y 0.052 por contexto._\n")

        A("---\n\n## 6. 🔴 La tabla de cuatro celdas, nunca en una sola direccion\n")
        A("Puerta: `p_adj` de los `de_real.parquet` / `de_pred.parquet` que escribe "
          "`run --write-degenes`, o sea la del marcador, con el gen propio fuera.\n\n")
        cc = [("real ∩ pred", "media_interseccion"), ("solo real", "media_solo_real"),
              ("solo pred (falsos positivos)", "media_solo_pred_FALSOS_POSITIVOS"),
              ("ambos vacios (de 60 blancos)", "ambos_vacios"),
              ("genes llamados, pred", "media_genes_llamados_pred"),
              ("genes llamados, real", "media_genes_llamados_real")]
        A(tabla(["media por blanco"] + [ETQ[n] for n in ok],
                [[lab] + [f"{B[n]['cuatro_celdas'][key]:,.3f}" for n in ok] for lab, key in cc]))
        # precision y cobertura, derivadas
        filas = []
        for lab, fn in (("precision = ∩ / llamados por pred",
                         lambda c: c["media_interseccion"] / c["media_genes_llamados_pred"]
                         if c["media_genes_llamados_pred"] > 0 else None),
                        ("cobertura = ∩ / llamados por real",
                         lambda c: c["media_interseccion"] / c["media_genes_llamados_real"]
                         if c["media_genes_llamados_real"] > 0 else None)):
            fila = [lab]
            for n in ok:
                v = fn(B[n]["cuatro_celdas"])
                fila.append("n/d" if v is None else f"{100*v:.1f}%")
            filas.append(fila)
        A(tabla([""] + [ETQ[n] for n in ok], filas))
        A("**Y la descomposicion del Jaccard**, para no confundir solape con vacio:\n\n")
        filas = [["`de_wilcoxon_sig_jaccard` medido"] + [f"{B[n]['u']['de_wilcoxon_sig_jaccard']:.4f}" for n in ok],
                 ["blancos con ambos conjuntos vacios"] + [f"{int(B[n]['cuatro_celdas']['ambos_vacios'])}/60" for n in ok],
                 ["lo que aportan esos vacios (1.0 cada uno)"] +
                 [f"{B[n]['cuatro_celdas']['ambos_vacios']/60:.4f}" for n in ok]]
        A(tabla([""] + [ETQ[n] for n in ok], filas))
        for n in ok:
            c = B[n]["cuatro_celdas"]
            pb = sorted(c["por_blanco"], key=lambda x: -x["interseccion"])
            idx = np.linspace(0, len(pb) - 1, min(20, len(pb))).astype(int)
            A(f"### 6.{ok.index(n)+1} {ETQ[n]} — 20 blancos crudos, de mayor a menor interseccion\n")
            A(tabla(["#", "blanco", "n real", "n pred", "∩", "solo real", "solo pred"],
                    [[i + 1, f"`{pb[j]['blanco']}`", f"{pb[j]['n_real']:,}", pb[j]["n_pred"],
                      pb[j]["interseccion"], f"{pb[j]['solo_real']:,}", pb[j]["solo_pred"]]
                     for i, j in enumerate(idx)]))
    else:
        A("_Aun no hay brazos cerrados._\n")

    A("---\n\n## 7. Defectos propios\n")
    A("1. **Es una sub-superficie de 2025, no la de entrega.** 60 blancos de 300, "
      "control submuestreado a 9,200 de 114,528, panel de 18,080 genes. La eleccion "
      "del submuestreo fue por coste del bundle (39.6 min medidos), no por biologia.\n")
    A("2. **Un solo pliegue.** Los seis brazos se puntuan sobre los mismos 60 blancos "
      "retenidos; no hay repeticion por pliegue, asi que no hay barra de error entre "
      "pliegues para las seis metricas, solo la de la CV de la seccion 2.\n")
    A("3. **Una sola semilla de emision por brazo.** La variabilidad de "
      "`_reescala_entero` no esta medida: brazos distintos usan semillas distintas, asi "
      "que parte de la diferencia entre columnas es ruido de emision y no efecto de k.\n")
    A("4. **12 de los 60 blancos no entran en `lfc_nmae`** (1 por puerta vacia, 11 por "
      "menos de 10 genes con puerta), y el conjunto es el mismo para todos los brazos "
      "porque la puerta es del lado real.\n")
    A("5. **El aviso del issue #291** afecta a las metricas de direccion de todos los "
      "brazos por igual y su efecto no es de un solo signo.\n")
    A("6. **k es un parametro post-hoc.** Multiplicar el delta no cambia la funcion ni "
      "su R²; mueve la escala de la emision. Nada de lo medido aqui dice que el k que "
      "mejor puntue en 2025 sea el que mejor puntue en D/E/F — son lineas distintas.\n")
    A(f"\n_Generado {time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime())} UTC desde "
      f"`barrido_k.json`. Reejecutable: `reporte_bk.py`._\n")

    open(MD, "w").write("".join(L))
    print(f"escrito {MD}: {len(''.join(L).splitlines())} lineas, {len(ok)} brazos")

if __name__ == "__main__":
    main()
