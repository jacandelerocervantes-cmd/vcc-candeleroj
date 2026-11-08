#!/usr/bin/env python
"""Escribe modelo/out/entrega.md de lo que haya en disco."""
import glob
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, "/home/ubuntu/vcc/modelo/src/noche")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import comun as C  # noqa: E402
import datos as D  # noqa: E402

OUT = "/home/ubuntu/vcc/modelo/out/entrega.md"
AJ = "/home/ubuntu/vcc/modelo/out/entrega_ajuste.json"
CT = "/home/ubuntu/vcc/modelo/out/entrega_contrato.json"
PR = "/datos/modelo/out/entrega/progreso.json"
LOGS = "/home/ubuntu/vcc/modelo/out"


def j(p):
    return json.load(open(p)) if os.path.exists(p) else None


def main():
    aj, ct, pr = j(AJ), j(CT), j(PR)
    W = ["# La entrega: la funcion, la emision y el contrato", "",
         "## 🔴 Que es esto y que no", "",
         "Un solo `.vcc` empaquetado de un solo `.h5ad`: 300 perturbaciones x 400 "
         "celulas x 3 contextos (D, E, F), **sin ninguna fila de control** — las "
         "celulas de control son la ENTRADA del modelo, no su salida.", "",
         "🔴 **La verdad del gen no viaja.** De VCC 2025 sale **solo la funcion** (los "
         "coeficientes de la seccion 1): ningun delta de 2025 se copia a la emision, "
         "tampoco el de los 28 blancos que 2025 y el panel 2026 comparten. Son lineas "
         "distintas; esos 28 pasan por la misma funcion que los otros 272.", "",
         "🔴 **No se ha enviado nada.** No se ha usado `vcc sample`, ni `vcc submit`, "
         "ni `~/vcc/benchmark/`, ni `/datos/benchmark_VEDADO_*`. No hay commits.", ""]
    if aj:
        W += ["---", "", "## 1. PASO 1 — la funcion aprendida en VCC 2025", "",
              "### 1.1 El objetivo y el modelo, declarados antes de ajustar", ""]
        W += C.tabla(["que", "valor"], [
            ["objetivo", "`y(g,h) = log2((cpm_pert + 1)/(cpm_ctrl + 1))`, cpm = media "
             "sobre celulas de `x/L*1e6`"],
            ["modelo", aj.get("modelo")],
            ["rasgos", ", ".join(f"`{r}`" for r in aj.get("rasgos", []))],
            ["rejilla de `alpha`", ", ".join(str(a) for a in aj.get("alphas", []))],
            ["blancos de 2025", aj.get("n_perts_2025")],
            ["genes de 2025", aj.get("n_genes_2025")],
            ["genes que pasan el filtro CPM>5 del control de 2025", aj.get("n_genes_filtro_2025")],
            ["forma del diseno", aj.get("forma_diseno")],
            ["blancos de 2025 con perfil gwps", aj.get("blancos_2025_con_gwps")],
            ["genes de 2025 con perfil gwps", aj.get("genes_2025_con_gwps")],
            ["celulas de 2025 leidas", f"{(aj.get('pseudobulk_2025') or {}).get('n_celulas', 0):,}"],
            ["celulas de control de 2025", f"{(aj.get('pseudobulk_2025') or {}).get('n_control_2025', 0):,}"]])
        W += ["", "### 1.2 🔴 La validacion cruzada PARTIDA POR BLANCO", "",
              f"5 pliegues, permutacion de los blancos con la semilla "
              f"{(aj.get('cv') or {}).get('semilla_baraja')}. Ningun blanco esta a la vez "
              "en ajuste y en prueba. Se elige `alpha` por **la mediana de la "
              "correlacion de Pearson POR BLANCO retenido**.", ""]
        fil = []
        for a, v in (aj.get("cv", {}).get("por_alpha") or {}).items():
            fil.append([f"`{a}`", C.f4(v["pearson_medio"]), C.f4(v["R2_medio"]),
                        f"**{C.f4(v['pearson_por_blanco_mediana'])}**",
                        "**elegida**" if float(a) == aj.get("alpha_elegida") else ""])
        W += C.tabla(["alpha", "Pearson (todas las filas)", "R²",
                      "Pearson por blanco (mediana)", ""], fil)
        a_el = str(aj.get("alpha_elegida"))
        pl = (aj.get("cv", {}).get("por_alpha") or {}).get(a_el, {}).get("pliegues", [])
        if pl:
            W += ["", f"**Los 5 pliegues de la `alpha` elegida ({a_el}):**", ""]
            W += C.tabla(["pliegue", "n blancos", "Pearson", "Spearman", "R²", "RMSE",
                          "Pearson por blanco: mediana", "[min, max]"],
                         [[p["pliegue"], p["n_blancos"], C.f4(p["pearson"]),
                           C.f4(p["spearman"]), C.f4(p["R2"]), C.f4a(p["rmse"]),
                           C.f4(p["pearson_por_blanco_mediana"]),
                           f"[{C.f4(p['pearson_por_blanco_rango'][0])}, "
                           f"{C.f4(p['pearson_por_blanco_rango'][1])}]"] for p in pl])
            todos = [x for p in pl for x in p["pearson_por_blanco"] if x is not None]
            if todos:
                t = np.asarray(todos)
                W += ["", "🔴 **La distribucion de la correlacion por blanco retenido** "
                      f"(los {t.size} blancos, cada uno predicho por un modelo que no lo "
                      "vio):", ""]
                W += C.tabla(["cantidad", "valor"], [
                    ["mediana", C.f4(float(np.median(t)))],
                    ["media", C.f4(float(t.mean()))],
                    ["[p10, p90]", f"[{C.f4(float(np.percentile(t,10)))}, "
                     f"{C.f4(float(np.percentile(t,90)))}]"],
                    ["[min, max]", f"[{C.f4(float(t.min()))}, {C.f4(float(t.max()))}]"],
                    ["fraccion > 0", C.f4a(float((t > 0).mean()))],
                    ["fraccion > 0.1", C.f4a(float((t > 0.1).mean()))]])
                sel = np.argsort(t)[np.linspace(0, t.size - 1, 20).round().astype(int)]
                nombres = [b for p, pr_ in zip(pl, aj["cv"]["blancos_de_cada_pliegue"])
                           for b in pr_]
                W += ["", "**20 blancos retenidos, equiespaciados por su correlacion** "
                      "(el intermedio crudo de la cifra de arriba):", ""]
                W += C.tabla(["#", "blanco", "Pearson del blanco retenido"],
                             [[i + 1, f"`{nombres[k]}`" if k < len(nombres) else "n/d",
                               C.f4(float(t[k]))] for i, k in enumerate(sel)])
        ctr = aj.get("controles") or {}
        W += ["", "### 1.3 Los dos controles declarados", ""]
        W += C.tabla(["control", "Pearson", "R²"], [
            ["(a) predecir 0 en todo gen", "n/d (constante)",
             C.f4(ctr.get("a_predecir_cero", {}).get("R2"))],
            ["(b) solo el perfil de gwps, una pendiente",
             C.f4(ctr.get("b_solo_perfil_gwps", {}).get("pearson_medio")),
             C.f4(ctr.get("b_solo_perfil_gwps", {}).get("R2_medio"))],
            ["el modelo de 10 rasgos",
             C.f4((aj.get("cv", {}).get("por_alpha") or {}).get(a_el, {}).get("pearson_medio")),
             C.f4((aj.get("cv", {}).get("por_alpha") or {}).get(a_el, {}).get("R2_medio"))]])
        co = aj.get("coeficientes") or {}
        if co:
            W += ["", "### 1.4 Los coeficientes (rasgos estandarizados)", "",
                  f"Intercepto `b0` = {co['b0']:+.6f}, `alpha` = {co['alpha']}.", ""]
            W += C.tabla(["rasgo", "peso (sobre rasgo estandarizado)", "media", "sd"],
                         [[f"`{r}`", C.f4(co["w"][r]), C.fe(co["mu"][r]), C.fe(co["sd"][r])]
                          for r in aj["rasgos"]])
            af = aj.get("ajuste_final_en_muestra") or {}
            W += ["", "**La escala de lo predicho contra lo observado** (en muestra, "
                  "sobre las filas del ajuste):", ""]
            W += C.tabla(["cantidad", "y observada", "y predicha"], [
                ["sd", C.f4(af.get("y_observada", {}).get("sd")),
                 C.f4(af.get("y_predicha", {}).get("sd"))],
                ["p1", C.f4(af.get("y_observada", {}).get("p1")),
                 C.f4(af.get("y_predicha", {}).get("p1"))],
                ["p50", C.f4(af.get("y_observada", {}).get("p50")),
                 C.f4(af.get("y_predicha", {}).get("p50"))],
                ["p99", C.f4(af.get("y_observada", {}).get("p99")),
                 C.f4(af.get("y_predicha", {}).get("p99"))]])
            W += ["", "🔴 La sd de lo predicho es mucho menor que la de lo observado: una "
                  "ridge encoge, y el encogimiento es parte del modelo declarado, no un "
                  "ajuste posterior.", ""]
    if pr:
        W += ["---", "", "## 2. PASO 2 — la emision sobre D/E/F", ""]
        W += C.tabla(["que", "valor"], [
            ["filas escritas", f"{pr.get('n_filas_escritas', 0):,} de "
             f"{pr.get('forma_objetivo', [0])[0]:,}"],
            ["nnz", f"{pr.get('nnz', 0):,}"],
            ["recorte de `y` en log2", f"±{pr.get('recorte_log2')}"],
            ["blancos sin perfil gwps (repliegue declarado)",
             f"{len(pr.get('sin_gwps', []))}: " + ", ".join(f"`{g}`" for g in pr.get("sin_gwps", []))],
            ["tamano del h5ad", f"{pr.get('tam_bytes', 0) / 2**30:.2f} GiB"
             if pr.get("tam_bytes") else "aun escribiendose"]])
        pc = pr.get("por_contexto") or {}
        if pc:
            W += ["", "**Por contexto** — los rasgos del gen se recalculan con el control "
                  "de ESE contexto:", ""]
            W += C.tabla(["ctx", "n control", "profundidad mediana", "genes que pasan CPM>5",
                          "sd de y", "y p1", "y p50", "y p99", "razon [min, max]",
                          "recortados", "semilla control", "semilla ruido"],
                         [[k, f"{v['n_control']:,}", f"{v['profundidad_mediana']:,.0f}",
                           f"{v['n_genes_filtro']:,}", C.f4(v["y_sd"]), C.f4(v["y_p1"]),
                           C.f4(v["y_p50"]), C.f4(v["y_p99"]),
                           f"[{v['razon_min']:.3f}, {v['razon_max']:.3f}]",
                           f"{v['n_recortados']:,} ({v['frac_recortados']:.2e})",
                           v["semilla_control"], v["semilla_ruido"]]
                          for k, v in pc.items()])
    if ct:
        W += ["", "---", "", "## 3. PASO 3 — el contrato, comprobado por mi antes de "
              "empaquetar", "",
              f"**{ct.get('n_ok')}/{ct.get('n_total')} comprobaciones OK.** "
              f"Archivo {ct.get('tam_GiB')} GiB, nnz {ct.get('nnz', 0):,} "
              f"(margen x{ct.get('margen_nnz')} contra el tope de 4,750,000,000), "
              f"densidad {ct.get('densidad')}, {ct.get('nnz_por_celula')} genes no nulos "
              f"por celula, suma media por celula {ct.get('suma_mediana_aprox')}.", ""]
        W += C.tabla(["#", "comprobacion", "medido", "esperado", "ok"],
                     [[c["#"], c["que"], str(c["medido"]), str(c["esperado"]),
                       "**si**" if c["ok"] else "**NO**"] for c in ct["comprobaciones"]])
    for nom, f in (("`vcc prep --dry-run`", os.path.join(LOGS, "entrega_prep_dryrun.log")),
                   ("`vcc prep`", os.path.join(LOGS, "entrega_prep.log"))):
        if os.path.exists(f):
            t = open(f).read().strip().splitlines()
            W += ["", f"### {nom}", "", "```", *t[-30:], "```", ""]
    W += ["", "---", "", "## 4. Defectos de la construccion", "",
          "1. **Una sola linea de entrenamiento.** 2025 es un contexto; los rasgos de "
          "CONTEXTO (profundidad, n de control) no son identificables ahi y por eso no "
          "entran en el ajuste. Lo que se transfiere a D/E/F es una funcion ajustada en "
          "un solo contexto.",
          "2. **El perfil de gwps es de K562 y cubre 7,681 de los 18,533 genes** y 277 de "
          "las 300 perturbaciones. Para el resto, el rasgo entra en su valor de "
          "no-informacion y la prediccion se queda con la parte comun.",
          "3. **La ridge encoge.** La sd de `y` predicha es una fraccion de la observada: "
          "la emision es conservadora por construccion.",
          "4. **El emisor desplaza medias, no dispersiones.** `_reescala_entero` aplica "
          "Binomial/Poisson sobre una celula de control; `sobredispersion.md` ya midio "
          "que ese dispositivo no reproduce la sobredispersion del delta real.",
          "5. **El objetivo es una razon de CPM medio, no un estadistico de DE.** Es lo "
          "que el emisor puede aplicar; las metricas de DE del marcador miran otra cosa "
          "(Wilcoxon por gen), y la relacion entre ambas no se mide aqui.",
          "6. **Las 400 celulas salen del propio pool de control con reemplazo**, asi que "
          "la variabilidad celula a celula de la emision es la del control mas el ruido "
          "de la inyeccion.", ""]
    W += [f"_Generado {time.strftime('%Y-%m-%d %H:%M:%S')} UTC._"]
    C.escribe_md(OUT, W)
    print(f"escrito {OUT}: {len(W)} lineas")


if __name__ == "__main__":
    main()
