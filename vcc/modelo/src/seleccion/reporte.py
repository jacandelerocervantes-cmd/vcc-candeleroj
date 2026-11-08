#!/usr/bin/env python
"""Ensambla modelo/out/seleccion_control.md y su json. Parte 1 y Parte 2 en
secciones separadas; ninguna tabla mezcla cifras de las dos."""
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import reporte_p1  # noqa: E402
import reporte_p2  # noqa: E402
import semillas as S  # noqa: E402
from reporte_p1 import tabla  # noqa: E402

OUT_MD = "/home/ubuntu/vcc/modelo/out/seleccion_control.md"
OUT_JSON = "/home/ubuntu/vcc/modelo/out/seleccion_control.json"


def main():
    W1, d1 = reporte_p1.main()
    W2, d2, m2 = reporte_p2.main()
    sem = S.comprueba()
    assert sem["sin_choques"]
    W = ["# PASO 5a — La seleccion: lo que se computa SOLO del control", "",
         "## 🔴 CÓMO SE LEE", "",
         "* **La Parte 1 es sobre la superficie de ENTREGA (D/E/F) y la de validacion (A/B/C).** Sus cifras son de "
         "ellas, computadas de sus 18,400 controles. El filtro del preset lo decide la referencia; en la puntuacion "
         "real la referencia es ese pool, asi que el filtro de cada contexto que se da aqui no es una aproximacion: "
         "es la misma cuenta (comprobada gen a gen contra `compute_de`, §1.1).",
         "* **La Parte 2 es sobre OTRAS superficies** (Jiang, 6 lineas; VCC 2025, un contexto), las que tienen verdad. "
         "De ella solo viaja el **MECANISMO** — que el filtro del contexto excluye genes de R(g), que la potencia "
         "selecciona — **nunca las fracciones**.",
         "* Las dos partes **no se mezclan en ninguna tabla**.",
         "* Formato: datos. Sin veredicto, sin conclusiones.", "",
         "## 0. Construccion, declarada antes de medir", ""]
    W += tabla(["# ", "cantidad", "valor"], [
        [1, "filtro", f"`filter_gene_min_cpm_cell = {d1['umbral']}` del preset `vcc2026`, leido de `EvalConfig`; pasa si CPM medio en la referencia `>` umbral"],
        [2, "referencia en Parte 1", "los 18,400 controles de cada contexto, todos"],
        [3, "comprobacion del filtro", "`compute_de` del preset, referencia = 18,400 controles, blanco = 400 controles; genes probados = genes de esta cuenta, gen a gen"],
        [4, "cantidades por gen (Parte 1)", "CPM medio, tasa de deteccion, Fano de cuentas crudas"],
        [5, "control de 1.4", "dos mitades de 9,200 del MISMO contexto: aleatoria (semilla) y por guia `ntc_id` (23/23, determinista)"],
        [6, "R(g,L) en Jiang", "`compute_de` + `prep_de_side` del preset, `p_adj < 0.05`, sin el propio gen; comprobado gen a gen contra `venn_ed.json`"],
        [7, "unidad de la Parte 2", "(g, gen, L1, L2), L1 ≠ L2, gen ∈ R(g,L1); y = 1 si gen ∉ R(g,L2)"],
        [8, "predictores (todos del control de L2), orden", "p1 pasa filtro → p2 CPM medio (deciles) → p3 deteccion (deciles) → p4 potencia (log n, log UMI mediana); fijado antes de medir"],
        [9, "medida A", "fraccion de las y = 1 cuyo gen no pasa el filtro de L2 (determinista)"],
        [10, "medida B", "D² de una logistica (devianza explicada), SOLO y ACUMULADO"],
        [11, "nulo", "20 permutaciones uniformes de las 6 lineas: las estadisticas de control de L2 se toman de π(L2); semillas 20281100 + k"],
        [12, "anclaje", "(obs − media nulo) / (1 − media nulo)"],
        [13, "igualacion en potencia (Jiang)", "control submuestreado a n_min y TODAS las celulas adelgazadas binomialmente a la UMI mediana minima; se re-corre la puerta"],
        [14, "VCC 2025", "la poblacion real de `vcc2025_suelo` (semillas 20280201/20280202 reusadas a proposito), verificada por sha256"],
    ])
    W += ["", "### 0.1 Semillas: familia nueva 20281000+, comprobadas EN CODIGO", ""]
    W += tabla(["uso", "semillas"], [[k, ", ".join(map(str, v))] for k, v in sem["semillas_enteras"].items()])
    W += ["", f"Choques internos: {len(sem['choques_internos'])}; con rangos gastados: {len(sem['choques_rangos_gastados'])}; "
          f"con la enumeracion de `noche/comun.py`: {len(sem['choques_enumeracion_comun'])}; con el barrido de "
          f"{sem['n_enteros_2026_2028_barridos']} enteros 2026xxxx–2028xxxx del codigo y las salidas (maximo previo "
          f"de la familia 2028: {sem['max_barrido_familia_2028']}): {len(sem['choques_barrido'])}. "
          f"**Sin choques: {'si' if sem['sin_choques'] else 'NO'}.** Guardadas en `{OUT_JSON}` y en "
          "`modelo/src/seleccion/semillas.py`.", ""]
    W += ["---", ""] + W1 + ["---", ""] + W2
    W += ["---", "", "## Defectos propios", "",
          "1. **El residuo de 2.3 incluye el ruido de muestreo de las DOS llamadas** (R(g,L1) y R(g,L2)), no solo "
          "biologia. En Jiang no hay replica dentro de linea para medir ese suelo, asi que el residuo es una cota "
          "superior de lo que el control no explica por razones biologicas.",
          "2. **El nulo usa permutaciones uniformes**, que dejan puntos fijos (columna `# puntos fijos`): en esos "
          "sorteos parte de las estadisticas son las verdaderas, lo que sube la media del nulo y hace conservador el anclaje.",
          "3. **p4 es de linea**: con 6 lineas, 2 coeficientes sobre 6 valores capturan cualquier diferencia de linea "
          "correlacionada con n y profundidad, no solo potencia.",
          "4. **D² depende de la codificacion** (deciles, logistica con L2 debil C = 1e4). Otra codificacion da otros numeros; "
          "esta se fijo antes de medir.",
          "5. **Las filas no son independientes** (el mismo gen entra en varios pares y perturbaciones): no se dan "
          "intervalos; la unica referencia de dispersion es el rango del nulo.",
          "6. **La igualacion solo iguala n de control y profundidad**. El lado perturbado ya tiene 400 celulas en todas "
          "las lineas; la eficiencia de la guia es del lado perturbado y queda fuera por diseño.",
          "7. **Jiang tiene pools de control de 1,379–3,228 celulas** y un espacio de 34,025 genes; el reto, 18,400 y "
          "18,533. El regimen de potencia es otro: por eso de la Parte 2 no viaja ninguna fraccion.",
          "8. **VCC 2025 es un solo contexto**: no da especificidad de contexto, solo 2.4.",
          "9. **El Fano es de cuentas crudas** y depende de la profundidad.",
          "10. **El control de 1.4 usa mitades** (9,200 celulas): su ruido es ~√2 mayor que el de un contexto entero; "
          "acota por arriba la fragilidad por muestreo, no la mide al tamaño real.",
          "11. **La medida A da a p1 todo el credito determinista** porque p1 va primero en el orden declarado; la "
          "columna SOLO de 2.2 es la que permite ver cuanto toma cada predictor por su cuenta.",
          "12. **Las cifras agregadas de 2.1–2.3 estan dominadas por MCF7 como L1** (Σ|R| = 12,867 de 24,410 en crudo): "
          "cada fila pesa igual, asi que las lineas con R(g,L) grande pesan mas. La tabla por par ordenado de 2.1 es la "
          "que lo deja ver.",
          "13. **El nulo de 2.1/2.2 ya es alto** (p. ej. ~12 % en la medida A) porque la expresion de un gen se parece "
          "entre lineas: las estadisticas de otra linea ya predicen parte de y. El anclaje descuenta eso; la cifra cruda no.", "",
          f"_Generado {time.strftime('%F %T')} UTC. Fuentes: `/datos/modelo/out/seleccion/p1_filtrables.json`, "
          "`p2_datos.json`, `p2_mecanismo.json`; codigo en `modelo/src/seleccion/`._"]
    with open(OUT_MD + ".tmp", "w") as fh:
        fh.write("\n".join(W) + "\n")
    os.replace(OUT_MD + ".tmp", OUT_MD)
    json.dump({"semillas": sem, "parte1_resumen": {L: {k: d1["ctx"][L][k] for k in ("n_pasa", "verif", "L_mediana")}
                                                   for L in d1["ctx"]},
               "parte2_mecanismo": {k: m2[k] for k in ("crudo", "igualado", "universo")}},
              open(OUT_JSON, "w"), default=str)


if __name__ == "__main__":
    main()
