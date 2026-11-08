#!/usr/bin/env python
"""B4b - PASO 4: los anclajes NATIVOS de la superficie VCC 2025, con la MISMA
maquinaria que en Jiang.

  * `b` de `cell-eval2 baseline` (baseline.py)   -- via `anclas_nativas.cmd_baseline`
  * `r` de `cell-eval2 prep-real-bundle`         -- via `anclas_nativas.cmd_bundle`
  * bundle de competencia, y `config_hash` VERIFICADO contra
    `competition.competition_payload()`

Ningun comando se reescribe: se llaman los armadores de `anclas_nativas.py` sin
modificarlos, con `--control` omitido para que el valor quede el del preset
(`non-targeting`) y el `config_hash` no se desvie -- el defecto que
`anclas_nativas.py` documenta.

🔴 ADVERTENCIA TRANSDUCTIVA, la misma de Jiang (baseline.py:6-10): el perfil del
baseline se promedia de las perturbaciones reales evaluadas, asi que `b` y `r`
son informacion de ORACULO. Son artefactos de PUNTUACION. Ninguna de sus cifras
puede entrar en un modelo.
"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import comun as C
import v25

sys.path.insert(0, C.SRC)

OUT_JSON = os.path.join(C.OUT, "vcc2025_anclajes.json")
OUT_MD = os.path.join(C.OUT, "vcc2025_anclajes.md")
BUNDLE = os.path.join(v25.WORK, "BUNDLE_2025_test")
BASE_DIR = os.path.join(v25.WORK, "BASELINE_2025_test")
BASE_PRED = os.path.join(v25.WORK, "baseline_pred_2025_test.h5ad")


def jiang_ref():
    """Los `b` y `r` de Jiang ya medidos, para la tabla de control."""
    p = os.path.join(C.OUT, "anclas_nativas_comp.json")
    if not os.path.exists(p):
        return []
    return [json.loads(l) for l in open(p) if l.strip()]


def md(res):
    o = []; A = o.append
    A("# VCC 2025 como superficie: los anclajes nativos (`b`, `r`, `config_hash`)")
    A("")
    A("## 🔴 0. Que es esto y que no")
    A("")
    A("Los dos extremos de la escala `vcc2026` construidos por el camino "
      "**NATIVO del marcador** sobre la superficie 2025, con los **mismos "
      "armadores de comando** que se usaron en Jiang "
      "(`anclas_nativas.cmd_baseline` y `anclas_nativas.cmd_bundle`, sin "
      "modificar).")
    A("")
    A("🔴 **`b` y `r` son informacion de ORACULO** (`baseline.py:6-10`): el "
      "perfil del baseline se promedia de las perturbaciones reales evaluadas. "
      "Son artefactos de PUNTUACION. **Ninguna de sus cifras puede entrar en un "
      "modelo.**")
    A("")
    A("## 1. La superficie, declarada antes de medir")
    A("")
    A("\n".join(v25.tabla_construccion(res["real_meta"])))
    A("")
    A("## 2. Las dos invocaciones y sus guardas")
    A("")
    fil = []
    for k in ("baseline", "bundle"):
        d = res.get(k, {})
        fil.append([f"`cell-eval2 {'baseline' if k=='baseline' else 'prep-real-bundle'}`",
                    "si" if d.get("ok") else "🔴 NO",
                    d.get("exitcode"), C.dur(d.get("duracion_s") or 0),
                    (d.get("error") or "")[:200]])
    A("\n".join(C.tabla(["comando", "ok", "exitcode", "duracion", "error (200 car.)"], fil)))
    A("")
    A("**Comandos literales:**")
    A("")
    for k in ("baseline", "bundle"):
        if res.get(k, {}).get("cmd"):
            A(f"```\n{' '.join(res[k]['cmd'])}\n```")
            A("")
    if not res.get("bundle", {}).get("ok"):
        A("🔴 **El bundle no se completo.** Lo que sigue esta incompleto por esa "
          "razon, no por eleccion. `vcc2025_inventario.md` queda en pie.")
        A("")
    A("## 3. 🔴 El `config_hash`, verificado contra `competition_payload()`")
    A("")
    h = res.get("hash_regla") or {}
    A("\n".join(C.tabla(["# ", "cantidad", "valor"], [
        [1, "`rule_config_hash` obtenido del `config.yaml` del bundle", f"`{h.get('obtenido')}`"],
        [2, "`config_hash` de referencia (`competition_payload()`)", f"`{h.get('referencia')}`"],
        [3, "**coincide**", "**si**" if h.get("coincide") else "🔴 **NO**"],
        [4, "`control` en el config", f"`{h.get('control_en_config')}`"],
        [5, "`pert_col` en el config", f"`{h.get('pert_col_en_config')}`"],
        [6, "bundle diagnostico", res.get("bundle_diagnostico")],
        [7, "`anchor_base_seed` / `anchor_splits`", f"{res.get('anchor_base_seed')} / {res.get('anchor_splits')} (convencion del marcador, `competition.py:38,41`)"],
    ])))
    A("")
    A("## 4. El `b` por metrica: los dos lados")
    A("")
    jr = res.get("jiang", [])
    fil = []
    for m in C.VCC2026_METRICS:
        v = (res.get("b") or {}).get(m)
        jv = [r["b"].get(m) for r in jr if r.get("b") and r["b"].get(m) is not None]
        fil.append([f"`{C.METRICA_CORTA[m]}`", m, C.f4(v) if v is not None else "n/d",
                    C.rango(jv), C.f4(sum(jv) / len(jv)) if jv else "n/d"])
    A("\n".join(C.tabla(
        ["corta", "metrica", "**`b` en VCC 2025 (test)**",
         "rango del `b` en las 6 lineas de Jiang", "media en Jiang"], fil)))
    A("")
    A("🔴 Los dos lados **no son comparables como numeros**: distinta superficie, "
      "distinto panel (18,080 contra 34,025 genes), distinta profundidad "
      "(~54,000 contra 3,319-6,203 UMI/celula). La columna de Jiang esta como "
      "**control de que la maquinaria corrio igual**, no como comparacion.")
    A("")
    A("## 5. El `r` por metrica (`replicate`, el anclaje de arriba)")
    A("")
    fil = []
    for m in C.VCC2026_METRICS:
        d = (res.get("r") or {}).get(m) or {}
        fil.append([f"`{C.METRICA_CORTA[m]}`",
                    C.f4(d.get("replicate")), C.f4(d.get("replicate_sd")),
                    C.f4(d.get("replicate_min")), C.f4(d.get("replicate_max")),
                    d.get("n_perturbations_min"), d.get("n_perturbations_max"),
                    d.get("estimator")])
    A("\n".join(C.tabla(
        ["metrica", "**`replicate`**", "sd", "min", "max",
         "n perts min", "n perts max", "estimador"], fil)))
    A("")
    A("### 5.1 🔴 Los intermedios crudos: los 5 splits de anclaje, uno por uno")
    A("")
    sp = res.get("splits") or {}
    fil = []
    n = 0
    for m in C.VCC2026_METRICS:
        for d in sp.get(m, []):
            n += 1
            fil.append([n, f"`{C.METRICA_CORTA[m]}`", d.get("split_index"),
                        d.get("seed"), C.f4(d.get("value")),
                        d.get("n_perturbations")])
    A("\n".join(C.tabla(
        ["# ", "metrica", "split", "semilla", "valor", "n perturbaciones"], fil))
      if fil else "(no se escribieron `anchor_splits`)")
    A("")
    A("## 6. 🔴 Lo que esta superficie NO puede dar")
    A("")
    A("\n".join(C.tabla(["# ", "lo que no da", "por que"], [
        [1, "**Dominancia entre contextos**", "Es UN solo contexto. `b` y `r` son de un punto, no de un rango."],
        [2, "**Rango entre contextos**", "En Jiang hay 6 `b` y 6 `r` y por eso hay rango. Aqui hay uno."],
        [3, "Un `b` comparable con el de Jiang", "Distinto panel y distinta profundidad; la comparacion es de maquinaria, no de valor."],
    ])))
    A("")
    A("## 7. Defectos de la construccion")
    A("")
    A("\n".join(C.tabla(["# ", "defecto", "que implica"], [
        [1, "`b` y `r` son transductivos.", "Artefactos de puntuacion. No entran en ningun modelo."],
        [2, "15 de 100 perturbaciones llevan < 400 celulas.", "La potencia por perturbacion no es uniforme; el `n_perturbations` de cada split de anclaje lo refleja."],
        [3, "El control se submuestreo a 18,400 de 38,176.", "Semilla declarada. El `b` depende del pool de control, asi que depende de ese submuestreo."],
        [4, "Una sola realizacion de la poblacion real.", "No hay barra de error sobre la eleccion de las celulas por perturbacion."],
        [5, "`--control` se omite a proposito.", "Para que el valor quede el del preset y el `config_hash` no se desvie -- el defecto que `anclas_nativas.py` documenta de Jiang."],
    ])))
    A("")
    A(f"*Corrida: {res['t_inicio']} -> {res['t_fin']} UTC, duracion "
      f"{C.dur(res['dt_total_s'])}.*")
    return o


def main():
    t0 = time.time()
    os.makedirs(v25.WORK, exist_ok=True)
    C.guarda_disco()
    import anclas_nativas as AN
    real_path, real_meta = v25.arma_real()
    C.guarda_disco()
    res = {"bloque": "B4b", "semillas": C.comprueba_semillas(),
           "real_meta": real_meta, "jiang": jiang_ref(),
           "anchor_base_seed": AN.ANCHOR_BASE_SEED,
           "anchor_splits": AN.ANCHOR_SPLITS,
           "t_inicio": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(t0))}

    os.makedirs(BASE_DIR, exist_ok=True)
    cmd = AN.cmd_baseline(real_path, BASE_DIR, BASE_PRED, control=None)
    C.log("B4b cell-eval2 baseline...")
    p, dt = C.corre(cmd, os.path.join(v25.WORK, "b4b_baseline.log"))
    res["baseline"] = {"cmd": cmd, "ok": p.returncode == 0,
                       "exitcode": p.returncode, "duracion_s": dt,
                       "error": (p.stderr[-3000:] if p.returncode else None)}
    C.log(f"B4b baseline exit={p.returncode} {dt:.0f}s")
    if p.returncode == 0:
        try:
            res["b"] = AN.leer_b(BASE_DIR)
        except Exception as e:
            res["b_error"] = repr(e)
        try:
            res["baseline_meta"] = AN.leer_meta(
                os.path.join(BASE_DIR, "run_meta.json"),
                ["config_hash", "n_perturbations", "control", "pert_col"])
        except Exception as e:
            res["baseline_meta_error"] = repr(e)
    res["t_fin"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
    res["dt_total_s"] = time.time() - t0
    C.escribe_json(OUT_JSON, res)

    if p.returncode == 0:
        C.guarda_disco()
        cmd2 = AN.cmd_bundle(real_path, BASE_PRED, BUNDLE, control=None)
        C.log("B4b cell-eval2 prep-real-bundle (5 splits de anclaje)...")
        p2, dt2 = C.corre(cmd2, os.path.join(v25.WORK, "b4b_bundle.log"))
        res["bundle"] = {"cmd": cmd2, "ok": p2.returncode == 0,
                         "exitcode": p2.returncode, "duracion_s": dt2,
                         "error": (p2.stderr[-3000:] if p2.returncode else None)}
        C.log(f"B4b bundle exit={p2.returncode} {dt2:.0f}s")
        if p2.returncode == 0:
            for nombre, f in (("r", lambda: AN.leer_r(BUNDLE)),
                              ("splits", lambda: AN.leer_splits(BUNDLE)),
                              ("hash_regla", lambda: AN.hash_de_la_regla(BUNDLE))):
                try:
                    res[nombre] = f()
                except Exception as e:
                    res[nombre + "_error"] = repr(e)
            md_path = os.path.join(BUNDLE, "metadata.json")
            if os.path.exists(md_path):
                try:
                    res["bundle_diagnostico"] = json.load(open(md_path)).get("diagnostic")
                except Exception:
                    pass
    else:
        res["bundle"] = {"ok": False, "exitcode": None, "duracion_s": 0,
                         "error": "no se intento: el baseline fallo"}

    # el pred del baseline es grande y ya no se necesita
    C.borra(BASE_PRED)
    res["t_fin"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
    res["dt_total_s"] = time.time() - t0
    C.escribe_json(OUT_JSON, res)
    C.escribe_md(OUT_MD, md(res))
    C.log(f"B4b listo: {OUT_MD}")


if __name__ == "__main__":
    main()
