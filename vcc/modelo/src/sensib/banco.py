#!/usr/bin/env python
"""El banco: una CELDA = una estructura inyectada en un contexto, con su propio
bundle, y uno o varios BRAZOS de prediccion puntuados contra ella.

🔴 La rejilla del BRAZO B esta DECLARADA AQUI, antes de correr nada. La del
BRAZO C no se declara: sale MEDIDA del brazo A (transferencia de ESTRUCTURA).

Coste de una celda (medido, no estimado, en `dt_*` de cada registro):
  baseline + prep-real-bundle una vez por celda (el lado real cambia con la
  estructura, asi que el bundle no se puede compartir entre celdas), y
  run + score por brazo.

Cada celda escribe su json EN CUANTO TERMINA y borra sus intermedios grandes.
Reentrante: si el json de la celda existe, se salta.
"""
import itertools
import json
import os
import shutil
import sys
import time

import numpy as np

sys.path.insert(0, "/home/ubuntu/vcc/modelo/src/noche")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import comun as C  # noqa: E402
import sem as S  # noqa: E402
import sup as U  # noqa: E402

CELDAS = os.path.join(U.WORK, "celdas")
TMP = os.path.join(U.WORK, "tmp")

# ───────── REJILLA DEL BRAZO B, DECLARADA ANTES DE CORRER ─────────
AMPLITUDES = [0.25, 0.50, 1.00, 2.00]        # |log2FC| inyectado
FRACCIONES = [0.002, 0.02, 0.20]             # fraccion de genes FILTRABLES movidos
COMUNES = [0.0, 0.5]                         # fraccion del movimiento que es comun
REJILLA_B = [{"k": i, "amplitud": a, "frac_genes": fg, "frac_comun": fc}
             for i, (a, fg, fc) in enumerate(itertools.product(AMPLITUDES, FRACCIONES, COMUNES))]
CENTRO_B = 9                                  # (0.50, 0.02, 0.5): la celda de calibracion
F_GRID = [0.0, 0.25, 0.50, 0.75, 1.0]         # barrido de calidad del brazo C


def cruz():
    """El subconjunto EN CRUZ: el centro y, desde el, un eje a la vez."""
    c = REJILLA_B[CENTRO_B]
    out = [CENTRO_B]
    for i, r in enumerate(REJILLA_B):
        dif = sum(1 for kk in ("amplitud", "frac_genes", "frac_comun") if r[kk] != c[kk])
        if dif == 1:
            out.append(i)
    return sorted(set(out))


def libre():
    return shutil.disk_usage("/datos")[2] / (1024 ** 3)


def brazos_B():
    return [{"id": "suelo", "modo": "suelo", "f": None},
            {"id": "oraculo", "modo": "magnitud", "f": 1.0}]


def brazos_C():
    b = [{"id": "suelo", "modo": "suelo", "f": None}]
    for i, f in enumerate(F_GRID):
        b.append({"id": f"mag_f{f:g}", "modo": "magnitud", "f": f, "slot": i})
        b.append({"id": f"ide_f{f:g}", "modo": "identidad", "f": f, "slot": 5 + i})
    return b


def nombre(ctx, etiqueta):
    return f"{ctx}_{etiqueta}"


def corre_celda(ctx, etiqueta, k_sem, params, brazos, cache=None, extra=None, fam="B"):
    """Una celda entera. `params`: amplitud, frac_genes, frac_comun.
    `k_sem`: indice de semilla de estructura (y de ruido de inyeccion)."""
    out_json = os.path.join(CELDAS, nombre(ctx, etiqueta) + ".json")
    if os.path.exists(out_json):
        C.log(f"{ctx}/{etiqueta}: ya hecho"); return json.load(open(out_json))
    g0 = libre()
    if g0 < C.MIN_FREE_GIB:
        raise SystemExit(f"PARAR: {g0:.1f} GiB libres en /datos")
    t0 = time.time()
    rec = {"ctx": ctx, "etiqueta": etiqueta, "k_semilla": k_sem, "params": params,
           "gib_libres_al_empezar": round(g0, 2), "brazos": {}}

    # ── superficie (cacheada en el proceso: leerla cuesta ~40 s y 1.3 GiB)
    if cache is not None and cache.get("ctx") == ctx:
        X, ntc, var_names, inv, par, filtr = (cache["X"], cache["ntc"], cache["var"],
                                              cache["inv"], cache["par"], cache["filtr"])
    else:
        t = time.time()
        X, ntc, var_names, inv = U.lee_contexto(ctx)
        par = U.particion(ctx, ntc, var_names)
        par["var_names"] = var_names
        filtr, cpm = U.filtrables(X[par["idx_control"]], var_names,
                                  fuera=par["etiquetas"])
        inv["dt_lectura_s"] = time.time() - t
        inv["n_filtrables"] = int(filtr.size)
        inv["cpm_mediana_filtrables"] = float(np.median(cpm[filtr]))
        if cache is not None:
            cache.clear()
            cache.update({"ctx": ctx, "X": X, "ntc": ntc, "var": var_names,
                          "inv": inv, "par": par, "filtr": filtr})
    rec["inventario"] = {k: v for k, v in inv.items() if k != "cpm"}
    rec["particion"] = {k: v for k, v in par.items()
                        if k in ("guias_control", "guias_pert", "semilla", "n_control",
                                 "n_pp", "guia_partida", "etiquetas",
                                 "guia_de_etiqueta", "diag_etiquetas")}
    Xc = X[par["idx_control"]]
    etiquetas = sorted(par["pp"])

    # ── estructura inyectada
    seeds = S.s_estructura(k_sem) if fam == "B" else S.s_estructura_C(k_sem)
    s_ruido = S.s_ruido(ctx, k_sem) if fam == "B" else S.s_ruido_C(ctx, k_sem)
    L, esp_idx, com, s_com, diag = U.estructura(
        k_sem, etiquetas, filtr, len(var_names),
        params["amplitud"], params["frac_genes"], params["frac_comun"], seeds=seeds)
    diag["familia_semillas"] = fam
    diag["semilla_ruido_inyeccion"] = int(s_ruido)
    rec["estructura"] = diag
    if extra:
        rec["estructura"].update(extra)

    d = os.path.join(TMP, nombre(ctx, etiqueta))
    os.makedirs(d, exist_ok=True)
    real = os.path.join(d, "real.h5ad")
    try:
        t = time.time()
        rec["real"] = U.escribe_real(X, par, L, s_ruido, real)
        rec["dt_real_s"] = time.time() - t

        # ── baseline + bundle (una vez por celda)
        bp = os.path.join(d, "base_pred.h5ad")
        bd = os.path.join(d, "BASE")
        p, dt = C.corre(U.cmd_baseline(real, bd, bp), os.path.join(d, "baseline.log"))
        rec["dt_baseline_s"] = dt
        rec["baseline_ok"] = p.returncode == 0
        if p.returncode != 0:
            raise RuntimeError(f"baseline fallo: {p.stderr[-1200:]}")
        bu = os.path.join(d, "BUNDLE")
        p, dt = C.corre(U.cmd_bundle(real, bp, bu), os.path.join(d, "bundle.log"))
        rec["dt_bundle_s"] = dt
        rec["bundle_ok"] = p.returncode == 0
        if p.returncode == 0:
            rec["hash_bundle"] = U.hash_ok(os.path.join(bu, "config.yaml")) \
                if os.path.exists(os.path.join(bu, "config.yaml")) else None
        else:
            C.log(f"{ctx}/{etiqueta}: BUNDLE FALLO, se puntua solo `u`")
        os.remove(bp)
        C.escribe_json(out_json, rec)          # salida incremental: ya hay bundle+tiempos

        # ── brazos
        for b in brazos:
            tb = time.time()
            rng = np.random.default_rng(S.s_parcial(ctx, b.get("slot", 0))
                                        if b["modo"] != "suelo" else S.s_suelo(ctx))
            raz = U.razones_pred(L, esp_idx, com, s_com, etiquetas, len(var_names),
                                 filtr, b["modo"], b.get("f") or 0.0,
                                 params["amplitud"], rng)
            sem_pred = (S.s_suelo(ctx) if b["modo"] == "suelo"
                        else S.s_parcial(ctx, b.get("slot", 0)))
            pp = os.path.join(d, f"pred_{b['id']}.h5ad")
            meta_p = U.escribe_pred(Xc, par, raz, ctx, sem_pred, pp)
            rd = os.path.join(d, f"run_{b['id']}")
            p, dt_run = C.corre(U.cmd_run(pp, real, rd), os.path.join(d, f"run_{b['id']}.log"))
            br = {"modo": b["modo"], "f": b.get("f"), "pred": meta_p,
                  "dt_run_s": dt_run, "run_ok": p.returncode == 0,
                  "razon_media_abs_log2": float(np.mean([np.abs(v).sum() for v in raz.values()])),
                  "n_genes_movidos_pred": float(np.mean([int((v != 0).sum()) for v in raz.values()]))}
            if p.returncode == 0:
                br["u"] = U.lee_u(rd)
                br["hash_run"] = U.hash_ok(os.path.join(rd, "run_params.yaml"))
                if rec.get("bundle_ok"):
                    sc = os.path.join(d, f"scored_{b['id']}.csv")
                    p2, dt_sc = C.corre(U.cmd_score(rd, bu, sc),
                                        os.path.join(d, f"score_{b['id']}.log"))
                    br["dt_score_s"] = dt_sc
                    br["score_ok"] = p2.returncode == 0
                    br["s"] = U.lee_s(sc) if p2.returncode == 0 else None
                else:
                    br["s"] = None
                    br["score_ok"] = False
            else:
                br["mensaje"] = p.stderr[-1200:]
            br["dt_brazo_s"] = time.time() - tb
            rec["brazos"][b["id"]] = br
            C.borra(pp, rd)
            C.escribe_json(out_json, rec)      # incremental: brazo a brazo
            u = br.get("u") or {}
            C.log(f"{ctx}/{etiqueta}/{b['id']}: run {dt_run:.0f}s "
                  f"pds={u.get('pds_cosine')} jac={u.get('de_wilcoxon_sig_jaccard')} "
                  f"s_avg={(br.get('s') or {}).get('avg_score', {}).get('from_replicate')}")
    finally:
        C.borra(d)
    rec["dt_celda_s"] = time.time() - t0
    rec["gib_libres_al_acabar"] = round(libre(), 2)
    C.escribe_json(out_json, rec)
    C.log(f"{ctx}/{etiqueta}: celda completa en {rec['dt_celda_s'] / 60:.1f} min")
    return rec


def main():
    os.makedirs(CELDAS, exist_ok=True)
    os.makedirs(TMP, exist_ok=True)
    sem = S.comprueba()
    assert sem["sin_choques"], sem
    C.escribe_json(os.path.join(U.WORK, "semillas.json"), sem)
    modo = sys.argv[1]
    if modo == "calibra":
        r = REJILLA_B[CENTRO_B]
        corre_celda("D", f"B{CENTRO_B:02d}", CENTRO_B,
                    {k: r[k] for k in ("amplitud", "frac_genes", "frac_comun")},
                    brazos_B(), cache={})
    elif modo == "cruz":
        print(json.dumps([REJILLA_B[i] for i in cruz()], indent=1))
    elif modo == "rejilla":
        print(json.dumps(REJILLA_B, indent=1))


if __name__ == "__main__":
    main()
