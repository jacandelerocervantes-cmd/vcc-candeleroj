#!/usr/bin/env python
"""BARRIDO DE k SOBRE VCC 2025, BLANCOS RETENIDOS.

Donde: VCC 2025, pliegue 0 de LA MISMA particion por blanco de la CV
(semilla 20284950, leida de entrega_ajuste.json, no reconstruida a ojo).
El modelo se REAJUSTA sobre los otros 4 pliegues: ningun blanco puntuado
estuvo en su ajuste.  NO se toca D/E/F: ahi no hay verdad.
Lo unico que viaja despues es k, un numero.

Brazos: k en {1.0, 3.4, 2.0, 1.5, 2.5} con los 10 rasgos, mas el modelo de
8 rasgos a k=1.0.  Orden declarado: extremos y ablacion primero, para que una
tabla parcial abarque el rango.
"""
import os, sys, json, time, glob, subprocess
import numpy as np, scipy.sparse as sp, h5py

for p in ("/home/ubuntu/vcc/modelo/src", "/home/ubuntu/vcc/modelo/src/entrega",
          "/home/ubuntu/vcc/modelo/src/sensib", "/home/ubuntu/vcc/modelo/src/ntc_nulo"):
    sys.path.insert(0, p)
import datos as D
import ajuste as A
import sup as S
import suelo_arnes as H

WORK = "/datos/modelo/out/barrido"
OUT = "/home/ubuntu/vcc/modelo/out/barrido_k.json"
AJ = "/home/ubuntu/vcc/modelo/out/entrega_ajuste.json"
N_POR = 400
N_CTRL = 9200          # regimen 400 vs 9,200, el que ntc_nulo declaro como el de entrega
RECORTE = 2.0
CTRL = "non-targeting"
PLIEGUE = 0
MIN_FREE_GIB = 10.0

# rasgos del BLANCO g, los que cae la ablacion de 8 (criterio declarado ANTES de ajustar)
FUERA_8 = ["g_es_TF", "log10_cpm_g"]

# ── semillas, familia 20284000+, comprobadas en codigo y guardadas enteras
SEM = {"muestra_real": 20284600, "control": 20284601}
for i, nom in enumerate(["k1.0", "k3.4", "r8_k1.0", "k2.0", "k1.5", "k2.5"]):
    SEM["pred_" + nom] = 20284700 + i

def comprueba_semillas():
    gastadas = set()
    for a, b in ((20260000, 20266010), (20270000, 20283000)):
        gastadas |= set(range(a, b))
    gastadas |= {20284100, 20284101, 20284200, 20284201, 20284300, 20284301, 20284950}
    gastadas |= set(range(20284900, 20284910))
    aj = json.load(open(AJ))
    def hoja(o):
        if isinstance(o, int) and 20000000 < o < 21000000: gastadas.add(o)
        elif isinstance(o, dict):
            for v in o.values(): hoja(v)
        elif isinstance(o, list):
            for v in o: hoja(v)
    hoja(aj.get("semillas", {}))
    nuevas = sorted(SEM.values())
    choques = sorted(set(nuevas) & gastadas)
    return {"semillas": dict(SEM), "n_nuevas": len(nuevas),
            "n_gastadas_enumeradas": len(gastadas),
            "choques": choques, "sin_choques": not choques}

def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

def libre():
    s = os.statvfs("/datos"); return s.f_bavail * s.f_frsize / 2**30

def guarda(res):
    tmp = OUT + ".tmp"
    json.dump(res, open(tmp, "w"), indent=1, default=str); os.replace(tmp, OUT)

def corre(cmd, logp):
    with open(logp, "w") as f:
        p = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    return p.returncode

# ───────────────────────────── el lado real de 2025 ─────────────────────────────
def lado_real(targets):
    """Celulas REALES de 2025: hasta 400 por blanco retenido + N_CTRL de control."""
    genes = D.genes_2025(); G = len(genes)
    rng = np.random.default_rng(SEM["muestra_real"])
    tset = set(targets)
    piezas, nreal, ctrl_por_split = {}, {}, {}
    for split, rel in D.SPLITS.items():
        path = os.path.join(D.V25, rel)
        with h5py.File(path, "r") as f:
            tg = np.asarray(D.decode_categorical(f["obs"]["target_gene"]),
                            dtype=object).astype(str)
        aqui = sorted(tset & set(np.unique(tg).tolist()))
        for t in aqui:
            idx = np.where(tg == t)[0]
            nreal[t] = int(len(idx))
            if len(idx) > N_POR:
                idx = np.sort(rng.choice(idx, N_POR, replace=False))
            piezas[t] = (path, idx)
        ctrl_por_split[split] = (path, np.where(tg == CTRL)[0])
        log(f"  real {split}: {len(aqui)} blancos del pliegue, {len(ctrl_por_split[split][1]):,} control")
    # control: N_CTRL repartido en proporcion al tamano de cada split
    tot = sum(len(v[1]) for v in ctrl_por_split.values())
    rc = np.random.default_rng(SEM["control"])
    bloques_ctrl, n_ctrl_real = [], 0
    for split, (path, ci) in ctrl_por_split.items():
        n = int(round(N_CTRL * len(ci) / tot))
        sel = np.sort(rc.choice(ci, min(n, len(ci)), replace=False))
        bloques_ctrl.append(H.read_rows_csr(path, sel, G)); n_ctrl_real += len(sel)
    Xc = sp.vstack(bloques_ctrl, format="csr")
    bl, lab = [], []
    for t in sorted(piezas):
        path, idx = piezas[t]
        bl.append(H.read_rows_csr(path, idx, G)); lab += [t] * len(idx)
    bl.append(Xc); lab += [CTRL] * Xc.shape[0]
    M = sp.vstack(bl, format="csr")
    return M, lab, Xc, genes, nreal, {"n_control": int(Xc.shape[0]),
                                      "n_control_disponible": int(tot)}

def escribe(M, lab, genes, path):
    import anndata as ad, pandas as pd
    a = ad.AnnData(X=M.astype(np.int64), obs=pd.DataFrame({"target_gene": lab}),
                   var=pd.DataFrame(index=list(genes)))
    a.obs_names_make_unique(); a.write_h5ad(path)
    return {"shape": list(a.shape), "nnz": int(M.nnz)}

def lado_pred(Xc, targets, Y, genes, semilla, path):
    """400 por blanco remuestreadas CON REEMPLAZO del pool de control + 400 de
    control.  _reescala_entero via sup.aplica, SIN MODIFICAR."""
    rng = np.random.default_rng(semilla)
    n = Xc.shape[0]; bl, lab = [], []
    for i, t in enumerate(targets):
        ch = rng.integers(0, n, size=N_POR)
        bl.append(S.aplica(Xc[ch], Y[i], rng)); lab += [t] * N_POR
    ch = rng.integers(0, n, size=N_POR)
    bl.append(Xc[ch].tocsr()); lab += [CTRL] * N_POR
    M = sp.vstack(bl, format="csr")
    return escribe(M, lab, genes, path)

# ─────────────────────────── la tabla de cuatro celdas ───────────────────────────
def cuatro_celdas(run_dir, targets):
    """De los de_*.parquet de `run --write-degenes`: real∩pred, solo real,
    solo pred (falsos positivos), ambos vacios.  La puerta es la del marcador."""
    import pandas as pd
    info = {"archivos": sorted(os.path.basename(p) for p in glob.glob(os.path.join(run_dir, "*.parquet")))}
    def carga(patron):
        g = [p for p in glob.glob(os.path.join(run_dir, "*.parquet")) if patron in os.path.basename(p)]
        if not g: return None, None
        df = pd.read_parquet(g[0])
        info[patron + "_columnas"] = list(df.columns)
        info[patron + "_filas"] = int(len(df))
        return df, os.path.basename(g[0])
    dr, nr = carga("de_real"); dp, npd = carga("de_pred")
    if dr is None or dp is None:
        info["error"] = "no aparecieron de_real/de_pred.parquet"; return info
    def conjuntos(df):
        cp = next((c for c in ("target", "perturbation", "target_gene", "pert") if c in df.columns), df.columns[0])
        cg = next((c for c in ("feature", "gene", "gene_name", "feature_name") if c in df.columns), df.columns[1])
        q = next((c for c in df.columns if "p_adj" in c or "padj" in c or "qval" in c), None)
        d = df if q is None else df.loc[df[q] < 0.05]
        info.setdefault("columnas_usadas", {})[cp] = cg
        if q: info.setdefault("filtro_q", q)
        return {str(k): set(map(str, v)) for k, v in d.groupby(cp)[cg].apply(set).items()}
    R, P = conjuntos(dr), conjuntos(dp)
    filas, inter, sr, sp_, vac = [], [], [], [], 0
    for t in targets:
        a = R.get(t, set()) - {t}; b = P.get(t, set()) - {t}
        if not a and not b: vac += 1
        filas.append({"blanco": t, "n_real": len(a), "n_pred": len(b),
                      "interseccion": len(a & b), "solo_real": len(a - b),
                      "solo_pred": len(b - a)})
        inter.append(len(a & b)); sr.append(len(a - b)); sp_.append(len(b - a))
    info.update({"n_blancos": len(targets),
                 "media_interseccion": float(np.mean(inter)),
                 "media_solo_real": float(np.mean(sr)),
                 "media_solo_pred_FALSOS_POSITIVOS": float(np.mean(sp_)),
                 "ambos_vacios": int(vac),
                 "media_genes_llamados_real": float(np.mean([f["n_real"] for f in filas])),
                 "media_genes_llamados_pred": float(np.mean([f["n_pred"] for f in filas])),
                 "por_blanco": filas})
    return info

# ────────────────────────────────── principal ──────────────────────────────────
def main():
    os.makedirs(WORK, exist_ok=True)
    t0 = time.time()
    res = {"que_es": "barrido de k sobre VCC 2025, blancos retenidos del pliegue "
           f"{PLIEGUE} de la MISMA particion por blanco de la CV; el modelo se "
           "reajusta sin ellos. No se toca D/E/F. Lo unico que viaja es k.",
           "semillas": comprueba_semillas(), "n_por_etiqueta": N_POR,
           "n_control_objetivo": N_CTRL, "recorte_log2": RECORTE,
           "orden_declarado": ["k1.0", "k3.4", "r8_k1.0", "k2.0", "k1.5", "k2.5"],
           "brazos": {}}
    assert res["semillas"]["sin_choques"], res["semillas"]
    guarda(res)

    aj = json.load(open(AJ))
    alpha = float(aj["alpha_elegida"])
    pliegues = aj["cv"]["blancos_de_cada_pliegue"]
    retenidos = sorted(pliegues[PLIEGUE])
    res["alpha"] = alpha; res["n_retenidos_pliegue"] = len(retenidos)
    log(f"pliegue {PLIEGUE}: {len(retenidos)} blancos retenidos, alpha={alpha}")

    perts, genes25, cpm_pert, ncel, cpm_ctrl, det, fano, meta = D.pseudobulk_2025()
    y_full = np.log2((cpm_pert + A.EPS_CPM) / (cpm_ctrl[None, :] + A.EPS_CPM)).astype(np.float32)
    st25 = {"cpm": cpm_ctrl, "det": det, "fano": np.where(fano < 0, 0.0, fano),
            "n": meta["n_control_2025"]}
    filtr = np.where(cpm_ctrl > D.UMBRAL_CPM)[0]
    tf = D.tf_set()
    gw_vn, gw_perf, _ = D.gwps(sorted(set(perts.tolist())))
    g_fit = genes25[filtr]
    st_fit = {k: (v[filtr] if isinstance(v, np.ndarray) else v) for k, v in st25.items()}
    X = D.matriz_rasgos(perts, g_fit, gw_vn, gw_perf, tf, st_fit)
    Y = y_full[:, filtr].reshape(-1)
    G = len(g_fit)
    log(f"diseno {X.shape[0]:,} x {X.shape[1]}")

    # ── 1-bis: los rasgos de CONTEXTO no estan en los 10. Se mide su varianza.
    prof = np.median(np.asarray([0.0]))  # marcador, se rellena abajo
    res["rasgos_de_contexto"] = {
        "en_el_modelo": False,
        "rasgos_del_modelo": list(D.RASGOS),
        "nota": "profundidad_mediana y n_control NO entran en matriz_rasgos; el "
                "modelo empaquetado ya es el que no los usa. Dentro de UNA linea "
                "ambos son escalares unicos: su sd sobre las filas del ajuste es 0 "
                "exactamente, asi que no se pueden estandarizar ni identificar.",
        "sd_profundidad_sobre_filas_del_ajuste": 0.0,
        "sd_n_control_sobre_filas_del_ajuste": 0.0,
        "n_control_2025": int(meta["n_control_2025"])}
    i_fuera = [D.RASGOS.index(r) for r in FUERA_8]
    i_8 = [i for i in range(len(D.RASGOS)) if i not in i_fuera]
    res["ablacion_8"] = {"criterio_declarado_antes_de_ajustar":
                         "los 10 menos los dos rasgos del BLANCO g",
                         "fuera": FUERA_8,
                         "dentro": [D.RASGOS[i] for i in i_8]}
    guarda(res)

    # ── CV completa de los dos modelos, por blanco: Pearson y R2
    def cv_modelo(cols, nombre):
        porb_p, porb_r = [], []
        for k, pr in enumerate(pliegues):
            ipr = [int(np.where(perts == t)[0][0]) for t in sorted(pr)]
            itr = [i for i in range(len(perts)) if i not in set(ipr)]
            fil = lambda ts: np.concatenate([np.arange(i * G, (i + 1) * G) for i in ts])
            m = A.ridge(X[fil(itr)][:, cols].astype(np.float64), Y[fil(itr)].astype(np.float64), alpha)
            ite = fil(ipr)
            yp = A.predice(m, X[ite][:, cols].astype(np.float64)); yo = Y[ite].astype(np.float64)
            for n_ in range(len(ipr)):
                s = slice(n_ * G, (n_ + 1) * G)
                porb_p.append(A.pearson(yp[s], yo[s]))
                d = yo[s] - yp[s]
                porb_r.append(float(1 - (d ** 2).sum() / ((yo[s] - yo[s].mean()) ** 2).sum()))
            log(f"  CV {nombre} pliegue {k} hecho")
        def dist(v):
            v = np.asarray([x for x in v if x is not None], dtype=float)
            return {"mediana": float(np.median(v)), "media": float(v.mean()),
                    "p10_p90": [float(np.percentile(v, 10)), float(np.percentile(v, 90))],
                    "min_max": [float(v.min()), float(v.max())],
                    "frac_mayor_0": float((v > 0).mean()), "frac_mayor_0.1": float((v > 0.1).mean()),
                    "n": int(len(v))}
        return {"pearson_por_blanco": dist(porb_p), "R2_por_blanco": dist(porb_r),
                "pearson_medio_por_blanco": float(np.mean([x for x in porb_p if x is not None]))}

    for nombre, cols in (("10_rasgos", list(range(len(D.RASGOS)))), ("8_rasgos", i_8)):
        res.setdefault("cv_por_modelo", {})[nombre] = cv_modelo(cols, nombre)
        guarda(res)
        log(f"CV {nombre}: Pearson por blanco mediana "
            f"{res['cv_por_modelo'][nombre]['pearson_por_blanco']['mediana']:+.4f}, "
            f"R2 mediana {res['cv_por_modelo'][nombre]['R2_por_blanco']['mediana']:+.4f}")

    # ── modelos reajustados SIN el pliegue 0, y la y predicha sobre TODOS los genes
    ipr = [int(np.where(perts == t)[0][0]) for t in retenidos]
    itr = [i for i in range(len(perts)) if i not in set(ipr)]
    fil = lambda ts: np.concatenate([np.arange(i * G, (i + 1) * G) for i in ts])
    Xtr, Ytr = X[fil(itr)].astype(np.float64), Y[fil(itr)].astype(np.float64)
    modelos = {"10_rasgos": (A.ridge(Xtr, Ytr, alpha), list(range(len(D.RASGOS)))),
               "8_rasgos": (A.ridge(Xtr[:, i_8], Ytr, alpha), i_8)}
    Xall = D.matriz_rasgos(np.asarray(retenidos), genes25, gw_vn, gw_perf, tf, st25)
    Ypred = {}
    for nom, (m, cols) in modelos.items():
        yp = A.predice(m, Xall[:, cols].astype(np.float64)).reshape(len(retenidos), len(genes25))
        Ypred[nom] = np.clip(yp, -RECORTE, RECORTE).astype(np.float32)
        res.setdefault("modelos", {})[nom] = {
            "alpha": alpha, "n_rasgos": len(cols), "b0": float(m["b0"]) if isinstance(m, dict) and "b0" in m else None,
            "sd_y_predicha_retenidos": float(Ypred[nom].std()),
            "sd_y_observada_retenidos": float(y_full[ipr, :].std())}
    guarda(res)
    log("modelos reajustados sin el pliegue 0")

    # ── el lado real, UNA vez
    real = os.path.join(WORK, "real.h5ad")
    if not os.path.exists(real):
        M, lab, Xc, genes, nreal, cmeta = lado_real(retenidos)
        np.save(os.path.join(WORK, "pool_ctrl_shape.npy"), np.asarray(Xc.shape))
        sp.save_npz(os.path.join(WORK, "pool_ctrl.npz"), Xc)
        res["lado_real"] = {**escribe(M, lab, genes, real), **cmeta,
                            "n_celulas_reales_por_blanco_disponibles":
                            {k: v for k, v in sorted(nreal.items())}}
        del M
    else:
        Xc = sp.load_npz(os.path.join(WORK, "pool_ctrl.npz")).tocsr()
        genes = D.genes_2025()
    guarda(res); log(f"lado real listo ({libre():.0f} GiB libres)")

    # ── baseline y bundle, UNA vez
    base = os.path.join(WORK, "base_pred.h5ad"); bdir = os.path.join(WORK, "base_out")
    bundle = os.path.join(WORK, "bundle")
    if not os.path.exists(base):
        ec = corre(S.cmd_baseline(real, bdir, base), os.path.join(WORK, "baseline.log"))
        res["baseline_ec"] = ec; guarda(res)
        assert ec == 0, open(os.path.join(WORK, "baseline.log")).read()[-3000:]
    log(f"baseline listo ({time.time()-t0:.0f}s)")
    if not os.path.exists(os.path.join(bundle, "config.yaml")) and not os.path.isdir(bundle):
        ec = corre(S.cmd_bundle(real, base, bundle), os.path.join(WORK, "bundle.log"))
        res["bundle_ec"] = ec; guarda(res)
        assert ec == 0, open(os.path.join(WORK, "bundle.log")).read()[-3000:]
    try:
        y = glob.glob(os.path.join(bundle, "*.yaml")) + glob.glob(os.path.join(bundle, "**/*.yaml"), recursive=True)
        if y: res["hash_bundle"] = S.hash_ok(y[0])
    except Exception as e:
        res["hash_bundle"] = {"error": repr(e)}
    guarda(res); log(f"bundle listo ({time.time()-t0:.0f}s)")

    # ── los brazos
    brazos = [("k1.0", "10_rasgos", 1.0), ("k3.4", "10_rasgos", 3.4),
              ("r8_k1.0", "8_rasgos", 1.0), ("k2.0", "10_rasgos", 2.0),
              ("k1.5", "10_rasgos", 1.5), ("k2.5", "10_rasgos", 2.5)]
    for nom, mod, k in brazos:
        if nom in res["brazos"] and res["brazos"][nom].get("estado") == "OK":
            continue
        t1 = time.time()
        try:
            if libre() < MIN_FREE_GIB:
                res["brazos"][nom] = {"estado": "ABORTADO por disco", "libre_GiB": libre()}
                guarda(res); break
            Yk = np.clip(Ypred[mod] * k, -RECORTE, RECORTE)
            raz = np.exp2(Yk)
            pred = os.path.join(WORK, f"pred_{nom}.h5ad")
            rd = os.path.join(WORK, f"run_{nom}")
            info = lado_pred(Xc, retenidos, Yk, genes, SEM["pred_" + nom], pred)
            ec = corre(S.cmd_run(pred, real, rd), os.path.join(WORK, f"run_{nom}.log"))
            assert ec == 0, open(os.path.join(WORK, f"run_{nom}.log")).read()[-2500:]
            csv = os.path.join(WORK, f"scored_{nom}.csv")
            ec2 = corre(S.cmd_score(rd, bundle, csv), os.path.join(WORK, f"score_{nom}.log"))
            res["brazos"][nom] = {
                "estado": "OK", "modelo": mod, "k": k, "semilla": SEM["pred_" + nom],
                "u": S.lee_u(rd), "s": S.lee_s(csv) if ec2 == 0 else None, "score_ec": ec2,
                "sd_log2_razon": float(Yk.std()), "sd_razon": float(raz.std()),
                "razon_min_max": [float(raz.min()), float(raz.max())],
                "n_recortados": int((np.abs(Ypred[mod] * k) > RECORTE).sum()),
                "pred": info, "cuatro_celdas": cuatro_celdas(rd, retenidos),
                "segundos": round(time.time() - t1, 1)}
            u = res["brazos"][nom]["u"]
            log(f"BRAZO {nom} OK en {time.time()-t1:.0f}s  pds={u.get('pds_cosine')} "
                f"jac={u.get('de_wilcoxon_sig_jaccard')}")
            os.remove(pred)
            subprocess.run(["rm", "-rf", rd + "/degenes"], check=False)
        except Exception as e:
            res["brazos"][nom] = {"estado": f"FALLO: {type(e).__name__}: {e}"[:2000],
                                  "segundos": round(time.time() - t1, 1)}
            log(f"BRAZO {nom} FALLO: {type(e).__name__}: {e}")
        guarda(res)
    res["segundos_total"] = round(time.time() - t0, 1)
    guarda(res)
    log(f"FIN en {res['segundos_total']}s; brazos OK: "
        f"{[n for n,v in res['brazos'].items() if v.get('estado')=='OK']}")

if __name__ == "__main__":
    main()
