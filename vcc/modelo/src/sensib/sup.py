#!/usr/bin/env python
"""LA SUPERFICIE DE ENTREGA como banco de sensibilidad. Solo controles.

🔴 ESTO NO ES BIOLOGIA DE D/E/F. Es la sensibilidad del MARCADOR sobre los
controles de la entrega. La biologia de D/E/F no esta en disco.

LA CONSTRUCCION, DECLARADA ANTES DE MEDIR
  * Cada contexto trae 18,400 celulas de control, 46 `ntc_id` de 400 exactas
    (verificado en paso1_contrato.md), 18,533 genes, ~20,000 UMI de mediana.
  * PARTICION POR GUIA, nunca por celula (como ntc_nulo.md §1.1): las 46 guias
    ordenadas alfabeticamente se permutan con `s_particion(ctx)`; las primeras
    23 son el POOL DE CONTROL (9,200 celulas, etiqueta `non-targeting`) y las
    otras 23 son PSEUDO-PERTURBACIONES de 400 celulas (`PP01..PP23`).
  * Regimen de la llamada: 400 contra 9,200. Es el brazo `n400_ref_completa`
    de ntc_nulo.md, el declarado ahi como el que se parece a la entrega, y el
    mismo en el que ese reporte midio que el nulo da CERO.
  * LADO REAL: a cada pseudo-perturbacion se le aplica su vector de razones con
    `cand01_covar._reescala_entero` SIN MODIFICAR (Binomial si r<=1, Poisson
    aditivo si r>1; esperanza x*r en las dos ramas). El pool de control va
    entero y sin tocar. 18,400 filas.
  * LADO PRED: 400 celulas por etiqueta remuestreadas CON REEMPLAZO del pool de
    control (el contrato de `suelo_arnes`), con el vector de razones que le
    toque al brazo. 23*400 + 400 = 9,600 filas.
  * 🔴 LAS ETIQUETAS SON NOMBRES DE GEN DEL PROPIO RETO. El marcador exige que
    cada perturbacion resuelva a un gen del panel (`exclude_target_gene=True`:
    si no resuelve ninguna, `baseline` se niega y lo dice). Asi que las 23
    pseudo-perturbaciones llevan 23 de las 300 etiquetas del reto, sorteadas con
    semilla declarada entre las que son gen del panel (las 300 lo son), las
    MISMAS en los tres contextos de una particion. Con eso la exclusion del gen
    propio y la del panel operan como en la entrega.
  * 🔴 El gen propio de cada etiqueta NO se mueve, y los 23 genes-etiqueta
    quedan FUERA del sorteo de genes movidos: la exclusion del panel los saca
    del espacio de rasgos, asi que mover ahi seria senal que ninguna metrica
    puede ver. No se modela eficiencia de knockdown: no es el asunto aqui.
  * Los genes movibles son los que pasan el filtro del preset
    (`filter_gene_min_cpm_cell = 5.0`) MEDIDO en el pool de control de ese
    contexto: un gen que la puerta descarta no puede detectarse nunca.
"""
import json
import os
import sys

import h5py
import numpy as np
import scipy.sparse as sp

sys.path.insert(0, "/home/ubuntu/vcc/modelo/src/noche")
sys.path.insert(0, "/home/ubuntu/vcc/modelo/src")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import comun as C  # noqa: E402
import sem as S  # noqa: E402
from cand01_covar import _reescala_entero  # noqa: E402  (sin modificar)
from suelo_arnes import decode_categorical  # noqa: E402

RETO = {"A": "/datos/reto/val/context_A.h5ad", "B": "/datos/reto/val/context_B.h5ad",
        "C": "/datos/reto/val/context_C.h5ad", "D": "/datos/reto/final/context_D.h5ad",
        "E": "/datos/reto/final/context_E.h5ad", "F": "/datos/reto/final/context_F.h5ad"}
WORK = "/datos/modelo/out/sensib"
CTRL = "non-targeting"
N_POR = 400
N_GUIAS = 46
N_CTRL_GUIAS = 23
UMBRAL_CPM = 5.0
MET = C.VCC2026_METRICS if hasattr(C, "VCC2026_METRICS") else [
    "pds_cosine", "expr_mse_unbiased_capped_norm", "de_wilcoxon_lfc_nmae",
    "de_wilcoxon_direction_fidelity_yield_raw", "de_wilcoxon_direction_reach_raw",
    "de_wilcoxon_sig_jaccard"]
CELLEVAL2 = os.path.expanduser("~/celleval2_venv/bin/cell-eval2")


def log(m):
    C.log(m)


# ───────────────────────────── la superficie ─────────────────────────────
def lee_contexto(ctx):
    """X (csr int64), ntc_id (str array), var_names. No densifica nunca."""
    p = RETO[ctx]
    with h5py.File(p, "r") as f:
        shape = tuple(int(x) for x in f["X"].attrs["shape"])
        data = f["X"]["data"][:]
        indices = f["X"]["indices"][:]
        indptr = f["X"]["indptr"][:]
        ntc = np.array(decode_categorical(f["obs"]["ntc_id"]), dtype=object)
        tg = np.array(decode_categorical(f["obs"]["target_gene"]), dtype=object)
        node = f["var"]["_index"]
        vals = node["values"][:] if hasattr(node, "keys") else node[:]
        var_names = [v.decode() if isinstance(v, bytes) else v for v in vals]
    assert shape == (18400, 18533), shape
    assert set(tg) == {CTRL}, sorted(set(tg))[:3]
    ent = np.all(data == np.rint(data))
    X = sp.csr_matrix((data.astype(np.int64), indices, indptr), shape=shape)
    return X, ntc.astype(str), np.asarray(var_names, dtype=str), {
        "ruta": p, "shape": list(shape), "nnz": int(X.nnz), "X_entero": bool(ent),
        "n_guias": int(len(set(ntc))), "celulas_por_guia": sorted(set(
            int((ntc == g).sum()) for g in set(ntc)))}


def etiquetas_reto(ctx, var_names):
    """23 etiquetas = 23 nombres de gen del reto, sorteadas con semilla
    declarada entre las 300 de la particion de ese contexto."""
    import pandas as pd
    part = "final" if ctx in ("D", "E", "F") else "val"
    base = "/datos/reto/" + part + "/pert_counts.csv"
    pc = pd.read_csv(base)["target_gene"].astype(str).tolist()
    cand = sorted(set(pc) & set(map(str, var_names)))
    rng = np.random.default_rng(S.S_ETIQUETAS[part])
    sel = sorted(rng.choice(np.asarray(cand, dtype=object), size=N_CTRL_GUIAS,
                            replace=False).tolist())
    return sel, {"particion_reto": part, "n_perts_reto": len(pc),
                 "n_candidatas_que_son_gen": len(cand),
                 "semilla_etiquetas": int(S.S_ETIQUETAS[part]), "etiquetas": sel}


def particion(ctx, ntc, var_names):
    """23 guias de control / 23 pseudo-perturbaciones. Nunca parte una guia."""
    guias = sorted(set(ntc))
    assert len(guias) == N_GUIAS, len(guias)
    rng = np.random.default_rng(S.s_particion(ctx))
    perm = rng.permutation(N_GUIAS)
    g_ctrl = sorted(guias[i] for i in perm[:N_CTRL_GUIAS])
    g_pert = sorted(guias[i] for i in perm[N_CTRL_GUIAS:])
    assert not (set(g_ctrl) & set(g_pert))
    i_ctrl = np.where(np.isin(ntc, g_ctrl))[0]
    etiq, diag_e = etiquetas_reto(ctx, var_names)
    pp, mapa = {}, {}
    for k, g in enumerate(g_pert):
        idx = np.where(ntc == g)[0]
        assert idx.size == N_POR, (g, idx.size)
        pp[etiq[k]] = idx
        mapa[etiq[k]] = g
    return {"guias_control": g_ctrl, "guias_pert": g_pert, "idx_control": i_ctrl,
            "pp": pp, "etiquetas": etiq, "guia_de_etiqueta": mapa,
            "diag_etiquetas": diag_e, "semilla": S.s_particion(ctx),
            "n_control": int(i_ctrl.size), "n_pp": len(pp),
            "guia_partida": False}


def filtrables(Xc, var_names, fuera=()):
    """CPM medio por gen del pool de control > 5 (el filtro del preset), menos
    los genes-etiqueta (`fuera`), que la exclusion del panel saca del espacio."""
    Xw = Xc.astype(np.float64).tocsr()
    n = Xw.shape[0]
    L = np.asarray(Xw.sum(axis=1)).ravel()
    L[L <= 0] = 1.0
    cpm = np.asarray(Xw.T @ ((1e6 / L) / n)).ravel()
    ok = np.where(cpm > UMBRAL_CPM)[0]
    if len(fuera):
        vn = np.asarray(var_names, dtype=str)
        mal = np.where(np.isin(vn, np.asarray(list(fuera), dtype=str)))[0]
        ok = np.setdiff1d(ok, mal, assume_unique=False)
    return ok, cpm


# ───────────────────────── la estructura inyectada ─────────────────────────
def estructura(k, etiquetas, filtr, n_genes, amplitud, frac_genes, frac_comun, seeds=None):
    """Devuelve (log2 ratios por etiqueta, diagnostico). Rejilla DECLARADA.

    amplitud   = |log2FC| inyectado (misma magnitud en todos los genes movidos)
    frac_genes = fraccion de los genes FILTRABLES que se mueven por pseudo-pert
    frac_comun = fraccion de esos genes que es COMUN a las 23 (mismo gen, mismo
                 signo, misma amplitud); el resto es ESPECIFICO, sorteado
                 aparte para cada pseudo-perturbacion y disjunto del comun.
    """
    sg, ss = seeds if seeds is not None else S.s_estructura(k)
    rg, rs = np.random.default_rng(sg), np.random.default_rng(ss)
    m = int(round(frac_genes * filtr.size))
    n_com = int(round(frac_comun * m))
    n_esp = m - n_com
    com = rg.choice(filtr, size=n_com, replace=False) if n_com else np.array([], dtype=int)
    resto = np.setdiff1d(filtr, com, assume_unique=False)
    s_com = rs.choice([-1.0, 1.0], size=n_com) if n_com else np.array([])
    L = {}
    esp_idx = {}
    for e in etiquetas:
        v = np.zeros(n_genes, dtype=np.float64)
        if n_com:
            v[com] = s_com * amplitud
        if n_esp:
            pick = rg.choice(resto, size=n_esp, replace=False)
            v[pick] = rs.choice([-1.0, 1.0], size=n_esp) * amplitud
            esp_idx[e] = pick
        else:
            esp_idx[e] = np.array([], dtype=int)
        L[e] = v
    diag = {"celda": k, "semilla_genes": sg, "semilla_signos": ss,
            "amplitud_log2": amplitud, "frac_genes": frac_genes, "frac_comun": frac_comun,
            "n_filtrables": int(filtr.size), "m_movidos": m, "n_comun": n_com,
            "n_especifico": n_esp,
            "solape_esp_entre_pp": int(np.mean([
                np.intersect1d(esp_idx[a], esp_idx[b]).size
                for a, b in zip(etiquetas[:-1], etiquetas[1:])])) if n_esp else 0,
            "genes_comun": sorted(int(x) for x in com[:50]),
            "n_genes_comun_guardados": min(50, n_com)}
    return L, esp_idx, com, s_com, diag


def razones_pred(L, esp_idx, com, s_com, etiquetas, n_genes, filtr, modo, f, amplitud, rng):
    """El lado pred. El eje COMUN se da por acertado en los dos modos: `f`
    modula SOLO la componente especifica, que es lo que pide C3.

      suelo     : sin estructura (razon 1 en todo gen)
      magnitud  : comun exacto + f * especifico (direccion correcta, encogida)
      identidad : comun exacto + especifico con floor(f*n_esp) genes CORRECTOS y
                  el resto en genes EQUIVOCADOS (sorteados fuera del verdadero),
                  misma amplitud y signo sorteado
    """
    out = {}
    base_com = np.zeros(n_genes)
    if com.size:
        base_com[com] = s_com * amplitud
    for e in etiquetas:
        v = np.zeros(n_genes)
        if modo == "suelo":
            pass
        elif modo == "magnitud":
            v = base_com + f * (L[e] - base_com)
        elif modo == "identidad":
            v = base_com.copy()
            tv = esp_idx[e]
            nok = int(np.floor(f * tv.size))
            if tv.size:
                ok = rng.choice(tv, size=nok, replace=False) if nok else np.array([], dtype=int)
                v[ok] = L[e][ok]
                nmal = tv.size - nok
                if nmal:
                    pool = np.setdiff1d(filtr, np.union1d(com, tv), assume_unique=False)
                    mal = rng.choice(pool, size=min(nmal, pool.size), replace=False)
                    v[mal] = rng.choice([-1.0, 1.0], size=mal.size) * amplitud
        else:
            raise ValueError(modo)
        out[e] = v
    return out


def aplica(bloque, log2r, rng):
    if not np.any(log2r):
        return bloque.tocsr()
    return _reescala_entero(bloque, np.exp2(log2r), rng)


def escribe_real(X, par, L, semilla_ruido, path):
    """18,400 filas: 23 pseudo-perts inyectadas + el pool de control entero."""
    import anndata as ad
    import pandas as pd
    rng = np.random.default_rng(semilla_ruido)
    bl, lab = [], []
    for e in sorted(par["pp"]):
        bl.append(aplica(X[par["pp"][e]], L[e], rng))
        lab += [e] * N_POR
    bl.append(X[par["idx_control"]])
    lab += [CTRL] * par["n_control"]
    M = sp.vstack(bl, format="csr")
    a = ad.AnnData(X=M.astype(np.int64), obs=pd.DataFrame({"target_gene": lab}),
                   var=pd.DataFrame(index=list(par["var_names"])))
    a.obs_names_make_unique()
    a.write_h5ad(path)
    return {"shape": list(a.shape), "nnz": int(M.nnz), "huella": C.huella_adata(a),
            "semilla_ruido": int(semilla_ruido)}


def escribe_pred(Xc, par, razones, ctx, slot_seed, path):
    """9,600 filas: 400 por pseudo-pert + 400 de control, remuestreadas CON
    REEMPLAZO del pool de control (contrato de suelo_arnes)."""
    import anndata as ad
    import pandas as pd
    rng = np.random.default_rng(slot_seed)
    n = Xc.shape[0]
    bl, lab = [], []
    for e in sorted(par["pp"]) + [CTRL]:
        ch = rng.integers(0, n, size=N_POR)
        b = Xc[ch]
        bl.append(b.tocsr() if e == CTRL else aplica(b, razones[e], rng))
        lab += [e] * N_POR
    M = sp.vstack(bl, format="csr")
    a = ad.AnnData(X=M.astype(np.int64), obs=pd.DataFrame({"target_gene": lab}),
                   var=pd.DataFrame(index=list(par["var_names"])))
    a.obs_names_make_unique()
    a.write_h5ad(path)
    return {"shape": list(a.shape), "nnz": int(M.nnz), "huella": C.huella_adata(a),
            "semilla": int(slot_seed)}


# ─────────────────────────── el marcador entero ───────────────────────────
def cmd_baseline(real, outdir, pred):
    return [CELLEVAL2, "baseline", "-ar", real, "--preset", "vcc2026", "--profile", "vcc2026",
            "--pert-col", "target_gene", "--emit", "dispersed", "--seed", "0",
            "--save-pred", pred, "-o", outdir]


def cmd_bundle(real, base_pred, bundle):
    return [CELLEVAL2, "prep-real-bundle", "--real", real, "--baseline", base_pred,
            "--preset", "vcc2026", "--profile", "vcc2026", "--pert-col", "target_gene",
            "--anchor-base-seed", "0", "--anchor-splits", "5", "-o", bundle]


def cmd_run(pred, real, outdir):
    return [CELLEVAL2, "run", "-ap", pred, "-ar", real, "--preset", "vcc2026",
            "--profile", "vcc2026", "--pert-col", "target_gene", "--write-degenes",
            "-o", outdir]


def cmd_score(run_dir, bundle, out_csv):
    return [CELLEVAL2, "score", "--user-agg", os.path.join(run_dir, "agg_results.csv"),
            "--user-meta", os.path.join(run_dir, "run_meta.json"),
            "--real-bundle", bundle, "-o", out_csv]


def hash_ok(path_yaml):
    from cell_eval2.config import EvalConfig
    from cell_eval2 import competition
    cfg = EvalConfig.from_yaml(path_yaml)
    o = competition.rule_config_hash(cfg)
    r = competition.competition_payload()["config_hash"]
    return {"obtenido": o, "referencia": r, "coincide": o == r}


def lee_u(run_dir):
    import pandas as pd
    df = pd.read_csv(os.path.join(run_dir, "agg_results.csv"))
    fila = df.loc[df["statistic"] == "mean"].iloc[0]
    return {m: (None if m not in df.columns or pd.isna(fila[m]) else float(fila[m])) for m in MET}


def lee_s(csv):
    import pandas as pd
    if not os.path.exists(csv):
        return None
    df = pd.read_csv(csv)
    cols = [c for c in df.columns if c != "metric"]
    out = {}
    for _, r in df.iterrows():
        out[str(r["metric"])] = {c: (None if pd.isna(r[c]) else
                                     (float(r[c]) if isinstance(r[c], (int, float, np.floating))
                                      else r[c])) for c in cols}
    return out
