#!/usr/bin/env python
"""PASO 3 — COMPROBAR EL CONTRATO del entregable ANTES de empaquetar.

Se comprueba sobre el h5ad escrito, leyendo X por bloques (nunca entero):
  1. forma EXACTA 360,000 x 18,533
  2. `var` = los 18,533 nombres en el ORDEN de gene_names.csv, uno a uno
  3. `obs.target_gene` = simbolos de gen; las 300 de pert_counts.csv; CERO
     filas `non-targeting`; 1,200 filas por perturbacion (400 x 3 contextos)
  4. `obs.context` = D, E, F exactamente como en los archivos de control;
     120,000 filas cada uno
  5. X entera, no negativa, finita
  6. ninguna celula con suma > 1,000,000
  7. nnz <= 4,750,000,000 y margen
  8. CSR
"""
import json
import os
import sys

import h5py
import numpy as np

sys.path.insert(0, "/home/ubuntu/vcc/modelo/src/noche")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import comun as C  # noqa: E402
import datos as D  # noqa: E402

H5 = "/datos/modelo/out/entrega/prediction_k3.4.h5ad"
OUT = "/home/ubuntu/vcc/modelo/out/entrega_contrato.json"
N_FILAS = 360000
N_GENES = 18533
MAX_SUMA = 1_000_000
MAX_NNZ = 4_750_000_000
CTX = ["D", "E", "F"]


def main():
    genes = D.genes_2026()
    perts = D.perts_2026()
    r = {"archivo": H5, "tam_bytes": os.path.getsize(H5),
         "tam_GiB": round(os.path.getsize(H5) / 2 ** 30, 3), "comprobaciones": []}

    def chk(n, ok, medido, esperado):
        r["comprobaciones"].append({"#": len(r["comprobaciones"]) + 1, "que": n,
                                    "medido": medido, "esperado": esperado,
                                    "ok": bool(ok)})
        C.log(f"{'OK ' if ok else 'MAL'} {n}: {medido}")

    with h5py.File(H5, "r") as f:
        X = f["X"]
        shape = [int(x) for x in X.attrs["shape"]]
        chk("encoding de X", X.attrs.get("encoding-type") == "csr_matrix",
            X.attrs.get("encoding-type"), "csr_matrix")
        chk("forma", shape == [N_FILAS, N_GENES], shape, [N_FILAS, N_GENES])
        vnode = f["var"]["_index"]
        vals = vnode["values"][:] if hasattr(vnode, "keys") else vnode[:]
        vn = [v.decode() if isinstance(v, bytes) else str(v) for v in vals]
        igual = len(vn) == len(genes) and all(a == b for a, b in zip(vn, genes))
        primer = next((i for i, (a, b) in enumerate(zip(vn, genes)) if a != b), None)
        chk("var = gene_names.csv, mismo orden elemento a elemento", igual,
            f"{len(vn)} nombres, primer desajuste: {primer}", f"{len(genes)}, None")
        tg = f["obs"]["target_gene"]
        tg = np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in tg[:]])
        cx = f["obs"]["context"]
        cx = np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in cx[:]])
        chk("filas de obs", len(tg) == N_FILAS and len(cx) == N_FILAS,
            [len(tg), len(cx)], [N_FILAS, N_FILAS])
        chk("CERO filas non-targeting", int((tg == "non-targeting").sum()) == 0,
            int((tg == "non-targeting").sum()), 0)
        chk("blancos = las 300 de pert_counts.csv", set(tg.tolist()) == set(perts.tolist()),
            f"{len(set(tg.tolist()))} distintos", "300")
        u, cnt = np.unique(tg, return_counts=True)
        chk("1,200 filas por blanco (400 x 3)", bool((cnt == 1200).all()),
            f"min {cnt.min()} max {cnt.max()}", 1200)
        chk("simbolos, no identificadores de constructo",
            not any("-" in x and x.rsplit("-", 1)[-1].isdigit() for x in u),
            "ninguno con sufijo -N", "ninguno")
        uc, cc = np.unique(cx, return_counts=True)
        chk("contextos", sorted(uc.tolist()) == CTX, dict(zip(uc.tolist(), cc.tolist())),
            {c: 120000 for c in CTX})
        chk("120,000 filas por contexto", bool((cc == 120000).all()),
            cc.tolist(), [120000] * 3)
        nnz = int(X["data"].shape[0])
        chk(f"nnz <= {MAX_NNZ:,}", nnz <= MAX_NNZ, f"{nnz:,}", f"<= {MAX_NNZ:,}")
        r["nnz"] = nnz
        r["margen_nnz"] = round(MAX_NNZ / max(nnz, 1), 3)
        r["densidad"] = round(nnz / (N_FILAS * N_GENES), 5)
        r["nnz_por_celula"] = round(nnz / N_FILAS, 1)
        # ── X por bloques
        ip = X["indptr"][:].astype(np.int64)
        chk("indptr consistente", ip[0] == 0 and int(ip[-1]) == nnz and len(ip) == N_FILAS + 1,
            [int(ip[0]), int(ip[-1]), len(ip)], [0, nnz, N_FILAS + 1])
        bl = 20000
        mn, mx, max_suma, no_ent, no_fin, neg = np.inf, -np.inf, 0.0, 0, 0, 0
        idx_min, idx_max = N_GENES, -1
        suma_total = 0.0
        for i0 in range(0, N_FILAS, bl):
            i1 = min(i0 + bl, N_FILAS)
            a, b = int(ip[i0]), int(ip[i1])
            d = X["data"][a:b]
            ii = X["indices"][a:b]
            if d.size:
                mn = min(mn, float(d.min())); mx = max(mx, float(d.max()))
                no_ent += int((d != np.rint(d)).sum())
                no_fin += int((~np.isfinite(d)).sum())
                neg += int((d < 0).sum())
                idx_min = min(idx_min, int(ii.min())); idx_max = max(idx_max, int(ii.max()))
                s = np.add.reduceat(d.astype(np.float64),
                                    np.r_[0, np.cumsum(np.diff(ip[i0:i1 + 1]))[:-1]])
                s[np.diff(ip[i0:i1 + 1]) == 0] = 0.0
                max_suma = max(max_suma, float(s.max()))
                suma_total += float(s.sum())
        chk("X entera", no_ent == 0, f"{no_ent} valores no enteros", 0)
        chk("X no negativa", neg == 0 and mn >= 0, f"min {mn}", ">= 0")
        chk("X finita", no_fin == 0, f"{no_fin} no finitos", 0)
        chk("indices dentro de rango", 0 <= idx_min and idx_max < N_GENES,
            [idx_min, idx_max], [0, N_GENES - 1])
        chk(f"ninguna celula suma > {MAX_SUMA:,}", max_suma <= MAX_SUMA,
            f"maxima {max_suma:,.0f}", f"<= {MAX_SUMA:,}")
        r["max_valor"] = mx
        r["suma_mediana_aprox"] = round(suma_total / N_FILAS, 1)
    r["todas_ok"] = all(c["ok"] for c in r["comprobaciones"])
    r["n_ok"] = sum(1 for c in r["comprobaciones"] if c["ok"])
    r["n_total"] = len(r["comprobaciones"])
    C.escribe_json(OUT, r)
    C.log(f"CONTRATO: {r['n_ok']}/{r['n_total']} comprobaciones OK | "
          f"nnz {r['nnz']:,} (margen x{r['margen_nnz']}) | {r['tam_GiB']} GiB")
    if not r["todas_ok"]:
        sys.exit(2)


if __name__ == "__main__":
    main()
