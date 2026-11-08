#!/usr/bin/env python
"""Desglose del discriminador: ¿lo sostiene la COEXPRESION o una propension por
gen recalibrada blanco a blanco?

POR QUE. El discriminador intra-blanco da +0.18 (GB) / +0.36 (ridge), y los pesos
mas grandes son de SOLO-h (`h_frac_cero`, `h_var`, `h_disp`), no de par. Eso
importa, porque las dos cosas no son la misma respuesta al encargo:

  rasgos de PAR (cov, Pearson, Spearman, parcial, modulo, embedding, familia):
      "los genes coexpresados con g se mueven" -- es lo que se pregunta.
  rasgos de SOLO-h (media, varianza, dispersion, fraccion de ceros, rango):
      "este gen tiende a moverse" -- NO depende del blanco. Dentro de UN blanco
      puede predecir; promediado sobre blancos vale cero exacto, porque
      `Σ_g E[h,g] = 0`. No es habilidad especifica del blanco: es una propension
      del gen con un coeficiente que cada blanco ajusta a su manera.

Asi que se mide el discriminador con los dos conjuntos separados. Si lo sostienen
los de solo-h, el +0.36 no es "otra senal del control que capta lo especifico".

El gen blanco va fuera de rasgos, filas y evaluacion en los tres casos.

🔴 TRANSDUCTIVO. Nada de aqui sale a un envio.
"""
import json
import os
import sys

import numpy as np

LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
DIR = "/home/ubuntu/vcc/modelo/out/control_o_no"
OUT = os.path.join(DIR, "desglose.json")
SEED_FAMILIA = 20268000
ALPHA = 10.0
N_PERM = 10
PAR = ["cov", "r_pearson", "r_spearman", "r_parcial_10", "r_parcial_50",
       "abs_r_pearson", "abs_r_spearman", "delta_r", "rango_cor",
       "mismo_modulo", "dist_emb", "cos_emb", "misma_familia"]
SOLO_H = ["h_media", "h_var", "h_disp", "h_frac_cero", "h_rango_expr"]


def main():
    sal = {}
    for line in LINES:
        fold = LINES.index(line)
        seed = SEED_FAMILIA + 100 * fold + 60
        d = np.load(os.path.join(DIR, "rasgos_%s.npz" % line), allow_pickle=True)
        F = d["F"].astype(np.float64)
        E1 = d["E1"].astype(np.float64); E2 = d["E2"].astype(np.float64)
        rasgos = [str(x) for x in d["rasgos"]]
        G, nb, P = F.shape
        auto = F[:, :, rasgos.index("es_el_blanco")] > 0.5
        conj = {"par": [i for i, r in enumerate(rasgos) if r in PAR],
                "solo_h": [i for i, r in enumerate(rasgos) if r in SOLO_H],
                "ambos": [i for i, r in enumerate(rasgos)
                          if r in PAR or r in SOLO_H]}
        ent = {}
        for nom, cols in conj.items():
          for modo in ("obs", "nulo"):
            reps = 1 if modo == "obs" else N_PERM
            acum = []
            for rep in range(reps):
              rp = np.random.default_rng(seed + 8000 + rep)
              N, D, por_b = 0.0, 0.0, []
              for j in range(nb):
                m = ~auto[:, j]
                X = F[m][:, j, :][:, cols]
                n_ = int(m.sum())
                o = np.arange(n_) if modo == "obs" else rp.permutation(n_)
                e1 = E1[m, j][o]; e2 = E2[m, j][o]
                rng = np.random.default_rng(seed + j)
                mitad = rng.permutation(n_) < n_ // 2
                rs = []
                for sel in (mitad, ~mitad):
                    tr, te = ~sel, sel
                    y = (e1[tr] + e2[tr]) / 2.0; y = y - y.mean()
                    mu = X[tr].mean(0); sg = X[tr].std(0); sg[sg == 0] = 1.0
                    Zt = np.column_stack([np.ones(int(tr.sum())),
                                          (X[tr] - mu) / sg])
                    Ze = np.column_stack([np.ones(int(te.sum())),
                                          (X[te] - mu) / sg])
                    A = Zt.T @ Zt + np.eye(Zt.shape[1]) * ALPHA
                    A[0, 0] -= ALPHA
                    p = Ze @ np.linalg.solve(A, Zt.T @ y)
                    pc = p - p.mean()
                    a = e1[te] - e1[te].mean(); b = e2[te] - e2[te].mean()
                    c = float(a @ b); nn = np.linalg.norm(pc)
                    if c > 0 and nn > 0:
                        num = float(pc @ ((a + b) / 2.0))
                        rs.append(num / (nn * np.sqrt(c)))
                        N += num; D += nn * np.sqrt(c)
                por_b.append(float(np.nanmean(rs)) if rs else np.nan)
              acum.append({"R_agr": (N / D) if D > 0 else np.nan,
                           "mediana": float(np.nanmedian(por_b)),
                           "n_positivos": int(np.nansum(np.array(por_b) > 0)),
                           "por_blanco": por_b})
            if modo == "obs":
                ent[nom] = acum[0]
            else:
                ent[nom]["nulo_R_agr"] = [a["R_agr"] for a in acum]
        sal[line] = ent
        sd = float(np.nanstd(ent["par"]["nulo_R_agr"], ddof=1))
        sys.stderr.write(
            "%-6s  par=%+.4f (pos %d/%d, nulo %+.4f sd %.4f z %+.1f)   "
            "solo_h=%+.4f (pos %d/%d)   ambos=%+.4f\n"
            % (line, ent["par"]["R_agr"], ent["par"]["n_positivos"], nb,
               float(np.nanmean(ent["par"]["nulo_R_agr"])), sd,
               (ent["par"]["R_agr"] - np.nanmean(ent["par"]["nulo_R_agr"])) / sd
               if sd else np.nan,
               ent["solo_h"]["R_agr"], ent["solo_h"]["n_positivos"], nb,
               ent["ambos"]["R_agr"]))
        sys.stderr.flush()
    with open(OUT, "w") as fh:
        json.dump(sal, fh)
    sys.stderr.write("escrito %s\n" % OUT)


if __name__ == "__main__":
    main()
