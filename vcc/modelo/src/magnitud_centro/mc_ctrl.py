#!/usr/bin/env python
"""PASO 2 / parte 3: el CONTROL NEGATIVO de la fraccion CENTRO.

🔴 POR QUE HACE FALTA. `V_sh = <m1, m2>/G` mide la norma reproducible de la
MEDIA EMPIRICA de los deltas sobre las lineas. Si las lineas no compartieran
NADA del delta especifico del blanco, esa media empirica **no seria cero**: la
media de `n` vectores mutuamente incorrelados conserva ~`1/n` de la norma
cuadrada tipica, y el producto cruzado de mitades NO lo quita, porque el efecto
especifico de linea es la misma cosa real en las dos mitades (son celulas
disjuntas de la misma linea, no ruido independiente). Asi que `V_sh/V_tot` tiene
un **suelo de no-comparticion** que hay que medir, no suponer.

EL CONTROL, declarado: para cada blanco `g` se toma, en cada linea, el delta de
un blanco DISTINTO elegido por una permutacion INDEPENDIENTE por linea. Eso
destruye cualquier comparticion del delta especifico del blanco entre lineas
**y deja intacto el eje comun al IFNB** (que es compartido entre blancos Y entre
lineas). Es la TRAMPA 2 aplicada a la descomposicion, no solo a las
correlaciones.

  * permutaciones sin punto fijo por linea (`pi_L(g) != g` para todo g), para
    que ningun blanco se compare consigo mismo
  * 20 sorteos, media y rango
  * semilla: `20278900 + sorteo`, de la familia nueva y fuera de las 48 de
    particion y de la del bootstrap

Se reporta junto a la fraccion observada; la DIFERENCIA de las dos es lo que
queda de comparticion atribuible al blanco.
"""
import json, os, sys
import numpy as np

sys.path.insert(0, "/home/ubuntu/vcc/modelo/src")
import gwps_via as GV

LINES = GV.LINES
WORK = "/datos/modelo/out/magnitud_centro"
OUTDIR = "/home/ubuntu/vcc/modelo/out/magnitud_centro"
SEED_PERM = 20278900
N_PERM = 20


def perm_sin_punto_fijo(rng, n):
    while True:
        p = rng.permutation(n)
        if not np.any(p == np.arange(n)):
            return p


def main():
    meta = json.load(open(os.path.join(OUTDIR, "deltas_meta.json")))
    panel = meta["panel"]; P = len(panel)
    D = {}
    for L in LINES:
        D[L] = np.load(os.path.join(WORK, f"deltas_{L}.npz"),
                       allow_pickle=True)["D"].astype(np.float64)
    J, G = D[LINES[0]].shape[0], D[LINES[0]].shape[-1]
    E = {L: D[L] - D[L].mean(axis=1, keepdims=True) for L in LINES}
    sys.stderr.write(f"cargado J={J} P={P} G={G}\n")

    pares = [(a, b) for i, a in enumerate(LINES) for b in LINES[i + 1:]]

    def descompone(Hs, S, perms=None):
        """perms: dict L -> array(P) de indices de blanco a usar en esa linea."""
        nS = len(S)
        Vt = np.zeros((J, P)); Vs = np.zeros((J, P))
        for j in range(J):
            h1, h2 = [], []
            for L in S:
                idx = np.arange(P) if perms is None else perms[L]
                h1.append(Hs[L][j, idx, 0]); h2.append(Hs[L][j, idx, 1])
            h1 = np.stack(h1); h2 = np.stack(h2)
            m1, m2 = h1.mean(axis=0), h2.mean(axis=0)
            Vt[j] = np.einsum("spg,spg->p", h1, h2) / nS / G
            Vs[j] = np.einsum("pg,pg->p", m1, m2) / G
        return Vt, Vs

    res = {"semilla_perm_base": SEED_PERM, "n_permutaciones": N_PERM,
           "n_lineas": len(LINES), "variantes": {}}
    etiquetas = [("las_6_lineas", LINES)] + [(f"{a}-{b}", [a, b]) for a, b in pares]

    for nombre, Hs in (("delta_crudo", D), ("delta_especifico_del_blanco", E)):
        out = {}
        for etiq, S in etiquetas:
            Vt0, Vs0 = descompone(Hs, S)
            obs = float(Vs0.mean() / Vt0.mean())
            fr = []
            for k in range(N_PERM):
                rng = np.random.default_rng(SEED_PERM + k)
                perms = {L: perm_sin_punto_fijo(rng, P) for L in S}
                Vt, Vs = descompone(Hs, S, perms)
                fr.append(float(Vs.mean() / Vt.mean()))
            fr = np.array(fr)
            n = len(S)
            out[etiq] = {
                "n_lineas": n,
                "frac_observada": obs,
                "frac_control_permutado_media": float(fr.mean()),
                "frac_control_permutado_sd": float(fr.std(ddof=1)),
                "frac_control_permutado_min": float(fr.min()),
                "frac_control_permutado_max": float(fr.max()),
                "suelo_analitico_1_sobre_n": 1.0 / n,
                "exceso_sobre_el_control": obs - float(fr.mean()),
                "frac_anclada_en_el_control": (
                    (obs - float(fr.mean())) / (1.0 - float(fr.mean()))
                    if fr.mean() < 1 else None),
                "frac_anclada_en_1_sobre_n": (obs - 1.0 / n) / (1.0 - 1.0 / n),
                "los_20_sorteos": [float(x) for x in fr],
            }
            sys.stderr.write(
                f"  [{nombre}][{etiq}] obs={obs:.4f} ctrl={fr.mean():.4f} "
                f"(sd {fr.std(ddof=1):.4f}) 1/n={1.0/n:.4f} "
                f"anclada={out[etiq]['frac_anclada_en_el_control']:.4f}\n")
        res["variantes"][nombre] = out

    json.dump(res, open(os.path.join(OUTDIR, "ctrl.json"), "w"),
              indent=1, sort_keys=True, default=str)
    sys.stderr.write(f"\nlisto: {OUTDIR}/ctrl.json\n")


if __name__ == "__main__":
    main()
