#!/usr/bin/env python
"""Semillas del MAPA DE SENSIBILIDAD DE LA ENTREGA. Familia NUEVA 20283000+.

Formula del encargo: 20283000 + 100*fold + s, fold DESDE 1.
  fold: D=1, E=2, F=3, A=4, B=5, C=6  (los tres de entrega primero)
  s = 0  particion de las 46 guias en pool de control / pseudo-perturbaciones
  s = 1  emision SUELO del lado pred (remuestreo del pool de control)
  s = 2  emision ORACULO del lado pred
  s = 10+i  emision PARCIAL i del barrido de calidad f (i = 0..9)

Estructura inyectada (la misma en todos los contextos, para que la comparacion
entre contextos no mezcle estructura con contexto):
  20283700 + 2*k      genes movidos de la celda k de la rejilla
  20283701 + 2*k      signos de la celda k
Ruido de inyeccion (Binomial/Poisson de `_reescala_entero`), por contexto y celda:
  20283750 + 25*(fold-1) + k
Brazo C (estructura MEDIDA, 3 puntos): genes/signos 20283960+2*i / 20283961+2*i;
ruido 20283970 + 3*(fold-1) + i. Mitades de 2025 del brazo A2b: 20283920 + j.
Remuestreos de la banda de confianza del brazo A: 20283900 + i
Etiquetas de las 23 pseudo-perturbaciones (nombres de gen del reto): 20283950
(particion `final`, para D/E/F) y 20283951 (`val`, para A/B/C).

COMPROBACION EN CODIGO contra todo lo gastado hasta 20282900:
  (1) rangos vetados [20260000, 20266009] y [20270000, 20282999];
  (2) la enumeracion de `noche/comun.py` y sus SEMILLAS_NUEVAS;
  (3) las semillas de PASO 5a (`seleccion/semillas.py`);
  (4) las del suelo de llamada (`suelo_llamada/sl_semillas.py`);
  (5) un barrido: regex de enteros 202[6-8]xxxx en `modelo/src/**/*.py` y en los
      JSON de `modelo/out/` PARSEANDO el json y recorriendo solo los enteros
      (no digitos de un float: defecto documentado en magnitud_centro.md §1.4).
"""
import json
import os
import re
import sys

sys.path.insert(0, "/home/ubuntu/vcc/modelo/src/noche")
sys.path.insert(0, "/home/ubuntu/vcc/modelo/src/seleccion")
sys.path.insert(0, "/home/ubuntu/vcc/modelo/src/suelo_llamada")
import comun as C  # noqa: E402
import semillas as S5a  # noqa: E402
import sl_semillas as SL  # noqa: E402

CTX = ["D", "E", "F", "A", "B", "C"]
ENTREGA = ["D", "E", "F"]
FOLD = {c: i + 1 for i, c in enumerate(CTX)}
N_CELDAS_MAX = 25
N_F_MAX = 10


def s_ctx(ctx, slot):
    assert 0 <= slot < 100, slot
    return 20283000 + 100 * FOLD[ctx] + slot


def s_particion(ctx):
    return s_ctx(ctx, 0)


def s_suelo(ctx):
    return s_ctx(ctx, 1)


def s_oraculo(ctx):
    return s_ctx(ctx, 2)


def s_parcial(ctx, i):
    assert 0 <= i < N_F_MAX, i
    return s_ctx(ctx, 10 + i)


def s_estructura(k):
    """(genes movidos, signos) de la celda k de la rejilla."""
    assert 0 <= k < N_CELDAS_MAX, k
    return 20283700 + 2 * k, 20283701 + 2 * k


def s_ruido(ctx, k):
    assert 0 <= k < N_CELDAS_MAX, k
    return 20283750 + N_CELDAS_MAX * (FOLD[ctx] - 1) + k


def s_boot(i):
    return 20283900 + i


S_ETIQUETAS = {"final": 20283950, "val": 20283951}


def s_part2025(j):
    """Particiones de mitades de VCC 2025 en el brazo A2b (j = 0..3)."""
    assert 0 <= j < 4, j
    return 20283920 + j


def s_estructura_C(i):
    """(genes, signos) del punto de estructura MEDIDA i del brazo C (i = 0..2)."""
    assert 0 <= i < 3, i
    return 20283960 + 2 * i, 20283961 + 2 * i


def s_ruido_C(ctx, i):
    assert 0 <= i < 3, i
    return 20283970 + 3 * (FOLD[ctx] - 1) + i


SEMILLAS = {
    **{f"particion_{c}": [s_particion(c)] for c in CTX},
    **{f"pred_{c}": [s_suelo(c), s_oraculo(c)] + [s_parcial(c, i) for i in range(N_F_MAX)]
       for c in CTX},
    "estructura_genes": [s_estructura(k)[0] for k in range(N_CELDAS_MAX)],
    "estructura_signos": [s_estructura(k)[1] for k in range(N_CELDAS_MAX)],
    **{f"ruido_{c}": [s_ruido(c, k) for k in range(N_CELDAS_MAX)] for c in CTX},
    "bootstrap_A": [s_boot(i) for i in range(20)],
    "etiquetas_pseudopert": [20283950, 20283951],
    "particion_2025_A2b": [20283920 + j for j in range(4)],
    "estructura_C_genes": [20283960 + 2 * i for i in range(3)],
    "estructura_C_signos": [20283961 + 2 * i for i in range(3)],
    **{f"ruido_C_{c}": [20283970 + 3 * (FOLD[c] - 1) + i for i in range(3)] for c in CTX},
}
VETADOS = [(20260000, 20266009), (20270000, 20282999)]
PAT = re.compile(rb"(?<![0-9.])(202[678][0-9]{4})(?![0-9])")


def _ints(o, out):
    if isinstance(o, bool):
        return
    if isinstance(o, int):
        if 20260000 <= o <= 20289999:
            out.add(o)
    elif isinstance(o, dict):
        for v in o.values():
            _ints(v, out)
    elif isinstance(o, list):
        for v in o:
            _ints(v, out)


def _barrido():
    vis = {}
    for r in ["/home/ubuntu/vcc/modelo/src", "/home/ubuntu/vcc/src"]:
        for d, _s, fs in os.walk(r):
            if "sensib" in d:
                continue
            for f in fs:
                if f.endswith(".py"):
                    p = os.path.join(d, f)
                    for m in set(PAT.findall(open(p, "rb").read())):
                        vis.setdefault(int(m), p)
    for d, _s, fs in os.walk("/home/ubuntu/vcc/modelo/out"):
        for f in fs:
            if not f.endswith(".json") or f.startswith("sensibilidad"):
                continue
            p = os.path.join(d, f)
            s = set()
            try:
                txt = open(p).read()
                try:
                    _ints(json.loads(txt), s)
                except Exception:
                    for ln in txt.splitlines():
                        try:
                            _ints(json.loads(ln), s)
                        except Exception:
                            pass
            except Exception:
                continue
            for k in s:
                vis.setdefault(k, p)
    return vis


def comprueba():
    usos = {}
    for fam, vals in SEMILLAS.items():
        for v in vals:
            usos.setdefault(int(v), []).append(fam)
    internos = {k: v for k, v in usos.items() if len(v) > 1}
    en_veto = sorted(k for k in usos if any(a <= k <= b for a, b in VETADOS))
    expl = C._gastadas_explicitas() | {int(v) for vs in C.SEMILLAS_NUEVAS.values() for v in vs}
    en_expl = sorted(set(usos) & expl)
    s5a = {int(v) for vs in S5a.SEMILLAS.values() for v in vs}
    sll = {int(v) for vs in SL.SEMILLAS.values() for v in vs}
    barr = _barrido()
    en_barr = {k: barr[k] for k in sorted(set(usos) & set(barr))}
    return {"familia": "20283000 + 100*fold + s, fold desde 1 (D=1,E=2,F=3,A=4,B=5,C=6)",
            "semillas_enteras": SEMILLAS, "lista_plana": sorted(usos), "n": len(usos),
            "rangos_vetados": VETADOS, "choques_internos": internos, "choques_veto": en_veto,
            "choques_comun": en_expl, "choques_5a": sorted(set(usos) & s5a),
            "choques_suelo_llamada": sorted(set(usos) & sll),
            "n_barrido": len(barr), "max_barrido": max(barr) if barr else None,
            "choques_barrido": en_barr,
            "sin_choques": not (internos or en_veto or en_expl or (set(usos) & s5a)
                                or (set(usos) & sll) or en_barr)}


if __name__ == "__main__":
    r = comprueba()
    print(json.dumps({k: v for k, v in r.items()
                      if k not in ("semillas_enteras", "lista_plana")}, indent=1))
    print("n semillas:", r["n"], "| min", min(r["lista_plana"]), "max", max(r["lista_plana"]))
