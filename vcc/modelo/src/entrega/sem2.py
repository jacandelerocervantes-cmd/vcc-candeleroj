#!/usr/bin/env python
"""Semillas de LA ENTREGA. Familia NUEVA 20284000+, fold DESDE 1.

  fold: D=1, E=2, F=3   (los tres contextos de la entrega)
  20284000 + 100*fold + 0   remuestreo del pool de control de ese contexto
                            (las 400 celulas de cada perturbacion)
  20284000 + 100*fold + 1   ruido de la inyeccion (`_reescala_entero`)
  20284900 + k              particion por BLANCO de la validacion cruzada en
                            2025 (k = 0..4, 5 pliegues)
  20284950                  orden de barajado de los blancos antes de partir

Comprobadas EN CODIGO contra todo lo gastado (rangos vetados hasta 20283999,
la enumeracion de `noche/comun.py`, PASO 5a, el suelo de llamada, la
sensibilidad, y un barrido de enteros 202[6-8]xxxx por `modelo/src/**/*.py` y
los json de `modelo/out/`). Se guardan ENTERAS.
"""
import json
import os
import re
import sys

sys.path.insert(0, "/home/ubuntu/vcc/modelo/src/noche")
sys.path.insert(0, "/home/ubuntu/vcc/modelo/src/seleccion")
sys.path.insert(0, "/home/ubuntu/vcc/modelo/src/suelo_llamada")
sys.path.insert(0, "/home/ubuntu/vcc/modelo/src/sensib")
import comun as C  # noqa: E402
import semillas as S5a  # noqa: E402
import sl_semillas as SL  # noqa: E402
import sem as SEN  # noqa: E402

CTX = ["D", "E", "F"]
FOLD = {c: i + 1 for i, c in enumerate(CTX)}
N_PLIEGUES = 5


def s_control(ctx):
    return 20284000 + 100 * FOLD[ctx] + 0


def s_ruido(ctx):
    return 20284000 + 100 * FOLD[ctx] + 1


def s_pliegue(k):
    return 20284900 + k


S_BARAJA = 20284950

SEMILLAS = {
    **{f"control_{c}": [s_control(c)] for c in CTX},
    **{f"ruido_{c}": [s_ruido(c)] for c in CTX},
    "cv_pliegues": [s_pliegue(k) for k in range(N_PLIEGUES)],
    "baraja_blancos": [S_BARAJA],
}
VETADOS = [(20260000, 20266009), (20270000, 20283999)]
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
            if "entrega" in d:
                continue
            for f in fs:
                if f.endswith(".py"):
                    p = os.path.join(d, f)
                    for m in set(PAT.findall(open(p, "rb").read())):
                        vis.setdefault(int(m), p)
    for d, _s, fs in os.walk("/home/ubuntu/vcc/modelo/out"):
        for f in fs:
            if not f.endswith(".json") or f.startswith("entrega"):
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
    gast = (C._gastadas_explicitas()
            | {int(v) for vs in C.SEMILLAS_NUEVAS.values() for v in vs}
            | {int(v) for vs in S5a.SEMILLAS.values() for v in vs}
            | {int(v) for vs in SL.SEMILLAS.values() for v in vs}
            | {int(v) for vs in SEN.SEMILLAS.values() for v in vs})
    barr = _barrido()
    en_barr = {k: barr[k] for k in sorted(set(usos) & set(barr))}
    return {"familia": "20284000 + 100*fold + s, fold desde 1 (D=1, E=2, F=3)",
            "semillas_enteras": SEMILLAS, "lista_plana": sorted(usos), "n": len(usos),
            "rangos_vetados": VETADOS, "choques_internos": internos,
            "choques_veto": en_veto, "choques_gastadas": sorted(set(usos) & gast),
            "n_gastadas_enumeradas": len(gast), "n_barrido": len(barr),
            "choques_barrido": en_barr,
            "sin_choques": not (internos or en_veto or (set(usos) & gast) or en_barr)}


if __name__ == "__main__":
    r = comprueba()
    print(json.dumps({k: v for k, v in r.items() if k != "semillas_enteras"}, indent=1))
