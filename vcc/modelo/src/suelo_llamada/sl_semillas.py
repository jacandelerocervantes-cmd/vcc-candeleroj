#!/usr/bin/env python
"""Semillas del SUELO DE REPRODUCIBILIDAD DE LA LLAMADA. Familia NUEVA 20282000+.

Formula del encargo: 20282000 + 100*fold + j, j = 0..7.
🔴 fold se numera DESDE 1 (Jiang: A549=1 .. MCF7=6; VCC 2025 = 7). Razon: el
rango vetado por el encargo, "2027000x-2028200x", llega hasta 20282009, asi que
fold = 0 (20282000..20282007) caeria DENTRO de el. Con fold desde 1 se cumplen a
la vez la formula y el veto literal.

Desempate del recorte igualado (clave aleatoria por gen): 20282800 + j (Jiang) y
20282810 + j (2025).

Comprobacion EN CODIGO contra:
  (1) los rangos vetados por el encargo: [20260000, 20266009] y [20270000, 20282009];
  (2) la enumeracion de `noche/comun.py` y sus SEMILLAS_NUEVAS;
  (3) las semillas de PASO 5a (`seleccion/semillas.py`);
  (4) un barrido: regex de enteros 202[6-8]xxxx en `modelo/src/**/*.py` y `src/**/*.py`,
      y en los JSON de `modelo/out/` PARSEANDO el json y recorriendo los valores
      (solo enteros), para no tomar digitos de un float (defecto que documento
      magnitud_centro.md §1.4).
"""
import json
import os
import re
import sys

sys.path.insert(0, "/home/ubuntu/vcc/modelo/src/noche")
sys.path.insert(0, "/home/ubuntu/vcc/modelo/src/seleccion")
import comun as C  # noqa: E402
import semillas as S5a  # noqa: E402

LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
FOLD = {L: i + 1 for i, L in enumerate(LINES)}
FOLD["v2025"] = 7
J = 8


def semilla(tag, j):
    return 20282000 + 100 * FOLD[tag] + j


SEMILLAS = {
    **{f"particion_{t}": [semilla(t, j) for j in range(J)] for t in LINES + ["v2025"]},
    "desempate_recorte_jiang": [20282800 + j for j in range(J)],
    "desempate_recorte_2025": [20282810 + j for j in range(J)],
}
VETADOS = [(20260000, 20266009), (20270000, 20282009)]
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
            if "suelo_llamada" in d:
                continue
            for f in fs:
                if f.endswith(".py"):
                    p = os.path.join(d, f)
                    for m in set(PAT.findall(open(p, "rb").read())):
                        vis.setdefault(int(m), p)
    for d, _s, fs in os.walk("/home/ubuntu/vcc/modelo/out"):
        for f in fs:
            if not f.endswith(".json") or f.startswith("suelo_llamada"):
                continue
            p = os.path.join(d, f)
            s = set()
            try:
                with open(p) as fh:
                    txt = fh.read()
                try:
                    _ints(json.loads(txt), s)
                except Exception:
                    for ln in txt.splitlines():          # json por lineas
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
    en_5a = sorted(set(usos) & s5a)
    barr = _barrido()
    en_barr = {k: barr[k] for k in sorted(set(usos) & set(barr))}
    return {"familia": "20282000 + 100*fold + j, fold desde 1 (A549=1..MCF7=6, v2025=7)",
            "semillas_enteras": SEMILLAS, "lista_plana": sorted(usos), "n": len(usos),
            "rangos_vetados": VETADOS, "choques_internos": internos, "choques_veto": en_veto,
            "choques_comun": en_expl, "choques_5a": en_5a,
            "n_barrido": len(barr), "max_barrido": max(barr) if barr else None,
            "choques_barrido": en_barr,
            "sin_choques": not (internos or en_veto or en_expl or en_5a or en_barr)}


if __name__ == "__main__":
    r = comprueba()
    print(json.dumps({k: v for k, v in r.items() if k not in ("semillas_enteras", "lista_plana")}, indent=1))
