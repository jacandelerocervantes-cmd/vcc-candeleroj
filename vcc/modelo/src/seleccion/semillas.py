#!/usr/bin/env python
"""Semillas de PASO 5a. Familia NUEVA 20281000+.

Se declaran ENTERAS antes de medir y se comprueban EN CODIGO contra:
  (1) los rangos gastados nombrados por encargos anteriores,
  (2) la enumeracion explicita de `noche/comun.py`,
  (3) un BARRIDO de todo entero de 8 cifras 2026xxxx-2028xxxx que aparezca en
      el codigo (`modelo/src/**/*.py`, `src/**/*.py`) y en las salidas
      (`modelo/out/**/*.json`), excluyendo este propio modulo y las salidas de
      este paso.
"""
import os
import re
import sys

sys.path.insert(0, "/home/ubuntu/vcc/modelo/src/noche")
import comun as C  # noqa: E402

SEMILLAS = {
    # P1: 400 controles para la comprobacion del filtro contra compute_de, uno por contexto
    "p1_verif_filtro": [20281001 + i for i in range(6)],
    # P1: control de fragilidad del filtro: mitades aleatorias de celulas, una por contexto
    "p1_mitades": [20281021 + i for i in range(6)],
    # P2: nulo de permutacion del contexto de las estadisticas de control (del encargo)
    "p2_nulo": [20281100 + k for k in range(20)],
    # P2: igualacion de potencia en Jiang: submuestreo del control y adelgazamiento binomial
    "p2_ig_control": [20281201 + i for i in range(6)],
    "p2_ig_adelgaza": [20281211 + i for i in range(6)],
}

PROPIOS = {"semillas.py", "seleccion_control.json", "seleccion_control.md",
           "p1_filtrables.json", "p2_mecanismo.json"}
RAICES_PY = ["/home/ubuntu/vcc/modelo/src", "/home/ubuntu/vcc/src"]
RAICES_JSON = ["/home/ubuntu/vcc/modelo/out"]
PAT = re.compile(rb"(?<![0-9])(202[678][0-9]{4})(?![0-9])")


def _barrido():
    vistos = {}
    for raices, ext in ((RAICES_PY, ".py"), (RAICES_JSON, ".json")):
        for r in raices:
            for d, _s, fs in os.walk(r):
                if "/seleccion" in d:
                    continue
                for f in fs:
                    if not f.endswith(ext) or f in PROPIOS:
                        continue
                    p = os.path.join(d, f)
                    try:
                        b = open(p, "rb").read()
                    except Exception:
                        continue
                    for m in set(PAT.findall(b)):
                        vistos.setdefault(int(m), p)
    return vistos


def comprueba():
    usos = {}
    for fam, vals in SEMILLAS.items():
        for v in vals:
            usos.setdefault(int(v), []).append(fam)
    internos = {k: v for k, v in usos.items() if len(v) > 1}
    en_rango = sorted(k for k in usos if any(a <= k <= b for a, b in C.GASTADAS_RANGOS))
    expl = C._gastadas_explicitas() | {int(v) for vs in C.SEMILLAS_NUEVAS.values() for v in vs}
    en_expl = sorted(set(usos) & expl)
    barr = _barrido()
    en_barr = {k: barr[k] for k in sorted(set(usos) & set(barr))}
    ok = not internos and not en_rango and not en_expl and not en_barr
    return {
        "familia_nueva": "20281000 +",
        "semillas_enteras": SEMILLAS,
        "lista_plana_ordenada": sorted(usos),
        "n": len(usos),
        "choques_internos": internos,
        "choques_rangos_gastados": en_rango,
        "choques_enumeracion_comun": en_expl,
        "n_enteros_2026_2028_barridos": len(barr),
        "max_barrido_familia_2028": max([k for k in barr if k >= 20280000], default=None),
        "choques_barrido": en_barr,
        "sin_choques": ok,
    }


if __name__ == "__main__":
    import json
    r = comprueba()
    print(json.dumps({k: v for k, v in r.items() if k != "semillas_enteras"}, indent=1))
