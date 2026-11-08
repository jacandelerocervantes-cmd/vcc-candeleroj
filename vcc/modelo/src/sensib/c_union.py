#!/usr/bin/env python
"""BRAZO C — LA UNION: la estructura MEDIDA en el brazo A, inyectada en los
controles de la entrega, y el mapa calidad -> puntuacion.

🔴 TRANSFERENCIA DE ESTRUCTURA, DECLARADA. Lo que viaja de Jiang/2025 a D/E/F
son TRES NUMEROS de forma (amplitud, fraccion de genes movidos, fraccion
comun), tomados como RANGO (min, mediana, max de `a_resumen.json`), no valores
de gen ni de perturbacion. No viaja ningun delta real: los contextos no son los
mismos y no se presupone que lo sean.

C2 EL TECHO: el brazo `mag_f1` es el ORACULO PERFECTO (la misma estructura, en
   celulas distintas del pool de control, con ruido nuevo). Su `u` y su `s` son
   lo que saca un predictor perfecto a n = 400 y a la profundidad de la
   entrega. No tiene por que ser 1.0.
C3 EL MAPA: `f` recorre {0, 0.25, 0.5, 0.75, 1} en DOS variantes declaradas
   distintas:
     magnitud  = acierta la direccion y encoge la magnitud (f * especifico)
     identidad = acierta floor(f*n) de los genes y pone el resto en genes
                 EQUIVOCADOS, con la misma amplitud
   En las dos el eje COMUN se da por acertado, para que `f` module solo la
   componente especifica, que es lo que pide el encargo.

Uso: c_union.py <CTX> [i ...]
"""
import json
import os
import sys

import banco as B
import comun as C
import sup as U

RES = os.path.join(U.WORK, "a_resumen.json")


def main():
    ctx = sys.argv[1]
    pts = json.load(open(RES))["puntos_C"]
    idx = [int(x) for x in sys.argv[2:]] or list(range(len(pts)))
    C.log(f"BRAZO C {ctx}: {len(idx)} puntos de estructura medida")
    cache = {}
    for i in idx:
        p = pts[i]
        try:
            B.corre_celda(ctx, f"C{i}_{p['nombre']}", i,
                          {k: p[k] for k in ("amplitud", "frac_genes", "frac_comun")},
                          B.brazos_C(), cache=cache, fam="C",
                          extra={"punto_C": p, "origen": "brazo A, medido"})
        except SystemExit:
            raise
        except Exception as e:
            C.log(f"BRAZO C {ctx}/C{i} FALLO: {e!r}")


if __name__ == "__main__":
    main()
