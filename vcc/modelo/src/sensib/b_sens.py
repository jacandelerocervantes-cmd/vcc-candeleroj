#!/usr/bin/env python
"""BRAZO B — LA SENSIBILIDAD: que registra el marcador a n=400 sobre los
controles de la entrega, con verdad SINTETICA de rejilla DECLARADA.

🔴 ENCABEZADO: esto NO es biologia de D/E/F. Es la sensibilidad del marcador
sobre sus controles. La rejilla la elegi yo y por eso es un peso INVENTADO: el
brazo C es el que la sustituye por estructura medida.

Dos brazos de prediccion por celda, con el MISMO lado real y el MISMO bundle:
  suelo   : 400 celulas remuestreadas del pool de control, sin estructura
            (`predict_floor` por construccion)
  oraculo : la MISMA estructura del lado real, aplicada a celulas DISTINTAS del
            pool de control y con ruido Binomial/Poisson nuevo
El par (suelo, oraculo) de una celda es el RECORRIDO del marcador para esa
estructura: lo que puede registrar entre no saber nada y saberlo todo.

Uso: b_sens.py <CTX> cruz|todo|<k> [<k> ...]
"""
import sys

import banco as B
import comun as C


def main():
    ctx = sys.argv[1]
    arg = sys.argv[2]
    if arg == "cruz":
        ks = B.cruz()
    elif arg == "todo":
        ks = [r["k"] for r in B.REJILLA_B]
    else:
        ks = [int(x) for x in sys.argv[2:]]
    C.log(f"BRAZO B {ctx}: {len(ks)} celdas {ks}")
    cache = {}
    fallos = []
    for k in ks:
        r = B.REJILLA_B[k]
        try:
            B.corre_celda(ctx, f"B{k:02d}", k,
                          {kk: r[kk] for kk in ("amplitud", "frac_genes", "frac_comun")},
                          B.brazos_B(), cache=cache, fam="B")
        except SystemExit:
            raise
        except Exception as e:
            fallos.append((ctx, k, repr(e)[:300]))
            C.log(f"BRAZO B {ctx}/B{k:02d} FALLO: {e!r}")
    if fallos:
        C.log(f"BRAZO B {ctx}: fallaron {len(fallos)} celdas: {fallos}")


if __name__ == "__main__":
    main()
