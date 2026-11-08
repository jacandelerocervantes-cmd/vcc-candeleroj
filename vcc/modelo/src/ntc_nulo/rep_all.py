#!/usr/bin/env python
"""Ensambla `modelo/out/ntc_nulo.md` de las seis secciones."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rep_a, rep_b, rep_c, rep_d, rep_e, rep_f, rep_g

OUT = "/home/ubuntu/vcc/modelo/out/ntc_nulo.md"

L = []
L += rep_a.main()
L += rep_b.main()
L += rep_c.main()[0]
L += rep_d.main()
L += rep_e.main()
L += rep_f.main()
L += rep_g.main()

txt = "\n".join(L).rstrip() + "\n"
open(OUT, "w").write(txt)
print("escrito", OUT)
print("lineas:", txt.count("\n"), "bytes:", len(txt))
print("tablas:", txt.count("|---"))
