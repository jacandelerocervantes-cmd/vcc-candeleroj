#!/usr/bin/env python
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mcrep_a, mcrep_b, mcrep_c, mcrep_d, mcrep_e, mcrep_f
OUT = "/home/ubuntu/vcc/modelo/out/magnitud_centro.md"
L = []
L += mcrep_a.main()
L += mcrep_b.main()[0]
L += mcrep_c.main()
L += mcrep_d.main()
L += mcrep_e.main()
L += mcrep_f.main()
txt = "\n".join(L).rstrip() + "\n"
open(OUT, "w").write(txt)
print("escrito", OUT, "| lineas:", txt.count("\n"), "| bytes:", len(txt))
