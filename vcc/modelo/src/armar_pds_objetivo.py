#!/usr/bin/env python
"""Ensambla `modelo/out/pds_objetivo.md`: prosa del scratchpad + tablas
generadas por `pds_objetivo_tablas.py`. Solo lee y escribe el .md; no corre
ninguna medicion.

🔴 Los cuatro brazos diagnosticos no son candidatos y ninguna de sus cifras
entra en una submission ni en un rasgo. El encabezado del reporte lo declara.
"""
import io
import json
import os
import re
import sys
from contextlib import redirect_stdout

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pds_objetivo_tablas as T

SP = ("/tmp/claude-1000/-home-ubuntu-vcc/"
      "677bba4a-90db-471e-9b97-94590ae26aac/scratchpad")
DEST = "/home/ubuntu/vcc/modelo/out/pds_objetivo.md"
COMP = "/home/ubuntu/vcc/modelo/out/pds_composicion.json"

# referencias internas: la prosa se escribio con una numeracion provisional
REF = {"§4.3": "la TABLA B3", "§4.4": "la TABLA B3b"}
REF_DEF = {"§4.1": "la TABLA B0", "§4.2": "la TABLA A1"}


def prosa(nombre, extra=None):
    s = open(os.path.join(SP, nombre)).read()
    for k, v in (extra or {}).items():
        s = s.replace(k, v)
    return s.strip("\n")


def cap(fn, *a):
    buf = io.StringIO()
    with redirect_stdout(buf):
        r = fn(*a)
    return buf.getvalue().strip("\n"), r


def placeholders():
    """Rellena los PLACEHOLDER_* de la prosa de defectos con lo medido."""
    P = [json.loads(l) for l in open(T.OUT + "/pds_objetivo.json") if l.strip()]
    P = [r for r in P if not r.get("fallo_en")]
    # el brazo de copia no emite por razon y no tiene `diag_emision`
    g = [r["diag_emision"]["_global"] for r in P
         if isinstance(r.get("diag_emision"), dict)
         and r["diag_emision"].get("_global")]
    cero = sorted({x["n_control_mean_cero"] for x in g
                   if x.get("n_control_mean_cero") is not None})
    fr = [x["fraccion_delta_real_en_genes_sin_control"] for x in g
          if x.get("fraccion_delta_real_en_genes_sin_control") is not None]
    comp = json.load(open(COMP))
    fp = []
    for l, v in comp.items():
        tot = sum(f["n_ed"] for f in v["por_pert"])
        den = sum(f["n_ed_en_espacio"] for f in v["por_pert"])
        if tot:
            fp.append((tot - den) / tot)
    s_cero = (f"{cero[0]}" if len(cero) == 1 else
              f"{min(cero)} a {max(cero)}")
    return {
        "PLACEHOLDER_GENES_CERO": f"**{s_cero}**",
        "PLACEHOLDER_FRAC_CERO": (f"**{100*min(fr):.2f} % a "
                                  f"{100*max(fr):.2f} %**" if fr else "n/d"),
        "PLACEHOLDER_FRAC_PANEL": (f"**{100*min(fp):.2f} % a "
                                   f"{100*max(fp):.2f} %**" if fp else "n/d"),
    }


def main():
    d, ya = T.carga()
    a, fil = cap(T.parte_a)
    b0, _ = cap(T.composicion)
    bb, _ = cap(T.parte_b, d, ya, fil)
    c, _ = cap(T.parte_c, ya)
    cr, _ = cap(T.crudas, d)
    ve, _ = cap(T.verificaciones, d)

    ph = placeholders()
    piezas = [
        prosa("prosa_cab.md"),
        "---",
        prosa("prosa_constr.md", REF).lstrip("-").lstrip("\n"),
        "---",
        "## 4. PARTE A: la escala de `pds_cosine` contra el azar",
        a,
        "---",
        "## 5. PARTE B: de que esta hecho `pds_cosine`, y si se alcanza "
        "desde la cola",
        b0,
        bb,
        "---",
        "## 6. PARTE C: las dos cuentas ya medidas, recomprobadas",
        c,
        "---",
        "## 7. Los intermedios crudos, corrida por corrida",
        cr,
        ve,
        "---",
        prosa("prosa_defectos.md", {**REF, **REF_DEF, **ph})
        .lstrip("-").lstrip("\n"),
    ]
    txt = "\n\n".join(p for p in piezas if p) + "\n"
    txt = re.sub(r"\n{4,}", "\n\n\n", txt)
    with open(DEST, "w") as fh:
        fh.write(txt)
    sys.stderr.write(f"escrito {DEST}: {len(txt.splitlines())} lineas\n")
    for k, v in ph.items():
        sys.stderr.write(f"  {k} = {v}\n")
    sobra = [m for m in re.findall(r"PLACEHOLDER_\w+", txt)]
    sys.stderr.write(f"  placeholders sin rellenar: {sorted(set(sobra))}\n")


if __name__ == "__main__":
    main()
