#!/usr/bin/env python
"""Comprobacion EN CODIGO de que la familia de semillas NUEVA no choca con
ninguna de las ya usadas, y volcado de LAS SEMILLAS (no solo su conteo).

La familia nueva: 20278000 + 100*fold + j, fold = 0..5, j = 0..7  -> 48 semillas.

Las familias ya usadas se expanden con los MISMOS generadores que usan los
scripts del proyecto (base + fold, base + 100*fold + j, base + replica), y
ademas se recogen las semillas literalmente REGISTRADAS en los JSON de
`modelo/out/`, que es la verdad de lo que se corrio.
"""
import glob, json, os, re, sys

NUEVA_BASE = 20279000
FOLDS, JS = range(6), range(33)   # j = indice de la perturbacion en el panel ordenado
NUEVAS = [NUEVA_BASE + 100 * f + j for f in FOLDS for j in JS]
# las otras dos de la misma familia nueva, declaradas: bootstrap y permutacion
SEM_BOOT = []
SEM_PERM = []
TODAS_NUEVAS = NUEVAS + SEM_BOOT + SEM_PERM

# bases halladas en los scripts (grep de constantes SEED_*/SEM_* y literales)
BASES = [20260000, 20260001, 20261000, 20261001, 20261029, 20262000, 20263000,
         20264000, 20265000, 20266000, 20267000, 20267999, 20268000, 20269000,
         20270000, 20271000, 20272000, 20272500, 20273000, 20273900, 20274000,
         20274500, 20274999, 20275000, 20276000, 20277000]

# 🔴 La familia del PASO 2 (`magnitud_centro`) se expande con SU generador real,
# no con el genérico. Expandir `20278900` con `base + 100*f + j` lo derramaria
# hasta 20279509 y daria un choque FALSO contra la familia de este encargo: sus
# semillas reales son `base + k`, k = 0..19. Queda declarado porque es la segunda
# vez que el chequeo de semillas produce un falso positivo por el barrido.
EXACTAS_PASO2 = ([20278000 + 100 * f + j for f in range(6) for j in range(8)]
                 + [20278800] + [20278900 + k for k in range(20)])

usadas = set()
detalle = {}
for b in BASES:
    s = set()
    s.add(b)
    s.update(b + f for f in range(7))                      # base + fold
    s.update(b + 100 * f + j for f in range(7) for j in range(10))   # base+100f+j
    s.update(b + r for r in range(200))                    # base + replica (bootstrap)
    detalle[b] = (min(s), max(s), len(s))
    usadas |= s
usadas |= set(EXACTAS_PASO2)
detalle["paso2_exactas"] = (min(EXACTAS_PASO2), max(EXACTAS_PASO2),
                            len(EXACTAS_PASO2))

# y las semillas REGISTRADAS en los JSON ya escritos.
# 🔴 Se PARSEA el JSON y se recorren los valores, no se hace regex sobre el texto:
# un regex `20\d{6}` tambien casa dentro de los digitos de un float
# (p. ej. `-0.05420278503638758` contiene "20278503") y daria un choque falso.
# Eso paso en la primera version de este chequeo y queda declarado.
registradas = set()
origen = {}


def camina(o, fuente, clave=None):
    if isinstance(o, bool):
        return
    if isinstance(o, int):
        if 20_200_000 <= o <= 20_400_000:
            registradas.add(o)
            origen.setdefault(o, []).append(f"{fuente}:{clave}")
    elif isinstance(o, dict):
        for k, v in o.items():
            camina(v, fuente, str(k))
            # una semilla puede estar en la CLAVE (p. ej. "20274300": {...})
            if isinstance(k, str) and k.isdigit():
                camina(int(k), fuente, f"clave:{k}")
    elif isinstance(o, (list, tuple)):
        for v in o:
            camina(v, fuente, clave)


# 🔴 Se excluye el PROPIO directorio de salida de este encargo: en cuanto
# `mc_deltas.py` y `mc_calc.py` escriben, sus JSON ya contienen estas semillas y
# el chequeo se encontraria a si mismo. Lo que se compara es contra lo ANTERIOR.
MIO_DIR = "/home/ubuntu/vcc/modelo/out/posicion_mezcla"
for p_ in sorted(glob.glob("/home/ubuntu/vcc/modelo/out/*.json")
                 + glob.glob("/home/ubuntu/vcc/modelo/out/*/*.json")):
    if os.path.abspath(os.path.dirname(p_)) == os.path.abspath(MIO_DIR):
        continue
    try:
        d = json.load(open(p_))
    except Exception:
        continue
    camina(d, os.path.relpath(p_, "/home/ubuntu/vcc/modelo/out"))

choque_gen = sorted(set(TODAS_NUEVAS) & usadas)
choque_reg = sorted(set(TODAS_NUEVAS) & registradas)
choque_interno = sorted(set(NUEVAS) & (set(SEM_BOOT) | set(SEM_PERM)))

res = {
    "familia_nueva": "20279000 + 100*fold + j, fold=0..5, j=0..32 (j = indice de la perturbacion)",
    "n_semillas_nuevas": len(NUEVAS),
    "semillas_nuevas": NUEVAS,
    "semilla_bootstrap": SEM_BOOT,
    "semillas_permutacion_control": SEM_PERM,
    "n_total_nuevas": len(TODAS_NUEVAS),
    "choque_interno_entre_las_tres_familias_nuevas": choque_interno,
    "rango_nuevas": [min(NUEVAS), max(NUEVAS)],
    "bases_ya_usadas": BASES,
    "semillas_exactas_del_paso2": sorted(EXACTAS_PASO2),
    "nota_generador": ("la familia del PASO 2 se expande con SU generador real "
                       "(base+k para la de permutacion), no con el genérico "
                       "base+100*f+j, que la derramaria hasta 20279509 y daria "
                       "un choque falso"),
    "n_semillas_ya_usadas_expandidas": len(usadas),
    "rango_ya_usadas": [min(usadas), max(usadas)],
    "n_semillas_registradas_en_los_json": len(registradas),
    "nota_exclusion": ("se excluye modelo/out/magnitud_centro/ del barrido: son "
                       "las salidas de ESTE encargo y ya contienen estas semillas"),
    "rango_registradas": [min(registradas), max(registradas)] if registradas else None,
    "choque_contra_expansion_de_generadores": choque_gen,
    "choque_contra_las_registradas": choque_reg,
    "registradas_de_donde": {str(k): sorted(set(v))[:3] for k, v in
                             sorted(origen.items())},
    "sin_choque": (not choque_gen) and (not choque_reg) and (not choque_interno),
}
out = "/home/ubuntu/vcc/modelo/out/posicion_mezcla/semillas.json"
json.dump(res, open(out, "w"), indent=1, sort_keys=True)
print(f"particion: {len(NUEVAS)} semillas, {min(NUEVAS)}..{max(NUEVAS)}")
print(f"otras familias nuevas en este encargo: ninguna (solo la de empates)")
print(f"choque interno entre las tres familias nuevas: {choque_interno}")
print(f"ya usadas (expandidas): {len(usadas)}, {min(usadas)}..{max(usadas)}")
print(f"registradas en JSON: {len(registradas)}")
print(f"choque con generadores: {choque_gen}")
print(f"choque con registradas: {choque_reg}")
print(f"SIN CHOQUE: {res['sin_choque']}")
print("escrito", out)
if not res["sin_choque"]:
    sys.exit(2)
