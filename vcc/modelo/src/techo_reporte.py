#!/usr/bin/env python
"""Arma modelo/out/bucle_techo.md a partir de modelo/out/bucle_techo.json
(JSON-lines, un objeto por corrida: linea x variante x direccion) y resta
columna a columna contra modelo/out/bucle_suelo.json."""
import json
import os

TECHO_JSONL = "/home/ubuntu/vcc/modelo/out/bucle_techo.json"
SUELO_JSONL = "/home/ubuntu/vcc/modelo/out/bucle_suelo.json"
PANEL_PATH = "/home/ubuntu/vcc/modelo/out/suelo_panel.json"
OUT_MD = "/home/ubuntu/vcc/modelo/out/bucle_techo.md"

LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
METRICS = [
    "pds_cosine",
    "expr_mse_unbiased_capped_norm",
    "de_wilcoxon_lfc_nmae",
    "de_wilcoxon_direction_fidelity_yield_raw",
    "de_wilcoxon_direction_reach_raw",
    "de_wilcoxon_sig_jaccard",
]
CORTO = {
    "pds_cosine": "pds_cosine",
    "expr_mse_unbiased_capped_norm": "mse_capped_norm",
    "de_wilcoxon_lfc_nmae": "lfc_nmae",
    "de_wilcoxon_direction_fidelity_yield_raw": "fid",
    "de_wilcoxon_direction_reach_raw": "reach",
    "de_wilcoxon_sig_jaccard": "jac",
}
# Sentido de cada metrica: True = mas alto es mejor.
MAS_ES_MEJOR = {
    "pds_cosine": True,
    "expr_mse_unbiased_capped_norm": False,
    "de_wilcoxon_lfc_nmae": False,
    "de_wilcoxon_direction_fidelity_yield_raw": True,
    "de_wilcoxon_direction_reach_raw": True,
    "de_wilcoxon_sig_jaccard": True,
}


def cargar(path):
    out = []
    with open(path) as fh:
        for ln in fh:
            ln = ln.strip()
            if ln:
                out.append(json.loads(ln))
    return out


def fmt(v, nd=6):
    if v is None:
        return "NaN"
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    if f != f:
        return "NaN"
    return f"{f:.{nd}f}"


def de_por_mitad(rec, pert):
    """Mapea las cuentas de DE (lado real / lado pred de esa corrida) a mitad A y
    mitad B, segun la direccion. direccion='A' => la mitad A fue la PREDICCION."""
    c = rec["de_counts"][pert]
    if rec["direccion"] == "A":
        return c["n_DE_pred"], c["n_DE_real"], c
    return c["n_DE_real"], c["n_DE_pred"], c


def main():
    techo = [r for r in cargar(TECHO_JSONL) if "error" not in r]
    suelo = {r["linea_excluida"]: r for r in cargar(SUELO_JSONL) if "error" not in r}
    panel_data = json.load(open(PANEL_PATH))
    panel = panel_data["panel_comun"]

    idx = {(r["linea_excluida"], r["variante"], r["direccion"]): r for r in techo}
    L = []

    L.append("# Bucle del techo sobre 6 pliegues -- condicion IFNB\n")
    L.append(
        "El techo es el MISMO arnes del suelo con la funcion de emision "
        "reemplazada: en vez de emitir una prediccion degenerada (400 celulas "
        "del control remuestreadas con reemplazo), emite LA OTRA MITAD DE LAS "
        "CELULAS REALES de cada perturbacion. Las seis metricas salen de "
        "`cell-eval2 run --preset vcc2026 --profile vcc2026` igual que en el "
        "suelo, asi que la tabla resta columna a columna contra "
        "`bucle_suelo.md`. No hay ninguna formula nueva.\n"
    )
    L.append(
        f"**{len(techo)} corridas**: 6 pliegues x 2 variantes de control "
        "(S compartido / I partido) x 2 direcciones (A como prediccion / B como "
        "prediccion). Ninguna se promedia con otra.\n"
    )

    # ---------- Lo que no se pudo reusar ----------
    L.append("## Punto 1 -- lo que se reuso y lo que no, con el motivo\n")
    L.append("| pieza de suelo_arnes.py | reusada sin modificar | motivo si no |")
    L.append("|---|---|---|")
    L.append("| `load_global_obs`, `load_var_names`, `load_line_block`, `read_rows_csr` | si | - |")
    L.append("| `score_vs_self`, `read_agg_metrics`, `read_avg_score` | si | - |")
    L.append(
        "| `build_fold_adatas` | **no** | Su contrato de emision es "
        "`predict_fn(control_pool_X, labels, rng)`: la funcion de emision recibe "
        "SOLO el pool de control, nunca las celulas perturbadas, asi que por "
        "construccion no puede emitir la otra mitad de las celulas reales. "
        "Ademas fija el lado VERDAD (400 celulas sorteadas con `seed_real`) y el "
        "techo necesita gobernar LOS DOS lados para correr las dos direcciones. "
        "`techo_arnes.build_techo_pair` ensambla los AnnData con las MISMAS "
        "reglas (columna obs `target_gene`, mismo var index, "
        "`obs_names_make_unique`). |"
    )
    L.append(
        "| `run_cell_eval2` | **no** | Su argv es fijo y no admite banderas "
        "extra. El techo necesita dos: `--write-degenes`, para contar n_DE y "
        "|union| con la definicion del propio marcador (punto 5), y "
        "`--set control_source=pred` para la variante I (punto 3). "
        "`techo_arnes.run_cell_eval2_techo` reproduce ese argv exacto "
        "(`argv_suelo`) y agrega solo esas banderas. |"
    )
    L.append(
        "\nLas tres piezas del suelo **no fueron modificadas**. No hay `git` en "
        "este arbol, asi que la evidencia es mtime + md5, leidos al escribir "
        "este reporte: siguen siendo los de la corrida del suelo (04:16-04:19 "
        "UTC), anteriores a la primera linea de `techo_arnes.py`. Las piezas "
        "nuevas viven solo en `modelo/src/techo_arnes.py`, "
        "`modelo/src/bucle_techo.py` y `modelo/src/techo_reporte.py`.\n"
    )
    L.append("| archivo | mtime | md5 |")
    L.append("|---|---|---|")
    import hashlib
    import datetime
    for nombre in ("suelo_arnes.py", "suelo_floor.py", "bucle_suelo.py",
                   "techo_arnes.py", "bucle_techo.py"):
        p = os.path.join("/home/ubuntu/vcc/modelo/src", nombre)
        if not os.path.exists(p):
            continue
        with open(p, "rb") as fh:
            h = hashlib.md5(fh.read()).hexdigest()
        ts = datetime.datetime.fromtimestamp(
            os.path.getmtime(p), datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
        L.append(f"| `{nombre}` | {ts} | `{h}` |")

    # ---------- Punto 3 / control_source ----------
    L.append("## Punto 3 -- las dos variantes se deciden por `control_source`, no por el h5ad\n")
    L.append(
        "El preset `vcc2026` trae `control_source: real` "
        "(`cell-eval2/src/cell_eval2/configs/vcc2026.yaml`). Bajo ese valor el DE "
        "del lado PREDICHO toma su referencia del bloque de control del lado "
        "REAL -- de modo que poner un control distinto en el h5ad de prediccion "
        "no lo independiza. Por eso:\n"
    )
    L.append("| variante | control en el h5ad | `control_source` | que comparte | techo |")
    L.append("|---|---|---|---|---|")
    L.append(
        "| S (compartido) | VERDAD = pool de control completo de la linea; "
        "PRED = 400 celulas de ese mismo pool con reemplazo (la regla exacta del "
        "suelo) | `real` (preset tal cual) | las dos mitades hacen su DE contra "
        "el MISMO control: errores correlacionados | **optimista** |"
    )
    L.append(
        "| I (partido) | el control tambien se parte en mitades; cada mitad usa "
        "la suya | `pred` (`--set control_source=pred`) | nada: cada mitad trae "
        "su propio control | **pesimista** |"
    )
    L.append(
        "\nConsecuencia que conviene tener presente: **el suelo ya estaba medido "
        "con control compartido** (`control_source: real`), asi que la variante S "
        "es la unica que comparte convencion de puntuacion con `bucle_suelo.md`. "
        "La variante I cambia esa convencion, de modo que su resta contra el "
        "suelo mezcla dos efectos (mitad de poblacion + control independiente) y "
        "no solo uno.\n"
    )
    L.append(
        "Esto no es una hipotesis de esta sesion: `cell_eval2/ceiling.py` "
        "implementa un techo por self-split disjunto y **fuerza** "
        "`control_source=\"pred\"`, documentando que bajo `\"real\"` ambas "
        "mitades calculan su log2FC contra el control de una sola de ellas, lo "
        "que comparte el ruido de muestreo entre las dos cantidades cuya "
        "concordancia se esta midiendo y **sesga el techo hacia arriba**. Trae la "
        "medicion: `lfc_spearman` 0.54 -> 0.74, y `nsig_spearman` paso de NaN "
        "(sin techo defendible) a un 0.94 de apariencia confiable. En sus propias "
        "palabras, el override \"no es opcional\".\n"
    )
    L.append(
        "Ese `ceiling.py` NO se uso como sustituto, por dos razones: aplica la "
        "correccion Spearman-Brown `2r/(1+r)` (que el encargo pide no hacer), y "
        "su lista verificada `SB_METRICS` contiene solo UNA de las seis metricas "
        "vcc2026 (`pds_cosine`): las otras cinco saldrian `NaN`. El arnes del "
        "suelo con la emision reemplazada si entrega las seis.\n"
    )

    # ---------- Punto 4 / contrato ----------
    L.append("## Punto 4 -- el marcador NO exige un numero fijo de celulas por perturbacion\n")
    L.append(
        "`cell_eval2.io.validate_pair` exige exactamente tres cosas: eje de "
        "genes identico en nombre y orden, conjuntos de etiquetas de "
        "perturbacion iguales entre PRED y REAL, y que la etiqueta de control "
        "este presente en REAL. **Nada sobre el numero de celulas por etiqueta.**\n"
    )
    L.append("Lo que hizo el arnes del suelo en ese punto, para copiarlo:\n")
    L.append("| lado | regla del suelo | profundidad |")
    L.append("|---|---|---|")
    L.append("| VERDAD, perturbaciones | 400 celulas sin reemplazo, con `assert pool >= 400` | 400 |")
    L.append("| VERDAD, control | TODO el pool de control, sin submuestrear | 1379-3228 |")
    L.append("| PREDICCION, toda etiqueta | 400 celulas con reemplazo | 400 |")
    L.append(
        "\nEs decir: el suelo **ya enviaba profundidades desiguales** entre PRED "
        "y REAL (su control real tenia 1379-3228 celulas contra 400 en pred; los "
        "totales fueron 14579 vs 13600 filas en A549). Las 400 eran convencion "
        "del arnes y cifra del reto, no requisito del marcador.\n"
    )
    L.append(
        "**Decision, declarada y no silenciosa:** las mitades se emiten a su "
        "tamano natural, sin remuestrear ni rellenar a 400. Emitir la mitad tal "
        "cual es lo que mide la reproducibilidad de la propia verdad; rellenar a "
        "400 con reemplazo habria inventado celulas en el lado VERDAD y habria "
        "sido, de hecho, una correccion por la perdida de potencia -- justo lo "
        "que el encargo prohibe.\n"
    )
    L.append(
        "**Consecuencia, declarada y NO corregida:** el techo se mide a ~la mitad "
        "de la profundidad por perturbacion del suelo. La resta columna a columna "
        "contra `bucle_suelo.md` mezcla el efecto techo con ese cambio de "
        "profundidad. Queda como **cota inferior**.\n"
    )

    # ---------- Punto 2: n de celulas ----------
    L.append("## Punto 2 -- n de celulas por perturbacion y sus mitades (no son 400)\n")
    L.append(
        "400 es la cifra del reto, no de Jiang. Los n reales de IFNB y las "
        "mitades efectivas que de ahi salen, por linea. `extra` marca a que "
        "mitad fue la celula sobrante de un n impar (desempate alternado por "
        "ordinal de grupo, no dirigido).\n"
    )
    ref = {}
    for line in LINES:
        r = idx.get((line, "S", "A"))
        if r:
            ref[line] = r
    if ref:
        L.append(
            "| perturbacion | " + " | ".join(f"{l} n (A/B)" for l in LINES if l in ref) + " |"
        )
        L.append("|---" * (1 + len([l for l in LINES if l in ref])) + "|")
        for pert in panel:
            celdas = []
            for line in LINES:
                if line not in ref:
                    continue
                r = ref[line]
                n = r["n_celulas_por_pert"][pert]
                a = r["n_mitadA_por_pert"][pert]
                b = r["n_mitadB_por_pert"][pert]
                celdas.append(f"{n} ({a}/{b})")
            L.append(f"| {pert} | " + " | ".join(celdas) + " |")
        celdas = []
        for line in LINES:
            if line not in ref:
                continue
            r = ref[line]
            n = r["n_celulas_por_pert"]["NT"]
            a = r["n_mitadA_por_pert"]["NT"]
            b = r["n_mitadB_por_pert"]["NT"]
            celdas.append(f"{n} ({a}/{b})")
        L.append("| **NT (control)** | " + " | ".join(celdas) + " |")

        todos_n = [r["n_celulas_por_pert"][p] for r in ref.values() for p in panel]
        L.append(
            f"\nRango de n por (linea, perturbacion): **{min(todos_n)} a "
            f"{max(todos_n)}** celulas, mediana "
            f"{sorted(todos_n)[len(todos_n)//2]}. Ninguna es 400. La mitad mas "
            f"delgada del estudio tiene **{min(todos_n)//2} celulas**, frente a "
            "las 400 por perturbacion del suelo.\n"
        )

    # ---------- Las tablas, una por variante y direccion ----------
    L.append("## Las tablas -- una por variante y por direccion (el entregable)\n")
    L.append(
        "Mismas seis columnas y mismo orden de lineas que `bucle_suelo.md`. "
        "`avg_score` se puntua contra si mismo igual que en el suelo, asi que "
        "vale 0.0 por construccion y no es informativo (ver nota al final).\n"
    )
    for variante in ("S", "I"):
        for direccion in ("A", "B"):
            etiqueta = {
                ("S", "A"): "Variante S (control COMPARTIDO) -- direccion: mitad A como PREDICCION",
                ("S", "B"): "Variante S (control COMPARTIDO) -- direccion: mitad B como PREDICCION",
                ("I", "A"): "Variante I (control PARTIDO) -- direccion: mitad A como PREDICCION",
                ("I", "B"): "Variante I (control PARTIDO) -- direccion: mitad B como PREDICCION",
            }[(variante, direccion)]
            L.append(f"\n### {etiqueta}\n")
            L.append(
                "| linea excluida | n control (real/pred) | "
                + " | ".join(CORTO[m] for m in METRICS) + " | avg_score |"
            )
            L.append("|---" * (3 + len(METRICS)) + "|")
            for line in LINES:
                r = idx.get((line, variante, direccion))
                if r is None:
                    L.append(f"| {line} | (pendiente) | " + " | ".join("-" for _ in METRICS) + " | - |")
                    continue
                L.append(
                    f"| {line} | {r['n_real_ctrl']}/{r['n_pred_ctrl']} | "
                    + " | ".join(fmt(r.get(m)) for m in METRICS)
                    + f" | {fmt(r.get('avg_score'))} |"
                )

    # ---------- Resta contra el suelo ----------
    L.append("\n## La resta columna a columna contra el suelo\n")
    L.append(
        "`techo - suelo` por metrica y por pliegue. El signo util depende del "
        "sentido de cada metrica: en `mse_capped_norm` y `lfc_nmae` mas bajo es "
        "mejor, en las otras cuatro mas alto es mejor. La columna `mejor` dice si "
        "el techo quedo del lado bueno del suelo.\n"
    )
    for variante in ("S", "I"):
        for direccion in ("A", "B"):
            L.append(f"\n### Variante {variante}, direccion {direccion}: techo - suelo\n")
            L.append("| linea | " + " | ".join(CORTO[m] for m in METRICS) + " | mejor en |")
            L.append("|---" * (2 + len(METRICS)) + "|")
            for line in LINES:
                r = idx.get((line, variante, direccion))
                s = suelo.get(line)
                if r is None or s is None:
                    continue
                celdas = []
                n_mejor = 0
                for m in METRICS:
                    tv, sv = r.get(m), s.get(m)
                    if tv is None or sv is None:
                        celdas.append("NaN")
                        continue
                    d = tv - sv
                    bueno = (d > 0) if MAS_ES_MEJOR[m] else (d < 0)
                    if bueno:
                        n_mejor += 1
                    celdas.append(f"{d:+.6f}")
                L.append(f"| {line} | " + " | ".join(celdas) + f" | {n_mejor}/6 |")

    # ---------- El denominador (r - b) ----------
    L.append("\n## El denominador `(r - b)`: donde la escala funciona y donde NO\n")
    L.append(
        "La escala `(u - b) / (r - b)` solo tiene sentido si el techo `r` esta "
        "del lado BUENO del suelo `b`. Aqui el `span` esta orientado por el "
        "sentido de cada metrica: `span = r - b` donde mas alto es mejor, "
        "`span = b - r` donde mas bajo es mejor. **Un `span` negativo significa "
        "que el techo quedo peor que el suelo y que esa metrica no admite esta "
        "escala.**\n"
    )
    L.append("| metrica | var | dir | spans por linea (A549,BXPC3,HAP1,HT29,K562,MCF7) | n negativos |")
    L.append("|---|---|---|---|---|")
    roto = {}
    for m in METRICS:
        for variante in ("S", "I"):
            for direccion in ("A", "B"):
                spans = []
                for line in LINES:
                    r = idx.get((line, variante, direccion))
                    s = suelo.get(line)
                    if r is None or s is None:
                        continue
                    tv, sv = r.get(m), s.get(m)
                    if tv is None or sv is None:
                        spans.append(None)
                        continue
                    spans.append((tv - sv) if MAS_ES_MEJOR[m] else (sv - tv))
                if not spans:
                    continue
                neg = sum(1 for v in spans if v is not None and v <= 0)
                roto.setdefault(m, []).append(neg)
                txt = ", ".join("NaN" if v is None else f"{v:+.4f}" for v in spans)
                L.append(f"| {CORTO[m]} | {variante} | {direccion} | {txt} | {neg}/{len(spans)} |")
    rotas = [m for m, negs in roto.items() if sum(negs) > 0]
    sanas = [m for m, negs in roto.items() if sum(negs) == 0]
    L.append("\n**Lectura.**\n")
    if sanas:
        L.append(
            "- Metricas donde el techo quedo del lado bueno del suelo en TODAS "
            "las lineas, variantes y direcciones (la escala funciona): "
            + ", ".join(f"`{CORTO[m]}`" for m in sanas) + ".\n"
        )
    if rotas:
        L.append(
            "- Metricas con al menos un `span` negativo o nulo (la escala NO "
            "funciona ahi): " + ", ".join(f"`{CORTO[m]}`" for m in rotas) + ".\n"
        )
        L.append(
            "- El caso de `jac` no es un fallo de medicion sino el artefacto "
            "D-023 operando en el SUELO: la prediccion degenerada del suelo no "
            "llama nada significativo, de modo que toda perturbacion cuya mitad "
            "real tampoco llama nada queda con union vacia y cobra un 1.0 "
            "regalado. Ese 1.0 inflaba el suelo. El techo, con senal real en los "
            "dos lados, tiene uniones NO vacias y jaccards genuinos menores que "
            "1 -- y sale por debajo. **El suelo de `jac` no es un suelo: es un "
            "artefacto de union vacia.** Mientras eso no se arregle, `jac` no "
            "admite la escala `(u - b)/(r - b)`, y la tabla de uniones vacias de "
            "mas abajo es el lugar donde se ve.\n"
        )

    # ---------- El intervalo S..I ----------
    L.append("\n## El intervalo del techo: entre S y I (no se elige una)\n")
    L.append(
        "El techo verdadero esta ENTRE la variante optimista (S, control "
        "compartido) y la pesimista (I, control partido). Aqui las dos, con el "
        "ancho del intervalo. Se promedia sobre las dos DIRECCIONES solo para "
        "dar el ancho; las direcciones crudas estan arriba, sin promediar.\n"
    )
    L.append("| linea | metrica | S dir A | S dir B | I dir A | I dir B | intervalo S..I | ancho |")
    L.append("|---|---|---|---|---|---|---|---|")
    for line in LINES:
        for m in METRICS:
            sa = idx.get((line, "S", "A"), {}).get(m)
            sb = idx.get((line, "S", "B"), {}).get(m)
            ia = idx.get((line, "I", "A"), {}).get(m)
            ib = idx.get((line, "I", "B"), {}).get(m)
            if None in (sa, sb, ia, ib):
                continue
            s_med = (sa + sb) / 2
            i_med = (ia + ib) / 2
            lo, hi = min(s_med, i_med), max(s_med, i_med)
            L.append(
                f"| {line} | {CORTO[m]} | {fmt(sa)} | {fmt(sb)} | {fmt(ia)} | "
                f"{fmt(ib)} | {fmt(lo)} .. {fmt(hi)} | {fmt(hi-lo)} |"
            )

    # ---------- El suelo de jac, con sus intermedios ----------
    L.append("\n### El suelo de `jac` medido por dentro: es la tasa de union vacia\n")
    L.append(
        "La afirmacion de arriba no se infiere, se cuenta. En `results.csv` de "
        "cada pliegue del suelo, el jaccard por perturbacion toma **solo dos "
        "valores: exactamente 0.0 o exactamente 1.0**, nunca uno intermedio (se "
        "verifico en los seis pliegues, 33 de 33 perturbaciones cada uno). Con "
        "una prediccion degenerada el conjunto predicho esta vacio, asi que un "
        "1.0 solo puede venir de `|R u P| = 0` y un 0.0 significa que el lado "
        "real si tenia senal. De modo que la columna `jac` del suelo **es la "
        "fraccion de perturbaciones sin ninguna senal de DE**, y no mide nada de "
        "la prediccion.\n"
    )
    L.append(
        "| linea | suelo: perts jac=1.0 (union vacia) | suelo: perts jac=0.0 | "
        "suelo: otros valores | suelo `jac` = 1.0s/33 | techo S: perts union "
        "vacia | techo S `jac` | span (techo-suelo) |"
    )
    L.append("|---" * 8 + "|")
    for line in LINES:
        s = suelo.get(line)
        r = idx.get((line, "S", "A"))
        if s is None or r is None:
            continue
        rc = os.path.join("/datos/modelo/out", f"run_{line}", "results.csv")
        n1 = n0 = notros = None
        if os.path.exists(rc):
            n1 = n0 = notros = 0
            with open(rc) as fh:
                next(fh, None)
                for ln in fh:
                    parts = ln.strip().split(",")
                    if len(parts) < 3 or parts[1] != "de_wilcoxon_sig_jaccard":
                        continue
                    v = float(parts[2])
                    if v == 1.0:
                        n1 += 1
                    elif v == 0.0:
                        n0 += 1
                    else:
                        notros += 1
        vacias_t = sum(
            1 for p in panel if r["de_counts"].get(p, {}).get("union_vacia")
        )
        sj, tj = s.get("de_wilcoxon_sig_jaccard"), r.get("de_wilcoxon_sig_jaccard")
        frac = "-" if n1 is None else f"{n1}/33 = {n1/33:.6f}"
        L.append(
            f"| {line} | {'-' if n1 is None else n1} | {'-' if n0 is None else n0} | "
            f"{'-' if notros is None else notros} | {frac} | {vacias_t} | "
            f"{fmt(tj)} | {tj - sj:+.6f} |"
        )
    L.append(
        "\nLa quinta columna reproduce exactamente la columna `jac` de "
        "`bucle_suelo.md` como un simple conteo de uniones vacias dividido por "
        "33. Donde el suelo tiene 0 uniones vacias (BXPC3 y MCF7) el techo queda "
        "limpiamente por encima; donde el suelo cobraba mas 1.0 regalados que el "
        "techo (A549 15 vs 10, HAP1 14 vs 4, HT29 6 vs 1) el span sale negativo. "
        "**Eso es todo lo que significa el span negativo de `jac`: no que la "
        "verdad sea menos reproducible que una prediccion degenerada, sino que el "
        "suelo de esa metrica estaba midiendo la tasa de perturbaciones sin "
        "senal.**\n"
    )

    L.append(
        "\n### ¿Es S realmente la optimista? (se comprueba, no se asume)\n"
    )
    L.append(
        "Por metrica, en cuantas de las 6 lineas la variante S quedo del lado "
        "BUENO de la variante I (promediando las dos direcciones solo para esta "
        "comparacion). El encargo esperaba S optimista en las seis.\n"
    )
    L.append("| metrica | S mejor que I en | S peor que I en | veredicto |")
    L.append("|---|---|---|---|")
    for m in METRICS:
        mejor = peor = 0
        for line in LINES:
            sa = idx.get((line, "S", "A"), {}).get(m)
            sb = idx.get((line, "S", "B"), {}).get(m)
            ia = idx.get((line, "I", "A"), {}).get(m)
            ib = idx.get((line, "I", "B"), {}).get(m)
            if None in (sa, sb, ia, ib):
                continue
            s_med, i_med = (sa + sb) / 2, (ia + ib) / 2
            s_gana = (s_med > i_med) if MAS_ES_MEJOR[m] else (s_med < i_med)
            if s_gana:
                mejor += 1
            else:
                peor += 1
        if mejor and not peor:
            ver = "S optimista, como se esperaba"
        elif peor and not mejor:
            ver = "**S NO es la optimista aqui**"
        else:
            ver = "mezclado, no domina ninguna"
        L.append(f"| {CORTO[m]} | {mejor}/6 | {peor}/6 | {ver} |")
    L.append(
        "\nDonde S no sale optimista, el sospechoso es otra vez la union vacia: "
        "la variante I parte el control, de modo que cada mitad hace su DE con "
        "la mitad de las celulas de referencia, detecta MENOS genes "
        "significativos, y produce MAS uniones vacias -- cada una un 1.0 "
        "regalado en `jac`. El efecto es un techo mas alto por perdida de "
        "potencia, no por mejor concordancia. La tabla de uniones vacias por "
        "variante, mas abajo, es la que lo distingue.\n"
    )

    # ---------- Punto 5: por perturbacion ----------
    L.append("\n## Punto 5 -- por perturbacion: n_DE mitad A, n_DE mitad B, |union|\n")
    L.append(
        "`de_wilcoxon_sig_jaccard` devuelve **1.0 cuando la union es vacia** "
        "(convencion `J(vacio, vacio) = 1`, D-023). Sin la columna `|union|` a la "
        "vista, una perturbacion sin ninguna senal aparece con techo perfecto y "
        "es invisible. Las cuentas se calculan con la definicion del propio "
        "marcador: `p_adj < 0.05`, `.unique()` sobre (target, feature) y el "
        "anti-join que excluye el gen propio de cada target en LOS DOS lados "
        "(issue #172), llamando a `cell_eval2.de.prepare_de` sobre las tablas que "
        "escribio `--write-degenes`.\n"
    )
    L.append(
        "`n_DE_mitadA` y `n_DE_mitadB` son invariantes a la direccion dentro de "
        "una variante (cada mitad hace su DE contra el mismo control en las dos "
        "direcciones), asi que se tabulan una vez por variante.\n"
    )

    # Resumen de uniones vacias
    L.append("### Resumen: perturbaciones con union vacia (jaccard=1.0 regalado)\n")
    L.append(
        "La ultima columna es un **diagnostico, no una metrica**: es la media de "
        "los jaccards por perturbacion dejando fuera las de union vacia, "
        "calculada aqui para medir cuanto del numero reportado es artefacto. NO "
        "es comparable contra las columnas de `bucle_suelo.md` y no sustituye a "
        "lo que agrega el marcador.\n"
    )
    L.append(
        "| variante | linea | n perts con \\|union\\|=0 | de 33 | jac reportado "
        "(marcador) | jac sin esas perts (diagnostico) |"
    )
    L.append("|---|---|---|---|---|---|")
    for variante in ("S", "I"):
        for line in LINES:
            r = idx.get((line, variante, "A"))
            if r is None:
                continue
            vacias = [p for p in panel if r["de_counts"].get(p, {}).get("union_vacia")]
            jacs = [r["de_counts"][p]["jaccard_recomputado"] for p in panel if p in r["de_counts"]]
            no_vacias = [
                r["de_counts"][p]["jaccard_recomputado"]
                for p in panel
                if p in r["de_counts"] and not r["de_counts"][p]["union_vacia"]
            ]
            jac_rep = sum(jacs) / len(jacs) if jacs else float("nan")
            jac_sin = sum(no_vacias) / len(no_vacias) if no_vacias else float("nan")
            L.append(
                f"| {variante} | {line} | {len(vacias)} | 33 | {fmt(jac_rep)} | {fmt(jac_sin)} |"
            )

    for variante in ("S", "I"):
        L.append(f"\n### Variante {variante} -- detalle por (linea, perturbacion)\n")
        L.append(
            "| linea | perturbacion | n celulas | n mitad A | n mitad B | "
            "n_DE_mitadA | n_DE_mitadB | interseccion | **\\|union\\|** | jac marcador |"
        )
        L.append("|---" * 10 + "|")
        for line in LINES:
            r = idx.get((line, variante, "A"))
            if r is None:
                continue
            for pert in panel:
                if pert not in r["de_counts"]:
                    continue
                n_a, n_b, c = de_por_mitad(r, pert)
                jm = r["per_pert_metrics"].get(pert, {}).get("de_wilcoxon_sig_jaccard")
                marca = " **(union vacia)**" if c["union_vacia"] else ""
                L.append(
                    f"| {line} | {pert} | {r['n_celulas_por_pert'][pert]} | "
                    f"{r['n_mitadA_por_pert'][pert]} | {r['n_mitadB_por_pert'][pert]} | "
                    f"{n_a} | {n_b} | {c['n_DE_interseccion']} | {c['n_DE_union']} | "
                    f"{fmt(jm, 4)}{marca} |"
                )

    # ---------- 20 filas trabajadas ----------
    L.append("\n## 20 filas trabajadas, con los intermedios completos\n")
    L.append(
        "Veinte filas (linea, perturbacion) elegidas para cubrir los dos "
        "extremos de n, las dos variantes y los casos de union vacia. Cada fila "
        "trae TODO lo que entra en sus numeros, de la celula cruda al jaccard: "
        "no hay ningun paso intermedio oculto.\n"
    )
    filas = []
    for variante in ("S", "I"):
        for line in LINES:
            r = idx.get((line, variante, "A"))
            if r is None:
                continue
            for pert in panel:
                if pert not in r["de_counts"]:
                    continue
                n_a, n_b, c = de_por_mitad(r, pert)
                filas.append({
                    "variante": variante, "linea": line, "pert": pert,
                    "n": r["n_celulas_por_pert"][pert],
                    "mA": r["n_mitadA_por_pert"][pert],
                    "mB": r["n_mitadB_por_pert"][pert],
                    "nA": n_a, "nB": n_b,
                    "inter": c["n_DE_interseccion"], "union": c["n_DE_union"],
                    "jac_rec": c["jaccard_recomputado"],
                    "jac_mar": r["per_pert_metrics"].get(pert, {}).get("de_wilcoxon_sig_jaccard"),
                    "pds": r["per_pert_metrics"].get(pert, {}).get("pds_cosine"),
                    "reach": r["per_pert_metrics"].get(pert, {}).get(
                        "de_wilcoxon_direction_reach_raw"),
                    "fid": r["per_pert_metrics"].get(pert, {}).get(
                        "de_wilcoxon_direction_fidelity_yield_raw"),
                    "vacia": c["union_vacia"],
                })
    if filas:
        elegidas = []
        por_n = sorted(filas, key=lambda f: f["n"])
        elegidas += por_n[:4]                      # los n mas delgados
        elegidas += por_n[-4:]                     # los n mas gruesos
        vacias = [f for f in filas if f["vacia"]]
        elegidas += vacias[:4]                     # union vacia
        grandes = sorted(filas, key=lambda f: -f["union"])
        elegidas += grandes[:4]                    # union mas grande
        vistos = set()
        unicas = []
        for f in elegidas:
            k = (f["variante"], f["linea"], f["pert"])
            if k not in vistos:
                vistos.add(k)
                unicas.append(f)
        for f in filas:
            if len(unicas) >= 20:
                break
            k = (f["variante"], f["linea"], f["pert"])
            if k not in vistos:
                vistos.add(k)
                unicas.append(f)
        unicas = unicas[:20]

        L.append(
            "| # | var | linea | pert | n celulas | mitad A | mitad B | n_DE_A | "
            "n_DE_B | inter | \\|union\\| | jac recomputado | jac marcador | dif | "
            "pds_cosine | reach | fid |"
        )
        L.append("|---" * 17 + "|")
        for i, f in enumerate(unicas, 1):
            jm = f["jac_mar"]
            dif = "-" if jm is None else f"{abs(jm - f['jac_rec']):.2e}"
            L.append(
                f"| {i} | {f['variante']} | {f['linea']} | {f['pert']} | {f['n']} | "
                f"{f['mA']} | {f['mB']} | {f['nA']} | {f['nB']} | {f['inter']} | "
                f"{f['union']} | {fmt(f['jac_rec'], 6)} | {fmt(jm, 6)} | {dif} | "
                f"{fmt(f['pds'], 6)} | {fmt(f['reach'], 6)} | {fmt(f['fid'], 6)} |"
            )
        L.append(
            "\nLa columna `dif` es |jac recomputado - jac marcador|: cierra el "
            "circuito entre las cuentas de DE de esta tabla y el numero que el "
            "marcador puntuo. Si fuera distinto de ~0, las columnas n_DE/|union| "
            "no describirian lo que se puntuo.\n"
        )

    # ---------- Verificacion: los per-pert componen el agregado ----------
    L.append("## Verificacion: las filas por perturbacion componen el agregado\n")
    L.append(
        "`agg_results.csv` reporta la fila `mean` sobre las 33 perturbaciones. "
        "Aqui se recalcula esa media desde `results.csv` (por perturbacion) y se "
        "compara con la que se reporto en las tablas de arriba, para cada "
        "corrida. `jac cruce` es el maximo |jac recomputado - jac marcador| sobre "
        "las 33 perturbaciones de esa corrida.\n"
    )
    L.append(
        "| linea | var | dir | jac agregado | media de las 33 | dif | "
        "pds agregado | media de las 33 | dif | jac cruce max |"
    )
    L.append("|---" * 10 + "|")
    for variante in ("S", "I"):
        for direccion in ("A", "B"):
            for line in LINES:
                r = idx.get((line, variante, direccion))
                if r is None:
                    continue
                pj = [r["per_pert_metrics"][p]["de_wilcoxon_sig_jaccard"]
                      for p in panel if p in r["per_pert_metrics"]
                      and "de_wilcoxon_sig_jaccard" in r["per_pert_metrics"][p]]
                pp = [r["per_pert_metrics"][p]["pds_cosine"]
                      for p in panel if p in r["per_pert_metrics"]
                      and "pds_cosine" in r["per_pert_metrics"][p]]
                mj = sum(pj) / len(pj) if pj else float("nan")
                mp = sum(pp) / len(pp) if pp else float("nan")
                aj = r.get("de_wilcoxon_sig_jaccard")
                ap = r.get("pds_cosine")
                L.append(
                    f"| {line} | {variante} | {direccion} | {fmt(aj)} | {fmt(mj)} | "
                    f"{abs(aj - mj):.2e} | {fmt(ap)} | {fmt(mp)} | "
                    f"{abs(ap - mp):.2e} | {r['jaccard_cruce_max_dif']:.2e} |"
                )

    # ---------- Intermedios de corrida ----------
    L.append("\n## Intermedios de corrida (tiempo, disco, formas, semillas)\n")
    L.append(
        "| linea | var | dir | control_source | real_shape | pred_shape | "
        "run_s | total_s | free_gib | seed_split | seed_ctrl_boot | score |"
    )
    L.append("|---" * 12 + "|")
    for variante in ("S", "I"):
        for direccion in ("A", "B"):
            for line in LINES:
                r = idx.get((line, variante, direccion))
                if r is None:
                    continue
                est = "nego (ver nota)" if r.get("score_error") else "ok"
                L.append(
                    f"| {line} | {variante} | {direccion} | {r['control_source']} | "
                    f"{r['real_shape']} | {r['pred_shape']} | {r['run_duration_s']} | "
                    f"{r['total_duration_s']} | {r['free_gib_before']} | "
                    f"{r['seed_split']} | {r['seed_ctrl_boot']} | {est} |"
                )

    # ---------- Semillas ----------
    L.append("\n## Semillas\n")
    L.append("| uso | familia | valor |")
    L.append("|---|---|---|")
    L.append("| VERDAD del suelo | 2026000x | 20260001..20260006 (no reusada) |")
    L.append("| PREDICCION del suelo | 2026100x | 20261001..20261006 (no reusada) |")
    L.append("| particion A/B del techo | 2026200x | 20262000..20262005 |")
    L.append("| bootstrap del control de PRED, variante S | 2026300x | 20263000..20263005 |")
    L.append("\nNunca 42.\n")

    # ---------- Lo que no se hizo ----------
    L.append("## Lo que NO se hizo, a proposito\n")
    L.append(
        "- **No se promedio ninguna metrica a mano.** Las seis columnas son la "
        "fila `mean` de `agg_results.csv` que agrega el marcador, igual que en el "
        "suelo. La seccion de verificacion de arriba recalcula esa media solo "
        "para comprobar que coincide, no para sustituirla.\n"
        "- **No se promediaron las dos direcciones** en las tablas del "
        "entregable: se reportan las dos. El unico promedio de direcciones esta "
        "en la tabla del intervalo S..I, y solo para dar su ancho.\n"
        "- **No se eligio entre S e I.** El techo es el intervalo.\n"
        "- **No se corrigio el techo por la perdida de potencia** de usar media "
        "poblacion (ni con Spearman-Brown ni de otro modo). Queda declarado como "
        "cota inferior.\n"
        "- **No se explica por que `reach` da no-nulo con conjunto predicho "
        "vacio.** Sigue declarado como no entendido, pendiente de medirse aparte.\n"
        "- No se envio nada ni se corrio `vcc sample`. No se hizo ningun commit.\n"
        "- No se leyo NADA bajo `~/vcc/benchmark/` ni los `.h5ad` de "
        "`/datos/benchmark_VEDADO_techo_metricas`. Ver la nota al final.\n"
    )

    L.append("## `avg_score`: 0.0 cuando corre, y donde `score` se niega a correr\n")
    L.append(
        "`score_vs_self` (reusado sin modificar) pasa el mismo `agg_results.csv` "
        "como `--user-agg` y `--baseline-agg`, asi que `from_baseline` cae en 0.0 "
        "por construccion. La columna se conserva solo por paridad de forma con "
        "`bucle_suelo.md`; lo informativo son las seis columnas de metricas "
        "crudas.\n"
    )
    con_error = [r for r in techo if r.get("score_error")]
    if con_error:
        L.append(
            f"En **{len(con_error)} de {len(techo)} corridas** `score` se nego a "
            "puntuar, y el motivo es un resultado del techo, no un fallo del "
            "arnes:\n"
        )
        L.append("| linea | var | dir | mse_capped_norm del techo | motivo |")
        L.append("|---|---|---|---|---|")
        for r in sorted(con_error, key=lambda x: (x["variante"], x["direccion"],
                                                  x["linea_excluida"])):
            L.append(
                f"| {r['linea_excluida']} | {r['variante']} | {r['direccion']} | "
                f"{fmt(r.get('expr_mse_unbiased_capped_norm'))} | "
                f"{r['score_error'][:150]} |"
            )
        L.append(
            "\n`expr_mse_unbiased_capped_norm` es un estimador **INSESGADO** y "
            "**con signo**. Cuando las dos mitades concuerdan -- que es "
            "exactamente lo que mide el techo -- el MSE verdadero tiende a 0, y "
            "un estimador insesgado de algo que vale ~0 sale negativo la mitad de "
            "las veces. `score` entonces no puede definir una escala "
            "(`denominator ... is not a finite positive number`) y se niega, "
            "incluso con `--allow-config-mismatch`. Consecuencia practica, que "
            "conviene no descubrir mas tarde: **el techo no sirve como "
            "`--baseline-agg` para esa metrica**; sirve como el `r` de la escala "
            "`(u - b)/(r - b)` calculada a mano, con `r` cerca de 0. Es tambien "
            "la razon por la que `cell_eval2/ceiling.py` excluye esa metrica de "
            "su lista `SB_METRICS`, citando que \"is additionally signed, so it "
            "is not a bounded reliability in the first place\".\n"
        )

    L.append("## La escala que esto habilita\n")
    L.append(
        "Con los dos anclajes medidos en la misma definicion del marcador, un "
        "candidato `u` se lee por metrica como `(u - b) / (r - b)`, con `b` el "
        "suelo de `bucle_suelo.md` y `r` el techo de aqui. Dos advertencias que "
        "vienen de este reporte y no deben perderse al usar esa escala: `r` esta "
        "medido a la mitad de la profundidad de `b` (cota inferior, sin "
        "corregir), y solo la variante S comparte `control_source` con `b`.\n"
    )

    L.append("## Nota sobre el limite duro de `~/vcc/benchmark/`\n")
    L.append(
        "No se abrio ningun archivo bajo `~/vcc/benchmark/`, ni "
        "`techo_metricas.md`, ni `techo_estado.md`, ni `b_techo_i6_ruido.md`, ni "
        "nada en `benchmark/work/techo_metricas_work/`, ni los `.h5ad` de "
        "`/datos/benchmark_VEDADO_techo_metricas`. Lo unico que esta sesion vio "
        "de ahi fue el **nombre** del directorio `benchmark/` en un `ls` del raiz "
        "del proyecto, al ubicar `modelo/`; no se listo su interior ni se leyo "
        "nada. Ese trabajo responde a otra pregunta (como responde el marcador a "
        "un envio transformado) sobre otra poblacion (170,846 celulas del split "
        "test 2025) y no es un techo, asi que no habria servido aqui ni para "
        "comparar.\n"
    )
    L.append(
        "Las unicas fuentes de datos de este reporte son: "
        "`/datos/ref/convertido/Seurat_object_IFNB_Perturb_seq.h5ad` (IFNB), "
        "`modelo/out/suelo_panel.json` (el panel comun de 33), "
        "`modelo/out/bucle_suelo.json` y los `results.csv` de "
        "`/datos/modelo/out/run_*` (el suelo), y el codigo de `cell_eval2` "
        "0.20.0 en `~/cell-eval2/src/cell_eval2`.\n"
    )

    with open(OUT_MD, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"escrito {OUT_MD} con {len(techo)} corridas")


if __name__ == "__main__":
    main()
