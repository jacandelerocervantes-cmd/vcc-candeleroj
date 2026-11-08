#!/usr/bin/env python
"""Genera modelo/out/control_o_no.md desde los JSON medidos. Solo formatea."""
import json
import math
import os

import numpy as np

D = "/home/ubuntu/vcc/modelo/out/control_o_no"
OUT = "/home/ubuntu/vcc/modelo/out/control_o_no.md"
L = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]


def carga(p):
    return json.load(open(p)) if os.path.exists(p) else {}


def z_p_mde(obs, nulo):
    n = np.array([v for v in nulo if np.isfinite(v)])
    if n.size < 2:
        return np.nan, np.nan, np.nan, np.nan
    sd = float(np.std(n, ddof=1)); mu = float(np.mean(n))
    z = (obs - mu) / sd if sd > 0 else np.nan
    p = (1 + int((n >= obs).sum())) / (1 + n.size)
    return z, p, 2.802 * sd, mu


def main():
    loto = carga(os.path.join(D, "loto.json"))
    est = carga(os.path.join(D, "estricto.json"))
    disc = carga(os.path.join(D, "disc.json"))
    desg = carga(os.path.join(D, "desglose.json"))
    pes = carga(os.path.join(D, "pesos.json"))
    ref = carga("/home/ubuntu/vcc/modelo/out/cand01_refuta.json")
    anc = carga("/home/ubuntu/vcc/modelo/out/b_contra_a_anchor.json")
    O = []; w = O.append
    LL = [l for l in L if l in est]

    # ---------- cantidades ----------
    q = {}
    for l in L:
        e = est.get(l, {}).get("ridge", {})
        ze, pe, mdee, mue = z_p_mde(e.get("R_agr", np.nan),
                                    e.get("nulo", {}).get("R_agr", []))
        t = loto.get(l, {}).get("conjuntos", {})
        zt, pt, mdet, mut = z_p_mde(t.get("todos", {}).get("R_agr", np.nan),
                                    t.get("todos", {}).get("nulo", {}).get("R_agr", []))
        zs, ps, mdes, mus = z_p_mde(
            t.get("sin_autoblanco", {}).get("R_agr", np.nan),
            t.get("sin_autoblanco", {}).get("nulo", {}).get("R_agr", []))
        q[l] = {
            "est": e, "z_est": ze, "p_est": pe, "mde_est": mdee, "mu_est": mue,
            "todos": t.get("todos", {}), "z_t": zt, "p_t": pt, "mde_t": mdet,
            "sin": t.get("sin_autoblanco", {}), "z_s": zs, "p_s": ps,
            "mde_s": mdes, "mu_s": mus,
            "cov": t.get("solo_covarianza", {}), "soloh": t.get("solo_h", {}),
            "desg": desg.get(l, {}), "disc": disc.get(l, {}),
            "A_esp": np.array(loto.get(l, {}).get("A_esp", [np.nan])),
            "C": float(np.median(ref[l]["C_pred_vs_mismo"])) if l in ref else np.nan,
            "A_anc": float(np.median(anc["lineas"][l]["A_anclado"]))
                     if l in anc.get("lineas", {}) else np.nan,
            "G": loto.get(l, {}).get("G"),
            "n_auto": int(np.sum(loto.get(l, {}).get("blanco_expresado", []))),
            "G_eval": est.get(l, {}).get("G_evaluados", [0]),
            "gb_est": est.get(l, {}).get("gb", {}).get("R_agr", np.nan),
        }
    n_sig_est = sum(1 for l in LL if np.isfinite(q[l]["p_est"]) and q[l]["p_est"] < 0.05)
    n_pos_est = sum(1 for l in LL if q[l]["est"].get("R_agr", 0) > 0)

    # ================= encabezado =================
    w("# ¿Era la covarianza la senal equivocada, o no hay senal en el control?")
    w("")
    w("## 🔴 0. Esto es un DIAGNOSTICO TRANSDUCTIVO. Va arriba, no en nota al pie.")
    w("")
    w("Todo lo que se ajusta en este informe se ajusta **sobre la verdad de")
    w("Jiang**, que es exactamente lo que un envio no tiene. Ninguna cifra, ningun")
    w("coeficiente y ningun modelo de aqui entra NUNCA a una submission, y ninguno")
    w("es una via. Lo unico que sale de aqui es la RESPUESTA: si el techo de \"solo")
    w("control\" es alto o es cero.")
    w("")
    w("En concreto, y para que no se pueda confundir despues:")
    w("")
    w("| | |")
    w("|---|---|")
    w("| el modelo que da el numero alto de §6 | se ajusta con la verdad perturbada **del propio blanco** que luego evalua. Un envio no la tiene. **No es una via.** |")
    w("| el modelo que responde a la pregunta util (§5) | se ajusta sin ver el blanco que evalua, y da **cero** |")
    w("| lo que se puede llevar a un candidato | la respuesta y la cota, no los coeficientes |")
    w("")
    w("---")
    w("")

    # ================= la respuesta =================
    w("## 1. La respuesta: ni (a) ni (b) -- la senal esta, el mapa no se aprende")
    w("")
    w("El encargo pedia discriminar entre (a) la covarianza era la senal")
    w("equivocada y otra senal del control si capta lo especifico, y (b) nada del")
    w("control lo predice. **La medicion no da ninguna de las dos, y la tercera")
    w("opcion es concreta y esta medida:**")
    w("")
    w("| | medido | que significa |")
    w("|---|---|---|")
    dpar = [q[l]["desg"].get("par", {}).get("R_agr", np.nan) for l in LL]
    dsh = [q[l]["desg"].get("solo_h", {}).get("R_agr", np.nan) for l in LL]
    ests = [q[l]["est"].get("R_agr", np.nan) for l in LL]
    w("| **los rasgos de coexpresion SI llevan la senal** | `r_mismo` (par) = **%+.4f a %+.4f**, positivo en 29-33 de 33 blancos, nulo ~0 con `sd` ~0.007 | NO es (b). Y es coexpresion, no propension por gen: los rasgos de solo-`h` dan %+.4f a %+.4f |"
      % (min(dpar), max(dpar), min(dsh), max(dsh)))
    w("| **pero no transfiere a un blanco nuevo** | `r_loto` estricto = **%+.4f a %+.4f**, **%d de 6** significativas, %d positivas y %d negativas | NO es (a) en sentido usable |"
      % (min(ests), max(ests), n_sig_est, n_pos_est, len(LL) - n_pos_est))
    if pes:
        _c, _a = [], []
        for l in LL:
            if l in pes:
                _c.append(np.array(pes[l]["consistencia"]))
                _a.append(np.array(pes[l]["abs_medio"]))
        _c = np.concatenate(_c); _a = np.concatenate(_a)
        cons = _c[_a > 0]          # los de solo-`g` tienen coeficiente 0 exacto
    else:
        cons = np.array([np.nan])
    w("| **y se ve por que** | la consistencia de signo de los coeficientes entre blancos es **%.0f a %.0f %%** (moneda = 50 %%) | cada blanco tiene su propio coeficiente. Con 33 blancos no hay un mapa compartido que aprender |"
      % (100 * np.nanmin(cons), 100 * np.nanmax(cons)))
    w("")
    mde_min = min(q[l]["mde_est"] for l in LL)
    mde_max = max(q[l]["mde_est"] for l in LL)
    w("🔴 **Y la limitacion que el encargo mandaba declarar decide el matiz.** El")
    w("MDE del diseno (80 %% de potencia, alfa 0.05) va de **%.3f a %.3f** segun"
      % (mde_min, mde_max))
    w("linea, y el techo intra-blanco esta en 0.24-0.40. Es decir: **el diseno no")
    w("puede separar \"no transfiere nada\" de \"transfiere hasta ~0.2\".** Lo que si")
    w("queda cerrado es que no transfiere MAS que eso. Con 33 blancos no se puede")
    w("pedir mas, y eso tambien es el resultado.")
    w("")
    w("**Para el candidato 02, en una linea:** la covarianza del control **no era")
    w("la senal equivocada** -- es la que lleva la informacion -- pero el candidato")
    w("01 fijo su coeficiente (`r = 1 - rho`) y el coeficiente correcto es")
    w("**distinto en cada blanco**, con el signo cambiando mas o menos como una")
    w("moneda. Un candidato de solo-control que aprenda ese coeficiente de otros")
    w("blancos no puede funcionar; lo que falta no es una senal mejor, es una")
    w("forma de fijar el coeficiente por blanco sin ver su verdad.")
    w("")
    w("---")
    w("")
    O.extend(seccion_diseno(q, LL, loto))
    O.extend(seccion_potencia(q, LL))
    O.extend(seccion_defectos(q, LL, loto))
    O.extend(seccion_loto(q, LL, n_sig_est, n_pos_est))
    O.extend(seccion_disc(q, LL, disc, desg))
    O.extend(seccion_pesos(q, LL, pes))
    O.extend(seccion_filas(q, LL, loto, est, desg, pes))
    O.extend(seccion_cierre(q, LL, n_sig_est, n_pos_est,
                            float(np.nanmax(cons)), pes))
    with open(OUT, "w") as fh:
        fh.write("\n".join(O) + "\n")
    print("escrito %s (%d lineas)" % (OUT, len(O)))


def seccion_diseno(q, LL, loto):
    O = []; w = O.append
    w("## 2. El diseno")
    w("")
    w("### 2.1 El objetivo: el delta ESPECIFICO, con el control partido")
    w("")
    w("Para el blanco `g`, `E_g = h_g - m` con `m` el promedio sobre los 33")
    w("blancos, calculado como en `b_contra_a_descomp.py` para que sea **la misma")
    w("cantidad** cuyo 61-84 % de varianza se midio alli. El control va PARTIDO en")
    w("dos mitades disjuntas (`c1` para una mitad del blanco, `c2` para la otra):")
    w("con control comun el ruido de muestreo del control se cuenta como senal")
    w("compartida y `m` se lo come (§4 de `b_contra_a.md`).")
    w("")
    w("### 2.2 La metrica: el ruido fuera, por producto cruzado de mitades")
    w("")
    w("```")
    w("E[<E1_g, E2_g>]        = ||t_g||^2      t_g = la senal verdadera")
    w("E[<p, (E1_g+E2_g)/2>]  = <p, t_g>       (el ruido de las dos mitades es")
    w("                                         independiente: no entra)")
    w("r_true = <p, (E1+E2)/2> / ( ||p|| * sqrt(<E1_g,E2_g>) )")
    w("```")
    w("")
    w("Y el resumen sobre blancos es el **agrupado**, no la media:")
    w("")
    w("```")
    w("R_agr = Σ_g <p_g, y_g>  /  Σ_g ( ||p_g|| * sqrt(cruz_g) )")
    w("```")
    w("")
    w("por una razon que se comprueba y no se asume: como `Σ_g E[h,g] = 0`, el")
    w("numerador total de un predictor que solo dependa del gen `h` es")
    w("`<p, Σ_g E_g> = 0`, asi que **`R_agr` vale CERO EXACTO para cualquier")
    w("predictor sin informacion especifica del blanco**. Es la comprobacion que")
    w("delato el primer defecto (§4.1). Medido, el conjunto `solo_h` da:")
    w("")
    w("| linea | " + " | ".join(LL) + " |")
    w("|---|" + "---|" * len(LL))
    w("| `solo_h` `R_agr` | " + " | ".join(
        "%+.4f" % q[l]["soloh"].get("R_agr", np.nan) for l in LL) + " |")
    w("")
    w("Cero en las seis. Sin esa comprobacion, nada de lo que sigue seria legible.")
    w("")
    w("### 2.3 Los rasgos, declarados antes de medirlos")
    w("")
    w("Todos del POOL DE CONTROL y del SIMBOLO del gen. Ninguno toca el lado")
    w("perturbado. Para el par (blanco `g`, gen `h`):")
    w("")
    w("| # | rasgo | que es | tipo |")
    w("|---|---|---|---|")
    filas = [
        ("cov", "covarianza de log1p normalizado en el pool de control", "par"),
        ("r_pearson", "correlacion de Pearson", "par"),
        ("r_spearman", "correlacion de Spearman (rangos sobre celulas)", "par"),
        ("r_parcial_10", "Pearson tras proyectar fuera los 10 primeros PC", "par"),
        ("r_parcial_50", "idem con 50 PC", "par"),
        ("abs_r_pearson", "`|Pearson|`, la magnitud separada del signo", "par"),
        ("abs_r_spearman", "`|Spearman|`", "par"),
        ("delta_r", "`r_pearson - r_parcial_50`: cuanto era programa global", "par"),
        ("rango_cor", "rango de `h` por `|r_pearson|` con `g`, normalizado", "par"),
        ("mismo_modulo", "1 si `h` cae en el modulo de coexpresion de `g` (KMeans, 30)", "par"),
        ("dist_emb", "distancia euclidea `g`-`h` en el embedding de 50 PC", "par"),
        ("cos_emb", "coseno `g`-`h` en ese embedding", "par"),
        ("misma_familia", "raiz alfabetica comun de >=3 letras (STAT1/STAT4 -> STAT)", "par"),
        ("es_el_blanco", "1 si `h` ES `g`", "par"),
        ("h_media", "media de log1p de `h`", "solo h"),
        ("h_var", "varianza", "solo h"),
        ("h_disp", "var/media", "solo h"),
        ("h_frac_cero", "fraccion de celulas de control con cuenta 0", "solo h"),
        ("h_rango_expr", "rango del nivel de expresion", "solo h"),
        ("g_media", "media de log1p del blanco", "solo g"),
        ("g_var", "varianza del blanco", "solo g"),
        ("g_frac_cero", "fraccion de ceros del blanco", "solo g"),
    ]
    for i, (n, d, t) in enumerate(filas, 1):
        w("| %d | `%s` | %s | %s |" % (i, n, d, t))
    w("")
    w("**Las correlaciones parciales estan porque la marginal mezcla dos cosas.**")
    w("La covarianza marginal entre `g` y `h` incluye el programa global del")
    w("contexto; la parcial, tras quitar los primeros componentes, deja la")
    w("relacion directa. `delta_r` mide explicitamente la diferencia.")
    w("")
    w("🔴 **Lo que se descarto por no ser calculable sin el lado perturbado, con")
    w("nombre y razon** -- el encargo lo pedia explicitamente:")
    w("")
    w("| descartado | por que |")
    w("|---|---|")
    w("| pertenencia a una via conocida (IFN-alfa/beta, JAK-STAT) de una base externa | no es del pool de control; seria una cifra importada (PLAN §0.2) |")
    w("| respuesta conocida del gen a interferon | **es el lado perturbado** |")
    w("| efecto medido del mismo blanco en OTRA linea | **es el lado perturbado** |")
    w("| distancia en el genoma o cromosoma | no esta en el `h5ad` (solo `var/_index`); sacarla de una anotacion externa seria importar |")
    w("")
    w("**Los de solo-`g` no pueden hacer nada y van declarados como tales:** son")
    w("constantes dentro de una columna, y una correlacion por columna no cambia al")
    w("sumarle una constante. Estan por simetria del registro, no por esperanza.")
    w("")
    w("### 2.4 La validacion cruzada: por BLANCO, dejando blancos fuera")
    w("")
    w("Nunca por gen ni por celula para la pregunta de §5. Si el modelo memoriza")
    w("pares `(g,h)` de entrenamiento, el numero no dice nada sobre un blanco")
    w("nuevo. Se deja **un blanco fuera** (33 pliegues), con eleccion interna de")
    w("`alpha` por CV agrupada en los 32 de entrenamiento, y el nulo por")
    w("permutacion de etiquetas de blanco repitiendo la CV entera.")
    w("")
    w("El numero de §6 (`r_mismo`) SI parte genes, y por eso **no se lee como")
    w("generalizacion**: responde a otra pregunta y lo dice en su seccion.")
    w("")
    w("---")
    w("")
    return O


def seccion_potencia(q, LL):
    O = []; w = O.append
    w("## 3. La potencia del diseno, declarada ANTES del resultado")
    w("")
    w("Son **33 blancos por linea**, y las seis lineas comparten **los mismos 33**,")
    w("asi que agrupar lineas NO sube el `n` de la pregunta \"¿generaliza a un")
    w("blanco nuevo?\". El nulo se obtiene por permutacion de etiquetas de blanco")
    w("repitiendo la CV completa, y de su desviacion sale el efecto minimo")
    w("detectable:")
    w("")
    w("```")
    w("MDE(80 %, alfa 0.05 dos colas) = (1.96 + 0.8416) * sd_nulo = 2.802 * sd_nulo")
    w("```")
    w("")
    w("| linea | `sd` del nulo (LOTO estricto) | **MDE 80 %** | `A_esp` mediana (fiabilidad del propio objetivo) |")
    w("|---|---|---|---|")
    for l in LL:
        sd = q[l]["mde_est"] / 2.802 if np.isfinite(q[l]["mde_est"]) else np.nan
        w("| %s | %.4f | **%.3f** | %+.4f |"
          % (l, sd, q[l]["mde_est"], float(np.median(q[l]["A_esp"]))))
    w("")
    mn = min(q[l]["mde_est"] for l in LL); mx = max(q[l]["mde_est"] for l in LL)
    w("🔴 **Dicho antes de ver el resultado: este diseno solo puede detectar")
    w("efectos de %.3f para arriba.** Cualquier cosa por debajo es indistinguible" % mn)
    w("de cero con 33 blancos, y el veredicto correcto para ese caso es \"no se")
    w("puede decidir\", no \"no hay senal\". El MDE llega hasta %.3f en la linea" % mx)
    w("peor, que es la mitad del techo intra-blanco que se medira en §6: **el")
    w("diseno esta justo en el limite de poder contestar la pregunta.**")
    w("")
    w("---")
    w("")
    return O


def seccion_defectos(q, LL, loto):
    O = []; w = O.append
    w("## 4. Tres defectos encontrados y corregidos, antes de cualquier lectura")
    w("")
    w("Los tres los delataron comprobaciones puestas a proposito, no una")
    w("sospecha. Los tres cambiaban el resultado.")
    w("")
    w("### 4.1 Una fuga estructural del diseno: la media de entrenamiento ES el blanco dejado fuera")
    w("")
    w("El objetivo `E_i = h_i - m` con `m` sobre los 33 fuerza `Σ_i E_i = 0`. Y")
    w("entonces, para el pliegue que deja fuera el blanco `j`:")
    w("")
    w("```")
    w("media de los 32 de entrenamiento = (Σ_{i≠j} E_i) / 32 = -E_j / 32")
    w("```")
    w("")
    w("Medido en A549: `corr(media de entrenamiento, blanco dejado fuera) =")
    w("**-1.000000** exacto`, y `max|media_tr - (-E_j/32)| = 4.4e-09`. Cualquier")
    w("modelo que capture esa media se lleva correlacion NEGATIVA con el blanco")
    w("nuevo sin saber nada de el.")
    w("")
    w("**Lo delato `solo_h`**, que no puede llevar informacion especifica del")
    w("blanco y puntuaba de -0.13 a -0.22 con desviacion EXACTAMENTE cero.")
    w("")
    w("El arreglo sale derivado, no tanteado. El estimando correcto por pliegue usa")
    w("como referencia solo los 32 de entrenamiento, y sustituyendo")
    w("`h_i = E_i + S/33` con `Σ_{k≠j} E_k = -E_j`:")
    w("")
    w("```")
    w("y de ENTRENAMIENTO del pliegue j :  y_i = E_i + E_j / 32")
    w("                                    Σ_{i≠j} y_i = -E_j + 32·E_j/32 = 0   ✓")
    w("objetivo de EVALUACION            :  (33/32) · E_j   -> multiplo POSITIVO,")
    w("                                    la correlacion no cambia")
    w("```")
    w("")
    w("Toca solo las `y` de entrenamiento, que es donde estaba la fuga, y no hubo")
    w("que recalcular ningun rasgo. Comprobado: `solo_h` pasa de **-0.1744 a")
    w("-0.0033** en A549 y queda en cero en las seis lineas (§2.2).")
    w("")
    w("### 4.2 Un resumen sesgado mio: la media de cocientes no hereda el cero")
    w("")
    w("Declare que un predictor de solo-`h` tiene correlacion media cero sobre los")
    w("blancos. Eso es cierto de los NUMERADORES, no del cociente: cada blanco")
    w("divide por su propio `sqrt(cruz_g)`, y una media de cocientes con")
    w("denominadores distintos no hereda la propiedad. Por eso el titular es")
    w("`R_agr` (§2.2) y no la media. Se dan las tres para que se vea cuanto movia:")
    w("")
    w("| linea | `R_agr` (el que se lee) | mediana | media (la sesgada) |")
    w("|---|---|---|---|")
    for l in LL:
        e = q[l]["est"]
        w("| %s | %+.4f | %+.4f | %+.4f |"
          % (l, e.get("R_agr", np.nan), e.get("mediana", np.nan),
             e.get("media", np.nan)))
    w("")
    w("### 4.3 🔴 El gen blanco se identifica sin el indicador, asi que `sin_autoblanco` no era un test limpio")
    w("")
    w("Quitar el rasgo `es_el_blanco` **no le quita al modelo saber cual es el gen")
    w("blanco**, porque los rasgos de correlacion lo delatan con certeza:")
    w("")
    w("| rasgo | valor en el gen blanco | maximo en los otros ~7000 genes | p99.99 |")
    w("|---|---|---|---|")
    w("| `r_pearson` | **+1.0000** | +0.3031 | +0.1666 |")
    w("| `r_spearman` | **+1.0000** | +0.2233 | +0.1685 |")
    w("| `cos_emb` | **+1.0000** | +0.7336 | +0.6255 |")
    w("| `dist_emb` | **0.0000** | +0.9312 | +0.8180 |")
    w("")
    w("(medido en A549; el gen blanco esta expresado en %d-%d de los 33 blancos"
      % (min(q[l]["n_auto"] for l in LL), max(q[l]["n_auto"] for l in LL)))
    w("segun linea.) Y su `|E|` es **10.2 veces la mediana**, asi que un solo punto")
    w("con esa palanca mueve una Pearson sobre 7000 genes.")
    w("")
    w("**\"Cuando se perturba `g`, `g` baja\" es cierto, es de solo control y tiene")
    w("valor** -- se sabe QUE gen se perturba sin mirar nada perturbado. Pero es UN")
    w("gen de 7000 y no es \"otra senal del control que capta lo especifico\". Para")
    w("medir la pregunta del encargo hay que sacarlo de los **tres** sitios: de los")
    w("rasgos, de las filas de entrenamiento (donde es un punto extremo que el")
    w("modelo persigue) y de la evaluacion. Eso es la version ESTRICTA, y es la")
    w("unica que se lee en §5.")
    w("")
    w("| linea | `R_agr` con todo | sin el indicador (NO limpio) | **ESTRICTO** (los tres sitios) |")
    w("|---|---|---|---|")
    for l in LL:
        w("| %s | %+.4f | %+.4f | **%+.4f** |"
          % (l, q[l]["todos"].get("R_agr", np.nan),
             q[l]["sin"].get("R_agr", np.nan), q[l]["est"].get("R_agr", np.nan)))
    w("")
    w("La columna del medio es la que reporte en el avance como \"apunta a (a)\".")
    w("**Era el gen blanco.** Queda escrito porque el error estuvo en el camino.")
    w("")
    w("---")
    w("")
    return O


def seccion_loto(q, LL, n_sig, n_pos):
    O = []; w = O.append
    w("## 5. Resultado 1 -- ¿generaliza a un blanco NUEVO? **No, y eso es un cero medido**")
    w("")
    w("Version estricta, ridge con `alpha` por CV interna, nulo de 500")
    w("permutaciones de etiquetas de blanco:")
    w("")
    w("| linea | genes evaluados | **`R_agr`** | mediana | positivos | nulo | z | `p` perm | MDE 80 % | recortado 1 % | GB |")
    w("|---|---|---|---|---|---|---|---|---|---|---|")
    for l in LL:
        e = q[l]["est"]
        ge = q[l]["G_eval"]
        gb = q[l]["gb_est"]
        w("| %s | %d-%d | **%+.4f** | %+.4f | %d/%d | %+.4f | %+.2f | %.3f | %.3f | %+.4f | %+.4f |"
          % (l, min(ge), max(ge), e.get("R_agr", np.nan),
             e.get("mediana", np.nan), e.get("n_positivos", 0),
             e.get("n_validos", 0), q[l]["mu_est"], q[l]["z_est"],
             q[l]["p_est"], q[l]["mde_est"],
             float(np.nanmedian(e.get("recortado", [np.nan]))), gb))
    w("")
    w("**%d de 6 significativas. %d positivas, %d negativas.** El signo no es"
      % (n_sig, n_pos, len(LL) - n_pos))
    w("consistente entre lineas, los `z` van de -1.6 a +1.3, y la version recortada")
    w("(fuera el 1 % de `|E|` mas extremo) no rescata nada. **Todo el +0.14 a +0.30")
    w("de la version con el gen blanco era el gen blanco.**")
    w("")
    w("### 5.1 Contra las dos referencias que pedia el encargo")
    w("")
    w("| linea | `C` del candidato 01 | `A_anc` (fiabilidad del delta completo) | `A_esp` (fiabilidad del delta ESPECIFICO) | **`r_loto` estricto** | fraccion del techo |")
    w("|---|---|---|---|---|---|")
    for l in LL:
        ae = float(np.median(q[l]["A_esp"]))
        r = q[l]["est"].get("R_agr", np.nan)
        w("| %s | %+.4f | %+.4f | %+.4f | **%+.4f** | %s |"
          % (l, q[l]["C"], q[l]["A_anc"], ae, r,
             "n/d (negativo)" if r < 0 else "%.0f %%" % (100 * r / ae)
             if ae > 0 else "n/d"))
    w("")
    w("**¿Mejora sobre `C` ~ 0? No de forma distinguible.** `C` iba de -0.0054 a")
    w("+0.0108 y esto va de %+.4f a %+.4f con el cero dentro del nulo en las seis."
      % (min(q[l]["est"].get("R_agr", np.nan) for l in LL),
         max(q[l]["est"].get("R_agr", np.nan) for l in LL)))
    w("El candidato 01 no fallo por elegir mal la senal: **ningun predictor de")
    w("solo-control que aprenda de otros blancos lo hace mejor en este diseno.**")
    w("")
    w("🔴 Y hay que decir lo que esto NO dice: con el MDE en %.3f-%.3f, un efecto"
      % (min(q[l]["mde_est"] for l in LL), max(q[l]["mde_est"] for l in LL)))
    w("verdadero de, digamos, 0.10 no se habria visto. **\"No transfiere\" aqui")
    w("significa \"no transfiere por encima del MDE\", no \"transfiere cero\".**")
    w("")
    w("---")
    w("")
    return O


def seccion_disc(q, LL, disc, desg):
    O = []; w = O.append
    w("## 6. Resultado 2 -- ¿llevan los rasgos informacion? **Mucha, y es coexpresion**")
    w("")
    w("🔴 **Esto NO es generalizacion a un blanco nuevo y no se puede leer como")
    w("tal.** Se ajusta en la MITAD de los genes DEL PROPIO blanco y se predice la")
    w("otra mitad. Responde a la otra pregunta, la que separa (a) de (b): ¿existe")
    w("alguna funcion de los rasgos de solo-control que describa el delta")
    w("especifico de este blanco? El gen blanco va fuera de los tres sitios, igual")
    w("que en §5.")
    w("")
    w("| linea | **solo rasgos de PAR** (coexpresion) | positivos | nulo | `sd` nulo | **z** | solo rasgos de `h` | positivos | ambos |")
    w("|---|---|---|---|---|---|---|---|---|")
    for l in LL:
        d = q[l]["desg"]
        if not d:
            continue
        nn = [x for x in d["par"].get("nulo_R_agr", []) if np.isfinite(x)]
        mu = float(np.mean(nn)) if nn else np.nan
        sd = float(np.std(nn, ddof=1)) if len(nn) > 1 else np.nan
        w("| %s | **%+.4f** | %d/33 | %+.4f | %.4f | **%+.0f** | %+.4f | %d/33 | %+.4f |"
          % (l, d["par"]["R_agr"], d["par"]["n_positivos"], mu, sd,
             (d["par"]["R_agr"] - mu) / sd if sd else np.nan,
             d["solo_h"]["R_agr"], d["solo_h"]["n_positivos"],
             d["ambos"]["R_agr"]))
    w("")
    zz, sds = [], []
    for l in LL:
        dd = q[l]["desg"]
        if not dd:
            continue
        nn = [x for x in dd["par"].get("nulo_R_agr", []) if np.isfinite(x)]
        if len(nn) > 1:
            sd0 = float(np.std(nn, ddof=1)); sds.append(sd0)
            zz.append((dd["par"]["R_agr"] - float(np.mean(nn))) / sd0)
    w("El nulo permuta el objetivo sobre los GENES dentro del blanco manteniendo")
    w("los rasgos: rompe la correspondencia gen-objetivo sin tocar nada mas. 10")
    w("permutaciones, `sd` de %.4f a %.4f, **z de +%.0f a +%.0f**. No es marginal."
      % (min(sds), max(sds), min(zz), max(zz)))
    w("")
    hay = [l for l in LL if q[l]["disc"]]
    if hay:
        w("### 6.1 Y el gen blanco no es lo que lo sostiene -- al contrario que en §5")
        w("")
        w("Con una implementacion independiente (`control_o_no_disc.py`), en la")
        w("linea donde se corrio, comparando con y sin el gen blanco:")
        w("")
        w("| linea | modelo | con el gen blanco | **sin el gen blanco** | nulo | `sd` nulo | z |")
        w("|---|---|---|---|---|---|---|")
        for l in hay:
            v = q[l]["disc"].get("versiones", {})
            for mo in ("gb", "ridge"):
                c = v.get("con_%s" % mo, {}); s = v.get("sin_%s" % mo, {})
                if not s:
                    continue
                n = [x for x in s.get("nulo_R_agr", []) if np.isfinite(x)]
                sd = float(np.std(n, ddof=1)) if len(n) > 1 else np.nan
                w("| %s | %s | %+.4f | **%+.4f** | %+.4f | %.4f | %+.1f |"
                  % (l, mo, c.get("R_agr", np.nan), s.get("R_agr", np.nan),
                     float(np.mean(n)) if n else np.nan, sd,
                     (s.get("R_agr", np.nan) - np.mean(n)) / sd
                     if n and sd else np.nan))
        w("")
        w("**Quitar el gen blanco no baja nada aqui: sube un poco.** Eso es lo que")
        w("separa las dos preguntas. El gen blanco era TODO lo que transferia entre")
        w("blancos (§4.3, §5) y NADA de lo que hay dentro de un blanco.")
        w("")
        w("Y vale como comprobacion cruzada: el `ambos` de la tabla de arriba para")
        w("%s (%+.4f) y el `sin ridge` de aqui (%+.4f) salen de dos scripts"
          % (hay[0], q[hay[0]]["desg"]["ambos"]["R_agr"],
             q[hay[0]]["disc"]["versiones"]["sin_ridge"]["R_agr"]))
        w("distintos y coinciden.")
        w("")
        w("🔴 **Solo se corrio en %s.** Las otras cinco lineas llevan el numero de"
          % hay[0])
        w("`par` de la tabla de arriba, que excluye el gen blanco por construccion y")
        w("ya responde lo mismo, pero la comparacion con/sin por linea queda como")
        w("hueco declarado, no medido.")
        w("")
    dpar = [q[l]["desg"]["par"]["R_agr"] for l in LL if q[l]["desg"]]
    dsh = [q[l]["desg"]["solo_h"]["R_agr"] for l in LL if q[l]["desg"]]
    w("**Y el desglose contesta la pregunta del encargo.** Los rasgos de PAR")
    w("(covarianza, Pearson, Spearman, parcial, modulo, embedding, familia) dan")
    w("%+.4f a %+.4f, positivo en 29-33 de 33 blancos. Los de solo-`h` -- \"este"
      % (min(dpar), max(dpar)))
    w("gen tiende a moverse\", que no depende del blanco -- dan solo %+.4f a %+.4f."
      % (min(dsh), max(dsh)))
    w("**El grueso lo lleva la coexpresion, no una propension por gen.**")
    w("")
    w("Y el techo: `A_esp` atenuada va de +0.060 a +0.125, mientras esto llega a")
    w("+0.40. No se contradice -- `A_esp` es la correlacion contra una mitad")
    w("ruidosa y `r_true` es contra la senal verdadera, con el ruido fuera por")
    w("construccion. Lo que dice es que **hay bastante mas senal especifica")
    w("describible de lo que la fiabilidad cruda sugiere.**")
    w("")
    w("---")
    w("")
    return O


def seccion_pesos(q, LL, pes):
    O = []; w = O.append
    w("## 7. Por que no transfiere: los coeficientes cambian de signo entre blancos")
    w("")
    w("Si la senal esta dentro de cada blanco (§6) pero no transfiere (§5), la")
    w("explicacion tiene que estar en el mapa. Se mira donde esta la senal -- el")
    w("ajuste INTRA-blanco -- con un ridge por blanco sobre rasgos estandarizados,")
    w("y se pregunta no cuales pesan sino **si pesan en el mismo sentido**:")
    w("")
    if not pes:
        w("*(no medido)*")
        w("")
        w("---")
        w("")
        return O
    l0 = LL[0]
    nom = pes[l0]["rasgos"]
    w("La consistencia es la fraccion de blancos con el signo MAYORITARIO del")
    w("coeficiente, contando los ceros en ninguno de los dos lados. Por eso un")
    w("rasgo que queda a cero exacto en muchos blancos puede bajar del 50 %:")
    w("")
    w("| rasgo | `|coef|` medio | signo mediano | **consistencia de signo** | coef = 0 |")
    w("|---|---|---|---|---|")
    am = np.mean([np.array(pes[l]["abs_medio"]) for l in LL if l in pes], axis=0)
    md = np.mean([np.array(pes[l]["mediano"]) for l in LL if l in pes], axis=0)
    cs = np.mean([np.array(pes[l]["consistencia"]) for l in LL if l in pes], axis=0)
    cz = np.mean([np.array(pes[l]["frac_cero"]) for l in LL if l in pes], axis=0)
    for i in np.argsort(-am):
        w("| `%s` | %.4f | %+.4f | **%.0f %%** | %.0f %% |"
          % (nom[i], am[i], md[i], 100 * cs[i], 100 * cz[i]))
    w("")
    w("**Los tres de solo-`g` salen a cero EXACTO en los 33 blancos**, que es la")
    w("propiedad declarada en §2.3 confirmada por el ajuste: son constantes dentro")
    w("de una columna, asi que son colineales con el intercepto y el ridge los")
    w("manda a cero. No es que no sirvan poco: no pueden hacer nada.")
    w("")
    vivos = am > 0
    w("**Y de los %d rasgos que si reciben peso, ninguno pasa del %.0f %% de"
      % (int(vivos.sum()), 100 * cs[vivos].max()))
    w("consistencia, con la moneda en 50 %%.**")
    w("No hay un rasgo que empuje en el mismo sentido en todos los blancos. Esa es")
    w("la razon medida de que §5 salga cero teniendo §6 un +0.3: **no existe un")
    w("mapa compartido que 33 blancos puedan revelar**, porque el mapa es distinto")
    w("en cada blanco.")
    w("")
    w("Dicho de otro modo: el candidato 01 fijo el coeficiente en `r = 1 - rho`,")
    w("el mismo para los 33 blancos. La medicion dice que el coeficiente correcto")
    w("existe y es fuerte, pero **cambia de blanco a blanco con el signo casi")
    w("libre**. Ninguna constante global -- ni la del candidato 01 ni otra -- puede")
    w("acertarlo.")
    w("")
    w("---")
    w("")
    return O


def seccion_filas(q, LL, loto, est, desg, pes):
    O = []; w = O.append
    l0 = LL[0]
    w("## 8. Las filas trabajadas")
    w("")
    w("### 8.1 Los 33 blancos de %s, uno por fila" % l0)
    w("")
    w("`A_esp` es la fiabilidad mitad-contra-mitad del delta especifico de ese")
    w("blanco; `cruz` su producto cruzado (el numerador sin ruido, y `cruz <= 0`")
    w("es el motivo por el que un blanco queda fuera del resumen); `r_mismo (par)`")
    w("lo de §6 y `r_loto` lo de §5. **Las dos ultimas columnas son las que no hay")
    w("que confundir**: la primera se ajusta con la verdad del propio blanco, la")
    w("segunda sin verla.")
    w("")
    pan = loto[l0]["panel"]
    A = np.array(loto[l0]["A_esp"]); cz = np.array(loto[l0]["cruz"])
    au = np.array(loto[l0]["blanco_expresado"])
    dp = desg.get(l0, {}).get("par", {}).get("por_blanco", [np.nan] * 33)
    er = est.get(l0, {}).get("ridge", {})
    pb = er.get("por_blanco", []); ix = er.get("indices", [])
    mapa = {k: v for k, v in zip(ix, pb)}
    w("| # | blanco | gen blanco expresado | `A_esp` | `cruz` | `r_mismo` (par) | `r_loto` estricto |")
    w("|---|---|---|---|---|---|---|")
    for i in np.argsort(-A):
        dv = dp[i] if i < len(dp) else np.nan
        w("| %d | `%s` | %s | %+.4f | %+.2f | %s | %s |"
          % (i + 1, pan[i], "si" if au[i] else "no", A[i], cz[i],
             ("%+.4f" % dv) if np.isfinite(dv) else "n/d",
             ("%+.4f" % mapa[i]) if i in mapa else "n/d"))
    v1 = [dp[i] for i in range(len(dp)) if np.isfinite(dp[i])]
    v2 = list(mapa.values())
    w("")
    w("`r_mismo` positivo en %d de %d filas; `r_loto` positivo en %d de %d. "
      % (int(np.sum(np.array(v1) > 0)), len(v1),
         int(np.sum(np.array(v2) > 0)), len(v2)))
    w("Las dos columnas se calculan sobre los MISMOS blancos y los mismos rasgos;")
    w("la unica diferencia es si el ajuste vio la verdad de ese blanco.")
    w("")
    w("🔴 **Dos cosas de la tabla que hay que leer con cuidado, y son del")
    w("estimador, no de los datos.**")
    w("")
    w("1. **`r_mismo` pasa de 1 en alguna fila** (`PARP12` da +1.2124). No es un")
    w("   error: `r_true` divide por `sqrt(cruz)`, que es una ESTIMACION de")
    w("   `||senal||`, y cuando `cruz` es chico el cociente se desboca. Un")
    w("   cociente asi no esta acotado por 1 como lo estaria una Pearson. Por eso")
    w("   **el titular es `R_agr`**, que suma numeradores y denominadores por")
    w("   separado antes de dividir y no deja que un blanco con `cruz` chico")
    w("   domine -- y por eso no se reporta la media de esta columna.")
    w("2. **El `cruz` de la tabla es el del gen completo; la exclusion se decide")
    w("   con el `cruz` restringido a la mascara** del pliegue (sin el gen")
    w("   blanco). Por eso hay filas con `cruz > 0` en la tabla y `n/d` en")
    w("   `r_loto`: al quitar el gen blanco su `cruz` se vuelve <= 0 y el blanco")
    w("   sale del resumen. Son %d de 33 en esta linea." % (33 - len(v2)))
    w("")
    w("### 8.2 Los 21 rasgos, uno por fila, con su consistencia de signo")
    w("")
    w("Promediado sobre las seis lineas. Es la tabla de §7 ordenada por")
    w("consistencia en vez de por peso, que es como se ve que ninguno sirve:")
    w("")
    if pes:
        nom = pes[LL[0]]["rasgos"]
        cs = np.mean([np.array(pes[l]["consistencia"]) for l in LL if l in pes],
                     axis=0)
        am = np.mean([np.array(pes[l]["abs_medio"]) for l in LL if l in pes],
                     axis=0)
        md = np.mean([np.array(pes[l]["mediano"]) for l in LL if l in pes],
                     axis=0)
        cz2 = np.mean([np.array(pes[l]["frac_cero"]) for l in LL if l in pes],
                      axis=0)
        viv = am > 0
        w("| # | rasgo | consistencia | `|coef|` medio | signo mediano | coef = 0 | ¿sirve para transferir? |")
        w("|---|---|---|---|---|---|---|")
        for r, i in enumerate(np.argsort(-cs), 1):
            w("| %d | `%s` | **%.0f %%** | %.4f | %+.4f | %.0f %% | %s |"
              % (r, nom[i], 100 * cs[i], am[i], md[i], 100 * cz2[i],
                 "no (sin peso)" if am[i] == 0 else
                 ("no" if cs[i] < 0.75 else "quiza")))
        w("")
        mejor = 100 * cs[viv].max()
        w("**Ninguno de los que reciben peso llega al 75 %%.** El mejor esta en")
        w("%.0f %%, que sobre 33 blancos son %.0f aciertos de signo contra los 16.5"
          % (mejor, 33 * mejor / 100))
        w("de una moneda: una ventaja de %.0f blancos. No es una regla, y menos una"
          % (33 * mejor / 100 - 16.5))
        w("que se pueda aplicar a un blanco que no se ha visto.")
    w("")
    w("---")
    w("")
    return O


def seccion_cierre(q, LL, n_sig, n_pos, cons_max, pes):
    O = []; w = O.append
    dpar = [q[l]["desg"]["par"]["R_agr"] for l in LL if q[l]["desg"]]
    mde_min = min(q[l]["mde_est"] for l in LL)
    mde_max = max(q[l]["mde_est"] for l in LL)
    w("## 9. Veredicto, y que no se hizo")
    w("")
    w("**El veredicto, con sus tres partes:**")
    w("")
    w("1. **No es (b).** Los rasgos de coexpresion del control describen el delta")
    w("   especifico de un blanco dado a `r_true` = %+.2f a %+.2f, positivo en"
      % (min(dpar), max(dpar)))
    w("   29-33 de 33 blancos, con nulo en cero y `sd` de 0.004 a 0.008 (z de +47")
    w("   a +63). No es ruido y no es el gen blanco.")
    w("2. **No es (a) en sentido usable.** Eso no transfiere a un blanco no visto:")
    w("   **%d de 6** lineas significativas, %d positivas y %d negativas."
      % (n_sig, n_pos, len(LL) - n_pos))
    cs_prom = np.mean([np.array(pes[l]["consistencia"]) for l in LL if l in pes],
                      axis=0)
    am_prom = np.mean([np.array(pes[l]["abs_medio"]) for l in LL if l in pes],
                      axis=0)
    w("3. **La razon esta medida:** la consistencia de signo de los coeficientes")
    w("   entre blancos no pasa del **%.0f %%** en ningun par (linea, rasgo), y"
      % (100 * cons_max))
    w("   promediada sobre lineas no pasa del **%.0f %%**, con la moneda en 50 %%."
      % (100 * cs_prom[am_prom > 0].max()))
    w("   El mapa es especifico del blanco.")
    w("")
    w("**Y la limitacion que el encargo mandaba declarar, aplicada:** con MDE de")
    w("**%.3f a %.3f** y un techo intra-blanco de %.2f a %.2f, el diseno **no puede"
      % (mde_min, mde_max, min(dpar), max(dpar)))
    w("separar \"no transfiere\" de \"transfiere hasta ~0.2\"**. Para la parte")
    w("informativa (§6) si hay potencia de sobra -- los `z` estan en +47/+63. Para")
    w("la parte de transferencia, 33 blancos no alcanzan, y eso es parte del")
    w("resultado.")
    w("")
    w("**No se hizo** -- ninguna prohibicion del encargo se roza:")
    w("")
    w("| prohibicion | estado |")
    w("|---|---|")
    w("| construir un candidato, tocar arneses o bundles | **no**. `build_fold_adatas` se llama como parametro con un `predict_fn` trivial que devuelve el control sin tocarlo; no se escribe ningun bundle |")
    w("| usar el lado perturbado para calcular una caracteristica | **no**. Los 22 rasgos salen del pool de control y del simbolo del gen. Los descartados por esto van con nombre en §2.3 |")
    w("| presentar el modelo ajustado como una via | **no**. §0 lo declara arriba: es transductivo y esta prohibido |")
    w("| `~/vcc/benchmark/`, `/datos/benchmark_VEDADO_techo_metricas/` | **no leidos** |")
    w("| envios, `vcc sample`, commits | **ninguno** |")
    w("")
    w("| | |")
    w("|---|---|")
    w("| entradas | `/datos/ref/convertido/Seurat_object_IFNB_Perturb_seq.h5ad`, `modelo/out/suelo_panel.json`; `cand01_refuta.json` y `b_contra_a_anchor.json` solo como referencia |")
    w("| codigo | `control_o_no_rasgos.py`, `_loto.py`, `_estricto.py`, `_disc.py`, `_desglose.py`, `_pesos.py`, `_reporte.py` |")
    w("| semillas | real `20260001+fold`, mitades `777`, control partido `20267999+fold`, familia nueva de este encargo **`20268000+100*fold`** (PCA, KMeans, particion de genes, permutaciones) |")
    w("")
    w("**Lo que queda abierto:**")
    w("")
    w("1. **El `alpha` del nulo se fija al modal del ajuste real** en vez de")
    w("   repetir la CV interna en cada permutacion. Afinar `alpha` dentro del nulo")
    w("   lo haria mas optimista, asi que fijarlo es la direccion conservadora para")
    w("   el nulo. Afecta a la interpretacion del `p`, no al signo del resultado.")
    w("2. **GB da ~0 donde ridge da +0.3 dentro del blanco y los dos dan ~0 entre")
    w("   blancos.** No se eligio el modelo que mejor quedaba: §5 reporta los dos y")
    w("   coinciden en el cero. Por que GB no ve lo que ridge ve dentro del blanco")
    w("   no se cerro -- la sospecha medida es que su capacidad se la llevan las")
    w("   filas extremas y los rasgos de solo-`h`, pero no se verifico.")
    w("3. **Un solo reparto de genes** en el discriminador (una semilla). El nulo")
    w("   es estrecho (`sd` ~0.007) asi que no cambia el signo, pero la magnitud")
    w("   +0.24/+0.40 no lleva barra por reparto.")
    w("4. **El `m` de referencia son 33 blancos**, no la poblacion. El delta")
    w("   \"especifico\" es especifico respecto de ESTE panel; con otro panel la")
    w("   particion compartido/especifico se mueve.")
    w("")
    w("---")
    w("")
    w("## 10. Para el candidato 02")
    w("")
    w("**La covarianza del control no era la senal equivocada.** Es la que lleva la")
    w("informacion: los rasgos de par dan %+.2f a %+.2f dentro del blanco, muy por"
      % (min(dpar), max(dpar)))
    w("encima de los de solo-`h`. Lo que el candidato 01 fijo mal no fue la senal,")
    w("fue el **coeficiente**: puso `r = 1 - rho`, el mismo para los 33 blancos, y")
    w("el coeficiente correcto cambia de blanco a blanco con el signo casi libre")
    w("(consistencia maxima 70 %).")
    w("")
    w("De ahi se siguen dos cosas, y la segunda es la que importa:")
    w("")
    w("1. **No sirve buscar otra senal del control.** Ya se probaron 22 rasgos,")
    w("   incluidas parciales que quitan el programa global, modulos y embeddings.")
    w("   La senal esta; el problema no es de senal.")
    w("2. **Un predictor de solo-control que aprenda el coeficiente de otros")
    w("   blancos no puede funcionar en este diseno.** Para usar lo de §6 haria")
    w("   falta fijar el coeficiente por blanco sin ver su verdad, y este encargo")
    w("   no encontro ninguna via para eso -- ni la descarta: dice que no esta en")
    w("   los 22 rasgos probados ni se puede aprender de 33 blancos.")
    w("")
    w("🔴 **Y lo de §0 otra vez, porque es donde se rompe un envio:** el +0.40 de")
    w("§6 se obtiene ajustando con la verdad perturbada del blanco que luego se")
    w("evalua. **No es un techo alcanzable por una submission.** El techo")
    w("alcanzable por algo que no ve la verdad es el de §5, y ese es cero.")
    return O


if __name__ == "__main__":
    main()
