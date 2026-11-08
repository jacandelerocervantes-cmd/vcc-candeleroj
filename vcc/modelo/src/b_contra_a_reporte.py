#!/usr/bin/env python
"""Genera modelo/out/b_contra_a.md desde los tres JSON medidos. Solo formatea:
no recalcula ninguna cantidad ni elige ningun umbral."""
import json
import math
import sys

import numpy as np

B = "/home/ubuntu/vcc/modelo/out/b_contra_a.json"
ANC = "/home/ubuntu/vcc/modelo/out/b_contra_a_anchor.json"
DES = "/home/ubuntu/vcc/modelo/out/b_contra_a_descomp.json"
REF = "/home/ubuntu/vcc/modelo/out/cand01_refuta.json"
OUT = "/home/ubuntu/vcc/modelo/out/b_contra_a.md"
LINEAS = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]


def binom_cola(k, n, p):
    """P(X >= k) exacta, n chico."""
    return sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i)
               for i in range(k, n + 1))


def ic(v):
    a = np.array(v)
    return float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))


def main():
    b = json.load(open(B)); ref = json.load(open(REF))
    anc = json.load(open(ANC)) if __import__("os").path.exists(ANC) else None
    des = json.load(open(DES)) if __import__("os").path.exists(DES) else None
    L = [l for l in LINEAS if l in b["lineas"]]
    nrep = b["meta"]["n_replicas"]
    nb = b["meta"]["n_blancos"]
    npar = b["meta"]["n_pares_completo"]
    O = []
    w = O.append

    # ---------- cantidades ----------
    q = {}
    for l in L:
        d = b["lineas"][l]
        A = np.array(d["A_mitades"]); Bc = np.array(d["B_completa"])
        Bm = np.array(d["B_mitades"])
        bo = d["bootstrap"]
        lo, hi = ic(bo["d_med"]); lom, him = ic(bo["d_mean"])
        loc, hic = ic(bo["dc_med"])
        q[l] = {
            "A": float(np.median(A)), "Bc": float(np.median(Bc)),
            "Bm": float(np.median(Bm)),
            "Am": float(A.mean()), "Bcm": float(Bc.mean()), "Bmm": float(Bm.mean()),
            "d_vieja": float(np.median(Bc) - np.median(A)),
            "d_nueva": float(np.median(Bm) - np.median(A)),
            "prof": float(np.median(Bc) - np.median(Bm)),
            "lo": lo, "hi": hi, "lom": lom, "him": him,
            "loc": loc, "hic": hic,
            "A_boot": float(np.mean(bo["A_med"])),
            "Bm_boot": float(np.mean(bo["Bm_med"])),
            "Bc_boot": float(np.mean(bo["Bc_med"])),
            "d_boot": float(np.mean(bo["d_med"])),
            "sd_d": float(np.std(bo["d_med"], ddof=1)),
            "excluye": (lo > 0 or hi < 0), "signo_neg": hi < 0,
            "ne": d["n_genes_expresados"], "nctrl": d["n_control_cells"],
            "d": d,
        }
    n_b_mayor_vieja = sum(1 for l in L if q[l]["d_vieja"] > 0)
    n_b_mayor_nueva = sum(1 for l in L if q[l]["d_nueva"] > 0)
    n_excluye = sum(1 for l in L if q[l]["excluye"])
    n_neg = sum(1 for l in L if q[l]["signo_neg"])
    n = len(L)

    # ---------- encabezado ----------
    w("# `B` contra `A` a profundidad igualada: la hipotesis queda REFUTADA")
    w("")
    w("`candidato_01.md` §1 midio `B > A` en cinco de las seis lineas y lo dejo")
    w("como **hipotesis**, no como resultado, por dos razones: 5 de 6 bajo una")
    w("moneda tiene `p = 7/64 = 0.109`, y ninguna de las dos cantidades llevaba")
    w("dispersion. Este encargo la cierra.")
    w("")
    w("🔴 **Se cierra en contra.** A profundidad igualada `B - A` **cambia de")
    w(f"signo en las {n} lineas**: `A > B`, no `B > A`. El IC 95 % de `B - A`")
    w(f"excluye el 0 en **{n_excluye} de {n}** lineas, y en las {n_neg} lo excluye por")
    w("debajo. La diferencia original era **profundidad**, no biologia: `A` se")
    w("medía entre mitades de 200 celulas y `B` entre poblaciones completas de")
    w("400, y la atenuacion por profundidad de `A` era **mayor** que el `B - A`")
    w("que se reporto.")
    w("")
    w("Y al declarar la construccion aparecio un **segundo** defecto que el encargo")
    w("no pedia: las dos columnas de cada correlacion se dividen por el **mismo**")
    w("pseudobulk de control, lo que `cell-eval2/src/cell_eval2/anchor.py:9-18`")
    w("declara sesgo **hacia arriba** y `anchor.py:768-779` **rehusa**")
    w("directamente. Se midio, y hay que decir exactamente que hace y que no")
    w("hace (§4): **no** mueve `B - A` -- infla `A` y `B` casi por el mismo")
    w("absoluto y se cancela en la diferencia, lo que vuelve esta refutacion mas")
    w("robusta, no menos. Lo que si rompe es el **nivel**: el `A` que")
    w("`candidato_01.md` presento como \"el techo alcanzable\" lee entre 0.05 y 0.12")
    w("de mas.")
    w("")
    w("**Para el candidato 02 esto cambia el diseno**: ya no hay apoyo medido para")
    w("\"el `b` captura la componente dominante de la varianza y lo especifico es")
    w("minoritario\". §5 mide esa fraccion directamente, que es lo que reemplaza a")
    w("la hipotesis.")
    w("")
    w("---")
    w("")

    # ---------- PASO 0 ----------
    w("## 0. PASO 0 - la construccion, escrita ANTES de anadir una sola barra")
    w("")
    w("La barra depende de la estructura, asi que primero la estructura. Leida de")
    w("`modelo/src/cand01_refuta.py`, no reconstruida de memoria:")
    w("")
    w("```")
    w("por cada blanco g del panel (33 blancos, 400 celulas reales cada uno):")
    w("  cr          = pseudobulk(TODAS las celulas de control de la linea)   # CP50k")
    w("  R[:, g]     = log2( pseudobulk(400 celulas de g) / cr )              # :99")
    w("  perm        = rng(777).permutation(indices de g)                     # :96")
    w("  H1[:, g]    = log2( pseudobulk(perm[:200]) / cr )                    # :100")
    w("  H2[:, g]    = log2( pseudobulk(perm[200:]) / cr )                    # :101")
    w("")
    w("mascara de genes (la MISMA para las dos cantidades):                   # :89")
    w("  exp = (cr >= 1) & (cp >= 1)      cp = pseudobulk del control PREDICHO")
    w("")
    w("rangos por columna via argsort estable, luego centrar y normalizar      # :48,:59")
    w("  nX = normaliza(rankea(X))        -> Spearman = producto escalar de columnas")
    w("")
    w("A = [ nH1[:,i] . nH2[:,i]  para i en 33 ]                              # :111")
    w("B = [ nR[:,i]  . nR[:,k]   para i<k en 33 ] -> 528 pares               # :116-118")
    w("")
    w("y lo que entra a la tabla de candidato_01.md es la MEDIANA de cada lista")
    w("  A = median(33 valores)   B = median(528 valores)   # cand01_reporte.py:125-126")
    w("```")
    w("")
    w("**Respuesta a las dos preguntas del PASO 0.** Spearman **por par y luego")
    w("mediana**, no un Spearman agrupado -- en las dos cantidades, con la misma")
    w("funcion de agregacion. Y **sobre los mismos genes**: la mascara `exp` se")
    w("calcula una vez por linea y la comparten `A` y `B`.")
    w("")
    w("**Veredicto de comparabilidad: SI son estructuralmente comparables en")
    w("forma** -- misma correlacion, mismo estimador por par, misma agregacion,")
    w("mismo conjunto de genes, mismo denominador. Por eso se sigue adelante. Lo")
    w("que **no** esta igualado son tres cosas, y las tres se miden aqui:")
    w("")
    w(f"| # | lo que no esta igualado | `A` | `B` | donde se trata |")
    w("|---|---|---|---|---|")
    w("| 1 | **profundidad** por lado | 200 celulas | 400 celulas | §1-§3, igualando |")
    w("| 2 | **pool de control** de las dos columnas | el mismo `cr` | el mismo `cr` | §4, partiendolo |")
    w(f"| 3 | **n de la mediana** | 33 valores | {npar} valores | §2, bootstrap pareado |")
    w("")
    w("La tercera no es un sesgo, es una diferencia de forma de la distribucion que")
    w("agregar con la mediana no trata igual. Se mide: la mediana y la media de")
    w("`A` se separan mucho mas que las de `B`.")
    w("")
    w("| linea | `A` mediana | `A` media | media - mediana | `B` mediana | `B` media | media - mediana |")
    w("|---|---|---|---|---|---|---|")
    for l in L:
        x = q[l]
        w("| %s | %+.4f | %+.4f | **%+.4f** | %+.4f | %+.4f | %+.4f |"
          % (l, x["A"], x["Am"], x["Am"] - x["A"], x["Bc"], x["Bcm"],
             x["Bcm"] - x["Bc"]))
    w("")
    w("**Y eso ya mueve el conteo antes de tocar la profundidad.** Con la mediana,")
    w(f"`B > A` en {n_b_mayor_vieja} de {n}; con la media, en")
    w("%d de %d. El 5 de 6 no era robusto al agregador."
      % (sum(1 for l in L if q[l]["Bcm"] > q[l]["Am"]), n))
    w("")
    w("### 0.1 Reproduccion exacta antes de cambiar nada")
    w("")
    w("`modelo/src/b_contra_a.py` recalcula `A` y `B` con las mismas semillas")
    w("(real `20260001 + fold`, emision `20266000 + 100*fold`, mitades `777`) y se")
    w("contrasta contra `cand01_refuta.json` valor por valor, no por el agregado:")
    w("")
    w("| linea | max dif absoluta en los 33 `A` | max dif absoluta en los %d `B` | genes expresados igual | celulas de control igual |" % npar)
    w("|---|---|---|---|---|")
    for l in L:
        v = q[l]["d"]["verificacion_candidato_01"]
        w("| %s | **%.1e** | **%.1e** | **%s** | **%s** |"
          % (l, v["A_max_abs_dif"], v["B_completa_max_abs_dif"],
             "si" if v["n_genes_expresados_igual"] else "NO",
             "si" if v["n_control_cells_igual"] else "NO"))
    w("")
    w("Cero exacto en las seis lineas y en las dos cantidades: lo que sigue se")
    w("mide sobre la misma construccion, no sobre una reimplementacion parecida.")
    w("")
    w("### 0.2 Una cosa que se temia y NO pasa")
    w("")
    w("A 200 celulas hay mas genes con cuenta cero que a 400, y")
    w("`cand01_refuta.py:107` manda los `log2` no finitos a `0.0`. Un bloque de")
    w("empates exactos que `argsort(kind='stable')` ordena por indice de gen daria")
    w("**el mismo orden en las dos columnas** y por tanto correlacion positiva")
    w("espuria, mayor en las mitades que en las completas. Se conto:")
    w("")
    w("| linea | no finitos en completa (400) | en mitad_1 (200) | en mitad_2 (200) |")
    w("|---|---|---|---|")
    for l in L:
        z = q[l]["d"]["n_cero_forzado"]
        w("| %s | %d | %d | %d |" % (l, z["completa"], z["mitad_1"], z["mitad_2"]))
    w("")
    tot = {l: q[l]["d"]["n_cero_forzado"] for l in L}
    con = [l for l in L if max(tot[l].values()) > 0]
    if not con:
        w("**Cero en todas.** La mascara `exp` ya basta para que ningun gen quede a")
        w("cuenta cero en un pseudobulk de 200 celulas: el artefacto no existe y no")
        w("hay nada que descontar.")
    else:
        w("**Cero en %d de las %d lineas, y no cero en %s.** Donde no es cero, es"
          % (len(L) - len(con), len(L), ", ".join("`%s`" % l for l in con)))
        for l in con:
            z = tot[l]
            w("mas en las mitades que en la completa (%d y %d contra %d), que es la"
              % (z["mitad_1"], z["mitad_2"], z["completa"]))
            w("direccion que inflaria `A`. El tamano es **%d genes de %d = %.2f %%**"
              % (max(z["mitad_1"], z["mitad_2"]), q[l]["ne"],
                 100.0 * max(z["mitad_1"], z["mitad_2"]) / q[l]["ne"]))
            w("de la columna.")
        w("")
        w("🔴 **Y no se descuenta, por una razon que hay que decir: a profundidad")
        w("igualada el artefacto esta en las DOS cantidades.** `A` y `B_mitades` se")
        w("calculan los dos sobre columnas de 200 celulas, con el mismo numero de")
        w("empates forzados y el mismo orden de indice de gen dentro del bloque de")
        w("empates. Inflar las dos por igual es justamente lo que igualar la")
        w("profundidad consigue; lo que **no** estaba igualado era `B_completa`, que")
        w("se calcula a 400 y tiene menos. Asi que el artefacto empuja `B_completa`")
        w("hacia ABAJO respecto de `A`, es decir en contra de `B > A`, y no puede")
        w("ser la causa del resultado de aqui. Queda contado, no corregido.")
    w("")
    w("---")
    w("")

    # ---------- 1 ----------
    w("## 1. Las tres cantidades, y el cambio de signo")
    w("")
    w("`A` no se corrige hacia arriba: ni Spearman-Brown ni ninguna otra")
    w("extrapolacion. `anchor.py:3-7` dice por que -- \"Deliberately NOT")
    w("Spearman-Brown corrected and not depth-corrected in any other way -- the")
    w("correction `ceiling.py` applies is a bounded-reliability extrapolation that")
    w("does not hold for five of the six scored `vcc2026` members\". Se iguala la")
    w("profundidad bajando `B`, no subiendo `A`:")
    w("")
    w("| | construccion | celulas por lado |")
    w("|---|---|---|")
    w("| `A_mitades` | `mitad_1(p)` contra `mitad_2(p)` | 200 vs 200 |")
    w("| `B_completa` | `poblacion(p)` contra `poblacion(q)` | 400 vs 400 |")
    w("| `B_mitades` | `mitad_1(p)` contra `mitad_1(q)` | **200 vs 200** |")
    w("")
    w("Pares de un blanco consigo mismo se excluyen de `B_mitades`: a profundidad")
    w("igualada `mitad_1(p)` contra `mitad_1(p)` es el mismo vector y daria 1.0")
    w("exacto. Son %d pares, los mismos %d del punto." % (npar, npar))
    w("")
    w("| linea | `A` | `B_completa` | `B_mitades` | `B_comp - A` (lo que se reporto) | **`B_mit - A`** (profundidad igualada) | signo |")
    w("|---|---|---|---|---|---|---|")
    for l in L:
        x = q[l]
        w("| %s | %+.4f | %+.4f | %+.4f | %+.4f | **%+.4f** | **%s** |"
          % (l, x["A"], x["Bc"], x["Bm"], x["d_vieja"], x["d_nueva"],
             "B>A" if x["d_nueva"] > 0 else "A>B"))
    w("")
    w(f"**`B > A` en {n_b_mayor_vieja} de {n} a profundidad desigual; en")
    w(f"{n_b_mayor_nueva} de {n} a profundidad igualada.** No es que la diferencia se")
    w("encoja: cambia de signo en todas.")
    w("")
    w("---")
    w("")

    # ---------- 2 ----------
    w("## 2. La barra: bootstrap PAREADO sobre blancos")
    w("")
    w(f"{nrep} replicas. En cada una: se remuestrean los {nb} blancos **con")
    w("reemplazo**, se vuelve a partir cada blanco en mitades con semilla propia, y")
    w("se recalculan `A` y `B` **sobre los mismos blancos remuestreados**. La")
    w("cantidad que se guarda es `B - A` por replica, **pareada**, asi que la barra")
    w("es la de la diferencia y no la resta de dos barras sueltas.")
    w("")
    w("**Las %d replicas pedidas se corrieron enteras: no hubo que bajarlas.**" % nrep)
    w("Coste medido, %s por linea (seis lineas, una pasada), porque el"
      % ("/".join("%.0f s" % q[l]["d"]["segundos"] for l in L[:1] )
         if False else "%.0f-%.0f s" % (
             min(q[l]["d"]["segundos"] for l in L),
             max(q[l]["d"]["segundos"] for l in L))))
    w("lado real se lee una vez y lo unico que se repite por replica es la")
    w("particion en mitades y el rankeo.")
    w("")
    w("**Sobre BLANCOS y no sobre los %d pares.** Dos pares que comparten un blanco" % npar)
    w("no son independientes; remuestrear pares daria una barra falsamente")
    w("estrecha. Semillas: `20267000 + 100*fold + replica`, familia nueva.")
    w("")
    w("| linea | `A` punto | `A` media boot | `B_mit` punto | `B_mit` media boot | `B_mit - A` punto | media boot | `sd` boot | **IC 95 %** | excluye 0 |")
    w("|---|---|---|---|---|---|---|---|---|---|")
    for l in L:
        x = q[l]
        w("| %s | %+.4f | %+.4f | %+.4f | %+.4f | **%+.4f** | %+.4f | %.4f | **[%+.4f, %+.4f]** | **%s** |"
          % (l, x["A"], x["A_boot"], x["Bm"], x["Bm_boot"], x["d_nueva"],
             x["d_boot"], x["sd_d"], x["lo"], x["hi"],
             "SI" if x["excluye"] else "no"))
    w("")
    w(f"**El IC excluye 0 en {n_excluye} de {n} lineas**, y en las {n_neg} que lo excluyen")
    w("lo excluye **por debajo de 0** (`A > B`).")
    w("")
    w("### 2.1 La expectativa por azar del conteo, que es lo que exige §7")
    w("")
    w("Un conteo sobre seis lineas no se lee sin su nulo. Hay dos nulos distintos y")
    w("se dan los dos, porque responden a preguntas distintas:")
    w("")
    w("| pregunta | nulo | esperado | observado | `p` exacta |")
    w("|---|---|---|---|---|")
    p1 = binom_cola(n_b_mayor_vieja, n, 0.5) if n_b_mayor_vieja else 1.0
    w("| ¿cuantas lineas con `B_comp > A`? (candidato_01) | moneda, `p = 0.5` | %.1f | %d de %d | **%.4f** = %d/%d |"
      % (n * 0.5, n_b_mayor_vieja, n, p1, round(p1 * 2 ** n), 2 ** n))
    p2 = binom_cola(n - n_b_mayor_nueva, n, 0.5)
    w("| ¿cuantas lineas con `A > B_mit`? (aqui) | moneda, `p = 0.5` | %.1f | %d de %d | **%.4f** = %d/%d |"
      % (n * 0.5, n - n_b_mayor_nueva, n, p2, round(p2 * 2 ** n), 2 ** n))
    p3 = binom_cola(n_excluye, n, 0.05)
    w("| ¿cuantas lineas con IC 95 %% que excluye 0? | `B - A = 0` en todas, `p = 0.05` | **%.2f** | %d de %d | **%.2e** |"
      % (n * 0.05, n_excluye, n, p3))
    p4 = binom_cola(n_neg, n, 0.025)
    w("| ¿cuantas con IC que excluye 0 **por debajo**? | `p = 0.025` (una cola) | **%.2f** | %d de %d | **%.2e** |"
      % (n * 0.025, n_neg, n, p4))
    w("")
    w("")
    w("🔴 **Las dos ultimas `p` suponen las seis lineas independientes, y hay que")
    w("decir en que sentido lo son y en que sentido no.** Lo son en el ruido: las")
    w("seis lineas son poblaciones de celulas **disjuntas**, asi que bajo el nulo")
    w("`B - A = 0` el error de muestreo de una linea no sabe nada del de otra y el")
    w("binomial vale. **No** lo son en el mecanismo: las seis comparten el mismo")
    w("panel de %d blancos y la misma causa -- una diferencia de profundidad de" % nb)
    w("construccion, que esta en las seis por igual. Asi que un 6 de 6 **no son")
    w("seis confirmaciones independientes de un efecto biologico**: es un")
    w("mecanismo medido seis veces. Para refutar una hipotesis eso es lo que hace")
    w("falta y sobra; para afirmar una, no bastaria.")
    w("")
    w("🔴 **Y hay que decir que el 5 de 6 original no era el problema principal.**")
    w("Su `p = 0.109` no excluye la moneda, pero el conteo tampoco era la medicion")
    w("debil: la medicion debil era que `A` y `B` no median a la misma")
    w("profundidad. Un conteo de 6 de 6 sobre cantidades mal emparejadas habria")
    w("dado `p = 0.016` y seguiria estando mal. **La barra y el conteo no arreglan")
    w("un estimando torcido; esto se cierra por el PASO 0, no por la `p`.**")
    w("")
    w("### 2.2 Con la media en vez de la mediana, por si el agregador manda")
    w("")
    w("| linea | `B_mit - A` con mediana | IC 95 % | `B_mit - A` con media | IC 95 % | coinciden en signo |")
    w("|---|---|---|---|---|---|")
    for l in L:
        x = q[l]
        dm = x["Bmm"] - x["Am"]
        w("| %s | %+.4f | [%+.4f, %+.4f] | %+.4f | [%+.4f, %+.4f] | **%s** |"
          % (l, x["d_nueva"], x["lo"], x["hi"], dm, x["lom"], x["him"],
             "si" if (dm < 0) == (x["d_nueva"] < 0) else "NO"))
    w("")
    w("---")
    w("")

    # ---------- 3 ----------
    w("## 3. Cuanto de `B - A` era profundidad")
    w("")
    w("`B_completa - B_mitades` es la atenuacion por bajar de 400 a 200 celulas por")
    w("lado, medida sobre la MISMA cantidad (los mismos %d pares, los mismos" % npar)
    w("genes, el mismo denominador). Es lo que `A` estaba pagando y `B` no:")
    w("")
    w("| linea | `B_comp - A` reportado | **profundidad** = `B_comp - B_mit` | resto = `B_mit - A` | profundidad / \\|reportado\\| |")
    w("|---|---|---|---|---|")
    for l in L:
        x = q[l]
        r = (x["prof"] / abs(x["d_vieja"])) if x["d_vieja"] != 0 else float("nan")
        w("| %s | %+.4f | **%+.4f** | %+.4f | %s |"
          % (l, x["d_vieja"], x["prof"], x["d_nueva"],
             ("%.1f x" % r) if np.isfinite(r) else "n/d"))
    w("")
    pr = [q[l]["prof"] for l in L]
    dv = [abs(q[l]["d_vieja"]) for l in L]
    w("La atenuacion por profundidad va de **%+.4f a %+.4f** y el `B - A` que se" % (min(pr), max(pr)))
    w("reporto iba de %+.4f a %+.4f en valor absoluto. **La profundidad es de %.1f a"
      % (min(dv), max(dv), min(p / d for p, d in zip(pr, dv))))
    w("%.1f veces el efecto que se reporto**, asi que no explica parte de la"
      % max(p / d for p, d in zip(pr, dv)))
    w("diferencia: la explica entera y sobra para cambiarle el signo.")
    w("")
    w("---")
    w("")
    O.extend(seccion_anchor(anc, L, nb))
    O.extend(seccion_descomp(des, L))
    O.extend(seccion_filas(b, L, nb, npar))
    O.extend(seccion_cierre(L, q, n_excluye, n))

    with open(OUT, "w") as fh:
        fh.write("\n".join(O) + "\n")
    sys.stderr.write("escrito %s (%d lineas)\n" % (OUT, len(O)))


def seccion_anchor(anc, L, nb):
    O = []; w = O.append
    w("## 4. El segundo confundidor, que no estaba en el encargo: el control compartido")
    w("")
    w("Apareció al escribir el PASO 0 y no se puede dejar fuera, porque el propio")
    w("marcador lo declara inadmisible. Todo delta es `log2(pb) - log2(cr)`, y en")
    w("`cand01_refuta.py` **las dos columnas de cada correlacion llevan el mismo")
    w("`cr`** -- las dos mitades en `A`, los dos blancos en `B`.")
    w("")
    w("`cell-eval2/src/cell_eval2/anchor.py:9-18`, literal:")
    w("")
    w("> under \"real\" the pred side's DE reference comes from the real side, so")
    w("> scoring half_b against half_a computes BOTH halves' log2FCs against")
    w("> half_a's control. That shares the control's sampling noise between the two")
    w("> quantities whose agreement is being measured and **biases the anchor")
    w("> upward** [...] Measured for de_lfc_nmae on a real panel: shared controls")
    w("> read 0.5-2.3% optimistic [...] `ceiling.py` measured the same mistake")
    w("> turning a **0.54 split-half reliability into 0.74**.")
    w("")
    w("Y `anchor.py:768-779` no lo avisa, lo **rehusa**: \"a replicate anchor is")
    w("only valid under per-half controls ('pred')\". El estimando no es")
    w("negociable (su `spec 5.5`).")
    if anc is None:
        w("")
        w("**No medido en este encargo** -- el JSON `b_contra_a_anchor.json` no")
        w("existe. Queda declarado como hueco abierto.")
        w("")
        w("---")
        w("")
        return O
    w("")
    w("Asi que la version conforme al ancla, con el pool de control partido en dos")
    w("mitades disjuntas `c1` y `c2`, y las dos cantidades con la MISMA forma")
    w("`mitad_1` contra `mitad_2` y 200 celulas por lado:")
    w("")
    w("```")
    w("A_anc = corr( log2(mitad_1(p)/c1), log2(mitad_2(p)/c2) )")
    w("B_anc = corr( log2(mitad_1(p)/c1), log2(mitad_2(q)/c2) )     p != q")
    w("```")
    w("")
    w("Misma barra: bootstrap pareado sobre blancos, misma familia de semillas, y")
    w("el control se reparte de nuevo en cada replica.")
    w("")
    w("| linea | celulas de control por mitad | `A_anc` | `B_anc` | `B_anc - A_anc` | media boot | **IC 95 %** | excluye 0 |")
    w("|---|---|---|---|---|---|---|---|")
    nex = 0
    for l in L:
        if l not in anc["lineas"]:
            continue
        d = anc["lineas"][l]
        A = float(np.median(d["A_anclado"])); Bv = float(np.median(d["B_anclado"]))
        lo, hi = ic(d["bootstrap"]["d_med"])
        ex = lo > 0 or hi < 0
        nex += int(ex)
        w("| %s | %d / %d | %+.4f | %+.4f | **%+.4f** | %+.4f | **[%+.4f, %+.4f]** | **%s** |"
          % (l, d["n_control_por_mitad"][0], d["n_control_por_mitad"][1], A, Bv,
             Bv - A, float(np.mean(d["bootstrap"]["d_med"])), lo, hi,
             "SI" if ex else "no"))
    w("")
    w("**Y el efecto del control compartido es mucho mayor que el 0.5-2.3 % que")
    w("`anchor.py` midio sobre `de_lfc_nmae`**, porque aqui la cantidad es un")
    w("Spearman sobre ~7000 genes y no un NMAE sobre genes DE:")
    w("")
    w("| linea | `A` control comun | `A_anc` control por mitad | `A` perdida | `B_mit` control comun | `B_anc` control por mitad | `B` perdida |")
    w("|---|---|---|---|---|---|---|")
    bb = json.load(open(B))["lineas"]
    for l in L:
        if l not in anc["lineas"]:
            continue
        d = anc["lineas"][l]
        Ac = float(np.median(bb[l]["A_mitades"]))
        Bcm = float(np.median(bb[l]["B_mitades"]))
        Aa = float(np.median(d["A_anclado"])); Ba = float(np.median(d["B_anclado"]))
        w("| %s | %+.4f | %+.4f | **%+.4f** | %+.4f | %+.4f | **%+.4f** |"
          % (l, Ac, Aa, Aa - Ac, Bcm, Ba, Ba - Bcm))
    w("")
    w(f"Conforme al ancla, el IC excluye 0 en **{nex} de {len(L)}** lineas.")
    w("")
    w("🔴 **Pero la tabla de arriba dice algo mas preciso que \"empuja en la misma")
    w("direccion\", y conviene no cobrarlo de mas: la perdida de `A` y la de `B`")
    w("son casi la misma en valor absoluto.** El ruido del control entraba")
    w("identico en las dos columnas de las dos cantidades, asi que al partirlo las")
    w("dos bajan a la par y la DIFERENCIA se queda donde estaba:")
    w("")
    w("| linea | `B_mit - A` con control comun | `B_anc - A_anc` con control por mitad | se mueve |")
    w("|---|---|---|---|")
    for l in L:
        if l not in anc["lineas"]:
            continue
        d = anc["lineas"][l]
        d1 = float(np.median(bb[l]["B_mitades"])) - float(np.median(bb[l]["A_mitades"]))
        d2 = float(np.median(d["B_anclado"])) - float(np.median(d["A_anclado"]))
        w("| %s | %+.4f | %+.4f | %+.4f |" % (l, d1, d2, d2 - d1))
    w("")
    w("**Asi que el control compartido NO es un confundidor de `B - A`: se cancela")
    w("en la resta.** Es un defecto del **nivel**, no de la comparacion, y tiene")
    w("dos consecuencias distintas que no hay que mezclar:")
    w("")
    w("1. **Refuerza la refutacion de §1-§2**, porque `B - A` sale igual por los")
    w("   dos caminos -- el resultado no cuelga de haber compartido el control.")
    w("2. **Rompe otra afirmacion de `candidato_01.md`, la del techo.** Alli se")
    w("   leyo \"`A` corre de +0.0820 a +0.2012 en las seis lineas: el test tiene")
    w("   potencia\". Conforme al ancla `A` corre de **%+.4f a %+.4f**: entre"
      % (min(float(np.median(anc["lineas"][l]["A_anclado"])) for l in L
             if l in anc["lineas"]),
         max(float(np.median(anc["lineas"][l]["A_anclado"])) for l in L
             if l in anc["lineas"])))
    w("   un tercio y dos tercios de lo que se reporto. El test sigue teniendo")
    w("   potencia -- `A` no es cero en ninguna linea -- pero el techo es")
    w("   bastante mas bajo, y eso **si** entra al diseno del candidato 02. No se")
    w("   corrige hacia arriba: se reporta medido.")
    w("")
    w("---")
    w("")
    return O


def seccion_descomp(des, L):
    O = []; w = O.append
    w("## 5. Lo que reemplaza a la hipotesis: la descomposicion, medida sin modelo")
    w("")
    w("El encargo pedia esta descomposicion **solo si la hipotesis sobrevivia**. No")
    w("sobrevive, asi que se mide por otra razon y se declara: `candidato_01.md`")
    w("uso `B > A` para afirmar que el `b` \"esta capturando la componente")
    w("dominante de la varianza\" y que lo especifico del blanco es \"la")
    w("minoritaria incluso en el lado real\". Esa afirmacion se queda sin apoyo")
    w("cuando `B > A` cae, y es exactamente la que decide el diseno del candidato")
    w("02. Asi que se mide directamente en vez de heredarla.")
    if des is None:
        w("")
        w("**No medido** -- falta `b_contra_a_descomp.json`. Hueco abierto.")
        w("")
        w("---")
        w("")
        return O
    w("")
    w("**Sin suponer ningun modelo.** Lo unico que se usa es que las dos mitades de")
    w("un blanco tienen ruido independiente -- son celulas disjuntas -- y media")
    w("cero. Con `h1_i`, `h2_i` los deltas de las dos mitades del blanco `i`,")
    w("centrados por columna sobre los genes, y `m1`, `m2` sus promedios sobre los")
    w("33 blancos:")
    w("")
    w("```")
    w("V_tot = (1/n) Σ_i <h1_i, h2_i> / G          senal reproducible total")
    w("V_sh  = <m1, m2> / G                        COMPARTIDA entre blancos")
    w("V_sp  = (1/n) Σ_i <h1_i-m1, h2_i-m2> / G    ESPECIFICA del blanco")
    w("```")
    w("")
    w("`V_tot = V_sh + V_sp` es una **identidad algebraica exacta**, no un ajuste:")
    w("los terminos cruzados se anulan porque `Σ_i (h1_i - m1) = 0`. Y el producto")
    w("**cruzado de mitades** es lo que quita el ruido sin modelarlo: su esperanza")
    w("es el producto de las senales verdaderas, porque el ruido de las dos mitades")
    w("es independiente. La varianza de **una** mitad si llevaria el ruido dentro y")
    w("habria que suponer un modelo para sacarlo. **Se puede medir sin modelo, y se")
    w("mide.**")
    w("")
    w("Y por lo mismo de §4, el control va **partido**: con un `cr` comun el")
    w("termino `-log2(cr)` es identico en todas las columnas, cae entero dentro de")
    w("`m1` y `m2` con la misma realizacion de ruido, y `<m1, m2>` contaria el")
    w("ruido del control como senal compartida. Se dan las dos para que se vea")
    w("cuanto era eso:")
    w("")
    w("| linea | | `V_tot` | `V_compartida` | `V_especifica` | **frac compartida** | **frac especifica** | residuo de la identidad |")
    w("|---|---|---|---|---|---|---|---|")
    for l in L:
        if l not in des["lineas"]:
            continue
        d = des["lineas"][l]
        for nom, k in (("control comun", "control_comun"),
                       ("**control por mitad**", "control_independiente")):
            v = d[k]
            w("| %s | %s | %.4f | %.4f | %.4f | **%.1f %%** | **%.1f %%** | %.1e |"
              % (l if k == "control_comun" else "", nom, v["V_tot"],
                 v["V_compartida"], v["V_especifica"],
                 100 * v["frac_compartida"], 100 * v["frac_especifica"],
                 abs(v["identidad_resid"])))
    w("")
    fc = [des["lineas"][l]["control_independiente"]["frac_compartida"]
          for l in L if l in des["lineas"]]
    fe = [des["lineas"][l]["control_independiente"]["frac_especifica"]
          for l in L if l in des["lineas"]]
    fcc = [des["lineas"][l]["control_comun"]["frac_compartida"]
           for l in L if l in des["lineas"]]
    w("**Con el control compartido la fraccion compartida sale de %.0f a %.0f %%**, que es"
      % (100 * min(fcc), 100 * max(fcc)))
    w("la lectura que respaldaba \"la componente dominante\". **Con el control")
    w("partido va de %.0f a %.0f %% y lo especifico de %.0f a %.0f %%.**"
      % (100 * min(fc), 100 * max(fc), 100 * min(fe), 100 * max(fe)))
    w("")
    if max(fe) > 0.5:
        nmaj = sum(1 for x in fe if x > 0.5)
        nmajc = sum(1 for x in [des["lineas"][l]["control_comun"]["frac_especifica"]
                                for l in L if l in des["lineas"]] if x > 0.5)
        w("🔴 **Da la vuelta a la premisa del candidato 02.** Con el control por")
        w("mitad lo especifico del blanco es la componente **mayoritaria en %d de las"
          % nmaj)
        w("%d lineas** (con el control comun lo era en %d de %d, que es por lo que la"
          % (len(fe), nmajc, len(fe)))
        w("lectura salia al contrario). Un candidato no tiene que acertar \"la")
        w("componente pequena\": la parte especifica del lado real es grande y")
        w("reproducible -- `V_especifica` no lleva ruido dentro, es un producto")
        w("cruzado de mitades -- y lo que la hacia parecer pequena era el ruido de")
        w("muestreo del control contado como senal compartida.")
        w("")
        w("Y notese que `V_especifica` es **la misma columna en las dos filas de cada")
        w("linea**: el control comun no la toca. Todo el movimiento de la fraccion")
        w("sale de `V_compartida`, que es donde cae el termino `-log2(cr)`. Eso es la")
        w("prueba de que el efecto es el que se dice y no una reescala global.")
    else:
        w("Lo compartido sigue siendo mayoritario incluso con el control partido, asi")
        w("que esa parte de la lectura de `candidato_01.md` se sostiene -- pero por")
        w("esta medicion, no por `B > A`.")
    w("")
    w("---")
    w("")
    return O


def seccion_filas(b, L, nb, npar):
    O = []; w = O.append
    panel = b["meta"]["panel"]
    lsel = L[0]
    d = b["lineas"][lsel]
    A = np.array(d["A_mitades"])
    Bm = np.array(d["B_mitades"]); Bc = np.array(d["B_completa"])
    # matriz de pares -> media por blanco
    M_m = np.zeros((nb, nb)); M_c = np.zeros((nb, nb))
    t = 0
    for i in range(nb):
        for k in range(i + 1, nb):
            M_m[i, k] = M_m[k, i] = Bm[t]
            M_c[i, k] = M_c[k, i] = Bc[t]
            t += 1
    off = ~np.eye(nb, dtype=bool)
    w("## 6. Las filas trabajadas")
    w("")
    w(f"### 6.1 Los {nb} blancos de {lsel}, uno por fila")
    w("")
    w("La columna `A_mitades` es la lista cuya **mediana** es el `A` de §1.")
    w("`B_mit medio` es el promedio de los 32 pares en que ese blanco entra a")
    w("profundidad igualada, y `B_comp medio` los mismos 32 a profundidad completa.")
    w("")
    w("⚠️ **La ultima columna NO es el estimador.** El numero de linea es")
    w("`mediana(33 valores de A) - mediana(%d pares de B)`, no la mediana de esta" % npar)
    w("columna: `B` se agrega sobre pares, no sobre blancos, y un blanco entra en")
    w("32 pares. La columna esta aqui como **intermedio y diagnostico de")
    w("concentracion** -- para ver si el resultado de linea lo arrastran dos o tres")
    w("filas o esta repartido -- no como la cantidad que se reporta.")
    w("")
    w("| # | blanco | celulas | `A_mitades` | `B_mit` medio (32 pares) | `B_comp` medio (32) | profundidad por blanco | **`B_mit` - `A`** | `A > B_mit` |")
    w("|---|---|---|---|---|---|---|---|---|")
    orden = np.argsort(-A)
    for r, i in enumerate(orden, 1):
        bm = M_m[i][off[i]].mean(); bc = M_c[i][off[i]].mean()
        w("| %d | `%s` | %d | %+.4f | %+.4f | %+.4f | %+.4f | **%+.4f** | %s |"
          % (r, panel[i], d["n_celulas_por_blanco"][i], A[i], bm, bc, bc - bm,
             bm - A[i], "**si**" if A[i] > bm else "no"))
    nmay = int(sum(1 for i in range(nb) if A[i] > M_m[i][off[i]].mean()))
    w("")
    w("**`A > B_mit` en %d de los %d blancos** de %s, y los tres que no son los tres"
      % (nmay, nb, lsel))
    w("de `A` mas baja de la tabla. El resultado de linea no lo arrastran dos o")
    w("tres filas: esta repartido. (Mediana de")
    w("la columna `A_mitades` = %+.4f, que es el `A` de §1; mediana de los %d pares"
      % (float(np.median(A)), npar))
    w("= %+.4f, que es el `B_mitades`.)" % float(np.median(Bm)))
    w("")
    w("### 6.2 Veinte replicas del bootstrap, una por fila")
    w("")
    w(f"De las {b['meta']['n_replicas']} de {lsel}. `blancos distintos` es cuantos de los {nb}")
    w("sobreviven al remuestreo con reemplazo (el resto son repeticiones);")
    w("`pares` es cuantos pares de blancos distintos quedan tras excluir los de un")
    w("blanco consigo mismo. El IC de §2 es el percentil 2.5/97.5 de la ultima")
    w("columna sobre las %d replicas." % b["meta"]["n_replicas"])
    w("")
    w("| replica | semilla | blancos distintos | pares | `A` mediana | `B_mit` mediana | **`B_mit - A`** |")
    w("|---|---|---|---|---|---|---|")
    bo = d["bootstrap"]
    base = b["meta"]["semilla_familia"] + 100 * d["fold_idx"]
    for r in range(20):
        w("| %d | %d | %d | %d | %+.4f | %+.4f | **%+.4f** |"
          % (r, base + r, bo["n_distintos"][r], bo["n_pares"][r],
             bo["A_med"][r], bo["Bm_med"][r], bo["d_med"][r]))
    dd = np.array(bo["d_med"])
    w("")
    w("Las 20 primeras replicas dan `B_mit - A` negativo en %d de 20; las %d"
      % (int((np.array(bo["d_med"][:20]) < 0).sum()), len(dd)))
    w("completas, en %d de %d. El signo no depende de que replicas se miren."
      % (int((dd < 0).sum()), len(dd)))
    w("")
    w("---")
    w("")
    return O


def seccion_cierre(L, q, n_excluye, n):
    O = []; w = O.append
    w("## 7. Que se hizo, que no, y que queda abierto")
    w("")
    w("**No se hizo** -- ninguna de las prohibiciones del encargo se roza:")
    w("")
    w("| prohibicion | estado |")
    w("|---|---|")
    w("| corregir `A` por profundidad (Spearman-Brown o cualquier otra) | **no**. Se bajo `B` a 200 celulas. `anchor.py:3-7` citado en §1 |")
    w("| construir un candidato, tocar un arnes, tocar un bundle | **no**. `build_fold_adatas` se llama como parametro, sin modificar; no se escribe ningun bundle |")
    w("| leer el lado real para CONSTRUIR algo | **no**. Todo esto es diagnostico del problema. Ninguna cifra de aqui entra a un modelo |")
    w("| `~/vcc/benchmark/`, `/datos/benchmark_VEDADO_techo_metricas/` | **no leidos**. Las unicas entradas son `IFNB`, `suelo_panel.json` y `cand01_refuta.json` |")
    w("| envios, `vcc sample`, commits | **ninguno** |")
    w("")
    w("**Lo que entro y lo que salio:**")
    w("")
    w("| | |")
    w("|---|---|")
    w("| entradas | `/datos/ref/convertido/Seurat_object_IFNB_Perturb_seq.h5ad`, `modelo/out/suelo_panel.json`, `modelo/out/cand01_refuta.json` (solo para verificar) |")
    w("| codigo | `modelo/src/b_contra_a.py`, `b_contra_a_anchor.py`, `b_contra_a_descomp.py`, `b_contra_a_reporte.py` |")
    w("| salidas | `b_contra_a.json`, `b_contra_a_anchor.json`, `b_contra_a_descomp.json`, este `.md` |")
    w("| semillas | real `20260001+fold`, emision `20266000+100*fold`, mitades del punto `777`, **bootstrap `20267000+100*fold+replica`**, control partido `20267999+fold` |")
    w("")
    w("**Lo que queda abierto, declarado el dia que aparece:**")
    w("")
    w("1. **La mascara de genes de `cand01_refuta.py` depende de la emision del")
    w("   candidato.** `exp = (cr >= 1) & (cp >= 1)` y `cp` es el pseudobulk del")
    w("   control **predicho**, que es un remuestreo de 400 celulas del control real.")
    w("   Para un diagnostico de la verdad eso no deberia entrar. Se midio la")
    w("   version de solo-verdad (`cr >= 1`) como verificacion y el resultado no")
    w("   cambia:")
    w("")
    w("| linea | genes con `exp` | genes con solo `cr >= 1` | `B_mit - A` con `exp` | con solo-verdad |")
    w("|---|---|---|---|---|")
    for l in L:
        d = q[l]["d"]; sv = d["sens_mascara_solo_verdad"]
        a2 = float(np.median(sv["A"])); b2 = float(np.median(sv["B_mitades"]))
        w("| %s | %d | %d | %+.4f | %+.4f |"
          % (l, d["n_genes_expresados"], d["n_genes_solo_verdad"],
             q[l]["d_nueva"], b2 - a2))
    w("")
    w("   No afecta al resultado; afecta a que `cand01_refuta.py` no es reutilizable")
    w("   tal cual como medidor de la verdad. Cerrarlo cuesta una linea.")
    w("")
    w("2. **El empate de agregador no se resolvio, se midio.** La mediana sobre 33")
    w("   valores y sobre 528 no son el mismo estimador aunque sea la misma funcion.")
    w("   Aqui la mediana y la media coinciden en signo en las seis lineas (§2.2),")
    w("   asi que no decide nada -- pero en `candidato_01.md` si decidia (5 de 6")
    w("   contra 4 de 6, §0). Afecta a la **interpretacion** de ese conteo, no a")
    w("   este resultado.")
    w("")
    w("3. **No se midio si `A` por si sola basta como techo.** `A` a 200 celulas es")
    w("   lo que es; lo que un candidato puede alcanzar a 400 no se extrapola desde")
    w("   aqui, por la misma razon por la que no se corrige con Spearman-Brown.")
    w("   Afecta a la interpretacion del techo, no a la comparacion `B` contra `A`.")
    w("")
    w("---")
    w("")
    w("## 8. Para el candidato 02, en una frase")
    w("")
    w("Tres cosas, y las tres son medidas, no lecturas:")
    w("")
    w("1. **`B > A` no existe.** Era la diferencia de profundidad entre 200 y 400")
    w("   celulas por lado. A profundidad igualada `A > B` en las seis lineas, con")
    w("   el IC de la diferencia excluyendo 0 en las seis. El control compartido")
    w("   **no** entra aqui: se cancela en la resta.")
    w("2. **El techo `A` es mas bajo de lo que se reporto**, porque el `A` de")
    w("   `candidato_01.md` comparte el control entre las dos mitades, que es")
    w("   exactamente lo que `anchor.py` rehusa. Sigue siendo distinto de cero en")
    w("   las seis lineas, asi que el test conserva potencia.")
    w("3. **Y lo especifico del blanco no es la componente minoritaria.** Medido")
    w("   sin modelo y con el control partido, es la mayoritaria. Asi que la frase")
    w("   \"un candidato que quiera superar al `b` tiene que acertar la componente")
    w("   pequena\" se queda sin su apoyo: la componente especifica es la grande.")
    w("")
    w("**Lo que esto le quita al candidato 02 es una excusa, no una via.** El")
    w("candidato 01 saco `C` ~ 0 y concordancia de signo 0.5; eso sigue en pie y")
    w("no lo explica que lo especifico fuera minoritario, porque no lo es.")
    return O


if __name__ == "__main__":
    main()
