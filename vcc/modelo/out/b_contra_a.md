# `B` contra `A` a profundidad igualada: la hipotesis queda REFUTADA

`candidato_01.md` §1 midio `B > A` en cinco de las seis lineas y lo dejo
como **hipotesis**, no como resultado, por dos razones: 5 de 6 bajo una
moneda tiene `p = 7/64 = 0.109`, y ninguna de las dos cantidades llevaba
dispersion. Este encargo la cierra.

🔴 **Se cierra en contra.** A profundidad igualada `B - A` **cambia de
signo en las 6 lineas**: `A > B`, no `B > A`. El IC 95 % de `B - A`
excluye el 0 en **6 de 6** lineas, y en las 6 lo excluye por
debajo. La diferencia original era **profundidad**, no biologia: `A` se
medía entre mitades de 200 celulas y `B` entre poblaciones completas de
400, y la atenuacion por profundidad de `A` era **mayor** que el `B - A`
que se reporto.

Y al declarar la construccion aparecio un **segundo** defecto que el encargo
no pedia: las dos columnas de cada correlacion se dividen por el **mismo**
pseudobulk de control, lo que `cell-eval2/src/cell_eval2/anchor.py:9-18`
declara sesgo **hacia arriba** y `anchor.py:768-779` **rehusa**
directamente. Se midio, y hay que decir exactamente que hace y que no
hace (§4): **no** mueve `B - A` -- infla `A` y `B` casi por el mismo
absoluto y se cancela en la diferencia, lo que vuelve esta refutacion mas
robusta, no menos. Lo que si rompe es el **nivel**: el `A` que
`candidato_01.md` presento como "el techo alcanzable" lee entre 0.05 y 0.12
de mas.

**Para el candidato 02 esto cambia el diseno**: ya no hay apoyo medido para
"el `b` captura la componente dominante de la varianza y lo especifico es
minoritario". §5 mide esa fraccion directamente, que es lo que reemplaza a
la hipotesis.

---

## 0. PASO 0 - la construccion, escrita ANTES de anadir una sola barra

La barra depende de la estructura, asi que primero la estructura. Leida de
`modelo/src/cand01_refuta.py`, no reconstruida de memoria:

```
por cada blanco g del panel (33 blancos, 400 celulas reales cada uno):
  cr          = pseudobulk(TODAS las celulas de control de la linea)   # CP50k
  R[:, g]     = log2( pseudobulk(400 celulas de g) / cr )              # :99
  perm        = rng(777).permutation(indices de g)                     # :96
  H1[:, g]    = log2( pseudobulk(perm[:200]) / cr )                    # :100
  H2[:, g]    = log2( pseudobulk(perm[200:]) / cr )                    # :101

mascara de genes (la MISMA para las dos cantidades):                   # :89
  exp = (cr >= 1) & (cp >= 1)      cp = pseudobulk del control PREDICHO

rangos por columna via argsort estable, luego centrar y normalizar      # :48,:59
  nX = normaliza(rankea(X))        -> Spearman = producto escalar de columnas

A = [ nH1[:,i] . nH2[:,i]  para i en 33 ]                              # :111
B = [ nR[:,i]  . nR[:,k]   para i<k en 33 ] -> 528 pares               # :116-118

y lo que entra a la tabla de candidato_01.md es la MEDIANA de cada lista
  A = median(33 valores)   B = median(528 valores)   # cand01_reporte.py:125-126
```

**Respuesta a las dos preguntas del PASO 0.** Spearman **por par y luego
mediana**, no un Spearman agrupado -- en las dos cantidades, con la misma
funcion de agregacion. Y **sobre los mismos genes**: la mascara `exp` se
calcula una vez por linea y la comparten `A` y `B`.

**Veredicto de comparabilidad: SI son estructuralmente comparables en
forma** -- misma correlacion, mismo estimador por par, misma agregacion,
mismo conjunto de genes, mismo denominador. Por eso se sigue adelante. Lo
que **no** esta igualado son tres cosas, y las tres se miden aqui:

| # | lo que no esta igualado | `A` | `B` | donde se trata |
|---|---|---|---|---|
| 1 | **profundidad** por lado | 200 celulas | 400 celulas | §1-§3, igualando |
| 2 | **pool de control** de las dos columnas | el mismo `cr` | el mismo `cr` | §4, partiendolo |
| 3 | **n de la mediana** | 33 valores | 528 valores | §2, bootstrap pareado |

La tercera no es un sesgo, es una diferencia de forma de la distribucion que
agregar con la mediana no trata igual. Se mide: la mediana y la media de
`A` se separan mucho mas que las de `B`.

| linea | `A` mediana | `A` media | media - mediana | `B` mediana | `B` media | media - mediana |
|---|---|---|---|---|---|---|
| A549 | +0.1575 | +0.1951 | **+0.0376** | +0.2004 | +0.1989 | -0.0014 |
| BXPC3 | +0.1544 | +0.1963 | **+0.0419** | +0.1198 | +0.1290 | +0.0092 |
| HAP1 | +0.0820 | +0.0954 | **+0.0134** | +0.1053 | +0.1076 | +0.0023 |
| HT29 | +0.1020 | +0.1327 | **+0.0308** | +0.1298 | +0.1341 | +0.0043 |
| K562 | +0.1609 | +0.1786 | **+0.0177** | +0.2184 | +0.2172 | -0.0012 |
| MCF7 | +0.2012 | +0.2200 | **+0.0187** | +0.2133 | +0.2106 | -0.0026 |

**Y eso ya mueve el conteo antes de tocar la profundidad.** Con la mediana,
`B > A` en 5 de 6; con la media, en
4 de 6. El 5 de 6 no era robusto al agregador.

### 0.1 Reproduccion exacta antes de cambiar nada

`modelo/src/b_contra_a.py` recalcula `A` y `B` con las mismas semillas
(real `20260001 + fold`, emision `20266000 + 100*fold`, mitades `777`) y se
contrasta contra `cand01_refuta.json` valor por valor, no por el agregado:

| linea | max dif absoluta en los 33 `A` | max dif absoluta en los 528 `B` | genes expresados igual | celulas de control igual |
|---|---|---|---|---|
| A549 | **0.0e+00** | **0.0e+00** | **si** | **si** |
| BXPC3 | **0.0e+00** | **0.0e+00** | **si** | **si** |
| HAP1 | **0.0e+00** | **0.0e+00** | **si** | **si** |
| HT29 | **0.0e+00** | **0.0e+00** | **si** | **si** |
| K562 | **0.0e+00** | **0.0e+00** | **si** | **si** |
| MCF7 | **0.0e+00** | **0.0e+00** | **si** | **si** |

Cero exacto en las seis lineas y en las dos cantidades: lo que sigue se
mide sobre la misma construccion, no sobre una reimplementacion parecida.

### 0.2 Una cosa que se temia y NO pasa

A 200 celulas hay mas genes con cuenta cero que a 400, y
`cand01_refuta.py:107` manda los `log2` no finitos a `0.0`. Un bloque de
empates exactos que `argsort(kind='stable')` ordena por indice de gen daria
**el mismo orden en las dos columnas** y por tanto correlacion positiva
espuria, mayor en las mitades que en las completas. Se conto:

| linea | no finitos en completa (400) | en mitad_1 (200) | en mitad_2 (200) |
|---|---|---|---|
| A549 | 0 | 0 | 0 |
| BXPC3 | 21 | 27 | 23 |
| HAP1 | 0 | 0 | 0 |
| HT29 | 0 | 0 | 0 |
| K562 | 0 | 0 | 0 |
| MCF7 | 0 | 0 | 0 |

**Cero en 5 de las 6 lineas, y no cero en `BXPC3`.** Donde no es cero, es
mas en las mitades que en la completa (27 y 23 contra 21), que es la
direccion que inflaria `A`. El tamano es **27 genes de 7475 = 0.36 %**
de la columna.

🔴 **Y no se descuenta, por una razon que hay que decir: a profundidad
igualada el artefacto esta en las DOS cantidades.** `A` y `B_mitades` se
calculan los dos sobre columnas de 200 celulas, con el mismo numero de
empates forzados y el mismo orden de indice de gen dentro del bloque de
empates. Inflar las dos por igual es justamente lo que igualar la
profundidad consigue; lo que **no** estaba igualado era `B_completa`, que
se calcula a 400 y tiene menos. Asi que el artefacto empuja `B_completa`
hacia ABAJO respecto de `A`, es decir en contra de `B > A`, y no puede
ser la causa del resultado de aqui. Queda contado, no corregido.

---

## 1. Las tres cantidades, y el cambio de signo

`A` no se corrige hacia arriba: ni Spearman-Brown ni ninguna otra
extrapolacion. `anchor.py:3-7` dice por que -- "Deliberately NOT
Spearman-Brown corrected and not depth-corrected in any other way -- the
correction `ceiling.py` applies is a bounded-reliability extrapolation that
does not hold for five of the six scored `vcc2026` members". Se iguala la
profundidad bajando `B`, no subiendo `A`:

| | construccion | celulas por lado |
|---|---|---|
| `A_mitades` | `mitad_1(p)` contra `mitad_2(p)` | 200 vs 200 |
| `B_completa` | `poblacion(p)` contra `poblacion(q)` | 400 vs 400 |
| `B_mitades` | `mitad_1(p)` contra `mitad_1(q)` | **200 vs 200** |

Pares de un blanco consigo mismo se excluyen de `B_mitades`: a profundidad
igualada `mitad_1(p)` contra `mitad_1(p)` es el mismo vector y daria 1.0
exacto. Son 528 pares, los mismos 528 del punto.

| linea | `A` | `B_completa` | `B_mitades` | `B_comp - A` (lo que se reporto) | **`B_mit - A`** (profundidad igualada) | signo |
|---|---|---|---|---|---|---|
| A549 | +0.1575 | +0.2004 | +0.1155 | +0.0429 | **-0.0420** | **A>B** |
| BXPC3 | +0.1544 | +0.1198 | +0.0742 | -0.0347 | **-0.0802** | **A>B** |
| HAP1 | +0.0820 | +0.1053 | +0.0578 | +0.0233 | **-0.0242** | **A>B** |
| HT29 | +0.1020 | +0.1298 | +0.0739 | +0.0278 | **-0.0281** | **A>B** |
| K562 | +0.1609 | +0.2184 | +0.1322 | +0.0575 | **-0.0287** | **A>B** |
| MCF7 | +0.2012 | +0.2133 | +0.1261 | +0.0120 | **-0.0752** | **A>B** |

**`B > A` en 5 de 6 a profundidad desigual; en
0 de 6 a profundidad igualada.** No es que la diferencia se
encoja: cambia de signo en todas.

---

## 2. La barra: bootstrap PAREADO sobre blancos

200 replicas. En cada una: se remuestrean los 33 blancos **con
reemplazo**, se vuelve a partir cada blanco en mitades con semilla propia, y
se recalculan `A` y `B` **sobre los mismos blancos remuestreados**. La
cantidad que se guarda es `B - A` por replica, **pareada**, asi que la barra
es la de la diferencia y no la resta de dos barras sueltas.

**Las 200 replicas pedidas se corrieron enteras: no hubo que bajarlas.**
Coste medido, 68-97 s por linea (seis lineas, una pasada), porque el
lado real se lee una vez y lo unico que se repite por replica es la
particion en mitades y el rankeo.

**Sobre BLANCOS y no sobre los 528 pares.** Dos pares que comparten un blanco
no son independientes; remuestrear pares daria una barra falsamente
estrecha. Semillas: `20267000 + 100*fold + replica`, familia nueva.

| linea | `A` punto | `A` media boot | `B_mit` punto | `B_mit` media boot | `B_mit - A` punto | media boot | `sd` boot | **IC 95 %** | excluye 0 |
|---|---|---|---|---|---|---|---|---|---|
| A549 | +0.1575 | +0.1577 | +0.1155 | +0.1165 | **-0.0420** | -0.0412 | 0.0112 | **[-0.0602, -0.0194]** | **SI** |
| BXPC3 | +0.1544 | +0.1578 | +0.0742 | +0.0737 | **-0.0802** | -0.0840 | 0.0107 | **[-0.1028, -0.0640]** | **SI** |
| HAP1 | +0.0820 | +0.0783 | +0.0578 | +0.0576 | **-0.0242** | -0.0207 | 0.0069 | **[-0.0337, -0.0069]** | **SI** |
| HT29 | +0.1020 | +0.1018 | +0.0739 | +0.0731 | **-0.0281** | -0.0287 | 0.0077 | **[-0.0441, -0.0135]** | **SI** |
| K562 | +0.1609 | +0.1649 | +0.1322 | +0.1281 | **-0.0287** | -0.0368 | 0.0133 | **[-0.0630, -0.0143]** | **SI** |
| MCF7 | +0.2012 | +0.2141 | +0.1261 | +0.1324 | **-0.0752** | -0.0817 | 0.0130 | **[-0.1067, -0.0551]** | **SI** |

**El IC excluye 0 en 6 de 6 lineas**, y en las 6 que lo excluyen
lo excluye **por debajo de 0** (`A > B`).

### 2.1 La expectativa por azar del conteo, que es lo que exige §7

Un conteo sobre seis lineas no se lee sin su nulo. Hay dos nulos distintos y
se dan los dos, porque responden a preguntas distintas:

| pregunta | nulo | esperado | observado | `p` exacta |
|---|---|---|---|---|
| ¿cuantas lineas con `B_comp > A`? (candidato_01) | moneda, `p = 0.5` | 3.0 | 5 de 6 | **0.1094** = 7/64 |
| ¿cuantas lineas con `A > B_mit`? (aqui) | moneda, `p = 0.5` | 3.0 | 6 de 6 | **0.0156** = 1/64 |
| ¿cuantas lineas con IC 95 % que excluye 0? | `B - A = 0` en todas, `p = 0.05` | **0.30** | 6 de 6 | **1.56e-08** |
| ¿cuantas con IC que excluye 0 **por debajo**? | `p = 0.025` (una cola) | **0.15** | 6 de 6 | **2.44e-10** |


🔴 **Las dos ultimas `p` suponen las seis lineas independientes, y hay que
decir en que sentido lo son y en que sentido no.** Lo son en el ruido: las
seis lineas son poblaciones de celulas **disjuntas**, asi que bajo el nulo
`B - A = 0` el error de muestreo de una linea no sabe nada del de otra y el
binomial vale. **No** lo son en el mecanismo: las seis comparten el mismo
panel de 33 blancos y la misma causa -- una diferencia de profundidad de
construccion, que esta en las seis por igual. Asi que un 6 de 6 **no son
seis confirmaciones independientes de un efecto biologico**: es un
mecanismo medido seis veces. Para refutar una hipotesis eso es lo que hace
falta y sobra; para afirmar una, no bastaria.

🔴 **Y hay que decir que el 5 de 6 original no era el problema principal.**
Su `p = 0.109` no excluye la moneda, pero el conteo tampoco era la medicion
debil: la medicion debil era que `A` y `B` no median a la misma
profundidad. Un conteo de 6 de 6 sobre cantidades mal emparejadas habria
dado `p = 0.016` y seguiria estando mal. **La barra y el conteo no arreglan
un estimando torcido; esto se cierra por el PASO 0, no por la `p`.**

### 2.2 Con la media en vez de la mediana, por si el agregador manda

| linea | `B_mit - A` con mediana | IC 95 % | `B_mit - A` con media | IC 95 % | coinciden en signo |
|---|---|---|---|---|---|
| A549 | -0.0420 | [-0.0602, -0.0194] | -0.0791 | [-0.1102, -0.0462] | **si** |
| BXPC3 | -0.0802 | [-0.1028, -0.0640] | -0.1155 | [-0.1469, -0.0848] | **si** |
| HAP1 | -0.0242 | [-0.0337, -0.0069] | -0.0368 | [-0.0492, -0.0175] | **si** |
| HT29 | -0.0281 | [-0.0441, -0.0135] | -0.0547 | [-0.0783, -0.0293] | **si** |
| K562 | -0.0287 | [-0.0630, -0.0143] | -0.0458 | [-0.0709, -0.0305] | **si** |
| MCF7 | -0.0752 | [-0.1067, -0.0551] | -0.0912 | [-0.1140, -0.0660] | **si** |

---

## 3. Cuanto de `B - A` era profundidad

`B_completa - B_mitades` es la atenuacion por bajar de 400 a 200 celulas por
lado, medida sobre la MISMA cantidad (los mismos 528 pares, los mismos
genes, el mismo denominador). Es lo que `A` estaba pagando y `B` no:

| linea | `B_comp - A` reportado | **profundidad** = `B_comp - B_mit` | resto = `B_mit - A` | profundidad / \|reportado\| |
|---|---|---|---|---|
| A549 | +0.0429 | **+0.0848** | -0.0420 | 2.0 x |
| BXPC3 | -0.0347 | **+0.0456** | -0.0802 | 1.3 x |
| HAP1 | +0.0233 | **+0.0475** | -0.0242 | 2.0 x |
| HT29 | +0.0278 | **+0.0559** | -0.0281 | 2.0 x |
| K562 | +0.0575 | **+0.0862** | -0.0287 | 1.5 x |
| MCF7 | +0.0120 | **+0.0872** | -0.0752 | 7.3 x |

La atenuacion por profundidad va de **+0.0456 a +0.0872** y el `B - A` que se
reporto iba de +0.0120 a +0.0575 en valor absoluto. **La profundidad es de 1.3 a
7.3 veces el efecto que se reporto**, asi que no explica parte de la
diferencia: la explica entera y sobra para cambiarle el signo.

---

## 4. El segundo confundidor, que no estaba en el encargo: el control compartido

Apareció al escribir el PASO 0 y no se puede dejar fuera, porque el propio
marcador lo declara inadmisible. Todo delta es `log2(pb) - log2(cr)`, y en
`cand01_refuta.py` **las dos columnas de cada correlacion llevan el mismo
`cr`** -- las dos mitades en `A`, los dos blancos en `B`.

`cell-eval2/src/cell_eval2/anchor.py:9-18`, literal:

> under "real" the pred side's DE reference comes from the real side, so
> scoring half_b against half_a computes BOTH halves' log2FCs against
> half_a's control. That shares the control's sampling noise between the two
> quantities whose agreement is being measured and **biases the anchor
> upward** [...] Measured for de_lfc_nmae on a real panel: shared controls
> read 0.5-2.3% optimistic [...] `ceiling.py` measured the same mistake
> turning a **0.54 split-half reliability into 0.74**.

Y `anchor.py:768-779` no lo avisa, lo **rehusa**: "a replicate anchor is
only valid under per-half controls ('pred')". El estimando no es
negociable (su `spec 5.5`).

Asi que la version conforme al ancla, con el pool de control partido en dos
mitades disjuntas `c1` y `c2`, y las dos cantidades con la MISMA forma
`mitad_1` contra `mitad_2` y 200 celulas por lado:

```
A_anc = corr( log2(mitad_1(p)/c1), log2(mitad_2(p)/c2) )
B_anc = corr( log2(mitad_1(p)/c1), log2(mitad_2(q)/c2) )     p != q
```

Misma barra: bootstrap pareado sobre blancos, misma familia de semillas, y
el control se reparte de nuevo en cada replica.

| linea | celulas de control por mitad | `A_anc` | `B_anc` | `B_anc - A_anc` | media boot | **IC 95 %** | excluye 0 |
|---|---|---|---|---|---|---|---|
| A549 | 689 / 690 | +0.0464 | +0.0038 | **-0.0426** | -0.0330 | **[-0.0567, -0.0141]** | **SI** |
| BXPC3 | 1614 / 1614 | +0.1040 | +0.0251 | **-0.0789** | -0.0793 | **[-0.0952, -0.0600]** | **SI** |
| HAP1 | 1586 / 1586 | +0.0240 | +0.0038 | **-0.0202** | -0.0186 | **[-0.0297, -0.0066]** | **SI** |
| HT29 | 1368 / 1368 | +0.0402 | +0.0056 | **-0.0346** | -0.0251 | **[-0.0418, -0.0087]** | **SI** |
| K562 | 699 / 699 | +0.0402 | +0.0099 | **-0.0303** | -0.0311 | **[-0.0582, -0.0107]** | **SI** |
| MCF7 | 1334 / 1335 | +0.1261 | +0.0590 | **-0.0671** | -0.0773 | **[-0.0972, -0.0588]** | **SI** |

**Y el efecto del control compartido es mucho mayor que el 0.5-2.3 % que
`anchor.py` midio sobre `de_lfc_nmae`**, porque aqui la cantidad es un
Spearman sobre ~7000 genes y no un NMAE sobre genes DE:

| linea | `A` control comun | `A_anc` control por mitad | `A` perdida | `B_mit` control comun | `B_anc` control por mitad | `B` perdida |
|---|---|---|---|---|---|---|
| A549 | +0.1575 | +0.0464 | **-0.1111** | +0.1155 | +0.0038 | **-0.1118** |
| BXPC3 | +0.1544 | +0.1040 | **-0.0504** | +0.0742 | +0.0251 | **-0.0491** |
| HAP1 | +0.0820 | +0.0240 | **-0.0580** | +0.0578 | +0.0038 | **-0.0540** |
| HT29 | +0.1020 | +0.0402 | **-0.0618** | +0.0739 | +0.0056 | **-0.0682** |
| K562 | +0.1609 | +0.0402 | **-0.1207** | +0.1322 | +0.0099 | **-0.1222** |
| MCF7 | +0.2012 | +0.1261 | **-0.0751** | +0.1261 | +0.0590 | **-0.0671** |

Conforme al ancla, el IC excluye 0 en **6 de 6** lineas.

🔴 **Pero la tabla de arriba dice algo mas preciso que "empuja en la misma
direccion", y conviene no cobrarlo de mas: la perdida de `A` y la de `B`
son casi la misma en valor absoluto.** El ruido del control entraba
identico en las dos columnas de las dos cantidades, asi que al partirlo las
dos bajan a la par y la DIFERENCIA se queda donde estaba:

| linea | `B_mit - A` con control comun | `B_anc - A_anc` con control por mitad | se mueve |
|---|---|---|---|
| A549 | -0.0420 | -0.0426 | -0.0006 |
| BXPC3 | -0.0802 | -0.0789 | +0.0013 |
| HAP1 | -0.0242 | -0.0202 | +0.0040 |
| HT29 | -0.0281 | -0.0346 | -0.0065 |
| K562 | -0.0287 | -0.0303 | -0.0015 |
| MCF7 | -0.0752 | -0.0671 | +0.0081 |

**Asi que el control compartido NO es un confundidor de `B - A`: se cancela
en la resta.** Es un defecto del **nivel**, no de la comparacion, y tiene
dos consecuencias distintas que no hay que mezclar:

1. **Refuerza la refutacion de §1-§2**, porque `B - A` sale igual por los
   dos caminos -- el resultado no cuelga de haber compartido el control.
2. **Rompe otra afirmacion de `candidato_01.md`, la del techo.** Alli se
   leyo "`A` corre de +0.0820 a +0.2012 en las seis lineas: el test tiene
   potencia". Conforme al ancla `A` corre de **+0.0240 a +0.1261**: entre
   un tercio y dos tercios de lo que se reporto. El test sigue teniendo
   potencia -- `A` no es cero en ninguna linea -- pero el techo es
   bastante mas bajo, y eso **si** entra al diseno del candidato 02. No se
   corrige hacia arriba: se reporta medido.

---

## 5. Lo que reemplaza a la hipotesis: la descomposicion, medida sin modelo

El encargo pedia esta descomposicion **solo si la hipotesis sobrevivia**. No
sobrevive, asi que se mide por otra razon y se declara: `candidato_01.md`
uso `B > A` para afirmar que el `b` "esta capturando la componente
dominante de la varianza" y que lo especifico del blanco es "la
minoritaria incluso en el lado real". Esa afirmacion se queda sin apoyo
cuando `B > A` cae, y es exactamente la que decide el diseno del candidato
02. Asi que se mide directamente en vez de heredarla.

**Sin suponer ningun modelo.** Lo unico que se usa es que las dos mitades de
un blanco tienen ruido independiente -- son celulas disjuntas -- y media
cero. Con `h1_i`, `h2_i` los deltas de las dos mitades del blanco `i`,
centrados por columna sobre los genes, y `m1`, `m2` sus promedios sobre los
33 blancos:

```
V_tot = (1/n) Σ_i <h1_i, h2_i> / G          senal reproducible total
V_sh  = <m1, m2> / G                        COMPARTIDA entre blancos
V_sp  = (1/n) Σ_i <h1_i-m1, h2_i-m2> / G    ESPECIFICA del blanco
```

`V_tot = V_sh + V_sp` es una **identidad algebraica exacta**, no un ajuste:
los terminos cruzados se anulan porque `Σ_i (h1_i - m1) = 0`. Y el producto
**cruzado de mitades** es lo que quita el ruido sin modelarlo: su esperanza
es el producto de las senales verdaderas, porque el ruido de las dos mitades
es independiente. La varianza de **una** mitad si llevaria el ruido dentro y
habria que suponer un modelo para sacarlo. **Se puede medir sin modelo, y se
mide.**

Y por lo mismo de §4, el control va **partido**: con un `cr` comun el
termino `-log2(cr)` es identico en todas las columnas, cae entero dentro de
`m1` y `m2` con la misma realizacion de ruido, y `<m1, m2>` contaria el
ruido del control como senal compartida. Se dan las dos para que se vea
cuanto era eso:

| linea | | `V_tot` | `V_compartida` | `V_especifica` | **frac compartida** | **frac especifica** | residuo de la identidad |
|---|---|---|---|---|---|---|---|
| A549 | control comun | 0.0177 | 0.0096 | 0.0081 | **54.4 %** | **45.6 %** | 6.9e-18 |
|  | **control por mitad** | 0.0098 | 0.0017 | 0.0081 | **17.3 %** | **82.7 %** | 1.7e-18 |
| BXPC3 | control comun | 0.0157 | 0.0054 | 0.0103 | **34.4 %** | **65.6 %** | 1.0e-17 |
|  | **control por mitad** | 0.0134 | 0.0031 | 0.0103 | **22.9 %** | **77.1 %** | 8.7e-18 |
| HAP1 | control comun | 0.0083 | 0.0031 | 0.0052 | **37.4 %** | **62.6 %** | 0.0e+00 |
|  | **control por mitad** | 0.0062 | 0.0010 | 0.0052 | **15.7 %** | **84.3 %** | 1.7e-18 |
| HT29 | control comun | 0.0112 | 0.0044 | 0.0068 | **39.3 %** | **60.7 %** | 0.0e+00 |
|  | **control por mitad** | 0.0082 | 0.0014 | 0.0068 | **16.6 %** | **83.4 %** | 1.7e-18 |
| K562 | control comun | 0.0130 | 0.0077 | 0.0053 | **59.3 %** | **40.7 %** | 1.7e-18 |
|  | **control por mitad** | 0.0065 | 0.0012 | 0.0053 | **18.8 %** | **81.2 %** | 5.2e-18 |
| MCF7 | control comun | 0.0093 | 0.0050 | 0.0043 | **53.8 %** | **46.2 %** | 1.7e-18 |
|  | **control por mitad** | 0.0071 | 0.0028 | 0.0043 | **39.1 %** | **60.9 %** | 8.7e-19 |

**Con el control compartido la fraccion compartida sale de 34 a 59 %**, que es
la lectura que respaldaba "la componente dominante". **Con el control
partido va de 16 a 39 % y lo especifico de 61 a 84 %.**

🔴 **Da la vuelta a la premisa del candidato 02.** Con el control por
mitad lo especifico del blanco es la componente **mayoritaria en 6 de las
6 lineas** (con el control comun lo era en 3 de 6, que es por lo que la
lectura salia al contrario). Un candidato no tiene que acertar "la
componente pequena": la parte especifica del lado real es grande y
reproducible -- `V_especifica` no lleva ruido dentro, es un producto
cruzado de mitades -- y lo que la hacia parecer pequena era el ruido de
muestreo del control contado como senal compartida.

Y notese que `V_especifica` es **la misma columna en las dos filas de cada
linea**: el control comun no la toca. Todo el movimiento de la fraccion
sale de `V_compartida`, que es donde cae el termino `-log2(cr)`. Eso es la
prueba de que el efecto es el que se dice y no una reescala global.

---

## 6. Las filas trabajadas

### 6.1 Los 33 blancos de A549, uno por fila

La columna `A_mitades` es la lista cuya **mediana** es el `A` de §1.
`B_mit medio` es el promedio de los 32 pares en que ese blanco entra a
profundidad igualada, y `B_comp medio` los mismos 32 a profundidad completa.

⚠️ **La ultima columna NO es el estimador.** El numero de linea es
`mediana(33 valores de A) - mediana(528 pares de B)`, no la mediana de esta
columna: `B` se agrega sobre pares, no sobre blancos, y un blanco entra en
32 pares. La columna esta aqui como **intermedio y diagnostico de
concentracion** -- para ver si el resultado de linea lo arrastran dos o tres
filas o esta repartido -- no como la cantidad que se reporta.

| # | blanco | celulas | `A_mitades` | `B_mit` medio (32 pares) | `B_comp` medio (32) | profundidad por blanco | **`B_mit` - `A`** | `A > B_mit` |
|---|---|---|---|---|---|---|---|---|
| 1 | `TYK2` | 400 | +0.5271 | +0.1410 | +0.2136 | +0.0726 | **-0.3861** | **si** |
| 2 | `IFNAR2` | 400 | +0.5241 | +0.1335 | +0.2060 | +0.0725 | **-0.3907** | **si** |
| 3 | `IFNAR1` | 400 | +0.4553 | +0.1460 | +0.2236 | +0.0776 | **-0.3094** | **si** |
| 4 | `JAK1` | 400 | +0.3205 | +0.1372 | +0.2187 | +0.0814 | **-0.1832** | **si** |
| 5 | `STAT3` | 400 | +0.2732 | +0.0988 | +0.1531 | +0.0543 | **-0.1743** | **si** |
| 6 | `JUN` | 400 | +0.2476 | +0.1138 | +0.1927 | +0.0789 | **-0.1338** | **si** |
| 7 | `MAP3K14` | 400 | +0.2239 | +0.1175 | +0.1906 | +0.0731 | **-0.1065** | **si** |
| 8 | `IRF1` | 400 | +0.2189 | +0.1200 | +0.2115 | +0.0915 | **-0.0989** | **si** |
| 9 | `STAT1` | 400 | +0.1918 | +0.1072 | +0.1874 | +0.0802 | **-0.0845** | **si** |
| 10 | `RAP1GAP` | 400 | +0.1806 | +0.1304 | +0.2197 | +0.0893 | **-0.0503** | **si** |
| 11 | `BATF2` | 400 | +0.1798 | +0.1045 | +0.1847 | +0.0802 | **-0.0753** | **si** |
| 12 | `TRAFD1` | 400 | +0.1766 | +0.0939 | +0.1712 | +0.0773 | **-0.0827** | **si** |
| 13 | `ZBP1` | 400 | +0.1736 | +0.1181 | +0.1888 | +0.0707 | **-0.0555** | **si** |
| 14 | `MAPK8` | 400 | +0.1699 | +0.1232 | +0.2134 | +0.0903 | **-0.0468** | **si** |
| 15 | `ETS2` | 400 | +0.1698 | +0.1073 | +0.1827 | +0.0755 | **-0.0626** | **si** |
| 16 | `CEBPG` | 400 | +0.1620 | +0.1115 | +0.2018 | +0.0903 | **-0.0505** | **si** |
| 17 | `ELK1` | 400 | +0.1575 | +0.1005 | +0.1591 | +0.0586 | **-0.0570** | **si** |
| 18 | `SP100` | 400 | +0.1573 | +0.1133 | +0.2060 | +0.0928 | **-0.0440** | **si** |
| 19 | `ETS1` | 400 | +0.1544 | +0.1317 | +0.2209 | +0.0892 | **-0.0226** | **si** |
| 20 | `IRF7` | 400 | +0.1530 | +0.1040 | +0.1843 | +0.0804 | **-0.0491** | **si** |
| 21 | `HERC6` | 400 | +0.1524 | +0.0862 | +0.1782 | +0.0920 | **-0.0662** | **si** |
| 22 | `DTX3L` | 400 | +0.1403 | +0.1149 | +0.2011 | +0.0862 | **-0.0254** | **si** |
| 23 | `PARP12` | 400 | +0.1394 | +0.1226 | +0.2218 | +0.0992 | **-0.0168** | **si** |
| 24 | `H1F0` | 400 | +0.1375 | +0.1155 | +0.2067 | +0.0912 | **-0.0220** | **si** |
| 25 | `STAT4` | 400 | +0.1321 | +0.1159 | +0.2001 | +0.0842 | **-0.0163** | **si** |
| 26 | `STAT5A` | 400 | +0.1311 | +0.1150 | +0.2044 | +0.0894 | **-0.0161** | **si** |
| 27 | `RNF114` | 400 | +0.1245 | +0.1183 | +0.2062 | +0.0879 | **-0.0062** | **si** |
| 28 | `AKT1` | 400 | +0.1237 | +0.1102 | +0.1971 | +0.0868 | **-0.0135** | **si** |
| 29 | `BRD9` | 400 | +0.1177 | +0.1095 | +0.2009 | +0.0914 | **-0.0082** | **si** |
| 30 | `SOCS1` | 400 | +0.1173 | +0.1053 | +0.1888 | +0.0834 | **-0.0119** | **si** |
| 31 | `ID2` | 400 | +0.1144 | +0.1251 | +0.2094 | +0.0844 | **+0.0107** | no |
| 32 | `HERC5` | 400 | +0.0984 | +0.1174 | +0.1999 | +0.0825 | **+0.0190** | no |
| 33 | `VAV1` | 400 | +0.0924 | +0.1191 | +0.2205 | +0.1014 | **+0.0267** | no |

**`A > B_mit` en 30 de los 33 blancos** de A549, y los tres que no son los tres
de `A` mas baja de la tabla. El resultado de linea no lo arrastran dos o
tres filas: esta repartido. (Mediana de
la columna `A_mitades` = +0.1575, que es el `A` de §1; mediana de los 528 pares
= +0.1155, que es el `B_mitades`.)

### 6.2 Veinte replicas del bootstrap, una por fila

De las 200 de A549. `blancos distintos` es cuantos de los 33
sobreviven al remuestreo con reemplazo (el resto son repeticiones);
`pares` es cuantos pares de blancos distintos quedan tras excluir los de un
blanco consigo mismo. El IC de §2 es el percentil 2.5/97.5 de la ultima
columna sobre las 200 replicas.

| replica | semilla | blancos distintos | pares | `A` mediana | `B_mit` mediana | **`B_mit - A`** |
|---|---|---|---|---|---|---|
| 0 | 20267000 | 22 | 512 | +0.1507 | +0.1108 | **-0.0399** |
| 1 | 20267001 | 17 | 505 | +0.1546 | +0.1097 | **-0.0449** |
| 2 | 20267002 | 23 | 513 | +0.1409 | +0.1165 | **-0.0245** |
| 3 | 20267003 | 21 | 510 | +0.1518 | +0.1167 | **-0.0351** |
| 4 | 20267004 | 19 | 507 | +0.1701 | +0.1171 | **-0.0530** |
| 5 | 20267005 | 22 | 515 | +0.1787 | +0.1172 | **-0.0615** |
| 6 | 20267006 | 22 | 515 | +0.1463 | +0.1139 | **-0.0324** |
| 7 | 20267007 | 20 | 511 | +0.1664 | +0.1113 | **-0.0551** |
| 8 | 20267008 | 18 | 507 | +0.1611 | +0.1109 | **-0.0501** |
| 9 | 20267009 | 22 | 513 | +0.1602 | +0.1079 | **-0.0523** |
| 10 | 20267010 | 21 | 514 | +0.1536 | +0.1199 | **-0.0337** |
| 11 | 20267011 | 22 | 512 | +0.1735 | +0.1164 | **-0.0571** |
| 12 | 20267012 | 17 | 501 | +0.1528 | +0.1236 | **-0.0291** |
| 13 | 20267013 | 21 | 514 | +0.1542 | +0.1204 | **-0.0339** |
| 14 | 20267014 | 22 | 515 | +0.1372 | +0.1203 | **-0.0169** |
| 15 | 20267015 | 21 | 514 | +0.1616 | +0.1163 | **-0.0453** |
| 16 | 20267016 | 23 | 515 | +0.1470 | +0.1190 | **-0.0279** |
| 17 | 20267017 | 20 | 511 | +0.1532 | +0.1226 | **-0.0307** |
| 18 | 20267018 | 20 | 512 | +0.1718 | +0.1172 | **-0.0546** |
| 19 | 20267019 | 19 | 508 | +0.1674 | +0.1169 | **-0.0505** |

Las 20 primeras replicas dan `B_mit - A` negativo en 20 de 20; las 200
completas, en 200 de 200. El signo no depende de que replicas se miren.

---

## 7. Que se hizo, que no, y que queda abierto

**No se hizo** -- ninguna de las prohibiciones del encargo se roza:

| prohibicion | estado |
|---|---|
| corregir `A` por profundidad (Spearman-Brown o cualquier otra) | **no**. Se bajo `B` a 200 celulas. `anchor.py:3-7` citado en §1 |
| construir un candidato, tocar un arnes, tocar un bundle | **no**. `build_fold_adatas` se llama como parametro, sin modificar; no se escribe ningun bundle |
| leer el lado real para CONSTRUIR algo | **no**. Todo esto es diagnostico del problema. Ninguna cifra de aqui entra a un modelo |
| `~/vcc/benchmark/`, `/datos/benchmark_VEDADO_techo_metricas/` | **no leidos**. Las unicas entradas son `IFNB`, `suelo_panel.json` y `cand01_refuta.json` |
| envios, `vcc sample`, commits | **ninguno** |

**Lo que entro y lo que salio:**

| | |
|---|---|
| entradas | `/datos/ref/convertido/Seurat_object_IFNB_Perturb_seq.h5ad`, `modelo/out/suelo_panel.json`, `modelo/out/cand01_refuta.json` (solo para verificar) |
| codigo | `modelo/src/b_contra_a.py`, `b_contra_a_anchor.py`, `b_contra_a_descomp.py`, `b_contra_a_reporte.py` |
| salidas | `b_contra_a.json`, `b_contra_a_anchor.json`, `b_contra_a_descomp.json`, este `.md` |
| semillas | real `20260001+fold`, emision `20266000+100*fold`, mitades del punto `777`, **bootstrap `20267000+100*fold+replica`**, control partido `20267999+fold` |

**Lo que queda abierto, declarado el dia que aparece:**

1. **La mascara de genes de `cand01_refuta.py` depende de la emision del
   candidato.** `exp = (cr >= 1) & (cp >= 1)` y `cp` es el pseudobulk del
   control **predicho**, que es un remuestreo de 400 celulas del control real.
   Para un diagnostico de la verdad eso no deberia entrar. Se midio la
   version de solo-verdad (`cr >= 1`) como verificacion y el resultado no
   cambia:

| linea | genes con `exp` | genes con solo `cr >= 1` | `B_mit - A` con `exp` | con solo-verdad |
|---|---|---|---|---|
| A549 | 6934 | 7169 | -0.0420 | -0.0386 |
| BXPC3 | 7475 | 7673 | -0.0802 | -0.0774 |
| HAP1 | 7857 | 8045 | -0.0242 | -0.0257 |
| HT29 | 6999 | 7224 | -0.0281 | -0.0322 |
| K562 | 7403 | 7616 | -0.0287 | -0.0280 |
| MCF7 | 7266 | 7444 | -0.0752 | -0.0771 |

   No afecta al resultado; afecta a que `cand01_refuta.py` no es reutilizable
   tal cual como medidor de la verdad. Cerrarlo cuesta una linea.

2. **El empate de agregador no se resolvio, se midio.** La mediana sobre 33
   valores y sobre 528 no son el mismo estimador aunque sea la misma funcion.
   Aqui la mediana y la media coinciden en signo en las seis lineas (§2.2),
   asi que no decide nada -- pero en `candidato_01.md` si decidia (5 de 6
   contra 4 de 6, §0). Afecta a la **interpretacion** de ese conteo, no a
   este resultado.

3. **No se midio si `A` por si sola basta como techo.** `A` a 200 celulas es
   lo que es; lo que un candidato puede alcanzar a 400 no se extrapola desde
   aqui, por la misma razon por la que no se corrige con Spearman-Brown.
   Afecta a la interpretacion del techo, no a la comparacion `B` contra `A`.

---

## 8. Para el candidato 02, en una frase

Tres cosas, y las tres son medidas, no lecturas:

1. **`B > A` no existe.** Era la diferencia de profundidad entre 200 y 400
   celulas por lado. A profundidad igualada `A > B` en las seis lineas, con
   el IC de la diferencia excluyendo 0 en las seis. El control compartido
   **no** entra aqui: se cancela en la resta.
2. **El techo `A` es mas bajo de lo que se reporto**, porque el `A` de
   `candidato_01.md` comparte el control entre las dos mitades, que es
   exactamente lo que `anchor.py` rehusa. Sigue siendo distinto de cero en
   las seis lineas, asi que el test conserva potencia.
3. **Y lo especifico del blanco no es la componente minoritaria.** Medido
   sin modelo y con el control partido, es la mayoritaria. Asi que la frase
   "un candidato que quiera superar al `b` tiene que acertar la componente
   pequena" se queda sin su apoyo: la componente especifica es la grande.

**Lo que esto le quita al candidato 02 es una excusa, no una via.** El
candidato 01 saco `C` ~ 0 y concordancia de signo 0.5; eso sigue en pie y
no lo explica que lo especifico fuera minoritario, porque no lo es.
