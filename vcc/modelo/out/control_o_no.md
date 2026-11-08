# ¿Era la covarianza la senal equivocada, o no hay senal en el control?

## 🔴 0. Esto es un DIAGNOSTICO TRANSDUCTIVO. Va arriba, no en nota al pie.

Todo lo que se ajusta en este informe se ajusta **sobre la verdad de
Jiang**, que es exactamente lo que un envio no tiene. Ninguna cifra, ningun
coeficiente y ningun modelo de aqui entra NUNCA a una submission, y ninguno
es una via. Lo unico que sale de aqui es la RESPUESTA: si el techo de "solo
control" es alto o es cero.

En concreto, y para que no se pueda confundir despues:

| | |
|---|---|
| el modelo que da el numero alto de §6 | se ajusta con la verdad perturbada **del propio blanco** que luego evalua. Un envio no la tiene. **No es una via.** |
| el modelo que responde a la pregunta util (§5) | se ajusta sin ver el blanco que evalua, y da **cero** |
| lo que se puede llevar a un candidato | la respuesta y la cota, no los coeficientes |

---

## 1. La respuesta: ni (a) ni (b) -- la senal esta, el mapa no se aprende

El encargo pedia discriminar entre (a) la covarianza era la senal
equivocada y otra senal del control si capta lo especifico, y (b) nada del
control lo predice. **La medicion no da ninguna de las dos, y la tercera
opcion es concreta y esta medida:**

| | medido | que significa |
|---|---|---|
| **los rasgos de coexpresion SI llevan la senal** | `r_mismo` (par) = **+0.2446 a +0.3998**, positivo en 29-33 de 33 blancos, nulo ~0 con `sd` ~0.007 | NO es (b). Y es coexpresion, no propension por gen: los rasgos de solo-`h` dan +0.0904 a +0.2021 |
| **pero no transfiere a un blanco nuevo** | `r_loto` estricto = **-0.2043 a +0.0685**, **0 de 6** significativas, 2 positivas y 4 negativas | NO es (a) en sentido usable |
| **y se ve por que** | la consistencia de signo de los coeficientes entre blancos es **39 a 73 %** (moneda = 50 %) | cada blanco tiene su propio coeficiente. Con 33 blancos no hay un mapa compartido que aprender |

🔴 **Y la limitacion que el encargo mandaba declarar decide el matiz.** El
MDE del diseno (80 % de potencia, alfa 0.05) va de **0.137 a 0.297** segun
linea, y el techo intra-blanco esta en 0.24-0.40. Es decir: **el diseno no
puede separar "no transfiere nada" de "transfiere hasta ~0.2".** Lo que si
queda cerrado es que no transfiere MAS que eso. Con 33 blancos no se puede
pedir mas, y eso tambien es el resultado.

**Para el candidato 02, en una linea:** la covarianza del control **no era
la senal equivocada** -- es la que lleva la informacion -- pero el candidato
01 fijo su coeficiente (`r = 1 - rho`) y el coeficiente correcto es
**distinto en cada blanco**, con el signo cambiando mas o menos como una
moneda. Un candidato de solo-control que aprenda ese coeficiente de otros
blancos no puede funcionar; lo que falta no es una senal mejor, es una
forma de fijar el coeficiente por blanco sin ver su verdad.

---

## 2. El diseno

### 2.1 El objetivo: el delta ESPECIFICO, con el control partido

Para el blanco `g`, `E_g = h_g - m` con `m` el promedio sobre los 33
blancos, calculado como en `b_contra_a_descomp.py` para que sea **la misma
cantidad** cuyo 61-84 % de varianza se midio alli. El control va PARTIDO en
dos mitades disjuntas (`c1` para una mitad del blanco, `c2` para la otra):
con control comun el ruido de muestreo del control se cuenta como senal
compartida y `m` se lo come (§4 de `b_contra_a.md`).

### 2.2 La metrica: el ruido fuera, por producto cruzado de mitades

```
E[<E1_g, E2_g>]        = ||t_g||^2      t_g = la senal verdadera
E[<p, (E1_g+E2_g)/2>]  = <p, t_g>       (el ruido de las dos mitades es
                                         independiente: no entra)
r_true = <p, (E1+E2)/2> / ( ||p|| * sqrt(<E1_g,E2_g>) )
```

Y el resumen sobre blancos es el **agrupado**, no la media:

```
R_agr = Σ_g <p_g, y_g>  /  Σ_g ( ||p_g|| * sqrt(cruz_g) )
```

por una razon que se comprueba y no se asume: como `Σ_g E[h,g] = 0`, el
numerador total de un predictor que solo dependa del gen `h` es
`<p, Σ_g E_g> = 0`, asi que **`R_agr` vale CERO EXACTO para cualquier
predictor sin informacion especifica del blanco**. Es la comprobacion que
delato el primer defecto (§4.1). Medido, el conjunto `solo_h` da:

| linea | A549 | BXPC3 | HAP1 | HT29 | K562 | MCF7 |
|---|---|---|---|---|---|---|
| `solo_h` `R_agr` | -0.0033 | +0.0000 | +0.0000 | +0.0000 | -0.0009 | -0.0000 |

Cero en las seis. Sin esa comprobacion, nada de lo que sigue seria legible.

### 2.3 Los rasgos, declarados antes de medirlos

Todos del POOL DE CONTROL y del SIMBOLO del gen. Ninguno toca el lado
perturbado. Para el par (blanco `g`, gen `h`):

| # | rasgo | que es | tipo |
|---|---|---|---|
| 1 | `cov` | covarianza de log1p normalizado en el pool de control | par |
| 2 | `r_pearson` | correlacion de Pearson | par |
| 3 | `r_spearman` | correlacion de Spearman (rangos sobre celulas) | par |
| 4 | `r_parcial_10` | Pearson tras proyectar fuera los 10 primeros PC | par |
| 5 | `r_parcial_50` | idem con 50 PC | par |
| 6 | `abs_r_pearson` | `|Pearson|`, la magnitud separada del signo | par |
| 7 | `abs_r_spearman` | `|Spearman|` | par |
| 8 | `delta_r` | `r_pearson - r_parcial_50`: cuanto era programa global | par |
| 9 | `rango_cor` | rango de `h` por `|r_pearson|` con `g`, normalizado | par |
| 10 | `mismo_modulo` | 1 si `h` cae en el modulo de coexpresion de `g` (KMeans, 30) | par |
| 11 | `dist_emb` | distancia euclidea `g`-`h` en el embedding de 50 PC | par |
| 12 | `cos_emb` | coseno `g`-`h` en ese embedding | par |
| 13 | `misma_familia` | raiz alfabetica comun de >=3 letras (STAT1/STAT4 -> STAT) | par |
| 14 | `es_el_blanco` | 1 si `h` ES `g` | par |
| 15 | `h_media` | media de log1p de `h` | solo h |
| 16 | `h_var` | varianza | solo h |
| 17 | `h_disp` | var/media | solo h |
| 18 | `h_frac_cero` | fraccion de celulas de control con cuenta 0 | solo h |
| 19 | `h_rango_expr` | rango del nivel de expresion | solo h |
| 20 | `g_media` | media de log1p del blanco | solo g |
| 21 | `g_var` | varianza del blanco | solo g |
| 22 | `g_frac_cero` | fraccion de ceros del blanco | solo g |

**Las correlaciones parciales estan porque la marginal mezcla dos cosas.**
La covarianza marginal entre `g` y `h` incluye el programa global del
contexto; la parcial, tras quitar los primeros componentes, deja la
relacion directa. `delta_r` mide explicitamente la diferencia.

🔴 **Lo que se descarto por no ser calculable sin el lado perturbado, con
nombre y razon** -- el encargo lo pedia explicitamente:

| descartado | por que |
|---|---|
| pertenencia a una via conocida (IFN-alfa/beta, JAK-STAT) de una base externa | no es del pool de control; seria una cifra importada (PLAN §0.2) |
| respuesta conocida del gen a interferon | **es el lado perturbado** |
| efecto medido del mismo blanco en OTRA linea | **es el lado perturbado** |
| distancia en el genoma o cromosoma | no esta en el `h5ad` (solo `var/_index`); sacarla de una anotacion externa seria importar |

**Los de solo-`g` no pueden hacer nada y van declarados como tales:** son
constantes dentro de una columna, y una correlacion por columna no cambia al
sumarle una constante. Estan por simetria del registro, no por esperanza.

### 2.4 La validacion cruzada: por BLANCO, dejando blancos fuera

Nunca por gen ni por celula para la pregunta de §5. Si el modelo memoriza
pares `(g,h)` de entrenamiento, el numero no dice nada sobre un blanco
nuevo. Se deja **un blanco fuera** (33 pliegues), con eleccion interna de
`alpha` por CV agrupada en los 32 de entrenamiento, y el nulo por
permutacion de etiquetas de blanco repitiendo la CV entera.

El numero de §6 (`r_mismo`) SI parte genes, y por eso **no se lee como
generalizacion**: responde a otra pregunta y lo dice en su seccion.

---

## 3. La potencia del diseno, declarada ANTES del resultado

Son **33 blancos por linea**, y las seis lineas comparten **los mismos 33**,
asi que agrupar lineas NO sube el `n` de la pregunta "¿generaliza a un
blanco nuevo?". El nulo se obtiene por permutacion de etiquetas de blanco
repitiendo la CV completa, y de su desviacion sale el efecto minimo
detectable:

```
MDE(80 %, alfa 0.05 dos colas) = (1.96 + 0.8416) * sd_nulo = 2.802 * sd_nulo
```

| linea | `sd` del nulo (LOTO estricto) | **MDE 80 %** | `A_esp` mediana (fiabilidad del propio objetivo) |
|---|---|---|---|
| A549 | 0.0815 | **0.228** | +0.0633 |
| BXPC3 | 0.0872 | **0.244** | +0.1245 |
| HAP1 | 0.1061 | **0.297** | +0.0603 |
| HT29 | 0.1015 | **0.284** | +0.0601 |
| K562 | 0.0673 | **0.189** | +0.0607 |
| MCF7 | 0.0490 | **0.137** | +0.1253 |

🔴 **Dicho antes de ver el resultado: este diseno solo puede detectar
efectos de 0.137 para arriba.** Cualquier cosa por debajo es indistinguible
de cero con 33 blancos, y el veredicto correcto para ese caso es "no se
puede decidir", no "no hay senal". El MDE llega hasta 0.297 en la linea
peor, que es la mitad del techo intra-blanco que se medira en §6: **el
diseno esta justo en el limite de poder contestar la pregunta.**

---

## 4. Tres defectos encontrados y corregidos, antes de cualquier lectura

Los tres los delataron comprobaciones puestas a proposito, no una
sospecha. Los tres cambiaban el resultado.

### 4.1 Una fuga estructural del diseno: la media de entrenamiento ES el blanco dejado fuera

El objetivo `E_i = h_i - m` con `m` sobre los 33 fuerza `Σ_i E_i = 0`. Y
entonces, para el pliegue que deja fuera el blanco `j`:

```
media de los 32 de entrenamiento = (Σ_{i≠j} E_i) / 32 = -E_j / 32
```

Medido en A549: `corr(media de entrenamiento, blanco dejado fuera) =
**-1.000000** exacto`, y `max|media_tr - (-E_j/32)| = 4.4e-09`. Cualquier
modelo que capture esa media se lleva correlacion NEGATIVA con el blanco
nuevo sin saber nada de el.

**Lo delato `solo_h`**, que no puede llevar informacion especifica del
blanco y puntuaba de -0.13 a -0.22 con desviacion EXACTAMENTE cero.

El arreglo sale derivado, no tanteado. El estimando correcto por pliegue usa
como referencia solo los 32 de entrenamiento, y sustituyendo
`h_i = E_i + S/33` con `Σ_{k≠j} E_k = -E_j`:

```
y de ENTRENAMIENTO del pliegue j :  y_i = E_i + E_j / 32
                                    Σ_{i≠j} y_i = -E_j + 32·E_j/32 = 0   ✓
objetivo de EVALUACION            :  (33/32) · E_j   -> multiplo POSITIVO,
                                    la correlacion no cambia
```

Toca solo las `y` de entrenamiento, que es donde estaba la fuga, y no hubo
que recalcular ningun rasgo. Comprobado: `solo_h` pasa de **-0.1744 a
-0.0033** en A549 y queda en cero en las seis lineas (§2.2).

### 4.2 Un resumen sesgado mio: la media de cocientes no hereda el cero

Declare que un predictor de solo-`h` tiene correlacion media cero sobre los
blancos. Eso es cierto de los NUMERADORES, no del cociente: cada blanco
divide por su propio `sqrt(cruz_g)`, y una media de cocientes con
denominadores distintos no hereda la propiedad. Por eso el titular es
`R_agr` (§2.2) y no la media. Se dan las tres para que se vea cuanto movia:

| linea | `R_agr` (el que se lee) | mediana | media (la sesgada) |
|---|---|---|---|
| A549 | -0.1232 | -0.0125 | -0.0600 |
| BXPC3 | -0.0091 | +0.0307 | -0.0002 |
| HAP1 | -0.2043 | -0.0636 | -0.0711 |
| HT29 | +0.0685 | +0.0457 | +0.0963 |
| K562 | -0.0560 | -0.0370 | -0.0674 |
| MCF7 | +0.0560 | +0.0526 | +0.0564 |

### 4.3 🔴 El gen blanco se identifica sin el indicador, asi que `sin_autoblanco` no era un test limpio

Quitar el rasgo `es_el_blanco` **no le quita al modelo saber cual es el gen
blanco**, porque los rasgos de correlacion lo delatan con certeza:

| rasgo | valor en el gen blanco | maximo en los otros ~7000 genes | p99.99 |
|---|---|---|---|
| `r_pearson` | **+1.0000** | +0.3031 | +0.1666 |
| `r_spearman` | **+1.0000** | +0.2233 | +0.1685 |
| `cos_emb` | **+1.0000** | +0.7336 | +0.6255 |
| `dist_emb` | **0.0000** | +0.9312 | +0.8180 |

(medido en A549; el gen blanco esta expresado en 21-25 de los 33 blancos
segun linea.) Y su `|E|` es **10.2 veces la mediana**, asi que un solo punto
con esa palanca mueve una Pearson sobre 7000 genes.

**"Cuando se perturba `g`, `g` baja" es cierto, es de solo control y tiene
valor** -- se sabe QUE gen se perturba sin mirar nada perturbado. Pero es UN
gen de 7000 y no es "otra senal del control que capta lo especifico". Para
medir la pregunta del encargo hay que sacarlo de los **tres** sitios: de los
rasgos, de las filas de entrenamiento (donde es un punto extremo que el
modelo persigue) y de la evaluacion. Eso es la version ESTRICTA, y es la
unica que se lee en §5.

| linea | `R_agr` con todo | sin el indicador (NO limpio) | **ESTRICTO** (los tres sitios) |
|---|---|---|---|
| A549 | +0.1334 | +0.0820 | **-0.1232** |
| BXPC3 | +0.1542 | +0.1167 | **-0.0091** |
| HAP1 | +0.2955 | +0.2491 | **-0.2043** |
| HT29 | +0.2256 | +0.1957 | **+0.0685** |
| K562 | +0.2717 | +0.2182 | **-0.0560** |
| MCF7 | +0.1535 | +0.1170 | **+0.0560** |

La columna del medio es la que reporte en el avance como "apunta a (a)".
**Era el gen blanco.** Queda escrito porque el error estuvo en el camino.

---

## 5. Resultado 1 -- ¿generaliza a un blanco NUEVO? **No, y eso es un cero medido**

Version estricta, ridge con `alpha` por CV interna, nulo de 500
permutaciones de etiquetas de blanco:

| linea | genes evaluados | **`R_agr`** | mediana | positivos | nulo | z | `p` perm | MDE 80 % | recortado 1 % | GB |
|---|---|---|---|---|---|---|---|---|---|---|
| A549 | 7168-7169 | **-0.1232** | -0.0125 | 13/29 | -0.0235 | -1.22 | 0.888 | 0.228 | -0.0401 | -0.0007 |
| BXPC3 | 7672-7673 | **-0.0091** | +0.0307 | 20/33 | -0.0201 | +0.13 | 0.463 | 0.244 | +0.0409 | -0.0324 |
| HAP1 | 8044-8045 | **-0.2043** | -0.0636 | 12/33 | -0.0334 | -1.61 | 0.940 | 0.297 | -0.0014 | -0.0068 |
| HT29 | 7223-7224 | **+0.0685** | +0.0457 | 21/33 | -0.0313 | +0.98 | 0.186 | 0.284 | +0.0583 | +0.0725 |
| K562 | 7615-7616 | **-0.0560** | -0.0370 | 10/31 | -0.0091 | -0.70 | 0.780 | 0.189 | -0.0713 | -0.0524 |
| MCF7 | 7443-7444 | **+0.0560** | +0.0526 | 22/33 | -0.0090 | +1.33 | 0.100 | 0.137 | +0.0801 | -0.0153 |

**0 de 6 significativas. 2 positivas, 4 negativas.** El signo no es
consistente entre lineas, los `z` van de -1.6 a +1.3, y la version recortada
(fuera el 1 % de `|E|` mas extremo) no rescata nada. **Todo el +0.14 a +0.30
de la version con el gen blanco era el gen blanco.**

### 5.1 Contra las dos referencias que pedia el encargo

| linea | `C` del candidato 01 | `A_anc` (fiabilidad del delta completo) | `A_esp` (fiabilidad del delta ESPECIFICO) | **`r_loto` estricto** | fraccion del techo |
|---|---|---|---|---|---|
| A549 | +0.0009 | +0.0464 | +0.0633 | **-0.1232** | n/d (negativo) |
| BXPC3 | +0.0009 | +0.1040 | +0.1245 | **-0.0091** | n/d (negativo) |
| HAP1 | +0.0046 | +0.0240 | +0.0603 | **-0.2043** | n/d (negativo) |
| HT29 | -0.0054 | +0.0402 | +0.0601 | **+0.0685** | 114 % |
| K562 | +0.0108 | +0.0402 | +0.0607 | **-0.0560** | n/d (negativo) |
| MCF7 | -0.0042 | +0.1261 | +0.1253 | **+0.0560** | 45 % |

**¿Mejora sobre `C` ~ 0? No de forma distinguible.** `C` iba de -0.0054 a
+0.0108 y esto va de -0.2043 a +0.0685 con el cero dentro del nulo en las seis.
El candidato 01 no fallo por elegir mal la senal: **ningun predictor de
solo-control que aprenda de otros blancos lo hace mejor en este diseno.**

🔴 Y hay que decir lo que esto NO dice: con el MDE en 0.137-0.297, un efecto
verdadero de, digamos, 0.10 no se habria visto. **"No transfiere" aqui
significa "no transfiere por encima del MDE", no "transfiere cero".**

---

## 6. Resultado 2 -- ¿llevan los rasgos informacion? **Mucha, y es coexpresion**

🔴 **Esto NO es generalizacion a un blanco nuevo y no se puede leer como
tal.** Se ajusta en la MITAD de los genes DEL PROPIO blanco y se predice la
otra mitad. Responde a la otra pregunta, la que separa (a) de (b): ¿existe
alguna funcion de los rasgos de solo-control que describa el delta
especifico de este blanco? El gen blanco va fuera de los tres sitios, igual
que en §5.

| linea | **solo rasgos de PAR** (coexpresion) | positivos | nulo | `sd` nulo | **z** | solo rasgos de `h` | positivos | ambos |
|---|---|---|---|---|---|---|---|---|
| A549 | **+0.3447** | 32/33 | -0.0001 | 0.0074 | **+47** | +0.1543 | 28/33 | +0.3626 |
| BXPC3 | **+0.3160** | 33/33 | -0.0013 | 0.0050 | **+63** | +0.1155 | 28/33 | +0.3405 |
| HAP1 | **+0.3998** | 31/33 | -0.0045 | 0.0081 | **+50** | +0.1399 | 30/33 | +0.4064 |
| HT29 | **+0.2956** | 33/33 | +0.0025 | 0.0058 | **+50** | +0.0904 | 31/33 | +0.3201 |
| K562 | **+0.2823** | 29/33 | -0.0008 | 0.0058 | **+49** | +0.2021 | 29/33 | +0.3133 |
| MCF7 | **+0.2446** | 33/33 | -0.0005 | 0.0039 | **+62** | +0.1352 | 25/33 | +0.2626 |

El nulo permuta el objetivo sobre los GENES dentro del blanco manteniendo
los rasgos: rompe la correspondencia gen-objetivo sin tocar nada mas. 10
permutaciones, `sd` de 0.0039 a 0.0081, **z de +47 a +63**. No es marginal.

### 6.1 Y el gen blanco no es lo que lo sostiene -- al contrario que en §5

Con una implementacion independiente (`control_o_no_disc.py`), en la
linea donde se corrio, comparando con y sin el gen blanco:

| linea | modelo | con el gen blanco | **sin el gen blanco** | nulo | `sd` nulo | z |
|---|---|---|---|---|---|---|
| A549 | gb | +0.1763 | **+0.1771** | -0.0005 | 0.0057 | +31.1 |
| A549 | ridge | +0.3493 | **+0.3626** | -0.0001 | 0.0072 | +50.6 |

**Quitar el gen blanco no baja nada aqui: sube un poco.** Eso es lo que
separa las dos preguntas. El gen blanco era TODO lo que transferia entre
blancos (§4.3, §5) y NADA de lo que hay dentro de un blanco.

Y vale como comprobacion cruzada: el `ambos` de la tabla de arriba para
A549 (+0.3626) y el `sin ridge` de aqui (+0.3626) salen de dos scripts
distintos y coinciden.

🔴 **Solo se corrio en A549.** Las otras cinco lineas llevan el numero de
`par` de la tabla de arriba, que excluye el gen blanco por construccion y
ya responde lo mismo, pero la comparacion con/sin por linea queda como
hueco declarado, no medido.

**Y el desglose contesta la pregunta del encargo.** Los rasgos de PAR
(covarianza, Pearson, Spearman, parcial, modulo, embedding, familia) dan
+0.2446 a +0.3998, positivo en 29-33 de 33 blancos. Los de solo-`h` -- "este
gen tiende a moverse", que no depende del blanco -- dan solo +0.0904 a +0.2021.
**El grueso lo lleva la coexpresion, no una propension por gen.**

Y el techo: `A_esp` atenuada va de +0.060 a +0.125, mientras esto llega a
+0.40. No se contradice -- `A_esp` es la correlacion contra una mitad
ruidosa y `r_true` es contra la senal verdadera, con el ruido fuera por
construccion. Lo que dice es que **hay bastante mas senal especifica
describible de lo que la fiabilidad cruda sugiere.**

---

## 7. Por que no transfiere: los coeficientes cambian de signo entre blancos

Si la senal esta dentro de cada blanco (§6) pero no transfiere (§5), la
explicacion tiene que estar en el mapa. Se mira donde esta la senal -- el
ajuste INTRA-blanco -- con un ridge por blanco sobre rasgos estandarizados,
y se pregunta no cuales pesan sino **si pesan en el mismo sentido**:

La consistencia es la fraccion de blancos con el signo MAYORITARIO del
coeficiente, contando los ceros en ninguno de los dos lados. Por eso un
rasgo que queda a cero exacto en muchos blancos puede bajar del 50 %:

| rasgo | `|coef|` medio | signo mediano | **consistencia de signo** | coef = 0 |
|---|---|---|---|---|
| `h_frac_cero` | 0.0412 | -0.0086 | **58 %** | 0 % |
| `r_parcial_10` | 0.0338 | +0.0011 | **58 %** | 0 % |
| `h_var` | 0.0330 | -0.0093 | **64 %** | 0 % |
| `h_disp` | 0.0210 | -0.0003 | **54 %** | 0 % |
| `r_pearson` | 0.0208 | +0.0000 | **61 %** | 0 % |
| `dist_emb` | 0.0205 | +0.0063 | **63 %** | 0 % |
| `r_parcial_50` | 0.0176 | -0.0015 | **60 %** | 0 % |
| `h_rango_expr` | 0.0129 | +0.0004 | **57 %** | 0 % |
| `delta_r` | 0.0126 | +0.0013 | **57 %** | 0 % |
| `h_media` | 0.0110 | -0.0022 | **59 %** | 0 % |
| `r_spearman` | 0.0108 | -0.0005 | **55 %** | 0 % |
| `cos_emb` | 0.0090 | -0.0000 | **52 %** | 0 % |
| `cov` | 0.0076 | +0.0001 | **58 %** | 0 % |
| `abs_r_pearson` | 0.0064 | +0.0004 | **58 %** | 0 % |
| `abs_r_spearman` | 0.0061 | +0.0003 | **57 %** | 0 % |
| `mismo_modulo` | 0.0058 | +0.0006 | **61 %** | 0 % |
| `rango_cor` | 0.0052 | -0.0002 | **54 %** | 0 % |
| `misma_familia` | 0.0011 | +0.0000 | **41 %** | 30 % |
| `g_media` | 0.0000 | +0.0000 | **0 %** | 100 % |
| `g_var` | 0.0000 | +0.0000 | **0 %** | 100 % |
| `g_frac_cero` | 0.0000 | +0.0000 | **0 %** | 100 % |

**Los tres de solo-`g` salen a cero EXACTO en los 33 blancos**, que es la
propiedad declarada en §2.3 confirmada por el ajuste: son constantes dentro
de una columna, asi que son colineales con el intercepto y el ridge los
manda a cero. No es que no sirvan poco: no pueden hacer nada.

**Y de los 18 rasgos que si reciben peso, ninguno pasa del 64 % de
consistencia, con la moneda en 50 %%.**
No hay un rasgo que empuje en el mismo sentido en todos los blancos. Esa es
la razon medida de que §5 salga cero teniendo §6 un +0.3: **no existe un
mapa compartido que 33 blancos puedan revelar**, porque el mapa es distinto
en cada blanco.

Dicho de otro modo: el candidato 01 fijo el coeficiente en `r = 1 - rho`,
el mismo para los 33 blancos. La medicion dice que el coeficiente correcto
existe y es fuerte, pero **cambia de blanco a blanco con el signo casi
libre**. Ninguna constante global -- ni la del candidato 01 ni otra -- puede
acertarlo.

---

## 8. Las filas trabajadas

### 8.1 Los 33 blancos de A549, uno por fila

`A_esp` es la fiabilidad mitad-contra-mitad del delta especifico de ese
blanco; `cruz` su producto cruzado (el numerador sin ruido, y `cruz <= 0`
es el motivo por el que un blanco queda fuera del resumen); `r_mismo (par)`
lo de §6 y `r_loto` lo de §5. **Las dos ultimas columnas son las que no hay
que confundir**: la primera se ajusta con la verdad del propio blanco, la
segunda sin verla.

| # | blanco | gen blanco expresado | `A_esp` | `cruz` | `r_mismo` (par) | `r_loto` estricto |
|---|---|---|---|---|---|---|
| 14 | `IFNAR2` | si | +0.4998 | +423.84 | +0.3505 | +0.0362 |
| 31 | `TYK2` | si | +0.4975 | +409.26 | +0.3084 | +0.1315 |
| 13 | `IFNAR1` | si | +0.3308 | +195.98 | +0.2693 | +0.0288 |
| 27 | `STAT3` | si | +0.2238 | +119.99 | +0.4331 | -0.1663 |
| 17 | `JAK1` | si | +0.1832 | +91.44 | +0.1678 | +0.0364 |
| 18 | `JUN` | si | +0.1476 | +73.39 | +0.0943 | -0.0296 |
| 19 | `MAP3K14` | si | +0.1210 | +56.27 | +0.4140 | -0.2656 |
| 15 | `IRF1` | si | +0.1091 | +50.00 | +0.2081 | -0.0125 |
| 30 | `TRAFD1` | si | +0.1061 | +52.86 | +0.5536 | -0.4709 |
| 8 | `ETS2` | si | +0.1015 | +46.93 | +0.4080 | -0.1221 |
| 2 | `BATF2` | si | +0.1011 | +44.47 | +0.4920 | -0.4753 |
| 6 | `ELK1` | no | +0.0978 | +46.12 | +0.5832 | +0.1532 |
| 26 | `STAT1` | si | +0.0912 | +42.25 | +0.4770 | -0.2547 |
| 11 | `HERC6` | si | +0.0756 | +36.04 | +0.3455 | -0.2078 |
| 33 | `ZBP1` | no | +0.0699 | +32.83 | +0.1461 | +0.1275 |
| 16 | `IRF7` | si | +0.0683 | +29.25 | +0.1139 | +0.0338 |
| 5 | `DTX3L` | si | +0.0633 | +27.03 | +0.2949 | -0.1765 |
| 4 | `CEBPG` | no | +0.0515 | +21.72 | +0.0636 | -0.0113 |
| 9 | `H1F0` | si | +0.0420 | +17.76 | +0.3592 | -0.0338 |
| 25 | `SP100` | si | +0.0399 | +15.36 | +0.1825 | -0.0577 |
| 20 | `MAPK8` | si | +0.0380 | +15.59 | +0.2784 | +0.0947 |
| 29 | `STAT5A` | no | +0.0368 | +15.31 | +0.1605 | -0.1029 |
| 1 | `AKT1` | si | +0.0351 | +14.95 | +0.2690 | +0.0870 |
| 24 | `SOCS1` | no | +0.0346 | +14.97 | +0.1579 | -0.2046 |
| 3 | `BRD9` | si | +0.0322 | +13.31 | +0.5135 | +0.1406 |
| 22 | `RAP1GAP` | no | +0.0290 | +11.57 | +0.2245 | +0.0151 |
| 23 | `RNF114` | no | +0.0252 | +10.10 | +0.1823 | -0.1213 |
| 28 | `STAT4` | no | +0.0218 | +9.42 | +0.0112 | +0.0510 |
| 7 | `ETS1` | si | +0.0087 | +3.36 | +1.1708 | +0.0366 |
| 10 | `HERC5` | si | +0.0025 | +1.08 | +0.4440 | n/d |
| 21 | `PARP12` | si | +0.0013 | +0.53 | +1.2124 | n/d |
| 12 | `ID2` | no | -0.0067 | -2.86 | +0.3048 | n/d |
| 32 | `VAV1` | no | -0.0493 | -21.13 | n/d | n/d |

`r_mismo` positivo en 32 de 32 filas; `r_loto` positivo en 13 de 29. 
Las dos columnas se calculan sobre los MISMOS blancos y los mismos rasgos;
la unica diferencia es si el ajuste vio la verdad de ese blanco.

🔴 **Dos cosas de la tabla que hay que leer con cuidado, y son del
estimador, no de los datos.**

1. **`r_mismo` pasa de 1 en alguna fila** (`PARP12` da +1.2124). No es un
   error: `r_true` divide por `sqrt(cruz)`, que es una ESTIMACION de
   `||senal||`, y cuando `cruz` es chico el cociente se desboca. Un
   cociente asi no esta acotado por 1 como lo estaria una Pearson. Por eso
   **el titular es `R_agr`**, que suma numeradores y denominadores por
   separado antes de dividir y no deja que un blanco con `cruz` chico
   domine -- y por eso no se reporta la media de esta columna.
2. **El `cruz` de la tabla es el del gen completo; la exclusion se decide
   con el `cruz` restringido a la mascara** del pliegue (sin el gen
   blanco). Por eso hay filas con `cruz > 0` en la tabla y `n/d` en
   `r_loto`: al quitar el gen blanco su `cruz` se vuelve <= 0 y el blanco
   sale del resumen. Son 4 de 33 en esta linea.

### 8.2 Los 21 rasgos, uno por fila, con su consistencia de signo

Promediado sobre las seis lineas. Es la tabla de §7 ordenada por
consistencia en vez de por peso, que es como se ve que ninguno sirve:

| # | rasgo | consistencia | `|coef|` medio | signo mediano | coef = 0 | ¿sirve para transferir? |
|---|---|---|---|---|---|---|
| 1 | `h_var` | **64 %** | 0.0330 | -0.0093 | 0 % | no |
| 2 | `dist_emb` | **63 %** | 0.0205 | +0.0063 | 0 % | no |
| 3 | `r_pearson` | **61 %** | 0.0208 | +0.0000 | 0 % | no |
| 4 | `mismo_modulo` | **61 %** | 0.0058 | +0.0006 | 0 % | no |
| 5 | `r_parcial_50` | **60 %** | 0.0176 | -0.0015 | 0 % | no |
| 6 | `h_media` | **59 %** | 0.0110 | -0.0022 | 0 % | no |
| 7 | `h_frac_cero` | **58 %** | 0.0412 | -0.0086 | 0 % | no |
| 8 | `r_parcial_10` | **58 %** | 0.0338 | +0.0011 | 0 % | no |
| 9 | `abs_r_pearson` | **58 %** | 0.0064 | +0.0004 | 0 % | no |
| 10 | `cov` | **58 %** | 0.0076 | +0.0001 | 0 % | no |
| 11 | `h_rango_expr` | **57 %** | 0.0129 | +0.0004 | 0 % | no |
| 12 | `abs_r_spearman` | **57 %** | 0.0061 | +0.0003 | 0 % | no |
| 13 | `delta_r` | **57 %** | 0.0126 | +0.0013 | 0 % | no |
| 14 | `r_spearman` | **55 %** | 0.0108 | -0.0005 | 0 % | no |
| 15 | `h_disp` | **54 %** | 0.0210 | -0.0003 | 0 % | no |
| 16 | `rango_cor` | **54 %** | 0.0052 | -0.0002 | 0 % | no |
| 17 | `cos_emb` | **52 %** | 0.0090 | -0.0000 | 0 % | no |
| 18 | `misma_familia` | **41 %** | 0.0011 | +0.0000 | 30 % | no |
| 19 | `g_media` | **0 %** | 0.0000 | +0.0000 | 100 % | no (sin peso) |
| 20 | `g_var` | **0 %** | 0.0000 | +0.0000 | 100 % | no (sin peso) |
| 21 | `g_frac_cero` | **0 %** | 0.0000 | +0.0000 | 100 % | no (sin peso) |

**Ninguno de los que reciben peso llega al 75 %%.** El mejor esta en
64 %, que sobre 33 blancos son 21 aciertos de signo contra los 16.5
de una moneda: una ventaja de 4 blancos. No es una regla, y menos una
que se pueda aplicar a un blanco que no se ha visto.

---

## 9. Veredicto, y que no se hizo

**El veredicto, con sus tres partes:**

1. **No es (b).** Los rasgos de coexpresion del control describen el delta
   especifico de un blanco dado a `r_true` = +0.24 a +0.40, positivo en
   29-33 de 33 blancos, con nulo en cero y `sd` de 0.004 a 0.008 (z de +47
   a +63). No es ruido y no es el gen blanco.
2. **No es (a) en sentido usable.** Eso no transfiere a un blanco no visto:
   **0 de 6** lineas significativas, 2 positivas y 4 negativas.
3. **La razon esta medida:** la consistencia de signo de los coeficientes
   entre blancos no pasa del **73 %** en ningun par (linea, rasgo), y
   promediada sobre lineas no pasa del **64 %**, con la moneda en 50 %.
   El mapa es especifico del blanco.

**Y la limitacion que el encargo mandaba declarar, aplicada:** con MDE de
**0.137 a 0.297** y un techo intra-blanco de 0.24 a 0.40, el diseno **no puede
separar "no transfiere" de "transfiere hasta ~0.2"**. Para la parte
informativa (§6) si hay potencia de sobra -- los `z` estan en +47/+63. Para
la parte de transferencia, 33 blancos no alcanzan, y eso es parte del
resultado.

**No se hizo** -- ninguna prohibicion del encargo se roza:

| prohibicion | estado |
|---|---|
| construir un candidato, tocar arneses o bundles | **no**. `build_fold_adatas` se llama como parametro con un `predict_fn` trivial que devuelve el control sin tocarlo; no se escribe ningun bundle |
| usar el lado perturbado para calcular una caracteristica | **no**. Los 22 rasgos salen del pool de control y del simbolo del gen. Los descartados por esto van con nombre en §2.3 |
| presentar el modelo ajustado como una via | **no**. §0 lo declara arriba: es transductivo y esta prohibido |
| `~/vcc/benchmark/`, `/datos/benchmark_VEDADO_techo_metricas/` | **no leidos** |
| envios, `vcc sample`, commits | **ninguno** |

| | |
|---|---|
| entradas | `/datos/ref/convertido/Seurat_object_IFNB_Perturb_seq.h5ad`, `modelo/out/suelo_panel.json`; `cand01_refuta.json` y `b_contra_a_anchor.json` solo como referencia |
| codigo | `control_o_no_rasgos.py`, `_loto.py`, `_estricto.py`, `_disc.py`, `_desglose.py`, `_pesos.py`, `_reporte.py` |
| semillas | real `20260001+fold`, mitades `777`, control partido `20267999+fold`, familia nueva de este encargo **`20268000+100*fold`** (PCA, KMeans, particion de genes, permutaciones) |

**Lo que queda abierto:**

1. **El `alpha` del nulo se fija al modal del ajuste real** en vez de
   repetir la CV interna en cada permutacion. Afinar `alpha` dentro del nulo
   lo haria mas optimista, asi que fijarlo es la direccion conservadora para
   el nulo. Afecta a la interpretacion del `p`, no al signo del resultado.
2. **GB da ~0 donde ridge da +0.3 dentro del blanco y los dos dan ~0 entre
   blancos.** No se eligio el modelo que mejor quedaba: §5 reporta los dos y
   coinciden en el cero. Por que GB no ve lo que ridge ve dentro del blanco
   no se cerro -- la sospecha medida es que su capacidad se la llevan las
   filas extremas y los rasgos de solo-`h`, pero no se verifico.
3. **Un solo reparto de genes** en el discriminador (una semilla). El nulo
   es estrecho (`sd` ~0.007) asi que no cambia el signo, pero la magnitud
   +0.24/+0.40 no lleva barra por reparto.
4. **El `m` de referencia son 33 blancos**, no la poblacion. El delta
   "especifico" es especifico respecto de ESTE panel; con otro panel la
   particion compartido/especifico se mueve.

---

## 10. Para el candidato 02

**La covarianza del control no era la senal equivocada.** Es la que lleva la
informacion: los rasgos de par dan +0.24 a +0.40 dentro del blanco, muy por
encima de los de solo-`h`. Lo que el candidato 01 fijo mal no fue la senal,
fue el **coeficiente**: puso `r = 1 - rho`, el mismo para los 33 blancos, y
el coeficiente correcto cambia de blanco a blanco con el signo casi libre
(consistencia maxima 70 %).

De ahi se siguen dos cosas, y la segunda es la que importa:

1. **No sirve buscar otra senal del control.** Ya se probaron 22 rasgos,
   incluidas parciales que quitan el programa global, modulos y embeddings.
   La senal esta; el problema no es de senal.
2. **Un predictor de solo-control que aprenda el coeficiente de otros
   blancos no puede funcionar en este diseno.** Para usar lo de §6 haria
   falta fijar el coeficiente por blanco sin ver su verdad, y este encargo
   no encontro ninguna via para eso -- ni la descarta: dice que no esta en
   los 22 rasgos probados ni se puede aprender de 33 blancos.

🔴 **Y lo de §0 otra vez, porque es donde se rompe un envio:** el +0.40 de
§6 se obtiene ajustando con la verdad perturbada del blanco que luego se
evalua. **No es un techo alcanzable por una submission.** El techo
alcanzable por algo que no ve la verdad es el de §5, y ese es cero.
