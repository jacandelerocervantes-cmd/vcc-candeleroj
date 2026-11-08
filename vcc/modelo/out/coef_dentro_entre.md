# ¿Se reproduce el COEFICIENTE dentro de un blanco? El control que le faltaba al 58 %

## 🔴 0. Esto es un DIAGNOSTICO TRANSDUCTIVO. Va arriba.

Todo lo que se ajusta aqui se ajusta **sobre la verdad de Jiang**, que es
exactamente lo que un envio no tiene. **Ninguna cifra, ningun coeficiente y
ningun modelo de este informe entra NUNCA a una submission, y ninguno es una
via.** Lo unico que sale de aqui es la respuesta a una pregunta de lectura:
cuanto del 42 % de desacuerdo de signo de `control_o_no.md` §7 era diferencia
real entre blancos y cuanto era ruido de estimacion.

No se construyo ningun candidato, no se toco ningun arnes ni bundle, no se uso
el lado perturbado para calcular ningun rasgo, y no se eligio modelo ni `alpha`
por el resultado. El detalle, en §8.

---

## 1. La respuesta, en una tabla

`W` = acuerdo de signo entre el coeficiente de la **mitad_A** y el de la
**mitad_B del mismo blanco**. `E` = acuerdo entre la **mitad_A de un blanco y
la mitad_A de otro**. Las dos a **media mitad de genes**, que es la leccion de
`b_contra_a.md`. Moneda = 50 % en las dos.

| | medido | que significa |
|---|---|---|
| **`W` esta claramente por encima de `E`** | `W` = **0.653 a 0.726**, `E` = **0.485 a 0.495**; `W - E` = **+0.16 a +0.21**, IC 95 % excluye 0 en **6 de 6** lineas | **la explicacion (ii) queda cerrada.** Los coeficientes NO son "los mismos con ruido": si lo fueran, `W` y `E` saldrian iguales |
| **y la correlacion dice lo mismo, mas fuerte** | `W_corr` = **+0.48 a +0.59**, `E_corr` = **-0.01 a +0.03**; `W - E` = **+0.54 a +0.61**, positiva en **200 de 200** replicas en las 6 lineas | las dos medidas no discrepan |
| **pero `W` no esta cerca de reproducible** | **0.67 de media**: una de cada tres comparaciones de signo DENTRO de un blanco no coincide | el coeficiente tampoco se estima bien con 3500 genes. El ruido de estimacion es grande, solo no lo es todo |
| **y el acuerdo verdadero entre blancos, quitada la atenuacion, esta en el azar** | `T` = **0.467 a 0.487** en las 6 lineas, por DEBAJO de 0.50 | **la explicacion (i) se sostiene.** El mapa es de verdad distinto en cada blanco |

**En una linea:** el 58 % de §7 no significaba lo que §7 dijo que significaba,
pero la conclusion de §7 era correcta **por una razon que §7 no podia ver**. El
coeficiente se reproduce dentro de un blanco bastante mas que entre blancos
(`W - E` > 0 en 6 de 6, IC excluye 0), asi que los coeficientes **difieren de
verdad**; y una vez descontado el ruido de estimacion, el acuerdo verdadero
entre blancos queda **en el azar o por debajo**. La direccion NO es promediar
entre mas blancos.

🔴 **Y hay que decir que el desenlace que salio no es ninguno de los tres que se
declararon.** Se declararon `W ≈ E` (~58/58), `W >> E` (~90/58) y `W` tambien
bajo (~60). Salio **`W` = 0.67, `E` = 0.49**: por encima de `E` sin ambiguedad
-- eso cierra el "los coeficientes son ruido" -- y a la vez lejos de
reproducible. Las dos mitades de la lectura declarada se cumplen a la vez y no
se puede elegir una: **los coeficientes difieren de verdad Y estan mal
estimados.** §3 lo lee sin reinterpretarlo.

---

## 2. El diseno

### 2.1 Lo que se parte, y a que profundidad

Los ~7000-8000 genes de cada blanco se parten en **dos mitades disjuntas** con
semilla propia por replica. Para el blanco `g` se ajusta un ridge sobre la
mitad_A y otro sobre la mitad_B, y se comparan los dos vectores de
coeficientes.

```
W (dentro) : coef(mitad_A del blanco g)  vs  coef(mitad_B del blanco g)
E (entre)  : coef(mitad_A del blanco g)  vs  coef(mitad_A del blanco q)   g != q
```

🔴 **Las DOS a media mitad de genes.** El 58 % de §7 salia de ajustes sobre
TODOS los genes; compararlo contra un `W` de media mitad repetiria exactamente
el confundidor de profundidad que `b_contra_a.md` destapo. El `E` a genes
completos se da **en su propia columna y declarado** (§4.3), no mezclado.

### 2.2 Las condiciones, las mismas de `control_o_no.md`

| condicion | como |
|---|---|
| version **ESTRICTA** (§4.3) | el gen blanco fuera de los rasgos (`es_el_blanco` no esta en el ajuste), fuera de las filas y fuera de la evaluacion. Se hereda la misma mascara `auto` de los rasgos |
| los mismos **22 rasgos** | los de `control_o_no.md` §2.3, leidos de `rasgos_<linea>.npz`. Estandarizados **dentro de cada mitad** (`mu`, `sg` de esa mitad), `y` centrada dentro de la mitad, sin intercepto -- igual que `control_o_no_pesos.py` |
| mismo **ridge** | el de §7. `alpha` por CV interna de 4 pliegues **sobre genes dentro de cada mitad**, de `(1, 10, 100, 1000, 10000)`, criterio = Pearson prediccion-contra-`y` fuera de pliegue |
| 🔴 los **3 rasgos de solo-`g`** | `g_media`, `g_var`, `g_frac_cero` siguen EN el ajuste (mismo ridge) y salen **del conteo**. Son constantes dentro de una columna, asi que su columna estandarizada es cero exacto y su coeficiente es cero exacto: contarlos meteria ceros que no son desacuerdos. **El conteo es sobre 18 rasgos, no 21** |

### 2.3 El `alpha` no se fijo a mano, y su desacuerdo ES el resultado

El encargo lo pedia asi: si el `alpha` elegido difiere entre las dos mitades,
eso es parte del ruido de estimacion. **Difiere en el 64-76 % de los blancos**
(§4.4), a menudo por tres o cuatro ordenes de magnitud. Es una medicion
independiente de lo mismo que dice `W`, y no habria aparecido con `alpha` fijo.

Se da **tambien** la version con `alpha = 10` fijo, que es el de §7, por
comparabilidad (§4.2). No se eligio entre las dos por el resultado: las dos
estan y coinciden.

### 2.4 Las dos medidas, las dos declaradas

* **consistencia de signo pareada**: fraccion de los 18 rasgos en que dos
  vectores de coeficientes coinciden en signo. Moneda = 50 % **exactamente**,
  tanto para `W` como para `E`, que es lo que las hace comparables.
* **correlacion de Pearson** entre los dos vectores de 18 coeficientes, que usa
  la magnitud y no solo el signo.

Las dos, por si discrepan. **No discrepan** (§4.1, §4.3).

### 2.5 Las barras: bootstrap pareado sobre blancos

200 replicas, como `b_contra_a.md` §2. En cada una: **nueva particion de genes**,
se remuestrean los 33 blancos **con reemplazo**, y `W` y `E` se recalculan
**sobre los mismos blancos remuestreados**. La cantidad guardada es `W - E` por
replica, **PAREADA**. Semillas: familia nueva **`20269000 + 100*fold + replica`**.

**Sobre BLANCOS y no sobre los 528 pares**, por la misma razon que alli: dos
pares que comparten un blanco no son independientes. Las 200 replicas se
corrieron enteras en las seis lineas; coste 118-138 s por linea.

🔴 **Un par del bootstrap cuyos dos miembros son el MISMO blanco original se
excluye de `E`.** Con remuestreo con reemplazo un blanco puede salir repetido, y
compararlo consigo mismo da 100 % por construccion, lo que inflaria `E` y
hundiria `W - E`. Sin esa exclusion la barra seria optimista en el sentido
contrario al que conviene al resultado.

---

## 3. 🔴 Lo que cada desenlace significaba, declarado antes de medir, y el que salio

Esto se escribio antes de ver los numeros:

| desenlace declarado | que significaria | ¿salio? |
|---|---|---|
| `W ≈ E` (las dos ~58 %) | los coeficientes son **RUIDO**. El relato del "mapa por blanco" se cae, el cero de §5 pide otra explicacion, y la direccion seria promediar entre MAS blancos | **NO.** `W - E` = +0.16 a +0.21, IC excluye 0 en 6/6, y en correlacion +0.54 a +0.61 con 200/200 replicas positivas |
| `W >> E` (~90 % contra ~58 %) | los coeficientes **DIFIEREN de verdad**, el relato se sostiene, y hace falta estructura externa que diga que regula a que por gen | **A MEDIAS.** `W > E` sin ambiguedad, pero `W` = 0.65-0.73, no ~0.90 |
| `W` tambien bajo (~60 %) | el coeficiente **no es estimable** ni dentro de un blanco con 3500 genes, y entonces NADA de §7 era interpretable | **A MEDIAS.** `W` = 0.67, por encima del ~60 declarado pero lejos de reproducible: 1 de cada 3 comparaciones falla |

**La lectura, sin reinterpretar nada.** Los tres desenlaces se declararon como
excluyentes y no lo eran. Lo que salio es la **conjuncion de los dos ultimos**,
y cada uno se lee con la consecuencia que se le asigno de antemano:

1. **De `W > E` (6/6, IC excluye 0):** los coeficientes difieren de verdad. Se
   sostiene el relato del mapa por blanco, y **no** la direccion de promediar
   entre mas blancos. El cero de §5 no pide otra explicacion.
2. **De `W` = 0.67 y no 0.90:** el coeficiente esta mal estimado con 3500 genes.
   Una parte grande del 42 % de desacuerdo de §7 era ruido, no diferencia -- y
   por eso §7 no era interpretable **como §7 lo interpretaba**, aunque su
   conclusion sobreviva.

Y la consecuencia que se le asigno al segundo desenlace -- "hace falta
estructura externa que diga que regula a que por gen" -- **es la que queda en
pie**, porque la que se le asigno al primero esta cerrada.

### 3.1 Cuanto del 42 % era ruido: la cuenta, con su supuesto declarado

El modelo tiene un solo parametro y se escribe antes de usarlo. Si cada
coeficiente estimado invierte su signo respecto al verdadero con probabilidad
`e`, independiente entre mitades, y dos blancos tienen coeficientes verdaderos
que coinciden en signo con probabilidad `T`:

```
W = e^2 + (1-e)^2                 (dos estimaciones del MISMO coeficiente verdadero)
E = T*W + (1-T)*(1-W)             (dos verdaderos que coinciden con prob T, cada
                                   uno estimado con el mismo ruido)
  =>  e = (1 - sqrt(2W-1)) / 2         T = (E - (1-W)) / (2W - 1)
```

| linea | `W` | `E` | **`e`** = flips de signo por ruido | **`T`** = acuerdo VERDADERO entre blancos |
|---|---|---|---|---|
| A549 | 0.672 | 0.489 | **0.207** | **0.467** |
| BXPC3 | 0.690 | 0.495 | **0.192** | **0.487** |
| HAP1 | 0.657 | 0.491 | **0.220** | **0.472** |
| HT29 | 0.726 | 0.485 | **0.164** | **0.467** |
| K562 | 0.653 | 0.494 | **0.223** | **0.479** |
| MCF7 | 0.665 | 0.490 | **0.213** | **0.468** |

**Uno de cada cinco signos de coeficiente es ruido** (`e` = 0.16 a 0.22), y una
vez descontado eso, **el acuerdo verdadero entre blancos es 0.467 a 0.487: el
azar, o un pelo por debajo.** En las seis lineas, por debajo de 0.50.

🔴 **El supuesto, declarado:** la misma tasa de flip `e` dentro y entre blancos,
e independencia entre las dos mitades. Si `e` fuera mayor entre blancos que
dentro (por ejemplo porque `||E_g||` varia entre blancos), `T` estaria sesgada.
Lo que NO depende de este modelo es el resultado de §1: `W - E` > 0 con IC que
excluye 0 es una comparacion directa, sin parametros.

**Y hay una razon estructural por la que `T` puede caer por debajo de 0.50, y no
es ruido.** `Σ_g E_g = 0` por construccion del delta especifico (comprobado:
`max|Σ_g E_g| / mediana|E|` = 1.3e-06). Esa restriccion fuerza una dependencia
NEGATIVA entre los deltas de blancos distintos, del orden de `-1/32`, y con ella
entre sus coeficientes. Un `T` de 0.47 es compatible con "independientes mas la
restriccion", no exige antagonismo biologico.

---

## 4. Las cifras, por linea

### 4.1 El resultado principal: `W`, `E` y la barra pareada

`alpha` por CV interna. 200 replicas bootstrap pareadas sobre blancos.

| linea | genes | genes/mitad | **`W` signo** | **`E` signo** | **`W - E`** | **IC 95 %** | excluye 0 | `W` corr | `E` corr |
|---|---|---|---|---|---|---|---|---|---|
| A549 | 7169 | 3584 | **0.672** | **0.489** | **+0.1924** | **[+0.1081, +0.2631]** | **SI** | +0.583 | +0.019 |
| BXPC3 | 7673 | 3835 | **0.690** | **0.495** | **+0.2138** | **[+0.1604, +0.2680]** | **SI** | +0.592 | +0.009 |
| HAP1 | 8045 | 4021 | **0.657** | **0.491** | **+0.1687** | **[+0.1120, +0.2342]** | **SI** | +0.492 | +0.029 |
| HT29 | 7224 | 3612 | **0.726** | **0.485** | **+0.2109** | **[+0.1525, +0.2687]** | **SI** | +0.569 | -0.003 |
| K562 | 7616 | 3807 | **0.653** | **0.494** | **+0.1865** | **[+0.1281, +0.2464]** | **SI** | +0.478 | -0.010 |
| MCF7 | 7444 | 3722 | **0.665** | **0.490** | **+0.1615** | **[+0.1034, +0.2195]** | **SI** | +0.526 | +0.004 |

**El IC excluye 0 en 6 de 6 lineas, y en las 6 lo excluye por ARRIBA.**

Bajo el nulo `W - E = 0` en todas, 6 de 6 con IC que excluye 0 por arriba tiene
`p = 0.025^6` = **2.4e-10** si las seis lineas fueran independientes. **Y hay
que decir en que sentido no lo son**, igual que en `b_contra_a.md` §2.1: lo son
en el ruido (seis poblaciones de celulas disjuntas), **no** en el mecanismo (el
mismo panel de 33 blancos y la misma causa en las seis). Un 6 de 6 aqui es **un
mecanismo medido seis veces**, no seis confirmaciones independientes.

Y el dato que cierra la explicacion (ii) sin ningun modelo: **`E` corr va de
-0.010 a +0.029**, es decir **cero**. Los vectores de coeficientes de dos
blancos distintos no estan correlacionados ni un poco. No es que apunten al
mismo sitio con ruido: no apuntan al mismo sitio.

### 4.2 Con `alpha = 10` fijo, el de §7, por comparabilidad

| linea | `W` signo | `E` signo | `W - E` | `W` corr | `E` corr |
|---|---|---|---|---|---|
| A549 | 0.672 | 0.484 | +0.1873 | +0.677 | +0.010 |
| BXPC3 | 0.739 | 0.493 | +0.2460 | +0.746 | +0.008 |
| HAP1 | 0.685 | 0.495 | +0.1900 | +0.594 | +0.023 |
| HT29 | 0.754 | 0.487 | +0.2675 | +0.637 | -0.007 |
| K562 | 0.652 | 0.498 | +0.1536 | +0.484 | -0.004 |
| MCF7 | 0.704 | 0.483 | +0.2202 | +0.569 | -0.018 |

**Fijar `alpha` SUBE `W`** (de 0.653-0.726 a 0.652-0.754) y sube `W` corr (de
0.478-0.592 a 0.484-0.746), y deja `E` igual. Tiene sentido y hay que decirlo:
parte del desacuerdo entre mitades con CV era que las dos mitades elegian
`alpha` distinto, lo que cambia la ESCALA del vector ademas de su direccion.

🔴 **Eso tambien descarta un confundidor en la direccion opuesta**, que es la
razon de dar esta columna. Se podia sospechar que `W` corr estaba inflada por
una escala compartida: si las dos mitades de un blanco eligen `alpha` parecido,
sus coeficientes se encogen parecido y correlacionan por eso y no por su perfil.
**Con `alpha` fijo esa via desaparece por completo, y `W` corr no baja: sube.**
El confundidor existia pero iba al contrario.

### 4.3 El `E` a genes COMPLETOS, en su propia columna

Lo que `E` vale si los dos ajustes ven los ~7000 genes en vez de 3500. **No se
compara con `W`**: `W` es de media mitad y esta es la columna de profundidad
distinta, declarada aparte.

| linea | `E` signo (media mitad) | **`E` signo (genes completos)** | `E` corr (completos) |
|---|---|---|---|
| A549 | 0.489 | **0.494** | +0.038 |
| BXPC3 | 0.495 | **0.486** | +0.006 |
| HAP1 | 0.491 | **0.493** | +0.003 |
| HT29 | 0.485 | **0.490** | +0.017 |
| K562 | 0.494 | **0.493** | -0.006 |
| MCF7 | 0.490 | **0.499** | -0.002 |

**Doblar los genes no mueve `E`:** 0.485-0.495 a media mitad, 0.486-0.499 a
genes completos. El azar en las dos. Asi que la igualacion de profundidad, que
era obligatoria por la leccion de `b_contra_a.md`, **resulto no cambiar este
numero** -- pero eso solo se sabe por haberla hecho, y es `W` quien la
necesitaba, no `E`.

### 4.4 El `alpha` elegido: su desacuerdo es parte del resultado

| linea | **blancos donde `alpha_A != alpha_B`** |
|---|---|
| A549 | **70 %** |
| BXPC3 | **76 %** |
| HAP1 | **73 %** |
| HT29 | **64 %** |
| K562 | **64 %** |
| MCF7 | **70 %** |

Reparto global de los `alpha` elegidos (396 ajustes de media mitad, 6 lineas x 33
blancos x 2 mitades): `1`: 140, `10`: 93, `100`: 88, `1000`: 45, `10000`: 30. A
genes completos (198 ajustes): `1`: 69, `10`: 58, `100`: 41, `1000`: 20, `10000`: 10.

**Dos mitades del mismo blanco eligen `alpha` distinto dos veces de cada tres**,
y la discrepancia llega a cuatro ordenes de magnitud (§6: `AKT1` elige `1` en una
mitad y `10000` en la otra; `ETS1` igual; `STAT5A` igual). Es una medicion
independiente de lo que dice `W`: **el ajuste esta mal determinado con 3500
genes, hasta en cuanto hay que regularizarlo.**

---

## 5. 🔴 Lo que esto corrige de §7, que no es su conclusion sino su argumento

§7 escribio: «ninguno pasa del 64 % de consistencia, con la moneda en 50 %».
**La moneda no esta en 50 % para ese estadistico.** §7 definio la consistencia
como la fraccion de blancos con el signo **MAYORITARIO** (`max(pos, neg)`), y el
maximo de dos conteos complementarios esta sesgado hacia arriba. Con 33 blancos
y coeficientes de signo independiente:

```
E[ max(k, 33-k) / 33 ]  con  k ~ Binomial(33, 0.5)  =  0.5700   (exacto)
```

| estadistico | nulo correcto | lo que §7 uso |
|---|---|---|
| mayoritario `max(pos,neg)`, n = 33 (el de §7) | **0.570** | 0.50 |
| acuerdo de signo pareado (el de aqui, `W` y `E`) | **0.500** | -- |

**Asi que los "58 %" de §7 estaban EN el azar**, no 8 puntos por encima, y el
64 % estaba **+0.07** por encima, no +0.14. Recalculado aqui con la misma
definicion de §7 sobre ajustes de genes completos, de los 18 rasgos con peso
**10 quedan en 0.570 o por encima y 8 por debajo** -- el reparto de una moneda.
El rasgo mas consistente es `dist_emb` con 0.631, **+0.061** sobre su nulo.

**La conclusion de §7 no se cae: se refuerza, y por otra via.** §7 concluyo "no
hay mapa compartido" de un numero que no estaba por encima de su azar -- el
argumento no sostenia la conclusion. Lo que la sostiene es `W - E` > 0 en 6 de 6
con IC que excluye 0, y `T` = 0.47 por debajo de 0.50: eso si es una comparacion
con su nulo en el sitio correcto.

Por eso el encargo tenia razon en el motivo, aunque la aritmetica fuera otra:
con n = 33, el 64 % no excluia una moneda. **No la excluia por un margen todavia
menor del que el encargo calculo**, porque el nulo no era 50 % sino 57 %.

---

## 6. Las filas trabajadas

### 6.1 Los 18 rasgos contados, ordenados por `W`

Promedio sobre las seis lineas. `T` desatenuado por la formula de §3.1, y **a
`n/d` cuando `W <= 0.60`**, porque `T` divide por `2W-1` y con `W` cerca de 0.5
el cociente se desboca -- el mismo cuidado que §8.1 de `control_o_no.md` tuvo
con `r_mismo`.

| # | rasgo | **`W`** | rango `W` (6 lineas) | **`E`** | **`W - E`** | **`T`** | `r` perfil | `r` perfil (a=10) | mayorit. §7 | vs nulo 0.570 | `|coef|` | lectura |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | `r_parcial_10` | **0.864** | 0.788 - 0.909 | 0.509 | +0.354 | **0.513** | +0.847 | +0.924 | 0.606 | +0.036 | 0.0322 | **se estima**, difiere de verdad |
| 2 | `delta_r` | **0.808** | 0.758 - 0.909 | 0.492 | +0.316 | **0.488** | +0.755 | +0.793 | 0.571 | +0.001 | 0.0114 | **se estima**, difiere de verdad |
| 3 | `dist_emb` | **0.798** | 0.667 - 0.879 | 0.536 | +0.262 | **0.560** | +0.835 | +0.878 | 0.631 | +0.061 | 0.0195 | **se estima**, algo compartido |
| 4 | `r_pearson` | **0.778** | 0.727 - 0.818 | 0.493 | +0.285 | **0.487** | +0.798 | +0.870 | 0.591 | +0.021 | 0.0194 | **se estima**, difiere de verdad |
| 5 | `r_parcial_50` | **0.732** | 0.636 - 0.879 | 0.502 | +0.230 | **0.504** | +0.795 | +0.869 | 0.606 | +0.036 | 0.0165 | **se estima**, difiere de verdad |
| 6 | `mismo_modulo` | **0.732** | 0.636 - 0.848 | 0.503 | +0.229 | **0.507** | +0.859 | +0.851 | 0.606 | +0.036 | 0.0057 | **se estima**, difiere de verdad |
| 7 | `h_disp` | **0.722** | 0.667 - 0.818 | 0.497 | +0.225 | **0.493** | +0.535 | +0.628 | 0.530 | -0.040 | 0.0201 | **se estima**, difiere de verdad |
| 8 | `h_media` | **0.692** | 0.515 - 0.758 | 0.492 | +0.200 | **0.479** | +0.653 | +0.760 | 0.571 | +0.001 | 0.0107 | **se estima**, difiere de verdad |
| 9 | `h_var` | **0.667** | 0.545 - 0.788 | 0.509 | +0.157 | **0.528** | +0.665 | +0.776 | 0.606 | +0.036 | 0.0324 | **se estima**, difiere de verdad |
| 10 | `h_frac_cero` | **0.657** | 0.545 - 0.788 | 0.492 | +0.165 | **0.474** | +0.546 | +0.617 | 0.571 | +0.001 | 0.0413 | **se estima**, difiere de verdad |
| 11 | `h_rango_expr` | **0.646** | 0.515 - 0.788 | 0.498 | +0.148 | **0.494** | +0.491 | +0.499 | 0.566 | -0.004 | 0.0124 | debil, difiere de verdad |
| 12 | `cos_emb` | **0.641** | 0.545 - 0.788 | 0.496 | +0.145 | **0.487** | +0.546 | +0.614 | 0.556 | -0.014 | 0.0081 | debil, difiere de verdad |
| 13 | `r_spearman` | **0.626** | 0.515 - 0.727 | 0.499 | +0.128 | **0.495** | +0.499 | +0.467 | 0.540 | -0.030 | 0.0101 | debil, difiere de verdad |
| 14 | `misma_familia` ⚠ | **0.625** | 0.515 - 0.758 | 0.514 | +0.111 | **0.556** | +0.135 | +0.141 | 0.409 | -0.161 | 0.0011 | debil, algo compartido |
| 15 | `abs_r_pearson` | **0.606** | 0.485 - 0.727 | 0.489 | +0.117 | **0.446** | +0.445 | +0.489 | 0.566 | -0.004 | 0.0057 | debil, difiere de verdad |
| 16 | `cov` | **0.576** | 0.424 - 0.697 | 0.492 | +0.083 | n/d | +0.270 | +0.289 | 0.571 | +0.001 | 0.0068 | debil |
| 17 | `abs_r_spearman` | **0.540** | 0.455 - 0.697 | 0.505 | +0.035 | n/d | +0.069 | +0.033 | 0.556 | -0.014 | 0.0056 | **no se estima** |
| 18 | `rango_cor` | **0.510** | 0.394 - 0.697 | 0.492 | +0.018 | n/d | +0.341 | +0.346 | 0.530 | -0.040 | 0.0047 | **no se estima** |
| 19 | `g_media` | n/d | n/d | n/d | n/d | n/d | n/d | n/d | n/d | n/d | 0.0000 | cero exacto: FUERA del conteo |
| 20 | `g_var` | n/d | n/d | n/d | n/d | n/d | n/d | n/d | n/d | n/d | 0.0000 | cero exacto: FUERA del conteo |
| 21 | `g_frac_cero` | n/d | n/d | n/d | n/d | n/d | n/d | n/d | n/d | n/d | 0.0000 | cero exacto: FUERA del conteo |

**Lo que contesta la pregunta por rasgo, que era el encargo:**

* **Con `W` alto y `E` en el azar -- difieren de verdad:** `r_parcial_10`
  (0.864 / 0.509), `delta_r` (0.808 / 0.492), `r_pearson` (0.778 / 0.493),
  `r_parcial_50`, `mismo_modulo`, `h_disp`. Son los rasgos de coexpresion
  **parcial** y de modulo, no la covarianza marginal.
* **Con `W` en el azar -- simplemente no se estiman:** `rango_cor` (0.510),
  `abs_r_spearman` (0.540), y `cov` al borde (0.576). 🔴 **De estos tres no se
  puede decir nada sobre si difieren entre blancos**, y §7 los listaba con su
  consistencia como si fuera comparable a la de los demas.
* **Y el que mas peso tenia en §7 no es el que mejor se estima.**
  `h_frac_cero` era el `|coef|` mayor (0.0413) y tiene `W` = 0.657;
  `r_parcial_10` tiene un tercio de su peso y `W` = 0.864. **Peso y
  estimabilidad no son lo mismo**, y §7 ordenaba por peso.

⚠ **`misma_familia` lleva marca y la razon es concreta.** Es el **unico** rasgo
con coeficientes exactamente cero: **41 %** de los ajustes de media mitad (cuando
el blanco no tiene ningun gen de su familia en esa mitad, la columna es constante
y el ridge la manda a cero, igual que los de solo-`g`). El acuerdo de signo
pareado cuenta "los dos exactamente cero" como coincidencia, lo que distorsiona
sobre todo `E`. Medido en crudo daba `E` = **0.333** y `T` = **-0.420**, fuera
del rango [0,1] -- el modelo declarandose roto. **Las cifras de su fila son las
calculadas excluyendo las celdas con cero** (`W` 0.591 -> 0.625, `E` 0.333 ->
0.514). Los otros 17 rasgos tienen **0 %** de ceros exactos, asi que la
correccion no les toca: el agregado pasa de `W` 0.677 / `E` 0.491 a **`W` 0.679
/ `E` 0.501**, y con `E` en 0.501 el azar queda clavado.

🔴 **Y la columna `r` perfil se da, pero NO sostiene nada.** Es la correlacion
**sobre los 33 blancos** entre el coeficiente de mitad_A y el de mitad_B para un
rasgo: mide si la variacion blanco-a-blanco de ese coeficiente se reproduce, y
llega a +0.85 en `r_parcial_10` (+0.92 con `alpha` fijo, asi que no es artefacto
de la escala compartida). Parece evidencia de (i) y **no lo es**: una reescala
positiva por blanco la haria alta sin que el mapa difiera, y una reescala
positiva **no rompe la transferencia**, porque la metrica de §5 es invariante a
escala. Lo que cierra (ii) es el signo y `E` corr ~ 0, no esta columna.

### 6.2 Los 33 blancos de A549, uno por fila

`W` signo y `W` corr son del propio blanco (mitad_A contra mitad_B). `A_esp` es
la fiabilidad del delta especifico de ese blanco, de `control_o_no.md` §8.1, para
ver si `W` sigue a la fiabilidad del objetivo. Ordenado por `W`.

| # | blanco | gen expr. | genes mit.A | genes mit.B | `alpha_A` | `alpha_B` | difiere | **`W` signo** | **`W` corr** | `A_esp` |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | `HERC5` | si | 3584 | 3584 | 10000 | 10000 | no | **0.944** | +0.837 | +0.0025 |
| 2 | `TRAFD1` | si | 3583 | 3585 | 1 | 1 | no | **0.889** | +0.960 | +0.1061 |
| 3 | `TYK2` | si | 3583 | 3585 | 10 | 1 | si | **0.889** | +0.980 | +0.4975 |
| 4 | `BATF2` | si | 3583 | 3585 | 1 | 10 | si | **0.833** | +0.754 | +0.1011 |
| 5 | `IFNAR2` | si | 3583 | 3585 | 10 | 1 | si | **0.833** | +0.962 | +0.4998 |
| 6 | `JUN` | si | 3583 | 3585 | 1 | 1 | no | **0.833** | +0.852 | +0.1476 |
| 7 | `H1F0` | si | 3583 | 3585 | 1000 | 100 | si | **0.833** | +0.311 | +0.0420 |
| 8 | `DTX3L` | si | 3583 | 3585 | 1 | 10 | si | **0.833** | +0.910 | +0.0633 |
| 9 | `STAT3` | si | 3584 | 3584 | 10 | 10 | no | **0.778** | +0.946 | +0.2238 |
| 10 | `IFNAR1` | si | 3583 | 3585 | 10 | 1 | si | **0.722** | +0.944 | +0.3308 |
| 11 | `JAK1` | si | 3583 | 3585 | 10 | 100 | si | **0.722** | +0.717 | +0.1832 |
| 12 | `PARP12` | si | 3583 | 3585 | 10 | 1000 | si | **0.722** | +0.448 | +0.0013 |
| 13 | `IRF7` | si | 3583 | 3585 | 10 | 10 | no | **0.722** | +0.807 | +0.0683 |
| 14 | `IRF1` | si | 3583 | 3585 | 100 | 10 | si | **0.722** | +0.706 | +0.1091 |
| 15 | `CEBPG` | no | 3584 | 3585 | 10000 | 10000 | no | **0.667** | -0.083 | +0.0515 |
| 16 | `ELK1` | no | 3584 | 3585 | 10 | 10 | no | **0.667** | +0.944 | +0.0978 |
| 17 | `RNF114` | no | 3584 | 3585 | 100 | 100 | no | **0.667** | +0.619 | +0.0252 |
| 18 | `HERC6` | si | 3583 | 3585 | 1 | 100 | si | **0.667** | +0.840 | +0.0756 |
| 19 | `MAP3K14` | si | 3583 | 3585 | 10 | 1 | si | **0.667** | +0.583 | +0.1210 |
| 20 | `MAPK8` | si | 3584 | 3584 | 10 | 100 | si | **0.611** | +0.794 | +0.0380 |
| 21 | `ETS1` | si | 3584 | 3584 | 10 | 10000 | si | **0.611** | +0.092 | +0.0087 |
| 22 | `SOCS1` | no | 3584 | 3585 | 1 | 1000 | si | **0.611** | +0.372 | +0.0346 |
| 23 | `STAT1` | si | 3584 | 3584 | 1 | 10 | si | **0.611** | +0.209 | +0.0912 |
| 24 | `ZBP1` | no | 3584 | 3585 | 10 | 1000 | si | **0.611** | +0.221 | +0.0699 |
| 25 | `STAT5A` | no | 3584 | 3585 | 1 | 10000 | si | **0.556** | +0.121 | +0.0368 |
| 26 | `RAP1GAP` | no | 3584 | 3585 | 1 | 1 | no | **0.556** | +0.661 | +0.0290 |
| 27 | `ID2` | no | 3584 | 3585 | 10000 | 10000 | no | **0.556** | +0.611 | -0.0067 |
| 28 | `ETS2` | si | 3583 | 3585 | 10 | 100 | si | **0.556** | +0.795 | +0.1015 |
| 29 | `STAT4` | no | 3584 | 3585 | 10000 | 10 | si | **0.500** | -0.013 | +0.0218 |
| 30 | `BRD9` | si | 3583 | 3585 | 10 | 100 | si | **0.500** | +0.623 | +0.0322 |
| 31 | `SP100` | si | 3584 | 3584 | 1000 | 10 | si | **0.444** | +0.302 | +0.0399 |
| 32 | `AKT1` | si | 3584 | 3584 | 1 | 10000 | si | **0.444** | +0.053 | +0.0351 |
| 33 | `VAV1` | no | 3584 | 3585 | 1 | 100 | si | **0.389** | +0.353 | -0.0493 |

`W` medio **0.672**, mediana 0.667, minimo **0.389** (`VAV1`), maximo **0.944**
(`HERC5`); `W` corr medio +0.583.

**Tres cosas de esta tabla, y las tres importan:**

1. **El rango es enorme: de 0.389 a 0.944.** "El coeficiente de un blanco" no es
   una sola cosa: en `TYK2` e `IFNAR2` se reproduce casi del todo (0.889, 0.833,
   con `W` corr +0.98 y +0.96), y en `VAV1`, `AKT1` y `SP100` esta por DEBAJO del
   azar. No tiene sentido hablar del coeficiente por blanco sin decir de que
   blanco.
2. **`W` sigue a `A_esp`, y eso es coherente.** Los tres `A_esp` mayores
   (`IFNAR2` +0.500, `TYK2` +0.498, `IFNAR1` +0.331) tienen `W` de 0.833, 0.889 y
   0.722, en los puestos 5, 3 y 10. `VAV1`, el unico con `A_esp` negativo
   (-0.0493), es el ultimo. **Donde el objetivo es fiable el coeficiente se
   reproduce; donde el objetivo es ruido el coeficiente es ruido.** La excepcion
   es `HERC5` (primero en `W` con `A_esp` = +0.0025), y se explica: sus dos
   mitades eligen `alpha` = 10000, el maximo, asi que los dos vectores estan
   encogidos casi a cero y sus signos los fija lo poco que queda -- un acuerdo
   alto sobre coeficientes practicamente nulos. `CEBPG` e `ID2` son el mismo
   caso.
3. **El `alpha` discrepa en 23 de 33**, y los desacuerdos grandes caen donde `W`
   corr se hunde: `AKT1` (`1` contra `10000`, `W` corr +0.053), `ETS1` (`10`
   contra `10000`, +0.092), `STAT5A` (`1` contra `10000`, +0.121), `STAT4`
   (`10000` contra `10`, -0.013). **No es que el `alpha` estropee la medida: es
   que donde el ajuste no esta determinado, ni la regularizacion ni el
   coeficiente se reproducen.**

---

## 7. Veredicto

1. **La explicacion (ii) queda cerrada.** Los coeficientes no son "los mismos
   estimados con ruido": `W - E` = **+0.16 a +0.21** con IC 95 % que excluye 0 en
   **6 de 6** lineas, y en correlacion **+0.54 a +0.61** con 200/200 replicas
   positivas en las seis. Y sin ningun modelo: **`E` corr esta en cero**
   (-0.010 a +0.029).
2. **La explicacion (i) se sostiene.** Descontada la atenuacion por ruido, el
   acuerdo de signo VERDADERO entre blancos es **`T` = 0.467 a 0.487**, en el
   azar o por debajo, en las seis lineas. El mapa es de verdad distinto en cada
   blanco, y el cero de §5 no pide otra explicacion.
3. **Pero el coeficiente tampoco se estima bien dentro de un blanco.** `W` =
   **0.653 a 0.726**: una de cada tres comparaciones de signo falla con 3500
   genes, el flip por ruido es `e` = **0.16 a 0.22**, y el `alpha` elegido difiere
   entre mitades en el **64-76 %** de los blancos. De tres rasgos de §7
   (`rango_cor`, `abs_r_spearman`, `cov`) **no se puede decir nada**: su `W` esta
   en el azar.
4. **Y el argumento de §7 era invalido aunque su conclusion no.** El nulo del
   estadistico mayoritario con n = 33 es **0.570**, no 0.50: los "58 %" de §7
   estaban EN el azar y el 64 % a +0.07 de el. Lo que sostiene la conclusion es
   esta medicion, no aquella.

**La direccion que se sigue, con lo que cierra cada parte:**

| | |
|---|---|
| promediar entre MAS blancos | **cerrado por el punto 1.** Los coeficientes difieren de verdad; promediarlos promedia cosas distintas |
| buscar otra senal del control | ya estaba cerrado por `control_o_no.md` §10 |
| **un coeficiente por blanco** | **sigue abierto, y es lo unico que queda** -- pero el punto 3 le pone precio: con 3500 genes el coeficiente de un blanco se reproduce a 0.67, asi que estimarlo por blanco **desde sus propios genes** hereda ese ruido. Hace falta estructura externa que diga que regula a que por gen, que es la consecuencia que el encargo le asigno de antemano a este desenlace |

🔴 **Lo que esto NO dice.** No dice que exista una forma de fijar el coeficiente
por blanco sin ver su verdad -- todo lo de aqui es transductivo y usa la verdad
de Jiang. Dice que si tal forma existiera, el objeto que buscaria (un
coeficiente propio de cada blanco) **es real y no un artefacto del ruido**, y
que su estimacion desde ~3500 genes del propio blanco es, por si sola,
insuficiente.

---

## 8. Lo que no se hizo, y las entradas

| prohibicion del encargo | estado |
|---|---|
| construir un candidato, tocar arneses o bundles | **no**. No se escribio ningun bundle ni se llamo a ningun arnes; el script solo lee `rasgos_<linea>.npz` y resuelve sistemas de 21x21 |
| usar el lado perturbado para calcular algun rasgo | **no**. Los 22 rasgos son los ya calculados en `control_o_no_rasgos.py`, del pool de CONTROL y del simbolo del gen. No se calculo ningun rasgo nuevo |
| elegir el modelo o el `alpha` por el resultado | **no**. El ridge es el de §7 y el `alpha` sale de CV interna sobre genes **dentro de cada mitad**, sin ver `W` ni `E`. La version con `alpha` = 10 fijo se reporta entera al lado (§4.2), no se elige entre las dos |
| presentar alguna cifra como via o meterla en una submission | **no**. §0 lo declara arriba |
| `~/vcc/benchmark/`, `/datos/benchmark_VEDADO_techo_metricas/` | **no leidos** |
| envios, `vcc sample`, commits | **ninguno** |

| | |
|---|---|
| entradas | `modelo/out/control_o_no/rasgos_{A549,BXPC3,HAP1,HT29,K562,MCF7}.npz` (los rasgos, `E1`, `E2`, `A_esp`, `cruz` ya calculados por `control_o_no_rasgos.py`). `control_o_no.md` §7 y §8.2 como referencia a corregir |
| codigo | `modelo/src/coef_dentro_entre.py` |
| salidas | `modelo/out/control_o_no/coef_dentro_entre.json`, `modelo/out/control_o_no/coef_dentro_entre_run.log` |
| semillas | familia nueva **`20269000 + 100*fold + replica`**, con `fold` = indice de linea (0-5) y `replica` = 0 para el punto, 1-200 para el bootstrap. Gobierna la particion de genes y los pliegues de la CV interna |
| coste | 118-138 s por linea, 200 replicas enteras en las seis. No hubo que bajarlas |

**Lo que queda abierto:**

1. **`T` depende de un modelo de un parametro** (§3.1): misma tasa de flip dentro
   y entre blancos, e independencia entre mitades. El resultado de §1 (`W - E` >
   0 con IC) **no** depende de el; la cifra 0.47 si.
2. **El bootstrap remuestrea blancos y reparte genes de nuevo, pero el reparto
   del punto es uno solo** (replica 0). Las medias bootstrap de `W` (0.659 a
   0.704) y de `E` (0.489 a 0.497) estan al lado de las del punto, asi que el
   reparto no mueve el titular.
3. **`W` se mide con el gen blanco fuera, igual que §7.** Lo que la version con
   el gen blanco daria no se midio, y por §4.3 de `control_o_no.md` no deberia
   leerse.
4. **No se midio `W` a genes completos.** `W` a profundidad completa seria mayor
   que 0.67 y diria cuanto de la irreproducibilidad es solo tamano de muestra;
   `E` a genes completos (§4.3) no se movio, pero de eso no se sigue que `W` no
   se mueva. Es el hueco declarado mas util de esta medicion.
