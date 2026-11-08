# PLAN — Virtual Cell Challenge 2026

**Cierre: 5 de noviembre, 23:59 UTC.** Hoy es 26 de octubre → **10 días.**
Documento maestro: `https://claude.ai/artifact/JjLaaEyzc9xB2FP3ab97y3` (rev 197)

---

## 0 · Reglas que no se negocian

1. **Nada se envía sin orden explícita.** El envío lo decide el autor, no el plan.
2. **No se usan métricas ni datos de modelos anteriores.** Ni el 84 %, ni el 7×, ni el 0.80. Si un
   benchmark publica un valor, ése es la referencia; lo que medimos es nuestra fracción de él (T5).
3. **No se explota la función de puntuación** (D-023). La unión vacía da `jac = 1.0`; no se usa.
4. **Lo que se elimine, se elimina su referencia.** Sin notas de hueco.
5. **El registro de limitaciones se llena el día en que la limitación aparece**, no al cierre (§6.9).
   Un registro escrito al final conserva solo lo que salió bien.

## 0.1 · Qué transfiere, y qué no

> A D, E y F pasa **la forma en que trabaja el modelo**.
> **No** los datos — no los tienen. **No** el error — es de cada contexto. **Los rangos, solo parcialmente.**

De ahí se sigue todo lo demás: el método **no necesita que los seis brazos sean buenos**. Necesita que
la forma sea la misma en todos y que **el benchmark la haya validado uno a uno**.

---

## 1 · La cadena crítica

```
T13  ->  X3  ->  X1        (no se puede paralelizar)
 |        |       |
 |        |       +-- cuál es el ciego: lo identifica el cribado, no una preferencia
 |        +---------- el arnés: identificar la configuración que reproduce lo publicado
 +------------------- contar réplicas por (línea, perturbación)
```

Todo lo demás cuelga de esto. **Si T13 no está, nada de lo que sigue tiene base.**

---

## 2 · Las cinco mediciones

| # | qué | cómo, en concreto |
|---|---|---|
| **T13** | réplicas por (línea, perturbación) | separar **(a)** pozo/transducción aparte — la que sirve para T1 — de **(b)** mismo pozo resecuenciado. **(c) guías distintas al mismo gen NO es réplica: es benchmark.** Más la columna: qué líneas aparecen en más de una fuente independiente |
| **T3** | cobertura | **el primitivo es la matriz de incidencia (fuente, línea, perturbación)**, al grano fino y sin colapsar. Toda cifra de cobertura es una **consulta** sobre ella → el error de «medido contra el panel equivocado» deja de ser posible |
| **T12** | tres niveles de varianza | entre células · entre réplicas · **entre preparaciones**, vía líneas repetidas en fuentes independientes. Confirmado desde afuera: la profundidad de silenciamiento la fija el nivel de sgRNA (MOI), no la biología del gen |
| **T11** | 11–14 % / 8,248 / 569 | **no son el mismo problema.** Vienen de contextos distintos → **tres validaciones independientes**, cada una en el suyo. No se reconcilian |
| **T5** | valores de referencia | los da **el benchmark publicado**. No se re-miden cifras con el pipeline viejo |

**Y los cálculos van después, sobre una base de regresión:** cuando no se conoce la ecuación pero sí los
dos lados, se recupera ajustando. Es un solo método aplicado tres veces — la descomposición
`μ + α_c + β_p + γ`, la regla de adopción de D5, y la estimación de `n_real`.

---

## 3 · El orden de construcción

```
0. EL BENCHMARK   <- la puerta
   -> tiene que dar UNO A UNO para que sea real
   -> y ES ALCANZABLE: tenemos los datos de laboratorio, así que dar uno a uno
      = RECUPERAR la configuración con la que se produjo el resultado publicado
   -> eso es identificar el instrumento, no ajustar un modelo

1. A, con lo que se tiene       -> residuo grande aquí = falla el PROCEDIMIENTO
2. B, sobre lo ya construido    -> residuo = coste de SU vista
3. C, sobre lo ya construido    -> idem

4. A D, E y F pasa LA FORMA
```

### 3.1 · La frontera del ajuste del arnés

Si se mueven suficientes perillas se acierta **cualquier** número sin haber recuperado nada — y
**el caso malo se ve idéntico al bueno**: también da uno a uno. Cuatro guardas:

- [ ] Las perillas son **convenciones, no biología**. Normalización, pseudobulk, pool de control,
      filtros: sí. Tamaños de efecto, dispersión, parámetros por gen: **no**.
- [ ] **Una sola configuración para todas las cifras publicadas.** Si hace falta una convención
      distinta para cada una, no se recuperó el pipeline: se acertaron números por separado.
- [ ] **Capa 2 reservada del ajuste**: ajustar contra lo publicado, luego **predecir las réplicas
      reservadas**, que no entraron.
- [ ] **Contar los grados de libertad.** Si las perillas son tantas como las cifras que se igualan,
      el acierto no informa.

---

## 4 · Calendario

| día | fecha | qué | hito |
|---|---|---|---|
| 1 | 27 oct | entorno, descarga, **T13** | réplicas contadas |
| 2 | 28 oct | **T3** (matriz de incidencia), cribado de admisibilidad | inventario = instrumento de partición |
| 3 | 29 oct | **X3**: barrer convenciones hasta el uno a uno | — |
| 4 | 30 oct | **PUERTA: ¿da uno a uno?** + capa 2 sobre réplicas reservadas | **go / no-go** |
| 5–6 | 31 oct – 1 nov | **A** completo con el procedimiento validado | residuo de A |
| 7 | 2 nov | **B** sobre lo construido | — |
| 8 | 3 nov | **C** sobre lo construido | intervalo de los tres |
| 9 | 4 nov | **D, E, F**: aplicar la forma · `n_real` · calibración | predicción completa |
| 10 | 5 nov | verificación, formato, **envío solo con orden explícita** | — |

**Regla de deslizamiento:** si un día se pasa, **se recorta alcance, no se recorta la puerta del día 4.**
Un número producido con un arnés que no da uno a uno no significa nada.

---

## 5 · El piso: la póliza de seguro

**Aunque todo lo demás falle, esto tiene `s > 0` garantizado y el campo no lo está usando.**

`nmae` es una suma de valores absolutos, así que el óptimo de `Σ|a·x − y|` **no es la pendiente de
mínimos cuadrados**: es la **mediana ponderada de los cocientes y/x con pesos |x|**.

- Sin ruido, NMAE = |a−1|.
- **Con ruido, `a` = 1 sin calibrar puede quedar BAJO el baseline** (medido: NMAE 1.2035, s = −0.339).
- **El óptimo L1 da s > 0 siempre** (medido: a ≈ 0.374 → NMAE 0.829 → **s = +0.285**).

→ No es «cualquier encogimiento sirve». Es **hay que calibrar, y calibrar está garantizado que ayuda.**

**Y la segunda palanca, independiente:** emitir con dispersión, con diana **n_pred ≈ n_real**.
La emisión degenerada no toca la magnitud (el `lfc` es idéntico hasta el último decimal, verificado),
pero mueve el **conjunto de llamadas** — el p de Wilcoxon se vuelve ~15 órdenes más extremo — y con eso
`fid`, `reach` y `jac`. Son las tres que **leen el conjunto de llamadas**, no «las tres de mayor span».

**Orden de palanca** (span menor = más palanca): **fid > jac > pds > nmae > reach > mse.**

---

## 6 · Disciplina de la máquina

Spot `r5ad.2xlarge`, us-east-2. **Interrupción = 20 minutos**, no un día. Pero solo si:

- [ ] `/datos` (NVMe local) = **solo datasets**. Se rebajan. Nada irreemplazable aquí.
- [ ] `out/` → **`aws s3 sync` cada 2 min** a S3, con versionado activado.
- [ ] Vigilante del aviso de interrupción corriendo (metadata **V2, con token**).
- [ ] *Interruption behavior: **stop*** → el volumen raíz sobrevive, así que **conda y los scripts no
      se reinstalan**. Solo `/datos` se rehace.

**Lo que NO se regenera descargando** — y por tanto va a S3 en cuanto exista:

1. La **matriz de incidencia** de T13.
2. La tabla del **MDE** por gen y contexto.
3. **La configuración del arnés** una vez identificada. ← la más cara de todas: no se recupera
   bajando datos, se recupera repitiendo el trabajo.
4. La partición declarada **corpus / reservado**.

---

## 7 · Estructura

```
vcc/
  PLAN.md              <- este archivo
  registro/            <- §6.9: se llena el día que aparece la limitación
    limitaciones.md       (qué falta · por qué · qué costaría cerrarlo · afecta resultado o
                           interpretación · fecha)
  corpus/
    particion.md       <- declarado: corpus vs reservado. Previene la fuga por declaración,
                          no por auditoría posterior
    inventario.parquet <- la matriz de incidencia (fuente, línea, perturbación)
  src/
  workflow/            <- snakemake: reanuda por existencia de archivos
  out/                 <- sincronizado a S3
  /datos               <- NVMe local, efímero, solo datasets
```

---

## 8 · Lo que queda abierto por diseño, no por olvido

- **D5** — criterio de adopción: no se fija antes de ver resultados. Instrumento definido: regresión
  lineal sobre diferencias con y sin el componente, y **decide hacia dónde apunta la pendiente, no un corte.**
- **T9** — parámetros del gemelo: **cotas, no estimaciones**, tomadas en su extremo conservador porque
  el error es asimétrico. Una cota holgada no cuesta nada; una apretada de más hace que D7 proyecte
  **lejos** de la verdad.
- **E** — es el brazo con menos apoyo (soporte de ML más delgado, ninguna vista que cubra el eje glial,
  la peor de las tres vistas de p53, única fila con línea única sin candidatas). Eso dice **dónde
  esperar más residuo**, no un orden que haya que acertar.

## 8.1 · Dos cotejos contra el archivo, en cuanto baje

- [ ] Procedencia de `ref/proteinturnover/` — si es marcaje isotópico por masas, asciende a pieza del
      lazo RNA → proteína.
- [ ] **46 guías NTC (documento) contra 100 (el reto).** Por la regla de T11 **no se reconcilian**:
      se validan por separado. Puede ser lo liberado en el archivo contra lo del experimento.
