#!/usr/bin/env python
"""CANDIDATO 01 -- respuesta lineal a lo largo de la covarianza del control.

LA PREGUNTA QUE INTENTA CONTESTAR. El extremo 0 de la escala es un comparador
oraculo que acota la TRIVIALIDAD: `baseline.py:6-10` dice que un puntaje igual o
por debajo "is not evidence of target-specific skill". El `b` emite el control
escalado a UN perfil de respuesta promedio, el MISMO para las 33
perturbaciones. El paso minimo por encima es un perfil POR BLANCO.

LA ESTRUCTURA. Para el blanco `g`, los genes que covarian con `g` EN EL POOL DE
CONTROL de ese contexto son los candidatos a moverse cuando `g` se perturba. Si
`g` baja y `h` covaria positivamente con `g`, `h` baja; si covaria
negativamente, `h` sube. La pendiente de esa respuesta se estima por regresion
DENTRO del mismo pool.

    y     = log1p(cuentas normalizadas por celula)    <- quita el confusor de
                                                         tamano de libreria
    rho_h = corr(y_h, y_g)  en el pool de control     <- orden, signo y magnitud
    r_h   = 1 - rho_h                                 <- factor multiplicativo

`r_h = 1 - rho_h` es lo unico que se aplica, y es exactamente "el orden y el
signo, nunca el valor": la unica cantidad que entra es la correlacion en
contexto. Es adimensional y **acotada por construccion** -- `rho` vive en
`[-1, 1]`, asi que `r` vive en `[0, 2]`. Un gen perfectamente coexpresado con el
blanco se apaga del todo; uno incorrelado se queda igual; uno anticorrelado
sube. No hay ninguna constante elegida aqui.

🔴 DOS VERSIONES ANTERIORES, MEDIDAS Y DESCARTADAS ANTES DE PUNTUAR NADA. Las
dos se descartaron por una propiedad medida de la EMISION, nunca por su `s`:
ninguna de las dos se puntuo.

  (v1) pendiente de regresion extrapolada a apagar `g` del todo:
       `beta_h = cov/var(y_g)`, `Delta_y_h = -beta_h * mean(y_g)`.
       Medido en A549: manda a CERO EXACTO a casi todos los genes que responden
       -- `r_min = 0.0000` con `r_max = 1.0000` en 24 de los 33 blancos, con
       cientos de genes por blanco (RAP1GAP 381, ZBP1 335, STAT4 332) a cero.
       `beta_h = rho * sd_h/sd_g` crece sin cota cuando `sd_h >> sd_g`.

  (v2) esperanza condicional a UNA `sd` de bajada del blanco:
       `Delta_y_h = -rho_h * sd(y_h)`.
       Medido en A549: sigue saturando -- la MEDIANA de los genes que responden
       va a cero exacto, 18 de 33 blancos quedan con `r_min = 0, r_max = 1`, y
       el pseudobulk muestra `|log2FC|` de hasta **17.4**, que es el artefacto de
       dividir por un cero. La causa es medible: para un gen disperso
       `sd(y_h) > mean(y_h)` -- con expresion en una fraccion `p` de celulas,
       `mean/sd ~ sqrt(p/(1-p))`, que para `p = 0.05` vale 0.23 -- asi que
       cualquier `|rho| > 0.23` lo cruza por debajo de cero. Un paso de una `sd`
       es la unidad natural para una conjunta normal, y log1p de cuentas
       dispersas no lo es: esta inflado de ceros.

LA MAGNITUD SIGUE SIENDO EL HUECO DECLARADO. Lo que esta estructura transfiere
con fundamento es el ORDEN y el SIGNO. Que la fraccion de cambio sea `rho` es
una eleccion de escala adimensional y acotada, no una cifra medida de cuanto
baja un gen cuando se perturba su regulador -- eso exigiria el lado real
perturbado, que esta prohibido. Queda como hueco, no como resultado, y es la
primera cosa que una iteracion deberia atacar (con su propia medicion, como
candidato nuevo).

POR QUE NORMALIZAR ANTES DE LA COVARIANZA, y no es un detalle. Sobre cuentas
crudas la covarianza entre dos genes cualesquiera esta dominada por el tamano de
libreria de la celula: una celula con libreria grande tiene cuentas altas de
TODO, asi que todos los genes correlacionan positivamente con todos. Usar eso
daria una respuesta casi uniforme y positiva -- que es exactamente el perfil
generico del `b` que se quiere superar. La normalizacion por celula es lo que
hace que la senal medida sea coexpresion y no profundidad.

LAS DOS RESTRICCIONES DURAS, y como se cumplen:
  1. NO TRANSDUCTIVO. El arnes lo impone estructuralmente: `predict_fn` recibe
     `(control_pool_X, labels, rng)` y nada mas (`suelo_arnes.py:176`). No ve la
     identidad de la linea ni una sola celula perturbada. `var_names` se pasa por
     clausura porque sin el no se puede saber QUE columna es el gen `g` -- son
     los nombres de las columnas del dataset, no el lado real perturbado. Se
     declara explicitamente.
  2. EL GEN BLANCO ESTA EXCLUIDO de la comparacion DE en los dos lados (#172).
     `beta_g = 1` por construccion, asi que `g` se apaga; eso NO cobra nada y se
     deja solo por coherencia biologica. Toda la senal puntuable esta en OTROS
     genes.

LO QUE SE ESTIMA EN CONTEXTO Y LO QUE ES CONVENCION, declarado sin esconderlo:
  * la correlacion `rho_h` y la desviacion `sd(y_h)`: estimadas del propio pool
    de control. No importadas.
  * el total de normalizacion: la MEDIANA de las librerias del propio pool.
  * el umbral de seleccion: se pone `beta_h = 0` donde la correlacion en
    contexto no pasa Benjamini-Hochberg a 0.05. **Ese 0.05 no se eligio aqui**:
    es el `de.p_adj_threshold` del propio preset vcc2026
    (`configs/vcc2026.yaml`), la misma convencion con la que el marcador decide
    que un gen es DE. No se probo ningun otro valor.
  * NO hay ninguna cifra importada: nada de Replogle, nada del benchmark,
    ninguna tabla de pesos, y ninguna cifra del `b` ni del anchor (que son
    artefactos de PUNTUACION).

EL CONTRATO DE EMISION -- conteos ENTEROS muestreados, nunca un perfil esperado
suave (`paso1_contrato.md`, HARD_CAP_NNZ). Se parte de celulas de control
remuestreadas CON REEMPLAZO, igual que el suelo, y se reescala gen a gen
manteniendo los enteros:

    r_h <= 1:  nuevo_h ~ Binomial(x_h, r_h)            (adelgazamiento)
    r_h  > 1:  nuevo_h = x_h + Poisson(x_h * (r_h-1))  (engorde)

Las dos ramas tienen esperanza exactamente `x_h * r_h` y las dos devuelven
enteros. El adelgazamiento puede crear ceros; el engorde actua solo donde la
celula elegida ya tiene cuentas, asi que un gen que vale 0 en esa celula se
queda en 0. Eso acota el `nnz` por arriba por el del suelo -- es una propiedad,
no un ajuste, y queda declarada.
"""
import numpy as np
import scipy.sparse as sp

N_PER_LABEL = 400
CONTROL_LABEL = "NT"

# el 0.05 del propio preset vcc2026 (`de.p_adj_threshold`), no una eleccion de aqui
FDR_Q = 0.05


def _bh_mascara(p, q=FDR_Q):
    """Benjamini-Hochberg: devuelve la mascara de los que pasan a nivel `q`."""
    n = p.size
    orden = np.argsort(p, kind="stable")
    p_ord = p[orden]
    umbral = q * np.arange(1, n + 1) / n
    pasa = p_ord <= umbral
    if not pasa.any():
        return np.zeros(n, dtype=bool)
    k = np.max(np.nonzero(pasa)[0])
    corte = p_ord[k]
    return p <= corte


def construye_respuesta(control_pool_X, var_names, targets):
    """Devuelve (ratios, diag). `ratios[g]` es un vector (n_genes,) de factores
    multiplicativos `r_h` para el blanco `g`. NO usa rng: la estructura del
    candidato es deterministica dado el pool de control, y solo el muestreo de
    la emision depende de la semilla."""
    from scipy.stats import norm as _norm

    X = control_pool_X.tocsr().astype(np.float64)
    n_cel, n_genes = X.shape
    libs = np.asarray(X.sum(axis=1)).ravel()
    total = float(np.median(libs[libs > 0])) if (libs > 0).any() else 1.0

    # normalizar por celula al total mediano del PROPIO pool, y log1p
    escala = np.zeros(n_cel)
    nz = libs > 0
    escala[nz] = total / libs[nz]
    Y = sp.diags(escala) @ X
    Y = Y.toarray()
    np.log1p(Y, out=Y)

    ybar = Y.mean(axis=0)
    Yc = Y - ybar
    var = (Yc * Yc).sum(axis=0) / max(n_cel - 1, 1)
    sd = np.sqrt(var)

    idx = {g: i for i, g in enumerate(var_names)}
    ratios = {}
    diag = {}
    for g in targets:
        ig = idx[g]
        vg = var[ig]
        if not np.isfinite(vg) or vg <= 0:
            # el blanco no se expresa en este pool: no hay nada que propagar
            ratios[g] = np.ones(n_genes)
            diag[g] = {"var_blanco": float(vg), "n_responde": 0,
                       "motivo": "blanco sin varianza en el pool de control"}
            continue
        cov = (Yc * Yc[:, ig:ig + 1]).sum(axis=0) / max(n_cel - 1, 1)

        # correlacion y su p en contexto (Fisher z), luego BH al 0.05 del preset
        with np.errstate(divide="ignore", invalid="ignore"):
            r = cov / (sd * sd[ig])
        r = np.nan_to_num(r, nan=0.0, posinf=0.0, neginf=0.0)
        r = np.clip(r, -0.999999, 0.999999)
        if n_cel > 3:
            z = np.arctanh(r) * np.sqrt(n_cel - 3.0)
            p = 2.0 * _norm.sf(np.abs(z))
        else:
            p = np.ones(n_genes)
        p[~np.isfinite(p)] = 1.0
        p[var <= 0] = 1.0
        pasa = _bh_mascara(p)
        pasa[ig] = True                      # el propio blanco siempre responde

        # r_h = 1 - rho_h, acotado en [0, 2] por construccion.
        rho_sel = np.where(pasa, r, 0.0)
        rat = 1.0 - rho_sel
        rat = np.clip(rat, 0.0, 2.0)
        rat[var <= 0] = 1.0
        ratios[g] = rat

        suben = int(np.sum((rat > 1.0) & pasa))
        bajan = int(np.sum((rat < 1.0) & pasa))
        diag[g] = {
            "var_blanco": float(vg),
            "media_log_blanco": float(ybar[ig]),
            "n_responde": int(pasa.sum()) - 1,
            "n_suben": suben, "n_bajan": bajan,
            "r_min": float(rat.min()), "r_max": float(rat.max()),
            "r_mediano_de_los_que_responden": float(np.median(rat[pasa]))
                if pasa.any() else 1.0,
            "corr_abs_max": float(np.abs(r[pasa & (np.arange(n_genes) != ig)]).max())
                if (pasa & (np.arange(n_genes) != ig)).any() else 0.0,
            "p_corte_bh": float(p[pasa].max()) if pasa.any() else 0.0,
        }
    diag["_pool"] = {"n_celulas_control": int(n_cel), "n_genes": int(n_genes),
                     "total_normalizacion_mediano": total,
                     "n_genes_con_varianza": int((var > 0).sum())}
    return ratios, diag


def _reescala_entero(bloque, rat, rng):
    """Reescala un bloque CSR de cuentas enteras por el factor por gen `rat`,
    devolviendo ENTEROS. Esperanza = x * rat en las dos ramas."""
    b = bloque.tocsr().astype(np.int64)
    x = b.data
    col = b.indices
    f = rat[col]
    out = np.zeros_like(x)

    baja = f <= 1.0
    if baja.any():
        out[baja] = rng.binomial(x[baja], np.clip(f[baja], 0.0, 1.0))
    sube = ~baja
    if sube.any():
        extra = rng.poisson(x[sube] * (f[sube] - 1.0))
        out[sube] = x[sube] + extra

    nb = sp.csr_matrix((out, b.indices.copy(), b.indptr.copy()), shape=b.shape)
    nb.eliminate_zeros()
    return nb


def hace_predict_fn(var_names, registro=None):
    """Devuelve el `predict_fn` que el arnes acepta sin modificarlo.
    `var_names` entra por clausura porque la firma del arnes no lo lleva y sin
    el no se puede localizar la columna del gen blanco. No es informacion del
    lado real perturbado: son los nombres de las columnas."""

    cache = {}

    def predict_cand01(control_pool_X, labels, rng):
        clave = (control_pool_X.shape, int(control_pool_X.sum()))
        if clave not in cache:
            targets = [l for l in labels if l != CONTROL_LABEL]
            cache[clave] = construye_respuesta(control_pool_X, var_names, targets)
        ratios, diag = cache[clave]
        if registro is not None:
            registro["diag"] = diag

        n_control = control_pool_X.shape[0]
        out = {}
        for label in labels:
            chosen = rng.integers(0, n_control, size=N_PER_LABEL)
            bloque = control_pool_X[chosen]
            if label == CONTROL_LABEL:
                out[label] = bloque.tocsr()      # el control predicho = el suelo
            else:
                out[label] = _reescala_entero(bloque, ratios[label], rng)
        return out

    return predict_cand01
