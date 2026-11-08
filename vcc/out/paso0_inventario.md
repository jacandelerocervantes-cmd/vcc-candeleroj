# Paso 0 — Inventario

**Fecha:** 2026-10-27 (día 1 del calendario de `PLAN.md`: «entorno, descarga, T13»).
**Alcance:** solo lectura — nombres, tamaños y metadatos de archivo. No se abrió ningún
`.h5ad`/`.rds`, no se escribió modelo, no se envió nada.

---

## 1 · Entorno

| ítem | estado |
|---|---|
| `~/cell-eval2` | clon git, commit `d86be3f4` (23 oct 2026) — **Import cell_eval2 0.20.0** desde `ArcInstitute/cell_eval2@a63231bb` |
| `~/vcc/mf.sh` | instalador de Miniforge (shell + payload binario embebido, 128 MB) — presente, no ejecutado en este paso |
| `~/vcc/src/` | `cpm_stats.py`, `pseudobulk.py` — los dos módulos ya copiados y auditados (cero valores calibrados; ver `registro/CONSTANTES_A_VERIFICAR.md` §«Lo que sí se copió») |
| `~/vcc/workflow/` | solo `vigilantes/sync_out.sh` y `vigilantes/spot_watch.sh` — snakemake todavía no tiene reglas |
| `~/vcc/corpus/` | vacío — `particion.md` e `inventario.parquet` (§7 del plan) **no existen aún** |
| `~/vcc/out/` | vacío hasta este archivo |

⚠️ **Desfase de versión:** el clon trae `cell-eval2` 0.20.0; los valores de referencia
publicados (`docs/vcc2026_metrics/vcc2026-metrics-brief.md`) se midieron con 0.15.0.
Registrado en `registro/limitaciones.md` §1 — no afecta a los números publicados, sí a
si este clon los reproduce bit a bit.

---

## 2 · Datos en `/datos` (NVMe local, efímero)

### 2.1 · `/datos/reto` — el panel de la competencia

| partición | contextos | archivos | tamaño |
|---|---|---|---|
| `val` | A, B, C | `context_{A,B,C}.h5ad`, `pert_counts.csv`, `gene_names.csv`, `manifest.json`, `vcc_2026_controls.zip` | 1.3 G |
| `final` | D, E, F | `context_{D,E,F}.h5ad`, `pert_counts.csv`, `gene_names.csv`, `manifest.json`, `vcc_2026_controls_final.zip` | 1.2 G |

De `manifest.json` (ambas particiones, sin abrir ningún `h5ad`):

- `n_genes = 18533`, `n_constructs = 300`, `cells_per_pert = 400`.
- Por contexto: `n_perturbations = 300`, `control_cells = 18400`,
  `ground_truth_cells = 138400`, **`n_ntc_ids = 46`** — igual en los seis contextos
  (A–F).
- `pert_col = target_gene`, `context_col = context`, `control_label = non-targeting`.

📌 **Dato nuevo para `PLAN.md` §8.1** («46 guías NTC (documento) contra 100 (el
reto)»): el manifest **del propio reto** declara `n_ntc_ids = 46` en los seis
contextos — coincide con el lado «documento» del cotejo, no con el «100». Falta el
cotejo contra el archivo GEO/SRA; no se cierra aquí.

### 2.2 · `/datos/ref` — fuentes de referencia (143 G total)

| fuente | tamaño | contenido | formato |
|---|---|---|---|
| `replogle/` | 80 G | K562 GWPS, K562 essential, RPE1 — sc y bulk | `h5ad` (6 archivos) |
| `otros/` | 41 G | Zenodo 14518762 (ver abajo) + GSE225775 + GSE139944 | `rds` ×9, `tar` ×2, `zip`, `pdf`, `txt` |
| `nadig/` | 14 G | Jurkat, HepG2 — sc crudo | `h5ad` |
| `orion/` | 4.9 G | HCT116, HEK293T | `parquet` |
| `feng/` | 2.7 G | Targeted/genome-wide screen, LFC por gen (y por línea) | `tsv.gz` |
| `tfact/` | 82 M | `TFactprofiler_all.csv` | `csv` |

`dl.log` muestra `rc=0` (descarga exitosa) para **las 31 entradas** de `/datos/ref`.
Ninguna descarga fallida registrada.

**`otros/` identificado:** `A_readme.txt` confirma que todo el directorio es el
depósito de Zenodo 14518762 del manuscrito «Systematic reconstruction of molecular
pathway signatures using scalable single-cell perturbation screens» — 5 objetos
Seurat por pathway (IFNB, IFNG, TNFA, TGFB1, INS), 3 genelists, un zip de resultados
DE, el protocolo de captura de guías, y `Bulk_RNAseq_Seurat_object_IFNG_and_TGFB_stim.rds`
(3.3 M, **confirmado descargado** — ver `registro/limitaciones.md` §5, cerrada hoy).

### 2.3 · Lo que requiere conversión antes de entrar al pipeline

`GSE225775_RAW.tar`, `GSE139944_RAW.tar` y los 9 `.rds` de `otros/` (Zenodo 14518762)
no son `h5ad` — requieren desempaquetado y/o un puente Seurat→AnnData antes de que
`cell-eval2` los pueda leer. Registrado en `registro/limitaciones.md` §4 (abierta,
no cerrada: es trabajo pendiente, no un bloqueo de disponibilidad).

---

## 3 · Estado de `registro/`

| archivo | estado |
|---|---|
| `CONSTANTES_A_VERIFICAR.md` | 14 constantes listadas; **#6 (span de `pds`) y #11 (desempate de `reach`) cerradas hoy contra el código** (`src/cell_eval2/scales.py:428`, `metrics/direction.py:508-511`, `docs/vcc2026_metrics/vcc2026-metrics-brief.md:341-346`) — el resto sigue pendiente |
| `limitaciones.md` | 5 entradas (desfase de versión de `cell-eval2` · potencia 80 % del MDE · α del intervalo conforme · datasets Seurat/RDS+tar · el `.rds` de 3.4 MB) — **§5 resuelta hoy mismo** en este paso 0 |

---

## 4 · Lo que sigue (calendario de `PLAN.md`, día 1→2)

- [ ] **T13** — contar réplicas por (línea, perturbación): separar pozo/transducción
  aparte de mismo-pozo-resecuenciado; guías distintas al mismo gen **no** son réplica.
  Aún no iniciado — requiere leer las capas de los `h5ad`/`rds`, fuera del alcance de
  este paso 0 (solo lectura de metadatos de archivo).
- [ ] **T3** — matriz de incidencia (fuente, línea, perturbación) → `corpus/inventario.parquet`.
  No existe todavía; este reporte es el inventario de **archivos**, no la matriz de
  incidencia fina que T3 pide.
- [ ] `corpus/particion.md` (corpus vs. reservado) — no existe todavía.
- [ ] Cotejo de `PLAN.md` §8.1: procedencia de `ref/otros/` ya resuelta arriba (es
  Zenodo 14518762); falta el cotejo de 46 vs. 100 NTC contra el archivo GEO/SRA.
