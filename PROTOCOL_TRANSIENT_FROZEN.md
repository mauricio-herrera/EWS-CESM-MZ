# PROTOCOLO CONGELADO — Experimento T: ramas a forzamiento fijo y escenarios transitorios (Saddle-Node)

**Fecha de congelamiento:** 2026-07-20. **Estado de ceguera:** solo se listaron
nombres de archivos del repositorio `RenevanWesten/AMOC-Saddle-Node` (Zenodo del
paper JGR-Oceans 2025); ningún valor de ninguna serie ha sido descargado ni
inspeccionado. Timestamp público (tag + Zenodo) requerido antes de ejecutar
`run_transient_experiment.py`.

## 1. Lecciones aplicadas (del caso de estudio previo)

(i) El instrumento es un medidor de proximidad; la hipótesis primaria es de
**discriminación por proximidad**, no de cronómetro. (ii) El umbral supercrítico no
transfiere entre regímenes: aquí se calibra en un **control del mismo experimento**
(rama CESM_0600). (iii) La verdad de terreno primaria es una **convención de onset**
prespecificada (no la regresión de quiebre, que fecha la fase rápida). (iv) La
certificación se define a nivel de refit desde el inicio.

## 2. Datos (fijos)

Repositorio `RenevanWesten/AMOC-Saddle-Node` (CC-BY; van Westen et al., 2025, JGR
Oceans 130(8):e2025JC022651). Series escalares anuales
`AMOC_transport_depth_0-1000m*.nc` (y `FOV_section_34S*.nc` como descriptivo):

- **Control:** `CESM_0600/Ocean/AMOC_transport_depth_0-1000m_branch_QE.nc`
  (rama a forzamiento fijo lejos del fold).
- **Runs de prueba a forzamiento fijo:** ramas `CESM_1500, CESM_1550, CESM_1600,
  CESM_1650, CESM_1700` (archivo `branch_QE` de cada una), ancladas crecientemente
  cerca del fold del QE (año de anclaje = nombre).
- **Runs transitorios (velocidad finita):** `CESM_0600_climate_change` y
  `CESM_1500_climate_change`, escenarios RCP4.5 y RCP8.5 (4 series).

Longitudes y desenlaces (colapsa o no) son desconocidos al congelar; el protocolo
los trata algorítmicamente.

## 3. Operador y estadístico (fijos; variante escalar d=1)

Estado: índice AMOC estandarizado causalmente dentro de la ventana [inicio del run,
t₀]. Ensemble: las 24 configuraciones exactas de `ews_analysis.py` (TAU_SETS ×
OSC_GRID × ALPHA_GRID) en su versión d=1 (mismas ecuaciones con q ∈ R¹). Covariable
exógena u: en ramas de forzamiento fijo, el tiempo normalizado causalmente; en runs
RCP, la rampa nominal (año − año_inicial)·F_esc/largo con F_esc = 4.5/8.5
normalizada por 8.5. ψ(t₀) = fracción supercrítica del ensemble; estadístico de
alarma = media móvil causal de ψ sobre 5 orígenes.

## 4. Calibración (fija; orden registrado en el log)

1. Umbral supercrítico por configuración: percentil mínimo de la serie de
   multiplicadores del control (orígenes cada 5 años desde el año 60 del control,
   configuración seleccionada por CV interna, regla a nivel de refit) con FPR
   sostenida ≤ 5 %/siglo — el mismo procedimiento del protocolo ascendente,
   recalculado en CESM_0600.
2. Umbral de alarma sobre la media móvil de ψ: máximo del control + 0.01.
3. Orígenes en runs de prueba: cada 5 años desde el año 100 del run (o desde un
   tercio del largo si el run es más corto que 300 años, para los RCP);
   refit en cada origen (la variante d=1 es barata).

## 5. Verdad de terreno (fija, se computa al final)

**Onset (primaria):** referencia = media y sd de los primeros 100 años del run
(o primer tercio si <300 años); onset = primer año de una racha de 5 años
consecutivos bajo (media − 3·sd). Si no existe, el run se clasifica
**no-colapsante**. **Secundaria (solo descriptiva):** regresión de quiebre de dos
segmentos.

## 6. Hipótesis y gates (fijos)

- **T1 — Paridad:** loader vs series archivadas, tolerancia 10⁻⁵ Sv (racional de la
  Enmienda 3 del caso previo).
- **T2 — Discriminación por proximidad (primaria):** (a) cero alarmas certificadas
  en el control; (b) en cada run colapsante, el máximo pre-onset de la media móvil
  de ψ supera el umbral; (c) el máximo pre-onset medio de los runs colapsantes
  supera al máximo de todos los runs no-colapsantes. Gate aprobado si (a) y (b) se
  cumplen y (c) se cumple; parcialmente aprobado si (a) y (b) sí, (c) no.
- **T3 — Timing (secundaria, por run):** alarma certificada (2 orígenes
  consecutivos ≥ umbral; tolerancia de retorno 1) viva en el onset con lead > 0, y
  conteo de alarmas certificadas-muertas pre-onset. En runs RCP (velocidad finita)
  este es el test de cronómetro; se reporta por run sin agregación.
- **T4 — Referencia F_ovS (descriptiva):** regla causal de giro con retractación
  (idéntica a la del caso previo) sobre las series F_ovS disponibles.

## 7. Interpretación prespecificada

T2 aprobado (con o sin T3) confirma el claim de medidor de proximidad en runs
independientes con desenlaces desconocidos al congelar. T3 aprobado en los runs
RCP habilita el claim de anticipación a velocidad finita. T2 fallado se reporta
como límite del método en régimen de forzamiento fijo cerca del fold. Ningún
umbral, configuración ni regla se ajusta tras ver los datos; análisis posteriores
se rotulan exploratorios.

## 8. Registro

`run_transient_experiment.py` registra timestamps del orden (control →
calibración → runs de prueba → onset → gates), hashes SHA-256 de código y de este
protocolo, y escribe `TRANSIENT_EXPERIMENT_RESULTS.json`.
