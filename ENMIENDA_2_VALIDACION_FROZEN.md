# ENMIENDA 2 (CONGELADA) — Validación ciega fuera de muestra del estadístico de ensemble

**Fecha de congelamiento:** 2026-07-19.
**Estado de ceguera:** al congelar esta enmienda NO se ha descargado ni inspeccionado
ningún valor de datos de la rama descendente. Las únicas operaciones realizadas fueron
listar nombres de archivos del repositorio para confirmar formato y cobertura
(años 2201–4400, misma estructura que la rama ascendente). Esta enmienda debe recibir
timestamp público antes de ejecutar `run_validation_downbranch.py`.

## 1. Motivación y estatus del hallazgo a validar

En la corrida ascendente (colapso en el año-modelo 1758), el análisis exploratorio
posterior a las reglas prespecificadas identificó un estadístico robusto: la fracción
supercrítica del ensemble de configuraciones del operador MZ, promediada en 50 años,
que produce alarma en 1681 (lead 77 años) con cero falsas alarmas respecto del
control. Por haber sido definido tras observar los datos, este estadístico requiere
confirmación en una transición independiente con todos sus componentes congelados.
Esta enmienda congela esos componentes y define el test.

## 2. Datos de validación (fijos)

- Fuente: van Westen & Dijkstra (2023), experimento de histéresis quasi-equilibrio,
  release GitHub `RenevanWesten/GRL-AMOC-Hysteresis` (Zenodo 10.5281/zenodo.8262424,
  CC-BY 4.0).
- Segmento: **rama descendente, años-modelo 2201–4400** (44 archivos de
  `AMOC_structure`, `AMOC_section_26N`, `FOV_section_34S`). La rama ascendente
  (1–2200) es la misma simulación ya analizada y queda excluida.
- Forzamiento: F_H(t) = 0.66 − 3·10⁻⁴·(t − 2200) Sv para t ∈ (2200, 4400]
  (rampa descendente simétrica; verificar tasa contra el paper antes de correr y
  documentar cualquier discrepancia como nota, sin alterar el resto).
- Transición objetivo: la recuperación (off→on) de la AMOC en la rama descendente.

## 3. Verdad de terreno (procedimiento fijo, se computa AL FINAL)

El año de transición t\* se estima con el mismo espíritu que el tipping publicado de
la rama ascendente: regresión lineal por tramos (dos segmentos, punto de quiebre
libre) sobre el transporte AMOC 0–1000 m en 26°N de la rama descendente, eligiendo el
quiebre que minimiza el error cuadrático total. El script lo computa después de emitir
las alarmas, y el orden de ejecución queda registrado en el log.

## 4. Estadístico congelado (sin grados de libertad restantes)

1. **Contrato de estado:** idéntico al protocolo v1.0 (5 coordenadas físicas en la
   latitud de grilla más cercana a 26.0°N + 3 EOFs del perfil 500–4500 m en la grilla
   de 33 niveles), con estandarización y EOFs causales dentro de la ventana
   [2201, t₀].
2. **Ensemble:** las 24 configuraciones exactas de `ews_analysis.py`
   (TAU_SETS × OSC_GRID × ALPHA_GRID), ajustadas por `fit_operator` con la covariable
   u = F_H estandarizada causalmente.
3. **Umbral supercrítico por configuración (transferido, congelado):**
   multiplicador exacto ≥ **0.9994927450828757** (calibración de control de la rama
   ascendente; no se recalibra).
4. **ψ(t₀):** fracción de las 24 configuraciones supercríticas en el origen t₀.
   Orígenes: t₀ = 2601, 2611, …, 4391 (ventana mínima de entrenamiento: 400 años).
5. **Estadístico de alarma:** media móvil causal de ψ sobre 5 orígenes (50 años).
6. **Umbral de alarma (transferido, congelado):** **0.46** (máximo del control
   ascendente 0.45 más margen 0.01).
7. **Regla de certificación:** regla A2 sobre la serie de la media móvil
   (2 orígenes consecutivos ≥ umbral; tolerancia de retorno 1 origen).

## 5. Gates de validación (fijos)

- **V1 — Paridad del loader en la rama descendente:** |Δ| ≤ 10⁻¹⁰ Sv contra
  `Data/CESM/Ocean/AMOC_transport_depth_0-1000m.nc` del repositorio de histéresis en
  los años usados (y contra `FOV_index_section_34S.nc` si se usa F_ovS descriptivo).
- **V2 — Éxito primario:** existe una alarma certificada con año < t\* (lead > 0), y
  el número de alarmas certificadas "muertas" (que tras certificarse violan la
  tolerancia de retorno antes de t\*) es cero.
- **V3 — Descriptivos (sin gate):** lead en años; ψ y su media móvil completas;
  F_ovS de la rama descendente como serie de referencia (su regla causal de giro fue
  diseñada para descensos y no se reinterpreta aquí; solo se grafica); AR(1)/varianza
  causales como comparación informativa.

## 6. Interpretación prespecificada

- **V2 aprueba:** el estadístico de ensemble anticipa dos transiciones distintas
  (colapso forzado y recuperación) del mismo modelo con componentes idénticos y cero
  falsas alarmas; el lead de ~77 años de la rama ascendente pasa a ser un resultado
  confirmado fuera de muestra a nivel de regla (se reportan ambos leads).
- **V2 falla:** el estadístico queda como hallazgo exploratorio de la rama ascendente;
  el paper lo reporta como tal, junto con el fallo de validación, y la comparación
  principal queda entre F_ovS (29 años, 45 retractaciones) y la ausencia de señal de
  los métodos clásicos.
- En ningún caso se ajustan umbrales, configuraciones ni reglas tras ver la rama
  descendente; cualquier análisis posterior a los gates se rotula exploratorio.

## 7. Registro

El script `run_validation_downbranch.py` escribe un log con timestamps del orden de
ejecución (ingesta → paridad → ψ → alarma → verdad de terreno → veredicto) y un JSON
final `VALIDATION_DOWNBRANCH_RESULTS.json`. El hash SHA-256 de `ews_analysis.py` y de
esta enmienda al momento de la corrida se incluyen en el JSON.
