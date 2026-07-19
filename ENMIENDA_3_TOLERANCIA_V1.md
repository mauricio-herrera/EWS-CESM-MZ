# ENMIENDA 3 — Tolerancia numérica del gate V1 (rama descendente)

**Fecha:** 2026-07-19. **Alcance:** exclusivamente el umbral numérico del gate de
paridad V1. No modifica el estadístico congelado, sus umbrales (0.9994927…, 0.46), la
regla de certificación, los orígenes ni el procedimiento de verdad de terreno de la
Enmienda 2, que permanecen byte a byte como en el tag `amendment-2-frozen`.

## Hecho que motiva la enmienda

La corrida de validación se detuvo en V1 con discrepancia máxima 2.05·10⁻⁷ Sv frente
a la serie procesada `AMOC_transport_depth_0-1000m.nc` del repositorio
`GRL-AMOC-Hysteresis`, superando la tolerancia original de 10⁻¹⁰ Sv (que provenía de
la paridad a precisión de máquina lograda contra el repositorio `SA-AMOC-Collapse`).

## Diagnóstico (sin inspección de valores de la rama descendente)

Se compararon entre sí las dos series procesadas oficiales en los años comunes
1–2200 de la rama ascendente (la misma simulación, archivada en ambos repositorios):

- |referencia_histéresis − referencia_SA|: máx = 5.51·10⁻⁸ Sv,
  media = 1.30·10⁻⁸ Sv (escala relativa ~5·10⁻⁹).

Es decir, los propios archivos oficiales difieren entre sí al nivel de 10⁻⁸–10⁻⁷ Sv
por procedencia numérica (orden de operaciones, precisión intermedia o versión de
librerías de la cadena de procesamiento original). El loader del paquete reproduce el
archivo SA a 1.8·10⁻¹⁵ Sv, por lo que queda exonerado: la discrepancia observada en la
rama descendente es del mismo carácter y orden que la discrepancia entre archivos.

## Cambio

- Tolerancia del gate V1: **10⁻¹⁰ Sv → 10⁻⁵ Sv** (unas 50 veces por encima de la
  discrepancia observada y ~6 órdenes de magnitud por debajo de la señal física de
  10–20 Sv; sigue detectando cualquier error real de niveles, máscaras o alineación,
  cuyos efectos son ≥ 10⁻¹ Sv).
- El script registra en el log la discrepancia observada y esta enmienda.

## Orden de registro

Esta enmienda se commitea y publica ANTES de re-ejecutar la validación; la
re-ejecución y sus resultados van en commits posteriores.
