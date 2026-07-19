# PROTOCOLO CONGELADO — Early warning de colapso simulado de AMOC con operador MZ vectorial (Nivel 2)

**Versión:** 1.0 (congelada). **Fecha de congelamiento:** 2026-07-19.
**Estado:** post-diseño, pre-análisis. Ningún diagnóstico de early warning (MZ, AR(1),
varianza, reglas de alarma) ha sido computado sobre los datos CESM al momento de
congelar este protocolo. Las únicas operaciones realizadas fueron: (i) inspección de la
estructura de archivos y variables; (ii) construcción del contrato de estado; (iii)
verificación de paridad del loader contra las series procesadas archivadas por los
autores originales (evidencia en `PARITY_EVIDENCE.json`). Este protocolo debe recibir
timestamp público (commit firmado, OSF o Zenodo) antes de ejecutar cualquier análisis.

## 1. Pregunta e hipótesis

**Pregunta.** En un colapso cuasi-estático simulado de la AMOC con ground truth
conocido, ¿el multiplicador lento del operador MZ vectorial aprendido causalmente
anticipa el colapso antes y con menos falsas alarmas que los indicadores estadísticos
clásicos?

**H1 (primaria).** El lead time de la alarma basada en el multiplicador MZ es mayor o
igual que el de la alarma AR(1)/varianza, con tasa de falsas alarmas equivalente
calibrada en el segmento de control.

**H2 (secundaria, descriptiva).** El balance de Schur δ_n separa la contribución del
feedback retardado de la del amortiguamiento instantáneo durante la aproximación al
colapso; se reporta la descomposición sin gate asociado.

**H3 (comparativa, sin gate).** El lead de MZ se compara con el del indicador físico
F_ovS (mínimo por splines publicado: año-modelo 1732; mínimo crudo: 1726; tipping:
1758; lead físico de referencia: 26–32 años según convención). F_ovS incorpora
información física privilegiada; la comparación es informativa, no un gate.

## 2. Datos (fijos)

- Fuente: van Westen, Kliphuis & Dijkstra (2024), *Sci. Adv.* 10:eadk1189. Release
  `SA-AMOC-Collapse_v1.0`, Zenodo DOI 10.5281/zenodo.10461549, licencia CC-BY 4.0.
- Corrida: hosing cuasi-estático CESM, años-modelo 1–2200, rampa
  F_H(t) = 3·10⁻⁴·t Sv (0.66 Sv en el año 2200). Resolución anual.
- Insumos: `AMOC_structure` (función de corriente año×profundidad×latitud),
  `AMOC_section_26N` (transporte 0–1000 m de referencia), `FOV_section_34S` (F_ovS).
- Ground truth del tipping: año-modelo **1758** (regresión con quiebre, publicado).
- Hashes SHA-256 de cada archivo descargado se archivan junto al caché
  (`cesm_ingest.download`).

## 3. Contrato de estado (fijo)

Estado vectorial anual en la latitud de grilla más cercana a 26.0°N (25.796°N):
d = 5 coordenadas físicas — máximo de celda superior, profundidad del máximo, mínimo
de retorno profundo (banda 1500–4500 m), transporte medio superior (0–1000 m),
transporte medio profundo (1500–4500 m) — más k = 3 componentes EOF del perfil
500–4500 m interpolado a la grilla de 33 niveles. Las EOFs y toda normalización se
ajustan exclusivamente dentro de la ventana de entrenamiento causal de cada origen.
Covariable exógena: u_n = F_H(año), conocida en forma cerrada. El índice escalar de
evento es el transporte 0–1000 m en 26°N (definición de referencia del release).

## 4. Diseño causal (fijo)

- **Segmento de control (calibración de umbrales):** años 1–400 (F_H ≤ 0.12 Sv).
  Todos los umbrales de alarma, para todos los métodos, se calibran únicamente aquí.
- **Orígenes de evaluación:** cada año desde 401 hasta 1757.
- **Entrenamiento causal:** ventana expansiva [1, t₀]. El operador MZ (y el AR de
  comparación) se reajustan cada 25 años; entre reajustes, el diagnóstico anual usa el
  último operador entrenado con datos ≤ t₀. Ningún dato posterior a t₀ entra en
  representación, normalización, EOFs, hiperparámetros ni residuos.
- **Hiperparámetros:** grilla idéntica a la del protocolo confirmatorio v3
  (kernels con bloques firmados/oscilatorios incluidos), seleccionados por validación
  interna dentro de [1, t₀] con bloques temporales contiguos.

## 5. Métodos comparados (fijos)

- **C1 — AR(1) + varianza (clásico):** coeficiente lag-1 y varianza en ventana móvil
  de 100 años sobre el índice escalar detrendado (spline como en la literatura EWS);
  alarma por cruce de umbral calibrado en control.
- **C2 — F_ovS:** serie del release; alarma con la misma regla de umbral y también
  lead publicado (mínimo por splines) como referencia.
- **C3 — Multiplicador MZ vectorial (primario):** multiplicador lento exacto del
  Jacobiano del operador aprendido y su reconstrucción de segundo orden; alarma por
  cruce de umbral calibrado en control.
- **C4 — Rollout de primer paso (secundario):** probabilidad a horizonte 50 años de
  caer bajo el índice de 10 Sv; solo descriptivo en este nivel.

## 6. Regla de alarma (fija, común a todos los métodos)

Umbral: percentil del diagnóstico en el control tal que la tasa de falsas alarmas
sostenidas sea ≤ 5% por siglo de control. Alarma sostenida: diagnóstico más allá del
umbral en ≥ 3 evaluaciones anuales consecutivas. Lead time: 1758 − (primer año de
alarma sostenida tras el cual el diagnóstico no retorna bajo el umbral por más de 10
años consecutivos). Se reporta además la curva completa alarma-vs-tiempo y la curva
lead/falsas-alarmas al variar el umbral (análoga a AUC dependiente del tiempo).

## 7. Gates de decisión (fijos)

- **E1 — Paridad del loader:** |Δ| ≤ 10⁻¹⁰ Sv contra `AMOC_transport_depth_0-1000m.nc`
  y `FOV_index_section_34S.nc` del release, en todos los años usados.
  *Estado al congelar: SATISFECHO (máx 1.8·10⁻¹⁵ Sv; evidencia archivada).*
- **E2 — Control:** todos los métodos alcanzan FPR ≤ 5% por siglo en el control con
  sus umbrales propios; si algún método no lo logra, se reporta y se excluye de H1.
- **E3 — Primario:** lead(C3) ≥ lead(C1). Si falla, el resultado se publica como
  negativo del diagnóstico local en este régimen.
- **E4 — Sensibilidad:** el signo de E3 no cambia bajo: latitud 26.5°N, banda
  profunda 1500–4000 m, ventana AR de 70 años, y reajuste cada 10 años.

## 8. Límites de inferencia (declarados)

Un solo evento de colapso: todo resultado es un estudio de caso, sin estimación de
tasas ni claims de generalización entre modelos. El segmento de control es un
pseudo-control (F_H pequeño pero no nulo); como robustez opcional se puede replicar la
calibración en el dataset de histéresis (Zenodo 10.5281/zenodo.8262424). Este
protocolo pertenece a un estudio separado del confirmatorio CMIP6 v3 y no comparte
claims con él.

## 9. Enmiendas

Cualquier cambio posterior al timestamp se documenta como enmienda numerada con fecha
y justificación, sin sobrescribir esta versión.

---

## ENMIENDA 1 (v1.1) — 2026-07-19

**Contexto.** Redactada tras la corrida primaria (post-datos, se declara). La corrida
primaria reveló tres defectos de especificación de la regla de alarma y de la
evaluación de gates; ninguna enmienda modifica datos, contrato de estado, métodos,
umbral de FPR ni la definición de lead. La corrida primaria se archiva sin cambios.

**A1 — Semántica del gate E3.** Cuando ambos leads son nulos, E3 es *no evaluable*
(no "aprobado"). Si solo el lead MZ existe, E3 aprueba; si solo existe el de AR(1),
falla. La evaluación `(lead or -1) >= (lead or -1)` de la corrida primaria queda
invalidada; el `E3: true` del RESULTS_SUMMARY primario es vacuo.

**A2 — Regla de alarma a nivel de refit para diagnósticos constantes por bloques.**
El diagnóstico MZ es constante entre refits, lo que vacía el sostenimiento anual y
endurece arbitrariamente el no-retorno. La regla pasa a operar sobre la serie a nivel
de refit: sostenimiento = 2 refits consecutivos sobre umbral; no-retorno = a lo más 1
refit consecutivo bajo umbral después de la alarma. Calibración de umbral y FPR en
control con la misma regla a nivel de refit. **Nota de honestidad:** bajo esta regla el
episodio supercrítico 1676–1700 de la corrida primaria (un solo refit) sigue sin
calificar como alarma; la regla no se ajustó para certificarlo. El test decisivo es la
variante refit_every=10 del gate E4, prevista en el protocolo original.

**A3 — Regla causal de mínimo para F_ovS.** El umbral de nivel es inadecuado para un
indicador monótonamente tendencial (disparo trivial en el año 405). F_ovS pasa a la
versión causal del indicador publicado: alarma en el primer año del período de
evaluación en que (i) el mínimo corriente está por debajo del mínimo del control
(excursión sin precedentes en 400 años) y (ii) la serie se sitúa al menos una
desviación estándar del control por sobre ese mínimo corriente durante ≥ 3 años
consecutivos. Toda alarma se retracta si ocurre un nuevo mínimo posterior; se reporta
la última alarma no retractada antes del tipping (lead) y el número de alarmas
retractadas como costo de falsas alarmas del método. Referencias retrospectivas: mínimo crudo 1726, spline
1732.

**A4 — Trazabilidad.** El lazo C3 registra por refit los hiperparámetros
seleccionados (taus, componente oscilatoria, α). Se agrega un análisis de robustez del
episodio supercrítico: fracción de configuraciones (t₀ ∈ vecindad, α × taus × osc)
con multiplicador ≥ umbral.

---

## ENMIENDA 2 (v1.2) — 2026-07-19 — Validación fuera de muestra del estadístico de ensemble

**Estado.** Congelada ANTES de computar cualquier diagnóstico sobre la rama de bajada.
Sobre esos datos solo se ha ejecutado una verificación mecánica del loader (ventana
4051–4100). Requiere timestamp público antes de ejecutar.

**Origen y honestidad.** La corrida primaria y su enmienda 1 dieron E3 = no evaluable.
El análisis exploratorio posterior (documentado en los derivados de la ronda 2)
identificó como estadístico robusto la fracción supercrítica del ensemble y mostró que
la separación frente a la variabilidad interna del control exige duración, no nivel.
El estadístico resultante se congela aquí en su forma funcional y se somete a un
evento independiente. Ningún resultado exploratorio de la rama de subida se reporta
como confirmatorio.

**Estadístico congelado (forma funcional).**
ψ(t₀) = fracción de configuraciones con multiplicador lento exacto ≥ umbral de
multiplicador, sobre la grilla fija de 24 configuraciones (3 conjuntos τ × 2 opciones
oscilatorias × 4 valores α del protocolo v1.0), con operador ajustado causalmente en
[inicio de rama, t₀]. Variable de alarma: media móvil causal de ψ sobre 50 años, con
orígenes cada 10 años. Regla de alarma: la de la Enmienda 1-A2 (2 orígenes
consecutivos sobre umbral; no-retorno = a lo más 1 origen bajo umbral).

**Diseño en la rama de bajada (GRL-AMOC-Hysteresis, Zenodo 10.5281/zenodo.10034589).**
- Datos: años 2201–4400, F_H decreciente a 3·10⁻⁴ Sv/año desde 0.66 Sv. El
  entrenamiento comienza en 2201; NUNCA cruza a la rama de subida (regímenes
  dinámicos distintos y datos ya explorados).
- Ground truth: inicio de la recuperación en el año-modelo 4091 (F_H = 0.093 Sv;
  recuperación completa 4090–4170, publicado).
- Pseudo-control: orígenes 2301–2700 (estado colapsado profundo). Calibraciones:
  (i) umbral de multiplicador = percentil de la Enmienda 1-A2 sobre el diagnóstico MZ
  puntual del pseudo-control; (ii) umbral de ψ̄₅₀ = máximo del pseudo-control. Los
  NÚMEROS se recalibran aquí; las REGLAS son las congeladas.
- Evaluación: orígenes 2701–4085.
- Comparadores con reglas espejadas: AR(1)/varianza (regla anual v1.0) y F_ovS con la
  regla causal de giro de la Enmienda 1-A3 con orientación invertida (detección de
  máximo), dado que en la rama de bajada F_ovS es positivo y el giro precede a la
  recuperación.

**Gates de validación.**
- V1: paridad del loader contra las series procesadas archivadas del release de
  histéresis (tolerancia 10⁻¹⁰ Sv).
- V2: FPR de ψ̄₅₀ en pseudo-control = 0 por construcción; se reporta la FPR de los
  comparadores con sus umbrales propios.
- V3 (primario): la alarma de ψ̄₅₀ ocurre con lead > 0 respecto de 4091 y sin
  retractación hasta la transición. Éxito ⇒ el lead exploratorio de la rama de subida
  (77 años) se reporta como replicado-en-forma en un evento independiente. Fracaso ⇒
  se reporta que el estadístico no transfiere y el único indicador con anticipación
  transferible es el físico.
- V4: sensibilidad mínima: media móvil de 30 y 70 años; pseudo-control 2301–2600; el
  signo de V3 no debe cambiar.

**Limitación declarada.** Mismo modelo climático y misma familia de forzamiento que la
rama de subida; la independencia es del EVENTO y del mecanismo de transición
(recuperación mediada por hielo marino, ~6× más rápida), no del modelo.
