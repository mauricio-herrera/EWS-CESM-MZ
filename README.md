# EWS_CESM_MZ_experiment — Nivel 2 (paper de seguimiento)

Early warning de colapso simulado de AMOC (CESM, hosing cuasi-estático, van Westen
et al. 2024) con el operador MZ vectorial v3. Estudio separado del confirmatorio CMIP6.

Contenido:
- `PROTOCOL_EWS_CESM_FROZEN.md` — protocolo congelado v1.0 (timestampear ANTES de analizar).
- `cesm_ingest.py` — ingestión Zenodo/GitHub -> contrato de estado v3; paridad exacta
  con las definiciones de referencia del release (AMOC 26N y FovS).
- `PARITY_EVIDENCE.json` — gate E1 satisfecho: diferencias <= 1.8e-15 Sv.
- `demo_outputs/` — tabla de contrato y perfiles para las ventanas 1-50, 1701-1750,
  1751-1800 (validación de extremo a extremo; el colapso de 1758 es visible).
- `cache_cesm/` NO se incluye en el zip (datos CC-BY descargables con el módulo).

Uso mínimo:
    from pathlib import Path
    import cesm_ingest as ci
    ci.build_contract_table(ci.PERIODS, Path("cache_cesm"), Path("outputs"))

Pasos siguientes (post-timestamp): implementar C1-C4 y la regla de alarma del
protocolo; correr sobre los 44 periodos (~460 MB de descarga total).
