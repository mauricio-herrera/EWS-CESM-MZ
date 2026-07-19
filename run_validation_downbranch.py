"""Validación ciega en la rama descendente de la histéresis (ENMIENDA 2, congelada).

Ejecutar SOLO después de dar timestamp público a ENMIENDA_2_VALIDACION_FROZEN.md.
Orden de ejecución registrado en el log: ingesta -> paridad (V1) -> psi -> alarma ->
verdad de terreno -> veredicto (V2). Ningún parámetro es ajustable por línea de
comandos: todos los valores están congelados en la enmienda.
"""
from __future__ import annotations

import hashlib
import json
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from netCDF4 import Dataset

import cesm_ingest as ci
import ews_analysis as ea

# ----------------------------- valores congelados ---------------------------
REPO_BASE = ("https://raw.githubusercontent.com/RenevanWesten/GRL-AMOC-Hysteresis/"
             "HEAD/Data/CESM/Data")
REF_BASE = ("https://raw.githubusercontent.com/RenevanWesten/GRL-AMOC-Hysteresis/"
            "HEAD/Data/CESM/Ocean")
PERIODS_DOWN = [(a, a + 49) for a in range(2201, 4400, 50)]
BRANCH_START = 2200.0
F_H_TOP = 0.66
F_RATE = 3.0e-4
THR_CONFIG = 0.9994927450828757      # umbral supercrítico por configuración (ascendente)
V1_TOL_SV = 1.0e-5                   # ENMIENDA 3: tolerancia entre archivos (ver ENMIENDA_3_TOLERANCIA_V1.md)
THR_PSI = 0.46                       # umbral de alarma sobre media móvil de psi
TRAIL_K = 5                          # 5 orígenes = 50 años
ORIGIN_START, ORIGIN_STEP = 2601.0, 10.0
MIN_ORIGIN_END = 4391.0

OUT = Path("outputs_validation"); OUT.mkdir(exist_ok=True)
LOG = OUT / "validation_log.txt"


def log(msg: str):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line)
    with open(LOG, "a") as f:
        f.write(line + "\n")


def forcing_down(years: np.ndarray) -> np.ndarray:
    return np.clip(F_H_TOP - F_RATE * (np.asarray(years, float) - BRANCH_START), 0.0, F_H_TOP)


def ingest_downbranch(cache: Path) -> tuple[pd.DataFrame, np.ndarray]:
    ci.RAW_BASE = REPO_BASE                       # redirige el loader al repo de histéresis
    st = ci.vector_state_26n(PERIODS_DOWN, cache)
    yrs, tr = ci.amoc_transport_26n(PERIODS_DOWN, cache)
    assert np.allclose(st["year"], yrs)
    df = pd.DataFrame({
        "year": st["year"], "F_H_Sv": forcing_down(st["year"]),
        "amoc_transport_0_1000m_26N_Sv": tr,
        "amoc_max_Sv": st["amoc_max"], "depth_of_max_m": st["depth_of_max_m"],
        "deep_return_min_Sv": st["deep_return_min"],
        "upper_mean_Sv": st["upper_mean"], "deep_mean_Sv": st["deep_mean"],
    })
    return df, st["profile_500_4500"]


def parity_v1(df: pd.DataFrame, cache: Path) -> dict:
    ref = cache / "ref_hyst_AMOC_transport.nc"
    if not ref.exists():
        urllib.request.urlretrieve(f"{REF_BASE}/AMOC_transport_depth_0-1000m.nc", ref)
    with Dataset(ref) as r:
        t = np.asarray(r["time"][:]); a = np.asarray(r["Transport"][:])
    sel = (t >= df.year.min()) & (t <= df.year.max())
    d = float(np.abs(df.amoc_transport_0_1000m_26N_Sv.values - a[sel]).max())
    return {"max_abs_diff_AMOC_Sv": d, "tolerance_Sv": V1_TOL_SV,
            "amendment": "ENMIENDA_3 (cross-archive numerical provenance)",
            "passed": d <= V1_TOL_SV}


def psi_series(df: pd.DataFrame, profiles: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    origins = np.arange(ORIGIN_START, MIN_ORIGIN_END + 1, ORIGIN_STEP)
    psis = []
    for t0 in origins:
        Q, _ = ea.causal_state(df, profiles, t0)
        trn = df.year.values <= t0
        u = df.F_H_Sv.values
        u = (u - u[trn].mean()) / (u[trn].std() + 1e-12)
        s = 0
        for taus in ea.TAU_SETS:
            for osc in ea.OSC_GRID:
                for alpha in ea.ALPHA_GRID:
                    op = ea.fit_operator(Q[trn], u[trn], taus, osc, alpha)
                    s += float(np.real(ea.slow_multiplier_exact(op))) >= THR_CONFIG
        psis.append(s / 24.0)
    return origins, np.asarray(psis)


def trailing_mean(x: np.ndarray, k: int = TRAIL_K) -> np.ndarray:
    out = np.full(len(x), np.nan)
    for i in range(k - 1, len(x)):
        out[i] = x[i - k + 1:i + 1].mean()
    return out


def certified_alarms(origins, stat, thr):
    """Regla A2 con contabilidad de alarmas muertas: lista de (año, murió_en|None)."""
    over = stat >= thr
    alarms, i, n = [], 0, len(over)
    while i <= n - 2:
        if over[i] and over[i + 1]:
            start = origins[i]
            run_below, died = 0, None
            j = i + 2
            while j < n:
                run_below = run_below + 1 if not over[j] else 0
                if run_below > 1:
                    died = origins[j]
                    break
                j += 1
            alarms.append((float(start), None if died is None else float(died)))
            i = j if died is not None else n
        else:
            i += 1
    return alarms


def transition_year(df: pd.DataFrame) -> float:
    """Verdad de terreno t*: quiebre de regresión lineal por tramos (2 segmentos)."""
    y = df.amoc_transport_0_1000m_26N_Sv.values
    t = df.year.values
    best, bk = np.inf, None
    for k in range(50, len(t) - 50):
        e = 0.0
        for a, b in [(0, k), (k, len(t))]:
            tt, yy = t[a:b], y[a:b]
            c = np.polyfit(tt, yy, 1)
            e += float(np.sum((yy - np.polyval(c, tt)) ** 2))
        if e < best:
            best, bk = e, t[k]
    return float(bk)


def main():
    cache = Path("cache_hysteresis")
    for name in ["ews_analysis.py", "ENMIENDA_2_VALIDACION_FROZEN.md"]:
        h = hashlib.sha256(open(name, "rb").read()).hexdigest()
        log(f"hash {name}: {h}")
    log("1/5 ingesta de la rama descendente (~460 MB la primera vez)")
    df, profiles = ingest_downbranch(cache)
    df.to_csv(OUT / "downbranch_state_contract.csv", index=False)
    np.savez_compressed(OUT / "downbranch_profiles.npz", year=df.year.values, profile=profiles)
    log(f"   {len(df)} años: {df.year.min():.0f}-{df.year.max():.0f}")

    log("2/5 gate V1: paridad del loader")
    v1 = parity_v1(df, cache)
    log(f"   V1 = {v1}")
    if not v1["passed"]:
        raise SystemExit("V1 FALLA: detener y revisar loader (no continuar la validación).")

    log("3/5 serie psi del ensemble (24 configs por origen; puede tardar ~1-2 min)")
    origins, psis = psi_series(df, profiles)
    stat = trailing_mean(psis)
    np.savez(OUT / "downbranch_psi.npz", origins=origins, psi=psis, trailing=stat)

    log("4/5 alarmas certificadas (regla A2, umbral congelado 0.46)")
    m = ~np.isnan(stat)
    alarms = certified_alarms(origins[m], stat[m], THR_PSI)
    log(f"   alarmas: {alarms}")

    log("5/5 verdad de terreno y veredicto V2")
    tstar = transition_year(df)
    live = [a for a, died in alarms if died is None]
    dead = [a for a, died in alarms if died is not None and a < tstar]
    lead = (tstar - live[0]) if live and live[0] < tstar else None
    v2 = bool(lead is not None and lead > 0 and len(dead) == 0)
    res = {"amendment": "ENMIENDA_2 v-congelada", "V1": v1,
           "transition_year_tstar": tstar,
           "certified_alarms": alarms, "dead_alarms_pre_tstar": dead,
           "alarm_year": (live[0] if live else None), "lead_years": lead,
           "V2_passed": v2,
           "frozen": {"thr_config": THR_CONFIG, "thr_psi": THR_PSI,
                      "trail_k": TRAIL_K, "origins": [ORIGIN_START, ORIGIN_STEP]}}
    json.dump(res, open(OUT / "VALIDATION_DOWNBRANCH_RESULTS.json", "w"), indent=2)
    log(f"   t*={tstar:.0f} | alarma={res['alarm_year']} | lead={lead} | V2={'APRUEBA' if v2 else 'FALLA'}")
    log("listo: outputs_validation/VALIDATION_DOWNBRANCH_RESULTS.json")


if __name__ == "__main__":
    main()
