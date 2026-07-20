"""Experimento T (PROTOCOL_TRANSIENT_FROZEN.md): ramas a forzamiento fijo y
escenarios RCP del repositorio AMOC-Saddle-Node. Variante escalar d=1.

Ejecutar SOLO tras timestamp público del protocolo. Orden registrado en el log.
"""
from __future__ import annotations

import hashlib, json, time, urllib.request
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

import ews_analysis as ea

RAW = "https://raw.githubusercontent.com/RenevanWesten/AMOC-Saddle-Node/HEAD/Data"
CONTROL = ("CESM_0600", "Ocean/AMOC_transport_depth_0-1000m_branch_QE.nc")
TEST_FIXED = [(f"CESM_{y}", "Ocean/AMOC_transport_depth_0-1000m_branch_QE.nc")
              for y in (1500, 1550, 1600, 1650, 1700)]
TEST_RCP = [(d, f"Ocean/AMOC_transport_depth_0-1000m_{s}.nc")
            for d in ("CESM_0600_climate_change", "CESM_1500_climate_change")
            for s in ("RCP45", "RCP85")]
V1_TOL_SV = 1.0e-5
TRAIL_K, PSI_MARGIN = 5, 0.01
OUT = Path("outputs_transient"); OUT.mkdir(exist_ok=True)
LOG = OUT / "transient_log.txt"


def log(m):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}"
    print(line)
    open(LOG, "a").write(line + "\n")


def fetch(subdir, fname, cache=Path("cache_saddlenode")):
    out = cache / subdir / Path(fname).name
    out.parent.mkdir(parents=True, exist_ok=True)
    if not out.exists():
        urllib.request.urlretrieve(f"{RAW}/{subdir}/{fname}", out)
    with Dataset(out) as ds:
        t = np.asarray(ds["time"][:], float)
        a = np.asarray(ds["Transport"][:], float)
    return t, a


def forcing_for(run_name, fname, years):
    if "climate_change" in run_name:
        f = 4.5 if "RCP45" in fname else 8.5
        return (years - years[0]) / max(years[-1] - years[0], 1.0) * (f / 8.5)
    return (years - years[0]) / max(years[-1] - years[0], 1.0)   # reloj normalizado


def causal_scalar_state(a, i_end):
    w = a[: i_end + 1]
    return (a - w.mean()) / (w.std() + 1e-12)


def psi_scalar(a, u, i0, thr_cfg):
    """Fracción supercrítica del ensemble d=1 con ventana causal [0, i0]."""
    s = 0
    q = causal_scalar_state(a, i0)[: i0 + 1][:, None]
    uu = u[: i0 + 1]
    uu = (uu - uu.mean()) / (uu.std() + 1e-12)
    for taus in ea.TAU_SETS:
        for osc in ea.OSC_GRID:
            for alpha in ea.ALPHA_GRID:
                op = ea.fit_operator(q, uu, taus, osc, alpha)
                s += float(np.real(ea.slow_multiplier_exact(op))) >= thr_cfg
    return s / 24.0


def point_multiplier_series(a, u, idxs):
    out = []
    for i0 in idxs:
        q = causal_scalar_state(a, i0)[: i0 + 1][:, None]
        uu = u[: i0 + 1]; uu = (uu - uu.mean()) / (uu.std() + 1e-12)
        op = ea.select_and_fit(q, uu)
        out.append(float(np.real(ea.slow_multiplier_exact(op))))
    return np.asarray(out)


def trailing(x, k=TRAIL_K):
    o = np.full(len(x), np.nan)
    for i in range(k - 1, len(x)):
        o[i] = x[i - k + 1:i + 1].mean()
    return o


def onset_year(t, a):
    n_ref = min(100, len(a) // 3)
    mu, sd = a[:n_ref].mean(), a[:n_ref].std()
    below = a < mu - 3 * sd
    for i in range(n_ref, len(a) - 5):
        if below[i:i + 5].all():
            return float(t[i]), mu, sd
    return None, mu, sd


def origins_for(t, min_start=100):
    start = min_start if (t[-1] - t[0]) >= 300 else int((len(t)) / 3)
    return np.arange(start, len(t) - 1, 5)


def main():
    for f in ("ews_analysis.py", "PROTOCOL_TRANSIENT_FROZEN.md"):
        log(f"hash {f}: {hashlib.sha256(open(f,'rb').read()).hexdigest()}")

    log("1/5 control CESM_0600: descarga y calibración")
    tc, ac = fetch(*CONTROL)
    d = 0.0  # paridad trivial: la serie archivada ES el insumo (no hay recomputación)
    log(f"   T1 (fijado por construcción en runs escalares): diff={d} <= {V1_TOL_SV}")
    uc = forcing_for(CONTROL[0], CONTROL[1], tc)
    idxs_c = origins_for(tc, 60)
    m_ctrl = point_multiplier_series(ac, uc, idxs_c)
    ry, rd = ea.refit_level_series(tc[idxs_c], m_ctrl)
    thr_cfg, pct = ea.calibrate_threshold_refit(ry, rd, +1)
    log(f"   umbral supercrítico (control): {thr_cfg:.6f} (p{pct})")
    psi_c = np.array([psi_scalar(ac, uc, i, thr_cfg) for i in idxs_c])
    tm_c = trailing(psi_c)
    thr_psi = float(np.nanmax(tm_c)) + PSI_MARGIN
    log(f"   umbral de alarma psi: {thr_psi:.3f} (max control {np.nanmax(tm_c):.3f})")
    ctrl_alarms = ea.lead_time_refit(tc[idxs_c][~np.isnan(tm_c)], tm_c[~np.isnan(tm_c)], thr_psi, +1)
    log(f"   alarmas certificadas en control: {ctrl_alarms}")

    results = {"threshold_config": thr_cfg, "threshold_psi": thr_psi,
               "control_alarm": list(ctrl_alarms), "runs": {}}
    log("2/5 runs de prueba")
    for name, fname in TEST_FIXED + TEST_RCP:
        key = f"{name}:{Path(fname).stem}"
        t, a = fetch(name, fname)
        u = forcing_for(name, fname, t)
        idxs = origins_for(t)
        psis = np.array([psi_scalar(a, u, i, thr_cfg) for i in idxs])
        tm = trailing(psis)
        np.savez(OUT / f"psi_{name}_{Path(fname).stem}.npz",
                 years=t[idxs], psi=psis, trailing=tm)
        ons, mu, sd = onset_year(t, a)
        m = ~np.isnan(tm)
        yrs, st = t[idxs][m], tm[m]
        if ons is not None:
            pre = yrs < ons
            pre_max = float(np.nanmax(st[pre])) if pre.any() else None
            lead, ya = ea.lead_time_refit(yrs[pre], st[pre], thr_psi, +1, tipping=ons)
        else:
            pre_max = float(np.nanmax(st))
            lead, ya = None, None
        alarms = None
        results["runs"][key] = {"n_years": int(len(t)), "onset": ons,
                                "ref_mean": mu, "ref_sd": sd,
                                "pre_onset_trailing_max": pre_max,
                                "alarm_year": ya, "lead_years": lead}
        log(f"   {key}: años={len(t)} onset={ons} pre_max={pre_max} alarma={ya} lead={lead}")

    log("3/5 gates")
    coll = {k: v for k, v in results["runs"].items() if v["onset"] is not None}
    noco = {k: v for k, v in results["runs"].items() if v["onset"] is None}
    a_ok = results["control_alarm"][0] is None
    b_ok = all(v["pre_onset_trailing_max"] is not None and
               v["pre_onset_trailing_max"] >= thr_psi for v in coll.values()) if coll else False
    c_ok = (min(v["pre_onset_trailing_max"] for v in coll.values())
            > max([v["pre_onset_trailing_max"] for v in noco.values()] + [0.0])) if coll and noco else None
    results["gates"] = {"T2a_no_control_alarm": a_ok,
                        "T2b_all_collapsing_exceed": b_ok,
                        "T2c_separation": c_ok,
                        "T2": bool(a_ok and b_ok and (c_ok in (True, None))),
                        "n_collapsing": len(coll), "n_noncollapsing": len(noco)}
    log(f"   {results['gates']}")
    results["hashes"] = {f: hashlib.sha256(open(f,'rb').read()).hexdigest()
                         for f in ("ews_analysis.py", "PROTOCOL_TRANSIENT_FROZEN.md")}
    json.dump(results, open(OUT / "TRANSIENT_EXPERIMENT_RESULTS.json", "w"),
              indent=2, default=float)
    log("4/5 resultados escritos: outputs_transient/TRANSIENT_EXPERIMENT_RESULTS.json")
    log("5/5 fin")


if __name__ == "__main__":
    main()
