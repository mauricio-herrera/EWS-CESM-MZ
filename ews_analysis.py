"""Análisis EWS-CESM según PROTOCOL_EWS_CESM_FROZEN.md (v1.0).

Implementa los cuatro métodos comparados y la regla de alarma común:
  C1  AR(1) + varianza en ventana móvil causal de 100 años.
  C2  F_ovS (serie del release, regla de alarma común; orientación decreciente).
  C3  Operador MZ vectorial lineal con kernels firmados/oscilatorios, entrenado
      causalmente con ventana expansiva; diagnóstico = multiplicador lento exacto
      del Jacobiano por bloques, con reconstrucción de segundo orden
      (Lyapunov–Schmidt) como verificación interna.
  C4  Primer paso por rollout estocástico (secundario, descriptivo).

Guardas de causalidad: para cada origen t0, estandarización, EOFs, selección de
hiperparámetros, residuos y ajuste usan exclusivamente años <= t0.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

TIPPING_YEAR = 1758.0
CONTROL_YEARS = (1.0, 400.0)
EVAL_YEARS = (401.0, 1757.0)
REFIT_EVERY = 25
AR_WINDOW = 100
N_EOF = 3
FPR_TARGET_PER_CENTURY = 0.05
SUSTAIN_YEARS = 3
NO_RETURN_YEARS = 10

TAU_SETS = ((1.0, 5.0, 20.0), (1.0, 6.0, 28.0), (2.0, 8.0, 30.0))
ALPHA_GRID = (1e-3, 1e-2, 1e-1, 1.0)
OSC_GRID = (None, (40.0, 0.90))          # (periodo [años], amortiguamiento rho)


# ---------------------------------------------------------------------------
# Estado causal: estandarización + EOFs dentro de la ventana de entrenamiento
# ---------------------------------------------------------------------------

PHYS_COLS = ["amoc_max_Sv", "depth_of_max_m", "deep_return_min_Sv",
             "upper_mean_Sv", "deep_mean_Sv"]


def causal_state(df: pd.DataFrame, profiles: np.ndarray, t0: float) -> tuple[np.ndarray, dict]:
    """Q (n_años x d) estandarizado y con EOFs ajustadas solo con años <= t0."""
    tr = df.year.values <= t0
    X = df[PHYS_COLS].values.astype(float)
    mu, sd = X[tr].mean(0), X[tr].std(0)
    sd[sd == 0] = 1.0
    Xs = (X - mu) / sd
    P = profiles.astype(float)
    pmu = P[tr].mean(0)
    Pc = P - pmu
    U, S, Vt = np.linalg.svd(Pc[tr], full_matrices=False)
    V = Vt[:N_EOF].T
    scale = Pc[tr] @ V
    ssd = scale.std(0)
    ssd[ssd == 0] = 1.0
    E = (Pc @ V) / ssd
    Q = np.concatenate([Xs, E], axis=1)
    meta = {"mu": mu, "sd": sd, "pmu": pmu, "V": V, "ssd": ssd}
    return Q, meta


# ---------------------------------------------------------------------------
# Operador MZ vectorial lineal (kernels firmados + oscilatorios)
# ---------------------------------------------------------------------------

def _memory_blocks(d: int, taus: Sequence[float], osc: tuple | None):
    """Bloques (D, C_drive) del lift; retorna también dimensión oculta."""
    blocks_D, blocks_C = [], []
    for tau in taus:
        rho = float(np.exp(-1.0 / tau))
        blocks_D.append(rho * np.eye(d))
        blocks_C.append((1.0 - rho) * np.eye(d))
    if osc is not None:
        period, rho = osc
        w = 2.0 * np.pi / period
        R = rho * np.array([[np.cos(w), np.sin(w)], [-np.sin(w), np.cos(w)]])
        blocks_D.append(np.kron(R, np.eye(d)))
        C = np.zeros((2 * d, d)); C[:d] = (1.0 - rho) * np.eye(d)
        blocks_C.append(C)
    from scipy.linalg import block_diag
    D = block_diag(*blocks_D)
    C = np.vstack(blocks_C)
    return D, C


def _run_memory(Q: np.ndarray, D: np.ndarray, C: np.ndarray) -> np.ndarray:
    """z_{n+1} = D z_n + C q_n, z_0 = 0. Devuelve Z alineado con Q (z_n)."""
    n, _ = Q.shape
    h = D.shape[0]
    Z = np.zeros((n, h))
    for t in range(n - 1):
        Z[t + 1] = D @ Z[t] + C @ Q[t]
    return Z


@dataclass
class MZOperator:
    A: np.ndarray          # d x d  (I + A_loc)
    B: np.ndarray          # d x h  (lectura de memoria, con signo libre)
    C: np.ndarray          # h x d
    D: np.ndarray          # h x h
    Eu: np.ndarray         # d      (acople al forzamiento)
    b0: np.ndarray         # d
    resid: np.ndarray      # residuos de entrenamiento (n_tr-1 x d)
    taus: tuple
    osc: tuple | None
    alpha: float

    def moments(self, k_max: int = 2) -> list[np.ndarray]:
        ImD_inv = np.linalg.inv(np.eye(self.D.shape[0]) - self.D)
        out, M = [], np.eye(self.D.shape[0])
        for _ in range(k_max + 1):
            M = M @ ImD_inv
            out.append(self.B @ M @ self.C)
        return out  # [N0, N1, N2]

    def full_map(self) -> np.ndarray:
        top = np.hstack([self.A, self.B])
        bot = np.hstack([self.C, self.D])
        return np.vstack([top, bot])


def fit_operator(Q: np.ndarray, u: np.ndarray, taus, osc, alpha: float) -> MZOperator:
    d = Q.shape[1]
    D, C = _memory_blocks(d, taus, osc)
    Z = _run_memory(Q, D, C)
    X = np.hstack([Q[:-1], Z[:-1], u[:-1, None], np.ones((len(Q) - 1, 1))])
    Y = Q[1:] - Q[:-1]
    G = X.T @ X + alpha * np.eye(X.shape[1])
    W = np.linalg.solve(G, X.T @ Y)              # (d+h+2) x d
    A_loc = W[:d].T
    B = W[d:d + Z.shape[1]].T
    Eu = W[-2]
    b0 = W[-1]
    resid = Y - X @ W
    return MZOperator(np.eye(d) + A_loc, B, C, D, Eu, b0, resid, tuple(taus), osc, alpha)


def select_and_fit(Q: np.ndarray, u: np.ndarray) -> MZOperator:
    """Selección causal: último 20% de [0,t0] como bloque de validación."""
    n = len(Q)
    n_val = max(20, int(0.2 * n))
    best, best_cfg = np.inf, None
    for taus in TAU_SETS:
        for osc in OSC_GRID:
            for alpha in ALPHA_GRID:
                op = fit_operator(Q[: n - n_val], u[: n - n_val], taus, osc, alpha)
                Z = _run_memory(Q, op.D, op.C)
                X = np.hstack([Q[:-1], Z[:-1], u[:-1, None], np.ones((n - 1, 1))])
                W = np.vstack([(op.A - np.eye(Q.shape[1])).T, op.B.T,
                               op.Eu[None, :], op.b0[None, :]])
                pred = X @ W
                err = np.mean((Q[1:] - Q[:-1] - pred)[-n_val:] ** 2)
                if err < best:
                    best, best_cfg = err, (taus, osc, alpha)
    taus, osc, alpha = best_cfg
    return fit_operator(Q, u, taus, osc, alpha)


# ---------------------------------------------------------------------------
# Diagnósticos de Schur: multiplicador exacto y expansión de segundo orden
# ---------------------------------------------------------------------------

def slow_multiplier_exact(op: MZOperator) -> complex:
    ev = np.linalg.eigvals(op.full_map())
    return ev[np.argmin(np.abs(ev - 1.0))]


def slow_multiplier_expansion(op: MZOperator) -> tuple[float, float]:
    """(rho_2do_orden, delta) vía Lyapunov–Schmidt sobre S = I - A - N0."""
    N0, N1, N2 = op.moments(2)
    d = op.A.shape[0]
    S = np.eye(d) - op.A - N0
    H = np.eye(d) + N1
    ev, R = np.linalg.eig(S)
    k = int(np.argmin(np.abs(ev)))
    delta = -np.real(ev[k])
    r = np.real(R[:, k]); r /= np.linalg.norm(r)
    evl, L = np.linalg.eig(S.T)
    kl = int(np.argmin(np.abs(evl - ev[k])))
    l = np.real(L[:, kl]); l /= (l @ r)
    gamma = float(l @ H @ r)
    Qp = np.eye(d) - np.outer(r, l)
    Pm = np.outer(r, l)
    G = np.linalg.inv(Qp @ S @ Qp + Pm) - Pm
    corr = float(l @ (N2 + H @ G @ Qp @ H) @ r)
    rho2 = 1.0 + delta / gamma + corr * delta**2 / gamma**3
    return rho2, delta


# ---------------------------------------------------------------------------
# C3: lazo causal de diagnóstico
# ---------------------------------------------------------------------------

def mz_causal_diagnostics(df: pd.DataFrame, profiles: np.ndarray,
                          eval_years=EVAL_YEARS, refit_every=REFIT_EVERY,
                          progress=True) -> pd.DataFrame:
    years = df.year.values
    uraw = df.F_H_Sv.values
    rows = []
    op = None
    refit_marks = np.arange(eval_years[0], eval_years[1] + 1, refit_every)
    for t0 in np.arange(eval_years[0], eval_years[1] + 1):
        if op is None or t0 in refit_marks:
            Q, _ = causal_state(df, profiles, t0)
            tr = years <= t0
            u = (uraw - uraw[tr].mean()) / (uraw[tr].std() + 1e-12)
            op = select_and_fit(Q[tr], u[tr])
            Qf, _ = causal_state(df, profiles, t0)   # estado congelado hasta refit
            uf = u
            if progress:
                print(f"refit t0={int(t0)}  taus={op.taus} osc={op.osc} alpha={op.alpha}")
        zex = slow_multiplier_exact(op)
        rho2, delta = slow_multiplier_expansion(op)
        rows.append({"year": t0, "mz_exact_re": float(np.real(zex)),
                     "mz_exact_mod": float(np.abs(zex)),
                     "mz_rho2": rho2, "mz_delta": delta,
                     "spectral_radius": float(np.max(np.abs(np.linalg.eigvals(op.full_map())))),
                     "taus": str(op.taus), "osc": str(op.osc), "alpha": op.alpha})
    return pd.DataFrame(rows)


def mz_control_diagnostics(df, profiles, control=CONTROL_YEARS, step=5):
    """Diagnóstico MZ dentro del control (para calibrar umbral), causal."""
    rows = []
    years = df.year.values
    uraw = df.F_H_Sv.values
    for t0 in np.arange(control[0] + 60, control[1] + 1, step):
        Q, _ = causal_state(df, profiles, t0)
        tr = years <= t0
        u = (uraw - uraw[tr].mean()) / (uraw[tr].std() + 1e-12)
        op = select_and_fit(Q[tr], u[tr])
        zex = slow_multiplier_exact(op)
        rho2, delta = slow_multiplier_expansion(op)
        rows.append({"year": t0, "mz_exact_re": float(np.real(zex)),
                     "mz_exact_mod": float(np.abs(zex)), "mz_rho2": rho2, "mz_delta": delta})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# C1: AR(1) + varianza, ventana causal con detrend lineal intra-ventana
# ---------------------------------------------------------------------------

def ar1_variance(index: np.ndarray, years: np.ndarray, window=AR_WINDOW) -> pd.DataFrame:
    rows = []
    for i in range(window, len(index)):
        w = index[i - window:i].astype(float)
        t = np.arange(window, dtype=float)
        b = np.polyfit(t, w, 1)
        r = w - np.polyval(b, t)
        ar1 = float(np.corrcoef(r[:-1], r[1:])[0, 1])
        rows.append({"year": years[i], "ar1": ar1, "variance": float(r.var())})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Regla de alarma común (protocolo §6)
# ---------------------------------------------------------------------------

def sustained_alarm_years(years, diag, thr, orient=+1, sustain=SUSTAIN_YEARS):
    over = orient * np.asarray(diag) >= orient * thr
    runs, start = [], None
    for i, o in enumerate(over):
        if o and start is None:
            start = i
        if (not o or i == len(over) - 1) and start is not None:
            end = i if not o else i + 1
            if end - start >= sustain:
                runs.append((start, end))
            start = None
    return runs, over


def fpr_per_century(years, diag, thr, orient=+1):
    runs, _ = sustained_alarm_years(years, diag, thr, orient)
    span = (years[-1] - years[0]) / 100.0
    return len(runs) / max(span, 1e-9)


def calibrate_threshold(control_diag, orient=+1,
                        target=FPR_TARGET_PER_CENTURY, years=None):
    """Umbral = percentil mínimo (más sensible) del control con FPR sostenida <= target."""
    if years is None:
        years = np.arange(len(control_diag), dtype=float)
    pcts = np.arange(80.0, 99.95, 0.5) if orient > 0 else np.arange(20.0, 0.05, -0.5)
    for p in pcts:
        thr = float(np.percentile(control_diag, p))
        if fpr_per_century(years, control_diag, thr, orient) <= target:
            return thr, p
    p = 99.9 if orient > 0 else 0.1
    return float(np.percentile(control_diag, p)), p


def lead_time(years, diag, thr, orient=+1, tipping=TIPPING_YEAR,
              sustain=SUSTAIN_YEARS, no_return=NO_RETURN_YEARS):
    runs, over = sustained_alarm_years(years, diag, thr, orient, sustain)
    for (s, e) in runs:
        below = ~over[s:]
        run_len, ok = 0, True
        for b in below:
            run_len = run_len + 1 if b else 0
            if run_len > no_return:
                ok = False
                break
        if ok:
            return float(tipping - years[s]), float(years[s])
    return None, None


# ---------------------------------------------------------------------------
# C4: primer paso por rollout (descriptivo)
# ---------------------------------------------------------------------------

def first_passage_probability(op: MZOperator, Q, u, idx_amoc_max=0,
                              thr_std=None, horizon=50, n_traj=200, seed=0):
    rng = np.random.default_rng(seed)
    d = Q.shape[1]
    Z = _run_memory(Q, op.D, op.C)
    q0, z0 = Q[-1], Z[-1]
    u_future = np.full(horizon, u[-1])          # rampa congelada en el origen (conservador)
    hits = 0
    for _ in range(n_traj):
        q, z = q0.copy(), z0.copy()
        hit = False
        for h in range(horizon):
            eps = op.resid[rng.integers(0, len(op.resid))]
            dq = (op.A - np.eye(d)) @ q + op.B @ z + op.Eu * u_future[h] + op.b0 + eps
            z = op.D @ z + op.C @ q
            q = q + dq
            if thr_std is not None and q[idx_amoc_max] <= thr_std:
                hit = True
                break
        hits += hit
    return (hits + 0.5) / (n_traj + 1)


# ---------------------------------------------------------------------------
# ENMIENDA 1 (v1.1) — reglas de alarma corregidas y evaluación de E3
# ---------------------------------------------------------------------------

SUSTAIN_REFITS = 2          # A2: sostenimiento = 2 refits consecutivos sobre umbral
NO_RETURN_REFITS = 1        # A2: tolerancia de retorno = a lo más 1 refit bajo umbral
FOV_RISE_YEARS = 10         # A3: años consecutivos sobre el mínimo corriente


def refit_level_series(years: np.ndarray, diag: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Colapsa una serie constante por bloques a un valor por bloque de refit."""
    ry, rd = [years[0]], [diag[0]]
    for y, v in zip(years[1:], diag[1:]):
        if v != rd[-1]:
            ry.append(y); rd.append(v)
    return np.asarray(ry), np.asarray(rd)


def lead_time_refit(refit_years, refit_diag, thr, orient=+1, tipping=TIPPING_YEAR,
                    sustain=SUSTAIN_REFITS, no_return=NO_RETURN_REFITS):
    """Regla A2 sobre la serie a nivel de refit: alarma = primer refit de una racha de
    >= `sustain` refits consecutivos sobre umbral, tras la cual nunca hay más de
    `no_return` refits consecutivos bajo umbral antes del tipping."""
    over = orient * np.asarray(refit_diag) >= orient * thr
    n = len(over)
    for i in range(n - sustain + 1):
        if all(over[i:i + sustain]):
            run_below, ok = 0, True
            for b in ~over[i + sustain:]:
                run_below = run_below + 1 if b else 0
                if run_below > no_return:
                    ok = False; break
            if ok:
                return float(tipping - refit_years[i]), float(refit_years[i])
    return None, None


def fpr_per_century_refit(refit_years, refit_diag, thr, orient=+1, sustain=SUSTAIN_REFITS):
    over = orient * np.asarray(refit_diag) >= orient * thr
    n, alarms, i = len(over), 0, 0
    while i <= n - sustain:
        if all(over[i:i + sustain]):
            alarms += 1
            while i < n and over[i]:
                i += 1
        else:
            i += 1
    span = (refit_years[-1] - refit_years[0]) / 100.0
    return alarms / max(span, 1e-9)


def calibrate_threshold_refit(ctrl_years, ctrl_diag, orient=+1,
                              target=FPR_TARGET_PER_CENTURY):
    pcts = np.arange(80.0, 99.95, 0.5) if orient > 0 else np.arange(20.0, 0.05, -0.5)
    for p in pcts:
        thr = float(np.percentile(ctrl_diag, p))
        if fpr_per_century_refit(ctrl_years, ctrl_diag, thr, orient) <= target:
            return thr, p
    p = 99.9 if orient > 0 else 0.1
    return float(np.percentile(ctrl_diag, p)), p


def fov_causal_min_alarm(years, fov, ctrl_min: float, ctrl_sd: float,
                         rise_sustain=3, tipping=TIPPING_YEAR,
                         eval_start=EVAL_YEARS[0]):
    """Regla A3 (versión causal del indicador publicado, con retractación):
    Una alarma se emite en el primer año t >= eval_start de cada episodio en que
    (i) el mínimo corriente de FovS es inferior al mínimo del control y (ii) FovS se
    sitúa >= mínimo corriente + 1 sd del control durante >= rise_sustain años
    consecutivos. La alarma queda RETRACTADA si posteriormente ocurre un nuevo mínimo
    (el supuesto punto de giro no lo era). Se reporta la última alarma no retractada
    antes del tipping (lead) y el número de alarmas retractadas (costo del método).

    Retorna: (lead, año_alarma_final, n_retractadas, [años de alarmas retractadas])."""
    run_min, above = np.inf, 0
    pending = None          # alarma vigente aún no retractada
    retracted = []
    for y, v in zip(years, fov):
        if y >= tipping:
            break
        if v < run_min:
            run_min, above = v, 0
            if pending is not None:
                retracted.append(pending)
                pending = None
        elif v >= run_min + ctrl_sd:
            above += 1
        else:
            above = 0
        if (pending is None and y >= eval_start and run_min < ctrl_min
                and above >= rise_sustain):
            pending = float(y)
    if pending is not None:
        return float(tipping - pending), pending, len(retracted), retracted
    return None, None, len(retracted), retracted


def evaluate_E3(lead_mz, lead_ar1) -> str:
    """A1: semántica correcta del gate primario."""
    if lead_mz is None and lead_ar1 is None:
        return "not_evaluable"
    if lead_mz is None:
        return "failed"
    if lead_ar1 is None:
        return "passed"
    return "passed" if lead_mz >= lead_ar1 else "failed"
