"""Ingestión CESM (SA-AMOC-Collapse, Zenodo 10.5281/zenodo.10461549) -> contrato de estado v3.

Fuentes (CC-BY 4.0, van Westen, Kliphuis & Dijkstra 2024, Sci. Adv. 10:eadk1189):
  GitHub  : RenevanWesten/SA-AMOC-Collapse, tag SA-AMOC-Collapse_v1.0
  Zenodo  : https://doi.org/10.5281/zenodo.10461549

Este módulo NO ejecuta ningún diagnóstico de early warning. Solo construye, de forma
estrictamente reproducible, las series que el protocolo congelado consume:

  1. amoc_transport_26n : transporte meridional 0–1000 m en 26N [Sv], espejo exacto de
     Program/CESM/Ocean/AMOC_transport.py (celdas parciales de fondo incluidas).
  2. estado vectorial v3 en ~26N desde AMOC_structure: máximo de celda superior,
     profundidad del máximo, mínimo de retorno profundo, transportes medios superior y
     profundo, y el perfil 500–4500 m interpolado a la grilla de 33 niveles del
     contrato (insumo de las EOFs, que se ajustan aguas abajo por fold).
  3. fov_34s : índice F_ovS [Sv], espejo de Program/CESM/Ocean/FOV_index_34S.py
     (velocidad baroclínica por remoción de la media de sección; salinidad zonal
     ponderada por DXT; referencia 35 g/kg).
  4. forzamiento F_H(año) = 3e-4 * año [Sv] (rampa publicada; 0.66 Sv en el año 2200).

La paridad numérica con los scripts de referencia es un gate del protocolo (E3):
`verify_parity()` debe ejecutarse y archivarse antes de cualquier análisis.
"""
from __future__ import annotations

import hashlib
import io
import json
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset

RAW_BASE = (
    "https://raw.githubusercontent.com/RenevanWesten/SA-AMOC-Collapse/"
    "SA-AMOC-Collapse_v1.0/Data/CESM/Data"
)
# Rama de bajada de la histeresis (van Westen & Dijkstra 2023, GRL;
# Zenodo 10.5281/zenodo.10034589): mismos formatos, anios 2201-4400.
RAW_BASE_HYSTERESIS = (
    "https://raw.githubusercontent.com/RenevanWesten/GRL-AMOC-Hysteresis/"
    "GRL-AMOC-Hysteresis/Data/CESM/Data"
)
PERIODS = [(a, a + 49) for a in range(1, 2200, 50)]        # 0001-0050 ... 2151-2200
PERIODS_DOWN = [(a, a + 49) for a in range(2201, 4400, 50)]  # 2201-2250 ... 4351-4400
FORCING_RATE_SV_PER_YR = 3.0e-4                      # rampa de hosing publicada
CONTRACT_DEPTHS_M = np.linspace(500.0, 4500.0, 33)   # grilla del contrato v3


def period_filename(a: int, b: int) -> str:
    return f"CESM_year_{a:04d}-{b:04d}.nc"


def download(subdir: str, a: int, b: int, cache: Path, base: str = None) -> Path:
    """Descarga (con caché local y hash) un archivo NetCDF del release."""
    base = base or RAW_BASE
    cache.mkdir(parents=True, exist_ok=True)
    name = period_filename(a, b)
    out = cache / subdir / name
    if out.exists():
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    url = f"{base}/{subdir}/{name}"
    with urllib.request.urlopen(url) as r:
        data = r.read()
    out.write_bytes(data)
    (out.with_suffix(".sha256")).write_text(hashlib.sha256(data).hexdigest())
    return out


# ---------------------------------------------------------------------------
# 1. Transporte AMOC 26N (espejo de AMOC_transport.py)
# ---------------------------------------------------------------------------

def _partial_bottom_layers(dz: np.ndarray, mask2d: np.ndarray, column_depth: np.ndarray,
                           depth_top: np.ndarray, depth_max: float) -> ma.MaskedArray:
    """Espesores de capa con ajuste de celdas parciales de fondo (lógica de referencia)."""
    nz, nx = mask2d.shape
    layer = ma.masked_array(np.repeat(dz[:, None], nx, axis=1), mask=mask2d)
    for k in range(nz):
        depth_diff = np.sum(layer[: k + 1], axis=0) - column_depth
        if k == nz - 1:
            depth_diff = layer[k] - (depth_max - depth_top[k])
        depth_diff = ma.masked_where(depth_diff < 0, depth_diff).filled(0.0)
        layer[k] = layer[k] - depth_diff
    return layer


def amoc_transport_26n(periods: Iterable[tuple[int, int]], cache: Path,
                       depth_min: float = 0.0, depth_max: float = 1000.0,
                       base: str = None) -> tuple[np.ndarray, np.ndarray]:
    years_all, tr_all = [], []
    layer_field = None
    for (a, b) in periods:
        p = download("AMOC_section_26N", a, b, cache, base)
        with Dataset(p) as ds:
            depth = ds["z_t"][:]
            i0 = int(np.abs(depth - depth_min).argmin())
            i1 = int(np.abs(depth - depth_max).argmin()) + 1
            years = np.asarray(ds["year"][:])
            dz = ds["dz"][i0:i1]
            dxu = ds["DXU"][:]
            hu = ds["HU"][:]
            v = ds["VVEL"][:, i0:i1, 0]          # (año, z, lon)
        if layer_field is None:
            depth_top = np.zeros(i1 - i0)
            depth_top[1:] = np.cumsum(np.asarray(dz)[:-1])
            layer_field = _partial_bottom_layers(np.asarray(dz), v[0].mask,
                                                 np.asarray(hu)[0], depth_top, depth_max)
        for t in range(len(years)):
            transport = v[t] * layer_field * np.asarray(dxu)[0]
            years_all.append(float(years[t]))
            tr_all.append(float(np.sum(transport)) / 1e6)
    return np.asarray(years_all), np.asarray(tr_all)


# ---------------------------------------------------------------------------
# 2. Estado vectorial v3 desde AMOC_structure (~26N)
# ---------------------------------------------------------------------------

@dataclass
class StateConfig:
    target_lat: float = 26.0
    upper_band_m: tuple[float, float] = (0.0, 1000.0)
    deep_band_m: tuple[float, float] = (1500.0, 4500.0)
    profile_depths_m: np.ndarray = field(default_factory=lambda: CONTRACT_DEPTHS_M.copy())


def vector_state_26n(periods: Iterable[tuple[int, int]], cache: Path,
                     cfg: StateConfig = StateConfig(), base: str = None) -> dict[str, np.ndarray]:
    years_all, rows, profiles = [], [], []
    for (a, b) in periods:
        p = download("AMOC_structure", a, b, cache, base)
        with Dataset(p) as ds:
            depth = np.asarray(ds["depth"][:])
            lat = np.asarray(ds["lat"][:])
            years = np.asarray(ds["time"][:])
            ilat = int(np.abs(lat - cfg.target_lat).argmin())
            prof = np.asarray(ds["AMOC"][:, :, ilat])   # (año, z) [Sv]
        up = (depth >= cfg.upper_band_m[0]) & (depth <= cfg.upper_band_m[1])
        dp = (depth >= cfg.deep_band_m[0]) & (depth <= cfg.deep_band_m[1])
        for t in range(len(years)):
            p_t = prof[t]
            k = int(np.argmax(p_t))
            rows.append([
                float(np.max(p_t)),                # amoc_max [Sv]
                float(depth[k]),                   # depth_of_max [m]
                float(np.min(p_t[dp])),            # deep_return_min [Sv]
                float(np.mean(p_t[up])),           # upper_mean [Sv]
                float(np.mean(p_t[dp])),           # deep_mean [Sv]
            ])
            profiles.append(np.interp(cfg.profile_depths_m, depth, p_t))
            years_all.append(float(years[t]))
    rows = np.asarray(rows)
    return {
        "year": np.asarray(years_all),
        "amoc_max": rows[:, 0],
        "depth_of_max_m": rows[:, 1],
        "deep_return_min": rows[:, 2],
        "upper_mean": rows[:, 3],
        "deep_mean": rows[:, 4],
        "profile_500_4500": np.asarray(profiles),   # (n_años, 33)
        "profile_depths_m": cfg.profile_depths_m,
        "lat_used": float(cfg.target_lat),
    }


# ---------------------------------------------------------------------------
# 3. F_ovS (espejo de FOV_index_34S.py)
# ---------------------------------------------------------------------------

def fov_34s(periods: Iterable[tuple[int, int]], cache: Path,
            depth_max: float = 6000.0, s_ref: float = 35.0,
            base: str = None) -> tuple[np.ndarray, np.ndarray]:
    years_all, fov_all = [], []
    geo = None
    for (a, b) in periods:
        p = download("FOV_section_34S", a, b, cache, base)
        with Dataset(p) as ds:
            depth = np.asarray(ds["z_t"][:])
            i1 = int(np.abs(depth - depth_max).argmin()) + 1
            years = np.asarray(ds["year"][:])
            dz = np.asarray(ds["dz"][:i1])
            dxu = np.asarray(ds["DXU"][:])[0]
            dxt = np.asarray(ds["DXT"][:])
            hu = np.asarray(ds["HU"][:])[0]
            ht = np.asarray(ds["HT"][:])
            v = ds["VVEL"][:, :i1, 0]               # (año, z, lon_u)
            s = ds["SALT"][:, :i1]                  # (año, z, 2, lon_t)
        if geo is None:
            depth_top = np.zeros(i1); depth_top[1:] = np.cumsum(dz[:-1])
            layer_u = _partial_bottom_layers(dz, v[0].mask, hu, depth_top, depth_max)
            # normalización zonal en la grilla t, capa a capa
            nz = i1
            gxt_norm = ma.masked_all((nz, 2, dxt.shape[1]))
            lat_w = ma.masked_all((nz, 2))
            for k in range(nz):
                gxt_k = ma.masked_array(dxt, mask=s[0, k].mask)
                lat_w[k] = np.sum(gxt_k, axis=1) / np.sum(gxt_k)
                for j in range(2):
                    gxt_norm[k, j] = gxt_k[j] / np.sum(gxt_k[j])
            geo = (layer_u, gxt_norm, lat_w)
        layer_u, gxt_norm, lat_w = geo
        area_u = layer_u * dxu
        for t in range(len(years)):
            transport = v[t] * area_u
            vel_barotropic = np.sum(transport) / np.sum(area_u)
            vel_clin = v[t] - vel_barotropic
            transport_clin = np.sum(vel_clin * area_u, axis=1)          # (z,)
            salt_zonal = np.sum(s[t] * gxt_norm, axis=2) - s_ref        # (z, 2)
            salt_zonal = np.sum(salt_zonal * lat_w, axis=1)             # (z,)
            fov = (-1.0 / s_ref) * np.sum(transport_clin * salt_zonal) / 1e6
            years_all.append(float(years[t]))
            fov_all.append(float(fov))
    return np.asarray(years_all), np.asarray(fov_all)


# ---------------------------------------------------------------------------
# 4. Forzamiento y tabla de contrato
# ---------------------------------------------------------------------------

def hosing_forcing(years: np.ndarray) -> np.ndarray:
    """F_H [Sv]: rampa lineal de 3e-4 Sv/año (0.66 Sv en el año-modelo 2200)."""
    return FORCING_RATE_SV_PER_YR * np.asarray(years, dtype=float)


def build_contract_table(periods: Sequence[tuple[int, int]], cache: Path,
                         out_dir: Path, with_fov: bool = True, base: str = None) -> Path:
    """Construye y archiva la tabla anual del contrato + perfiles para EOFs."""
    out_dir.mkdir(parents=True, exist_ok=True)
    st = vector_state_26n(periods, cache, base=base)
    yrs_t, tr = amoc_transport_26n(periods, cache, base=base)
    assert np.allclose(st["year"], yrs_t), "años inconsistentes entre structure y section"
    cols = {
        "year": st["year"],
        "F_H_Sv": hosing_forcing(st["year"]),
        "amoc_transport_0_1000m_26N_Sv": tr,
        "amoc_max_Sv": st["amoc_max"],
        "depth_of_max_m": st["depth_of_max_m"],
        "deep_return_min_Sv": st["deep_return_min"],
        "upper_mean_Sv": st["upper_mean"],
        "deep_mean_Sv": st["deep_mean"],
    }
    if with_fov:
        yrs_f, fov = fov_34s(periods, cache, base=base)
        assert np.allclose(st["year"], yrs_f)
        cols["FovS_Sv"] = fov
    import pandas as pd
    df = pd.DataFrame(cols)
    csv = out_dir / "cesm_state_contract.csv"
    df.to_csv(csv, index=False)
    np.savez_compressed(out_dir / "cesm_profiles_500_4500.npz",
                        year=st["year"], depths_m=st["profile_depths_m"],
                        profile=st["profile_500_4500"])
    meta = {
        "source_doi": "10.5281/zenodo.10461549",
        "license": "CC-BY 4.0 (van Westen, Kliphuis & Dijkstra 2024)",
        "periods": list(map(list, periods)),
        "lat_used_deg": st["lat_used"],
        "forcing_rate_Sv_per_yr": FORCING_RATE_SV_PER_YR,
        "contract_depths_m": CONTRACT_DEPTHS_M.tolist(),
    }
    (out_dir / "INGESTION_METADATA.json").write_text(json.dumps(meta, indent=2))
    return csv
