#!/usr/bin/env python3
"""Post-process OCIM CO2 perturbation outputs and replicate Atwood-style figures.

Outputs (under $DFP_DATA/output_V6/postproc_atwood_style by default):
  - timeseries_dredging.parquet (and CSV fallback)
  - maps_cumulative_1996_2020.nc
  - maps_cumulative_1996_2050.nc
  - map_dpH_upper1000m_2020.nc
  - Fig2_timeseries_atwood_style.{pdf,svg,png}
  - Fig3_maps_atwood_style.{pdf,svg,png}
  - Fig4_dpH_surface_2020_atwood_style.{pdf,svg,png}
  - qc_postproc.md

Design goals:
  - Memory-safe: process one scenario file at a time.
  - Robust variable discovery for MATLAB v7/v7.3 files.
  - No interactive display required (headless-safe plotting).
"""

from __future__ import annotations

import argparse
import dataclasses
import glob
import os
import re
from contextlib import suppress
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import scipy.io as sio
import xarray as xr

# Force non-interactive backend before importing pyplot.
import matplotlib

matplotlib.use("Agg")
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, TwoSlopeNorm

try:
    import cartopy.crs as ccrs
except Exception:  # pragma: no cover - handled by runtime checks
    ccrs = None

try:
    import h5py
except Exception:  # pragma: no cover - handled by runtime checks
    h5py = None


RHO_SW = 1025.0
MOLAR_MASS_C = 12.011   # g mol-1 (IUPAC 2021); matches step7_export_jtrawl.R M_C
CO2_PER_C = 44.010 / 12.011  # g CO2 per g C (exact molecular masses)

SCENARIOS = {
    "baseline": "ocim_baseline_*.mat",
    "dredge_1x": "ocim_dredge_1x_*.mat",
    "dredge_10x": "ocim_dredge_10x_*.mat",
    "dredge_100x": "ocim_dredge_100x_*.mat",
}

# Named dredging hotspot regions [lon_min, lon_max, lat_min, lat_max].
# Extents are tightly sized around each dredging cluster at OCIM 2° resolution.
# Coordinates from an earlier step5 projection diagnostic (not in this repository).
HOTSPOT_REGIONS: Dict[str, Dict] = {
    "manila": {
        "label": "Manila Bay / Philippines",
        "extent": [118.0, 124.0, 12.0, 17.0],   # step5 ref: lon [120,121.5], lat [14,15.5]
        "dredge_centroid": [120.9, 14.5],
    },
    "abu_dhabi": {
        "label": "Abu Dhabi / Persian Gulf",
        "extent": [51.0, 57.0, 22.5, 27.5],     # step5 ref: lon [53.5,55.5], lat [24,25.5]
        "dredge_centroid": [54.5, 24.5],
    },
    "singapore": {
        "label": "Singapore / Strait of Malacca",
        "extent": [101.5, 106.5, -1.5, 4.5],    # step5 ref: lon [103.5,104.2], lat [1,1.6]
        "dredge_centroid": [103.8, 1.3],
    },
}


@dataclasses.dataclass
class Geometry:
    lon2d: np.ndarray
    lat2d: np.ndarray
    m3d: np.ndarray
    dzt3d: np.ndarray
    dxt3d: np.ndarray
    dyt3d: np.ndarray
    dzt: np.ndarray
    iocn_1based: np.ndarray

    @property
    def ni(self) -> int:
        return int(self.m3d.shape[0])

    @property
    def nj(self) -> int:
        return int(self.m3d.shape[1])

    @property
    def nk(self) -> int:
        return int(self.m3d.shape[2])

    @property
    def surface_mask(self) -> np.ndarray:
        return self.m3d[:, :, 0] == 1

    @property
    def surface_linear_fortran(self) -> np.ndarray:
        return np.flatnonzero(self.surface_mask.ravel(order="F"))

    @property
    def area2d(self) -> np.ndarray:
        return self.dxt3d[:, :, 0] * self.dyt3d[:, :, 0]

    @property
    def area_surface_vec(self) -> np.ndarray:
        return self.area2d.ravel(order="F")[self.surface_linear_fortran]

    @property
    def dz_surface_vec(self) -> np.ndarray:
        return self.dzt3d[:, :, 0].ravel(order="F")[self.surface_linear_fortran]

    @property
    def ocean_cell_volume_vec(self) -> np.ndarray:
        v3d = self.dxt3d * self.dyt3d * self.dzt3d
        iocn0 = self.iocn_1based.astype(int).ravel() - 1
        return v3d.ravel(order="F")[iocn0]


@dataclasses.dataclass
class ScenarioData:
    name: str
    year: np.ndarray
    pco2: np.ndarray
    jco2: np.ndarray
    ph: Optional[np.ndarray]
    ph_kind: str
    path: str


def _expand(path: str) -> str:
    return os.path.abspath(os.path.expanduser(path))


def _latest_file(pattern: str) -> Optional[str]:
    files = glob.glob(pattern)
    if not files:
        return None
    return max(files, key=os.path.getmtime)


def detect_latest_inputs(output_root: str) -> Dict[str, str]:
    out = {}
    for scenario, patt in SCENARIOS.items():
        p = _latest_file(os.path.join(output_root, patt))
        if not p:
            raise FileNotFoundError(f"Missing required scenario file pattern: {patt}")
        out[scenario] = p

    forcing = _latest_file(os.path.join(output_root, "jdredge_ocim2_48l_*.mat"))
    if not forcing:
        raise FileNotFoundError("Missing forcing file pattern: jdredge_ocim2_48l_*.mat")
    out["forcing"] = forcing
    return out


def _is_hdf5_mat(path: str) -> bool:
    if h5py is None:
        return False
    try:
        with h5py.File(path, "r"):
            return True
    except Exception:
        return False


def _matlab_reorder(arr: np.ndarray) -> np.ndarray:
    """Convert h5py read order to MATLAB-like axis order."""
    if arr.ndim <= 1:
        return arr
    axes = tuple(range(arr.ndim - 1, -1, -1))
    return np.transpose(arr, axes=axes)


def _as_1d_float(arr: np.ndarray) -> np.ndarray:
    return np.asarray(arr).astype(float).squeeze().reshape(-1)


def _decode_char_array(arr: np.ndarray) -> str:
    a = np.asarray(arr).squeeze()
    if a.dtype.kind in ("U", "S"):
        return "".join(a.astype(str).reshape(-1)).strip()
    if a.dtype.kind in ("i", "u"):
        chars = [chr(int(c)) for c in a.reshape(-1) if int(c) > 0]
        return "".join(chars).strip()
    return str(a).strip()


class MatReader:
    """Read MATLAB v7/v7.3 with minimal memory usage for selected variables."""

    def __init__(self, path: str):
        self.path = path
        self._is_hdf5 = _is_hdf5_mat(path)
        self._h5 = None
        self._mat = None

    def __enter__(self) -> "MatReader":
        if self._is_hdf5:
            self._h5 = h5py.File(self.path, "r")
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        if self._h5 is not None:
            self._h5.close()
            self._h5 = None

    def has(self, key: str) -> bool:
        if self._is_hdf5:
            return key in self._h5
        if self._mat is None:
            self._mat = sio.loadmat(self.path, struct_as_record=False, squeeze_me=False)
        return key in self._mat

    def read(self, key: str) -> np.ndarray:
        if self._is_hdf5:
            node = self._h5[key]
            if isinstance(node, h5py.Dataset):
                if node.dtype.kind == "O":
                    if node.size == 1:
                        ref = node[tuple([0] * node.ndim)]
                        return self._read_ref(ref)
                    raise ValueError(f"Unsupported object dataset for key '{key}'")
                arr = np.array(node)
                return _matlab_reorder(arr)
            raise ValueError(f"Key '{key}' is not a dataset")

        if self._mat is None:
            self._mat = sio.loadmat(self.path, struct_as_record=False, squeeze_me=False)
        return np.array(self._mat[key])

    def _read_ref(self, ref) -> np.ndarray:
        arr = np.array(self._h5[ref])
        return _matlab_reorder(arr)

    def read_string(self, key: str) -> str:
        arr = self.read(key)
        return _decode_char_array(arr)

def _find_first(reader: MatReader, names: Sequence[str]) -> Optional[str]:
    for n in names:
        if reader.has(n):
            return n
    return None


def _orient_surface_time(arr: np.ndarray, nt: int) -> np.ndarray:
    a = np.asarray(arr).squeeze()
    if a.ndim == 1:
        if a.shape[0] == nt:
            return a.reshape(1, nt)
        raise ValueError(f"Cannot orient 1D array with length {a.shape[0]} to nt={nt}")
    if a.ndim != 2:
        raise ValueError(f"Expected 2D array for surface-time variable, got shape {a.shape}")
    if a.shape[1] == nt:
        return a
    if a.shape[0] == nt:
        return a.T
    raise ValueError(f"Cannot orient array shape {a.shape} with nt={nt}")


def _infer_year(nt: int, year_arr: Optional[np.ndarray], t_arr: Optional[np.ndarray]) -> np.ndarray:
    if year_arr is not None:
        y = _as_1d_float(year_arr)
        if y.size == nt:
            return y
    else:
        y = None

    if t_arr is not None:
        t = _as_1d_float(t_arr)
        if t.size == nt:
            if np.nanmax(t) < 1000:
                y0 = 1780.0
                if y is not None and y.size > 0:
                    y0 = float(y[0])
                return y0 + t
            return t

    if y is not None and y.size == nt - 1:
        return np.concatenate([y, [y[-1] + (y[-1] - y[-2] if y.size > 1 else 1.0)]])

    return np.arange(1780.0, 1780.0 + nt, 1.0)


def load_scenario(path: str, expected_year: Optional[np.ndarray] = None) -> ScenarioData:
    with MatReader(path) as reader:
        pco2_key = _find_first(reader, ["pco2a", "pCO2a", "pCO2atm", "pco2atm"])
        if pco2_key is None:
            raise KeyError(f"{path}: missing pCO2 variable")
        pco2 = _as_1d_float(reader.read(pco2_key))
        nt = pco2.size

        year_key = _find_first(reader, ["year", "years"])
        t_key = _find_first(reader, ["t_sim", "t"])
        year = _infer_year(nt, reader.read(year_key) if year_key else None, reader.read(t_key) if t_key else None)

        if expected_year is not None:
            if year.size != expected_year.size or not np.allclose(year, expected_year):
                raise ValueError(f"{path}: year axis mismatch with baseline")

        jco2_key = _find_first(reader, ["Jco2", "JCO2", "J_co2"])
        if jco2_key is None:
            raise KeyError(f"{path}: missing Jco2 variable")
        jco2 = _orient_surface_time(reader.read(jco2_key), nt)

        ph_key = _find_first(reader, ["pH_out", "pH", "ph", "pHsurf"])
        ph = None
        ph_kind = "missing"
        if ph_key is not None:
            ph_candidate = np.asarray(reader.read(ph_key)).squeeze()
            if ph_candidate.ndim == 2:
                ph = _orient_surface_time(ph_candidate, nt)
                ph_kind = "surface"
            elif ph_candidate.ndim == 1 and ph_candidate.size == nt:
                ph = ph_candidate.reshape(1, nt)
                ph_kind = "single_series"

        name_key = _find_first(reader, ["scenario_name", "scenario"])
        if name_key is not None:
            try:
                name = reader.read_string(name_key)
            except Exception:
                name = os.path.basename(path)
        else:
            m = re.search(r"ocim_(.*?)_\d{8}_\d{6}\.mat$", os.path.basename(path))
            name = m.group(1) if m else os.path.basename(path)

    return ScenarioData(
        name=name,
        year=year,
        pco2=pco2,
        jco2=jco2,
        ph=ph,
        ph_kind=ph_kind,
        path=path,
    )


def _load_ocim_ctl_geometry(ocim_ctl_path: str) -> Geometry:
    try:
        s = sio.loadmat(ocim_ctl_path, struct_as_record=False, squeeze_me=False)
        output = s["output"][0, 0]
        grid = output.grid[0, 0]
        m3d = np.array(output.M3d[0, 0], dtype=float)
        lon2d = np.array(grid.XT3d[0, 0][:, :, 0], dtype=float)
        lat2d = np.array(grid.YT3d[0, 0][:, :, 0], dtype=float)
        dzt3d = np.array(grid.DZT3d[0, 0], dtype=float)
        dxt3d = np.array(grid.DXT3d[0, 0], dtype=float)
        dyt3d = np.array(grid.DYT3d[0, 0], dtype=float)
        dzt = _as_1d_float(np.array(grid.dzt[0, 0], dtype=float))
    except NotImplementedError:
        if h5py is None:
            raise
        with h5py.File(ocim_ctl_path, "r") as f:
            output = f["output"]
            grid_ref = output["grid"]
            m3d = _matlab_reorder(np.array(output["M3d"]))

            def deref(name: str) -> np.ndarray:
                ref = grid_ref[name]
                if ref.shape == (1, 1):
                    return _matlab_reorder(np.array(f[ref[0, 0]]))
                return _matlab_reorder(np.array(ref))

            xt3d = deref("XT3d")
            yt3d = deref("YT3d")
            dzt3d = deref("DZT3d")
            dxt3d = deref("DXT3d")
            dyt3d = deref("DYT3d")
            dzt = _as_1d_float(deref("dzt"))
            lon2d = xt3d[:, :, 0]
            lat2d = yt3d[:, :, 0]

    iocn = np.flatnonzero(m3d.ravel(order="F") == 1) + 1
    return Geometry(
        lon2d=lon2d,
        lat2d=lat2d,
        m3d=m3d,
        dzt3d=dzt3d,
        dxt3d=dxt3d,
        dyt3d=dyt3d,
        dzt=dzt,
        iocn_1based=iocn.reshape(-1, 1),
    )


def _load_ocim_cache_geometry(cache_path: str) -> Geometry:
    with MatReader(cache_path) as reader:
        lon2d = np.asarray(reader.read("lon2d"), dtype=float)
        lat2d = np.asarray(reader.read("lat2d"), dtype=float)
        m3d = np.asarray(reader.read("M3d"), dtype=float)
        dzt3d = np.asarray(reader.read("DZT3d"), dtype=float)
        dxt3d = np.asarray(reader.read("DXT3d"), dtype=float)
        dyt3d = np.asarray(reader.read("DYT3d"), dtype=float)
        dzt = _as_1d_float(np.asarray(reader.read("dzt"), dtype=float))
        if reader.has("iocn"):
            iocn = np.asarray(reader.read("iocn"), dtype=float).astype(int)
        else:
            iocn = (np.flatnonzero(m3d.ravel(order="F") == 1) + 1).reshape(-1, 1)

    return Geometry(
        lon2d=lon2d,
        lat2d=lat2d,
        m3d=m3d,
        dzt3d=dzt3d,
        dxt3d=dxt3d,
        dyt3d=dyt3d,
        dzt=dzt,
        iocn_1based=iocn,
    )


def load_geometry(ocim_dir: str) -> Geometry:
    cache_path = os.path.join(ocim_dir, "ocim_cache.mat")
    if os.path.isfile(cache_path):
        print(f"[INFO] Geometry source: {cache_path}")
        return _load_ocim_cache_geometry(cache_path)

    ocim_ctl = os.path.join(ocim_dir, "OCIM2_48L_CTL.mat")
    if not os.path.isfile(ocim_ctl):
        raise FileNotFoundError(f"Missing geometry inputs: {cache_path} and {ocim_ctl}")
    print(f"[INFO] Geometry source: {ocim_ctl}")
    return _load_ocim_ctl_geometry(ocim_ctl)


def _dt_years(year: np.ndarray) -> np.ndarray:
    year = _as_1d_float(year)
    if year.size == 1:
        return np.array([1.0])
    d = np.diff(year)
    if np.allclose(d, d[0]):
        return np.full_like(year, float(d[0]))
    return np.concatenate([[d[0]], d])


def _surface_vec_to_map(vec: np.ndarray, geom: Geometry) -> np.ndarray:
    out = np.full(geom.ni * geom.nj, np.nan, dtype=float)
    out[geom.surface_linear_fortran] = vec
    return out.reshape((geom.ni, geom.nj), order="F")


def jco2_to_flux_mol_m2(jco2_ns_nt: np.ndarray, geom: Geometry) -> np.ndarray:
    """Convert Jco2 (umol kg-1 yr-1) to air-sea flux (mol m-2 yr-1), ocean->atm positive."""
    dz = geom.dz_surface_vec.reshape(-1, 1)
    return -jco2_ns_nt * RHO_SW * dz * 1e-6


def global_flux_gtc_yr(jco2_ns_nt: np.ndarray, geom: Geometry) -> np.ndarray:
    flux_mol_m2 = jco2_to_flux_mol_m2(jco2_ns_nt, geom)
    mol_yr = np.sum(flux_mol_m2 * geom.area_surface_vec.reshape(-1, 1), axis=0)
    return mol_yr * MOLAR_MASS_C / 1e15


def global_surface_mean(field_ns_nt: np.ndarray, geom: Geometry) -> np.ndarray:
    w = geom.area_surface_vec
    wsum = np.sum(w)
    if wsum <= 0:
        raise ValueError("Invalid surface area weights")
    return np.sum(field_ns_nt * w.reshape(-1, 1), axis=0) / wsum


def _window_mask(year: np.ndarray, y0: int, y1: int) -> np.ndarray:
    return (year >= y0) & (year <= y1)


def _hotspot_bbox(
    lon2d: np.ndarray,
    lat2d: np.ndarray,
    field2d: np.ndarray,
    frac_threshold: float = 0.05,
    padding_deg: float = 4.0,
) -> Optional[Tuple[float, float, float, float]]:
    """Return (lon_min, lon_max, lat_min, lat_max) bounding box of cells with signal
    exceeding frac_threshold * max(|field|), padded by padding_deg degrees."""
    vmax = float(np.nanmax(np.abs(field2d)))
    if vmax <= 0:
        return None
    mask = np.abs(field2d) > frac_threshold * vmax
    if not np.any(mask):
        return None
    lons = lon2d[mask]
    lats = lat2d[mask]
    lon_min = float(np.min(lons)) - padding_deg
    lon_max = float(np.max(lons)) + padding_deg
    lat_min = float(np.min(lats)) - padding_deg
    lat_max = float(np.max(lats)) + padding_deg
    return (
        max(lon_min, -180.0),
        min(lon_max, 180.0),
        max(lat_min, -90.0),
        min(lat_max, 90.0),
    )


def _region_has_signal(
    lon2d: np.ndarray,
    lat2d: np.ndarray,
    field2d: np.ndarray,
    extent: List[float],
    frac_threshold: float = 0.001,
) -> Tuple[bool, float]:
    """Return (has_signal, peak_fraction) for a geographic region.

    peak_fraction = max(|field|) inside region / max(|field|) globally.
    """
    lon_min, lon_max, lat_min, lat_max = extent
    j_mask = (lon2d[0, :] >= lon_min) & (lon2d[0, :] <= lon_max)
    i_mask = (lat2d[:, 0] >= lat_min) & (lat2d[:, 0] <= lat_max)
    region = field2d[np.ix_(i_mask, j_mask)]
    global_max = float(np.nanmax(np.abs(field2d)))
    if global_max <= 0:
        return False, 0.0
    region_max = float(np.nanmax(np.abs(region))) if region.size > 0 else 0.0
    frac = region_max / global_max
    return frac > frac_threshold, frac


def _sort_longitude_for_plot(lon2d: np.ndarray, field2d: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    # Assumes regular grid: all rows share the same longitude vector (lon2d[i,:] == lon2d[0,:]).
    lon = lon2d[0, :].copy()
    lon_wrapped = ((lon + 180.0) % 360.0) - 180.0
    idx = np.argsort(lon_wrapped)
    return lon_wrapped[idx], idx, field2d[:, idx]


def write_maps_netcdf(
    path: str,
    lon2d: np.ndarray,
    lat2d: np.ndarray,
    cum_emissions: np.ndarray,
    cum_df: np.ndarray,
    attrs: Dict[str, str],
) -> None:
    lon = lon2d[0, :]   # regular grid: all rows share same longitude vector
    lat = lat2d[:, 0]   # regular grid: all columns share same latitude vector
    ds = xr.Dataset(
        data_vars={
            "cum_emissions_pos_mol_m2": (("lat", "lon"), cum_emissions),
            "cum_dF_mol_m2": (("lat", "lon"), cum_df),
        },
        coords={"lat": lat, "lon": lon},
        attrs=attrs,
    )
    ds["cum_emissions_pos_mol_m2"].attrs["units"] = "mol m-2"
    ds["cum_dF_mol_m2"].attrs["units"] = "mol m-2"
    _safe_to_netcdf(ds, path)


def write_dph_netcdf(
    path: str,
    lon2d: np.ndarray,
    lat2d: np.ndarray,
    dph_map: np.ndarray,
    method_note: str,
    scenario: str,
) -> None:
    lon = lon2d[0, :]   # regular grid: all rows share same longitude vector
    lat = lat2d[:, 0]   # regular grid: all columns share same latitude vector
    ds = xr.Dataset(
        data_vars={"dpH_upper1000m_2020": (("lat", "lon"), dph_map)},
        coords={"lat": lat, "lon": lon},
        attrs={
            "scenario": scenario,
            "method": method_note,
            "note": (
                "Content is the SURFACE pH anomaly (OCIM layer k=1, ~50 m top layer). "
                "Variable named 'upper1000m' for consistency with Atwood et al. (2021) "
                "figure conventions; OCIM2-48L does not output depth-integrated pH."
            ),
        },
    )
    ds["dpH_upper1000m_2020"].attrs["units"] = "1"
    ds["dpH_upper1000m_2020"].attrs["long_name"] = "Surface pH anomaly due to dredging (OCIM layer k=1)"
    _safe_to_netcdf(ds, path)


def _safe_to_netcdf(ds: xr.Dataset, path: str) -> str:
    """Write NetCDF robustly, avoiding partial files and HDF backend issues.

    Strategy:
      1) Try default engine (typically netcdf4/h5netcdf if installed).
      2) Fallback to scipy engine (NetCDF3, no HDF5 dependency).
    """
    tmp_path = path + ".tmp"
    last_exc: Optional[Exception] = None
    engines: List[Optional[str]] = [None, "scipy"]

    for eng in engines:
        with suppress(FileNotFoundError):
            os.remove(tmp_path)
        try:
            kwargs = {}
            if eng is not None:
                kwargs["engine"] = eng
            ds.to_netcdf(tmp_path, **kwargs)
            os.replace(tmp_path, path)
            return eng or "default"
        except Exception as exc:
            last_exc = exc
            with suppress(FileNotFoundError):
                os.remove(tmp_path)
            continue

    raise RuntimeError(f"Failed to write NetCDF '{path}' with all engines") from last_exc


def _save_figure(fig: plt.Figure, out_base: str, dpi: int) -> None:
    fig.savefig(out_base + ".pdf", bbox_inches="tight")
    fig.savefig(out_base + ".svg", bbox_inches="tight")
    fig.savefig(out_base + ".png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def make_fig2(year: np.ndarray, ts_df: pd.DataFrame, out_dir: str, dpi: int) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True)
    scenarios = [("1x", "#1b9e77"), ("10x", "#d95f02"), ("100x", "#7570b3")]

    for tag, col in scenarios:
        axes[0].plot(year, ts_df[f"delta_cum_uptake_PgCO2_{tag}"], lw=2, color=col, label=f"dredge_{tag}")
    axes[0].axhline(0, color="0.2", lw=0.8)
    axes[0].set_ylabel("Delta cumulative ocean CO2 uptake (Pg CO2)")
    axes[0].set_title("A")

    for tag, col in scenarios:
        axes[1].plot(year, ts_df[f"delta_pco2_ppm_{tag}"], lw=2, color=col, label=f"dredge_{tag}")
    axes[1].axhline(0, color="0.2", lw=0.8)
    axes[1].set_ylabel("Delta atmospheric CO2 (ppm)")
    axes[1].set_title("B")

    for tag, col in scenarios:
        axes[2].plot(year, ts_df[f"delta_pH_global_{tag}"], lw=2, color=col, label=f"dredge_{tag}")
    axes[2].axhline(0, color="0.2", lw=0.8)
    axes[2].set_ylabel("Delta global ocean pH")
    axes[2].set_title("C")

    for ax in axes:
        for yy in (1996, 2020, 2050):
            ax.axvline(yy, color="0.5", ls="--", lw=0.7)
        ax.grid(alpha=0.25, ls=":", lw=0.6)

    axes[2].set_xlabel("Year")
    axes[0].legend(ncol=3, loc="best", frameon=False)
    fig.suptitle("Atwood-style Figure 2 (OCIM-only deltas vs baseline)", y=0.995)

    _save_figure(fig, os.path.join(out_dir, "Fig2_timeseries_atwood_style"), dpi=dpi)


def _require_cartopy() -> None:
    if ccrs is None:
        raise RuntimeError("cartopy is required for map figures but is not installed.")


def make_fig3(
    lon2d: np.ndarray,
    lat2d: np.ndarray,
    map_2020: Dict[str, np.ndarray],
    map_2050: Dict[str, np.ndarray],
    out_dir: str,
    dpi: int,
    scenario_label: str,
) -> None:
    _require_cartopy()

    proj = ccrs.Robinson()
    pc = ccrs.PlateCarree()
    fig, axes = plt.subplots(
        2, 2, figsize=(13, 7), subplot_kw={"projection": proj}, constrained_layout=True
    )

    left_fields = [map_2020["cum_emissions"], map_2050["cum_emissions"]]
    right_fields = [map_2020["cum_df"], map_2050["cum_df"]]
    vmax_left = np.nanmax([np.nanmax(f) for f in left_fields])
    vmax_right = np.nanmax([np.nanmax(np.abs(f)) for f in right_fields])
    vmax_left = max(vmax_left, 1e-12)
    vmax_right = max(vmax_right, 1e-12)

    rows = [
        ("1996-2020 cumulative", map_2020),
        ("1996-2050 cumulative", map_2050),
    ]
    left_mappable = None
    right_mappable = None

    for i, (row_title, maps) in enumerate(rows):
        lon_sorted, _, left_plot = _sort_longitude_for_plot(lon2d, maps["cum_emissions"])
        _, _, right_plot = _sort_longitude_for_plot(lon2d, maps["cum_df"])
        lat = lat2d[:, 0]
        lon_mesh, lat_mesh = np.meshgrid(lon_sorted, lat)

        ax_l = axes[i, 0]
        left_mappable = ax_l.pcolormesh(
            lon_mesh,
            lat_mesh,
            left_plot,
            transform=pc,
            cmap="YlOrRd",
            vmin=0.0,
            vmax=vmax_left,
            shading="auto",
        )
        ax_l.coastlines(linewidth=0.5)
        ax_l.set_title(f"{row_title}\nCumulative emissions (positive only)")

        ax_r = axes[i, 1]
        right_mappable = ax_r.pcolormesh(
            lon_mesh,
            lat_mesh,
            right_plot,
            transform=pc,
            cmap="RdBu_r",
            norm=TwoSlopeNorm(vmin=-vmax_right, vcenter=0.0, vmax=vmax_right),
            shading="auto",
        )
        ax_r.coastlines(linewidth=0.5)
        ax_r.set_title(f"{row_title}\nCumulative Delta air-sea CO2 flux (signed)")

    cbar_l = fig.colorbar(left_mappable, ax=axes[:, 0], orientation="horizontal", shrink=0.9, pad=0.04)
    cbar_l.set_label("mol m-2")
    cbar_r = fig.colorbar(right_mappable, ax=axes[:, 1], orientation="horizontal", shrink=0.9, pad=0.04)
    cbar_r.set_label("mol m-2")

    fig.suptitle(f"Atwood-style Figure 3 ({scenario_label})", y=1.01)
    _save_figure(fig, os.path.join(out_dir, "Fig3_maps_atwood_style"), dpi=dpi)


def make_fig4(
    lon2d: np.ndarray,
    lat2d: np.ndarray,
    dph_map: np.ndarray,
    out_dir: str,
    dpi: int,
    scenario_label: str,
    fig4_title: str,
    fig4_stem: str,
) -> None:
    _require_cartopy()

    proj = ccrs.Robinson()
    pc = ccrs.PlateCarree()
    fig = plt.figure(figsize=(11, 5.2))
    ax = plt.axes(projection=proj)

    lon_sorted, _, z_plot = _sort_longitude_for_plot(lon2d, dph_map)
    lat = lat2d[:, 0]
    lon_mesh, lat_mesh = np.meshgrid(lon_sorted, lat)

    vmax = np.nanmax(np.abs(z_plot))
    vmax = max(vmax, 1e-8)
    m = ax.pcolormesh(
        lon_mesh,
        lat_mesh,
        z_plot,
        transform=pc,
        cmap="RdBu_r",
        norm=TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax),
        shading="auto",
    )
    ax.coastlines(linewidth=0.6)
    cb = plt.colorbar(m, orientation="horizontal", pad=0.06, shrink=0.9)
    cb.set_label("Delta pH")
    ax.set_title(f"{fig4_title} ({scenario_label})")

    _save_figure(fig, os.path.join(out_dir, fig4_stem), dpi=dpi)


def make_fig_hotspot_zoom(
    lon2d: np.ndarray,
    lat2d: np.ndarray,
    map_2020: Dict[str, np.ndarray],
    map_2050: Dict[str, np.ndarray],
    dph_map: np.ndarray,
    out_dir: str,
    dpi: int,
    scenario_label: str,
) -> None:
    """Produce zoomed versions of Fig3 and Fig4 centered on the dredging hotspot.

    Bounding box is auto-detected from the 1996-2020 cumulative emissions signal
    (cells > 0.5% of peak), padded by 10 degrees.  Uses PlateCarree projection
    (no Robinson) for correct regional rendering.
    """
    _require_cartopy()

    bbox = _hotspot_bbox(lon2d, lat2d, map_2020["cum_emissions"], frac_threshold=0.005, padding_deg=10.0)
    if bbox is None:
        print("[WARN] make_fig_hotspot_zoom: no hotspot signal detected, skipping zoom figures")
        return

    lon_min, lon_max, lat_min, lat_max = bbox
    print(f"[INFO] hotspot bbox: lon [{lon_min:.0f}, {lon_max:.0f}]  lat [{lat_min:.0f}, {lat_max:.0f}]")

    pc = ccrs.PlateCarree()
    # Regular grid: all rows share the same lon/lat vectors
    lon_plot = lon2d[0, :]
    lat_plot = lat2d[:, 0]
    lon_mesh, lat_mesh = np.meshgrid(lon_plot, lat_plot)

    vmax_em = max(np.nanmax(map_2020["cum_emissions"]), np.nanmax(map_2050["cum_emissions"]), 1e-12)
    vmax_df = max(
        np.nanmax(np.abs(map_2020["cum_df"])),
        np.nanmax(np.abs(map_2050["cum_df"])),
        1e-12,
    )
    vmax_ph = max(float(np.nanmax(np.abs(dph_map))), 1e-8)

    # ── Fig3 zoom: 2×2 (emissions + delta-flux × 2 time windows) ──────────
    fig_a, axes_a = plt.subplots(
        2, 2, figsize=(12, 8),
        subplot_kw={"projection": pc},
        constrained_layout=True,
    )
    rows = [("1996\u20132020", map_2020), ("1996\u20132050", map_2050)]
    for i, (row_title, maps) in enumerate(rows):
        ax_l = axes_a[i, 0]
        ax_l.set_extent([lon_min, lon_max, lat_min, lat_max], crs=pc)
        m_l = ax_l.pcolormesh(
            lon_mesh, lat_mesh, maps["cum_emissions"],
            transform=pc, cmap="YlOrRd", vmin=0.0, vmax=vmax_em, shading="auto",
        )
        ax_l.coastlines(linewidth=0.8)
        gl = ax_l.gridlines(draw_labels=True, linewidth=0.4, color="gray", alpha=0.5, linestyle="--")
        gl.top_labels = False
        gl.right_labels = False
        ax_l.set_title(f"{row_title} — Cumulative emissions (mol m\u207b\u00b2)")
        fig_a.colorbar(m_l, ax=ax_l, orientation="horizontal", shrink=0.85, pad=0.12, label="mol m\u207b\u00b2")

        ax_r = axes_a[i, 1]
        ax_r.set_extent([lon_min, lon_max, lat_min, lat_max], crs=pc)
        m_r = ax_r.pcolormesh(
            lon_mesh, lat_mesh, maps["cum_df"],
            transform=pc, cmap="RdBu_r",
            norm=TwoSlopeNorm(vmin=-vmax_df, vcenter=0.0, vmax=vmax_df),
            shading="auto",
        )
        ax_r.coastlines(linewidth=0.8)
        gl2 = ax_r.gridlines(draw_labels=True, linewidth=0.4, color="gray", alpha=0.5, linestyle="--")
        gl2.top_labels = False
        gl2.right_labels = False
        ax_r.set_title(f"{row_title} — Cumulative \u0394F air-sea CO\u2082 (mol m\u207b\u00b2)")
        fig_a.colorbar(m_r, ax=ax_r, orientation="horizontal", shrink=0.85, pad=0.12, label="mol m\u207b\u00b2")

    fig_a.suptitle(
        f"Hotspot zoom — {scenario_label}  "
        f"(lon {lon_min:.0f}\u00b0\u2013{lon_max:.0f}\u00b0, lat {lat_min:.0f}\u00b0\u2013{lat_max:.0f}\u00b0)",
        y=1.01,
    )
    _save_figure(fig_a, os.path.join(out_dir, "Fig3_hotspot_zoom_atwood_style"), dpi=dpi)

    # ── Fig4 zoom: pH map ──────────────────────────────────────────────────
    fig_b = plt.figure(figsize=(9, 6))
    ax_b = plt.axes(projection=pc)
    ax_b.set_extent([lon_min, lon_max, lat_min, lat_max], crs=pc)
    m_b = ax_b.pcolormesh(
        lon_mesh, lat_mesh, dph_map,
        transform=pc, cmap="RdBu_r",
        norm=TwoSlopeNorm(vmin=-vmax_ph, vcenter=0.0, vmax=vmax_ph),
        shading="auto",
    )
    ax_b.coastlines(linewidth=0.8)
    gl_b = ax_b.gridlines(draw_labels=True, linewidth=0.4, color="gray", alpha=0.5, linestyle="--")
    gl_b.top_labels = False
    gl_b.right_labels = False
    cb_b = plt.colorbar(m_b, orientation="horizontal", pad=0.1, shrink=0.85)
    cb_b.set_label("\u0394pH")
    ax_b.set_title(
        f"\u0394pH surface 2020 — hotspot zoom — {scenario_label}  "
        f"(lon {lon_min:.0f}\u00b0\u2013{lon_max:.0f}\u00b0, lat {lat_min:.0f}\u00b0\u2013{lat_max:.0f}\u00b0)"
    )
    _save_figure(fig_b, os.path.join(out_dir, "Fig4_hotspot_zoom_dpH_surface_2020_atwood_style"), dpi=dpi)
    print("[INFO] hotspot zoom figures saved (PDF/SVG/PNG)")


def make_fig_named_hotspots(
    lon2d: np.ndarray,
    lat2d: np.ndarray,
    cum_em: np.ndarray,
    cum_df: np.ndarray,
    dph_map: np.ndarray,
    out_dir: str,
    dpi: int,
    scenario_label: str,
    regions: Optional[Dict[str, Dict]] = None,
) -> None:
    """3×3 figure: one row per named hotspot, three columns (emissions / delta-flux / delta-pH).

    Rows are ordered by peak signal fraction (strongest hotspot first).
    Regions with no detectable signal are skipped with a warning.
    """
    _require_cartopy()
    if regions is None:
        regions = HOTSPOT_REGIONS

    pc = ccrs.PlateCarree()
    lon_plot = lon2d[0, :]   # regular grid
    lat_plot = lat2d[:, 0]
    lon_mesh, lat_mesh = np.meshgrid(lon_plot, lat_plot)

    # Check signal in each region and rank by peak fraction.
    # All regions are always shown; inactive ones are labelled "(no signal at 1x)".
    active: List[Tuple[float, str, Dict]] = []
    for key, region in regions.items():
        has_sig, frac = _region_has_signal(lon2d, lat2d, cum_em, region["extent"])
        status = "ACTIVE" if has_sig else "no signal"
        print(f"[INFO] named hotspot '{key}': peak fraction = {frac:.4f}  ({status})")
        active.append((frac, key, region))
    active.sort(reverse=True, key=lambda x: x[0])

    n = len(active)
    vmax_em = max(float(np.nanmax(cum_em)), 1e-12)
    vmax_df = max(float(np.nanmax(np.abs(cum_df))), 1e-12)
    vmax_ph = max(float(np.nanmax(np.abs(dph_map))), 1e-8)

    col_cfgs = [
        ("Cumul. emissions 1996\u20132020\n(mol m\u207b\u00b2)", cum_em,  "YlOrRd", vmax_em, "mol m\u207b\u00b2", False),
        ("Cumul. \u0394F air-sea CO\u2082 1996\u20132020\n(mol m\u207b\u00b2)", cum_df, "RdBu_r", vmax_df, "mol m\u207b\u00b2", True),
        ("\u0394pH surface 2020",            dph_map, "RdBu_r", vmax_ph, "\u0394pH",       True),
    ]

    fig, axes = plt.subplots(
        n, 3, figsize=(15, 5 * n),
        subplot_kw={"projection": pc},
        constrained_layout=True,
    )
    if n == 1:
        axes = axes[np.newaxis, :]

    for i, (frac, key, region) in enumerate(active):
        extent = region["extent"]
        has_sig = frac > 0.001
        row_label = region["label"] + ("" if has_sig else "\n(no signal at 1x — dredging present)")
        for j, (col_title, field, cmap, vmax, cbar_label, diverging) in enumerate(col_cfgs):
            ax = axes[i, j]
            ax.set_extent(extent, crs=pc)
            if diverging:
                norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
                m = ax.pcolormesh(lon_mesh, lat_mesh, field, transform=pc,
                                  cmap=cmap, norm=norm, shading="auto")
            else:
                m = ax.pcolormesh(lon_mesh, lat_mesh, field, transform=pc,
                                  cmap=cmap, vmin=0.0, vmax=vmax, shading="auto")
            ax.coastlines(linewidth=0.8)
            gl = ax.gridlines(draw_labels=True, linewidth=0.4,
                              color="gray", alpha=0.5, linestyle="--")
            gl.top_labels = False
            gl.right_labels = False
            fig.colorbar(m, ax=ax, orientation="horizontal",
                         shrink=0.85, pad=0.12, label=cbar_label)
            if i == 0:
                ax.set_title(col_title, fontsize=10)
            if j == 0:
                ax.set_ylabel(f"{row_label}\n(peak {100*frac:.1f}% of global max)",
                              fontsize=9)

    fig.suptitle(
        f"Named dredging hotspots \u2014 {scenario_label}  "
        f"(ranked by signal intensity)",
        y=1.01, fontsize=12,
    )
    _save_figure(fig, os.path.join(out_dir, "Fig_hotspots_named_atwood_style"), dpi=dpi)
    print("[INFO] named hotspot figure saved (PDF/SVG/PNG)")


def make_fig_forcing_diagnostic(
    forcing_vec: np.ndarray,
    geom: Geometry,
    out_dir: str,
    dpi: int,
    regions: Optional[Dict[str, Dict]] = None,
) -> None:
    """Diagnostic map: WHERE are the Jdredge forcing cells on the OCIM grid?

    2x2 layout:
      (a) Top-left:     Global Robinson — land gray, ocean white, forcing cells colored
      (b) Top-right:    Singapore zoom (strongest signal)
      (c) Bottom-left:  Manila zoom
      (d) Bottom-right: Abu Dhabi zoom (0 forcing cells — all land at 2°)

    Each zoom panel shows gray-filled 2° land cells, white ocean, colored forcing
    cells, 2° grid lines, a red star at the step5 dredging centroid, and a text
    annotation explaining whether the centroid falls on an OCIM land or ocean cell.
    """
    _require_cartopy()
    if regions is None:
        regions = HOTSPOT_REGIONS

    # --- Map forcing vector indices to (i, j, k) on the OCIM grid ---
    ni, nj, nk = geom.ni, geom.nj, geom.nk
    iocn0 = geom.iocn_1based.ravel().astype(int) - 1  # 0-based Fortran-linear
    active_mask = np.abs(forcing_vec) > 0
    n_active = int(np.sum(active_mask))
    if n_active == 0:
        print("[WARN] forcing vector has zero active cells — skipping diagnostic figure")
        return

    active_linear = iocn0[active_mask]
    active_vals = forcing_vec[active_mask]
    abs_vals = np.abs(active_vals)
    # Unravel to (i, j, k) in Fortran order
    ijk = np.unravel_index(active_linear, (ni, nj, nk), order="F")
    cell_i, cell_j, cell_k = ijk[0], ijk[1], ijk[2]

    print(f"[INFO] forcing diagnostic: {n_active} active cells, "
          f"k range [{cell_k.min()}, {cell_k.max()}], "
          f"lon [{geom.lon2d[cell_i, cell_j].min():.1f}, {geom.lon2d[cell_i, cell_j].max():.1f}], "
          f"lat [{geom.lat2d[cell_i, cell_j].min():.1f}, {geom.lat2d[cell_i, cell_j].max():.1f}]")

    # --- Column-integrated 2D forcing map (sum |Jdredge| across depth layers) ---
    forcing_2d = np.zeros((ni, nj))
    for idx in range(len(cell_i)):
        forcing_2d[cell_i[idx], cell_j[idx]] += abs_vals[idx]
    forcing_2d[forcing_2d == 0] = np.nan

    unique_cols = len(set(zip(cell_i.tolist(), cell_j.tolist())))

    # --- Background masks ---
    land_2d = np.where(~geom.surface_mask, 1.0, np.nan)   # 1 on land, NaN on ocean
    ocean_2d = np.where(geom.surface_mask, 1.0, np.nan)    # 1 on ocean, NaN on land

    total_tgc_yr = float(np.sum(forcing_vec * geom.ocean_cell_volume_vec * RHO_SW * MOLAR_MASS_C * 1e-6) / 1e12)

    pc = ccrs.PlateCarree()
    rob = ccrs.Robinson()
    lon_plot = geom.lon2d[0, :]
    lat_plot = geom.lat2d[:, 0]
    lon_mesh, lat_mesh = np.meshgrid(lon_plot, lat_plot)

    vmax = float(np.nanmax(forcing_2d[np.isfinite(forcing_2d)])) if np.any(np.isfinite(forcing_2d)) else 1e-12
    vmax = max(vmax, 1e-12)

    # Gray and white colormaps for land/ocean background
    gray_cmap = ListedColormap(["#d9d9d9"])
    white_cmap = ListedColormap(["#ffffff"])

    # --- Layout: 2x2 using gridspec ---
    # Panel order: (a) global, (b) singapore, (c) manila, (d) abu_dhabi
    region_order = ["singapore", "manila", "abu_dhabi"]
    panel_labels = ["(a)", "(b)", "(c)", "(d)"]

    fig = plt.figure(figsize=(18, 12))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.28, wspace=0.15)

    # (a) Global panel — Robinson
    ax_global = fig.add_subplot(gs[0, 0], projection=rob)
    ax_global.set_global()

    # Layer 1: land as gray
    ax_global.pcolormesh(
        lon_mesh, lat_mesh, land_2d, transform=pc,
        cmap=gray_cmap, vmin=0.5, vmax=1.5, shading="auto",
    )
    # Layer 2: ocean as white
    ax_global.pcolormesh(
        lon_mesh, lat_mesh, ocean_2d, transform=pc,
        cmap=white_cmap, vmin=0.5, vmax=1.5, shading="auto",
    )
    # Layer 3: forcing cells colored
    m_global = ax_global.pcolormesh(
        lon_mesh, lat_mesh, forcing_2d, transform=pc,
        cmap="YlOrRd", vmin=0, vmax=vmax, shading="auto", zorder=3,
    )
    ax_global.coastlines(linewidth=0.5, color="0.3", zorder=4)

    # Red dashed rectangles for hotspot regions
    for key, region in regions.items():
        ext = region["extent"]
        rect_lons = [ext[0], ext[1], ext[1], ext[0], ext[0]]
        rect_lats = [ext[2], ext[2], ext[3], ext[3], ext[2]]
        ax_global.plot(rect_lons, rect_lats, transform=pc,
                       color="red", linewidth=1.5, linestyle="--", zorder=6)
    # Red star at each dredging centroid
    for key, region in regions.items():
        if "dredge_centroid" in region:
            cx, cy = region["dredge_centroid"]
            ax_global.plot(cx, cy, marker="*", color="red", markersize=10,
                           markeredgecolor="black", markeredgewidth=0.5,
                           transform=pc, zorder=7)

    cb_global = fig.colorbar(m_global, ax=ax_global, orientation="horizontal",
                             shrink=0.7, pad=0.05)
    cb_global.set_label("Column-integrated |Jdredge| (\u03bcmol C kg\u207b\u00b9 yr\u207b\u00b9)", fontsize=9)
    ax_global.set_title(
        f"{panel_labels[0]} Jdredge forcing on OCIM 2\u00b0\u00d72\u00b0 grid "
        f"({unique_cols} active columns, {total_tgc_yr:.4f} TgC/yr)",
        fontsize=10,
    )

    # --- Zoom panels ---
    zoom_positions = [gs[0, 1], gs[1, 0], gs[1, 1]]
    for panel_idx, region_key in enumerate(region_order):
        region = regions[region_key]
        ext = region["extent"]
        ax = fig.add_subplot(zoom_positions[panel_idx], projection=pc)
        ax.set_extent(ext, crs=pc)

        # Layer 1: land gray
        ax.pcolormesh(
            lon_mesh, lat_mesh, land_2d, transform=pc,
            cmap=gray_cmap, vmin=0.5, vmax=1.5, shading="auto",
        )
        # Layer 2: ocean white
        ax.pcolormesh(
            lon_mesh, lat_mesh, ocean_2d, transform=pc,
            cmap=white_cmap, vmin=0.5, vmax=1.5, shading="auto",
        )
        # Layer 3: forcing cells colored
        ax.pcolormesh(
            lon_mesh, lat_mesh, forcing_2d, transform=pc,
            cmap="YlOrRd", vmin=0, vmax=vmax, shading="auto", zorder=3,
        )
        ax.coastlines(linewidth=0.8, color="0.3", zorder=4)

        # Thin black grid lines at 2° cell boundaries
        for lon_val in lon_plot:
            if ext[0] <= lon_val <= ext[1]:
                ax.plot([lon_val, lon_val], [ext[2], ext[3]], color="black",
                        linewidth=0.3, alpha=0.5, transform=pc, zorder=2)
        for lat_val in lat_plot:
            if ext[2] <= lat_val <= ext[3]:
                ax.plot([ext[0], ext[1]], [lat_val, lat_val], color="black",
                        linewidth=0.3, alpha=0.5, transform=pc, zorder=2)

        in_region = np.zeros((ni, nj), dtype=bool)
        for ci, cj in zip(cell_i, cell_j):
            clon = geom.lon2d[ci, cj]
            clat = geom.lat2d[ci, cj]
            if ext[0] <= clon <= ext[1] and ext[2] <= clat <= ext[3]:
                in_region[ci, cj] = True
        n_cols_in = len(set(
            (ci, cj) for ci, cj in zip(cell_i.tolist(), cell_j.tolist())
            if ext[0] <= geom.lon2d[ci, cj] <= ext[1] and ext[2] <= geom.lat2d[ci, cj] <= ext[3]
        ))

        # Red star at dredging centroid + land/ocean annotation
        if "dredge_centroid" in region:
            cx, cy = region["dredge_centroid"]
            ax.plot(cx, cy, marker="*", color="red", markersize=14,
                    markeredgecolor="black", markeredgewidth=0.8,
                    transform=pc, zorder=8)

            # Check if centroid falls on land or ocean
            ci_near = int(np.argmin(np.abs(lat_plot - cy)))
            cj_near = int(np.argmin(np.abs(lon_plot - cx)))
            on_land = not geom.surface_mask[ci_near, cj_near]

            if on_land:
                annot_text = "Dredging centroid on\nOCIM land cell"
            else:
                annot_text = "Dredging centroid on\nOCIM ocean cell"

            # Place annotation with arrow pointing to centroid
            ax.annotate(
                annot_text,
                xy=(cx, cy), xycoords=pc._as_mpl_transform(ax),
                xytext=(25, 25), textcoords="offset points",
                fontsize=7.5, color="red", fontweight="bold",
                arrowprops=dict(arrowstyle="->", color="red", lw=1.2),
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="red", alpha=0.85),
                zorder=9,
            )

        gl = ax.gridlines(draw_labels=True, linewidth=0.0, alpha=0.0)
        gl.top_labels = False
        gl.right_labels = False

        if n_cols_in == 0:
            subtitle = f"{region['label']} \u2014 0 cells (all coastal = land at 2\u00b0)"
        else:
            subtitle = f"{region['label']} \u2014 {n_cols_in} forcing cell{'s' if n_cols_in != 1 else ''}"
        ax.set_title(f"{panel_labels[panel_idx + 1]} {subtitle}", fontsize=9)

    _save_figure(fig, os.path.join(out_dir, "Fig_forcing_diagnostic_atwood_style"), dpi=dpi)
    print(f"[INFO] forcing diagnostic figure saved ({n_active} cells in {unique_cols} columns)")


def _load_forcing_vector(forcing_path: str) -> Tuple[np.ndarray, str]:
    with MatReader(forcing_path) as reader:
        key = _find_first(reader, ["Jdredge", "Jtrawl"])
        if key is None:
            raise KeyError("forcing file missing Jdredge/Jtrawl")
        v = _as_1d_float(reader.read(key))
        return v, key


def write_qc_report(
    qc_path: str,
    files: Dict[str, str],
    forcing_total_tgc_yr: float,
    forcing_active: int,
    year: np.ndarray,
    ts: pd.DataFrame,
    map_stats: Dict[str, Dict[str, float]],
    ph_method: str,
    ocean_cells: int,
    map_ocean_active: Dict[str, float],
) -> None:
    d_p1 = ts["delta_pco2_ppm_1x"].to_numpy()
    d_p10 = ts["delta_pco2_ppm_10x"].to_numpy()
    d_p100 = ts["delta_pco2_ppm_100x"].to_numpy()
    u1 = ts["delta_cum_uptake_PgCO2_1x"].to_numpy()
    u10 = ts["delta_cum_uptake_PgCO2_10x"].to_numpy()
    u100 = ts["delta_cum_uptake_PgCO2_100x"].to_numpy()
    ph1 = ts["delta_pH_global_1x"].to_numpy()
    ph10 = ts["delta_pH_global_10x"].to_numpy()
    ph100 = ts["delta_pH_global_100x"].to_numpy()

    eps = 1e-30
    ratio10_p = d_p10[-1] / (d_p1[-1] + eps)
    ratio100_p = d_p100[-1] / (d_p1[-1] + eps)
    ratio10_u = u10[-1] / (u1[-1] + eps)
    ratio100_u = u100[-1] / (u1[-1] + eps)
    ratio10_ph = ph10[-1] / (ph1[-1] + eps)
    ratio100_ph = ph100[-1] / (ph1[-1] + eps)

    resid_p10 = float(np.nanmax(np.abs(d_p10 - 10.0 * d_p1)))
    resid_p100 = float(np.nanmax(np.abs(d_p100 - 100.0 * d_p1)))
    resid_u10 = float(np.nanmax(np.abs(u10 - 10.0 * u1)))
    resid_u100 = float(np.nanmax(np.abs(u100 - 100.0 * u1)))
    resid_ph10 = float(np.nanmax(np.abs(ph10 - 10.0 * ph1)))
    resid_ph100 = float(np.nanmax(np.abs(ph100 - 100.0 * ph1)))

    lines = []
    lines.append("# QC Post-processing (Atwood style)\n")
    lines.append("## Inputs")
    for k, v in files.items():
        lines.append(f"- `{k}`: `{v}`")

    lines.append("\n## Time axis")
    lines.append(f"- `n_time`: {len(year)}")
    lines.append(f"- `year_min`: {int(np.nanmin(year))}")
    lines.append(f"- `year_max`: {int(np.nanmax(year))}")
    lines.append("- Requested windows: 1996-2020 and 1996-2050")

    lines.append("\n## Forcing contract check")
    lines.append(f"- `forcing_total_TgC_yr`: {forcing_total_tgc_yr:.6f}")
    if not np.isnan(forcing_total_tgc_yr):
        ok = 0.001 < forcing_total_tgc_yr < 100
        lines.append(f"- Sanity range (0.001–100 TgC/yr): {'PASS' if ok else 'FAIL'}")
    lines.append(f"- `forcing_active_cells`: {forcing_active}")

    lines.append("\n## Linearity checks")
    lines.append("- Atmospheric CO2 deltas")
    lines.append(f"  - final ratio 10x/1x: {ratio10_p:.6f}")
    lines.append(f"  - final ratio 100x/1x: {ratio100_p:.6f}")
    lines.append(f"  - max residual vs exact scaling (10x): {resid_p10:.6e} ppm")
    lines.append(f"  - max residual vs exact scaling (100x): {resid_p100:.6e} ppm")
    lines.append("- Cumulative uptake deltas")
    lines.append(f"  - final ratio 10x/1x: {ratio10_u:.6f}")
    lines.append(f"  - final ratio 100x/1x: {ratio100_u:.6f}")
    lines.append(f"  - max residual vs exact scaling (10x): {resid_u10:.6e} PgCO2")
    lines.append(f"  - max residual vs exact scaling (100x): {resid_u100:.6e} PgCO2")
    lines.append("- Global pH deltas")
    lines.append(f"  - final ratio 10x/1x: {ratio10_ph:.6f}")
    lines.append(f"  - final ratio 100x/1x: {ratio100_ph:.6f}")
    lines.append(f"  - max residual vs exact scaling (10x): {resid_ph10:.6e}")
    lines.append(f"  - max residual vs exact scaling (100x): {resid_ph100:.6e}")

    lines.append("\n## pH method")
    lines.append(f"- Global pH anomaly method: `{ph_method}`")
    lines.append("- If upper-1000m pH is unavailable in OCIM outputs, surface pH fallback is used and documented.")

    lines.append("\n## Map sanity")
    lines.append(f"- Ocean surface cells: {ocean_cells}")
    for key, stats in map_stats.items():
        lines.append(f"- `{key}`:")
        lines.append(f"  - min: {stats['min']:.6e}")
        lines.append(f"  - max: {stats['max']:.6e}")
        lines.append(f"  - ocean active fraction (>0 abs): {map_ocean_active[key]:.6f}")
        if "emissions_nonneg" in stats:
            lines.append(f"  - emissions non-negative: {bool(stats['emissions_nonneg'])}")
        if "has_pos" in stats:
            lines.append(f"  - signed field has positive values: {bool(stats['has_pos'])}")
            lines.append(f"  - signed field has negative values: {bool(stats['has_neg'])}")

    with open(qc_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="OCIM post-processing and Atwood-style figure replication")
    parser.add_argument("--output-root", default=os.environ["DFP_DATA"] + "/output_V6", help="Directory containing ocim_*.mat outputs")
    parser.add_argument("--ocim-dir", default=os.path.join(os.environ.get("DFP_DATA", ""), "configuration/ocim"), help="OCIM configuration directory")
    parser.add_argument(
        "--postproc-dir",
        default=None,
        help="Output directory for derived datasets/figures (default: <output-root>/postproc_atwood_style)",
    )
    parser.add_argument(
        "--maps-scenario",
        default="dredge_1x",
        choices=["dredge_1x", "dredge_10x", "dredge_100x"],
        help="Scenario used for Fig3/Fig4 map products",
    )
    parser.add_argument("--dpi", type=int, default=300, help="PNG DPI")
    args = parser.parse_args(argv)

    output_root = _expand(args.output_root)
    ocim_dir = _expand(args.ocim_dir)
    postproc_dir = _expand(args.postproc_dir or os.path.join(output_root, "postproc_atwood_style"))
    os.makedirs(postproc_dir, exist_ok=True)

    print("=== POSTPROC ATWOOD STYLE (OCIM) ===")
    print(f"[INFO] output_root: {output_root}")
    print(f"[INFO] ocim_dir: {ocim_dir}")
    print(f"[INFO] postproc_dir: {postproc_dir}")

    files = detect_latest_inputs(output_root)
    for k, v in files.items():
        print(f"[INFO] latest {k}: {v}")

    geom = load_geometry(ocim_dir)
    print(f"[INFO] geometry dims: ni={geom.ni}, nj={geom.nj}, nk={geom.nk}, surface={geom.surface_mask.sum()}")

    baseline = load_scenario(files["baseline"])
    year = baseline.year
    nt = year.size
    dt = _dt_years(year)
    print(f"[INFO] baseline time: n={nt}, years=[{int(np.nanmin(year))}, {int(np.nanmax(year))}]")

    if baseline.jco2.shape[1] != nt:
        raise ValueError("Baseline Jco2 time dimension mismatch")
    if baseline.jco2.shape[0] != geom.surface_mask.sum():
        raise ValueError(
            f"Baseline Jco2 surface dimension mismatch: {baseline.jco2.shape[0]} vs {geom.surface_mask.sum()}"
        )

    mask_2020 = _window_mask(year, 1996, 2020)
    mask_2050 = _window_mask(year, 1996, 2050)
    if not np.any(mask_2020):
        raise ValueError("Window 1996-2020 is empty for this year axis")
    if not np.any(mask_2050):
        raise ValueError("Window 1996-2050 is empty for this year axis")
    print(f"[INFO] window 1996-2020 points: {int(mask_2020.sum())}")
    print(f"[INFO] window 1996-2050 points: {int(mask_2050.sum())}")

    flux_baseline_gtc = global_flux_gtc_yr(baseline.jco2, geom)

    ts = pd.DataFrame({"year": year.astype(int)})
    pH_method = "surface area-weighted mean from pH_out"

    scenario_data = {}
    for scen in ("dredge_1x", "dredge_10x", "dredge_100x"):
        print(f"[INFO] loading scenario: {scen}")
        sdata = load_scenario(files[scen], expected_year=year)
        scenario_data[scen] = sdata

        if sdata.jco2.shape[0] != geom.surface_mask.sum():
            raise ValueError(f"{scen}: Jco2 dimension mismatch with geometry surface mask")

        dpco2 = sdata.pco2 - baseline.pco2
        flux_gtc = global_flux_gtc_yr(sdata.jco2, geom)
        dflux_gtc = flux_gtc - flux_baseline_gtc
        uptake_gtc = -dflux_gtc
        cum_uptake_pgco2 = np.cumsum(uptake_gtc * dt) * CO2_PER_C

        if sdata.ph is not None and baseline.ph is not None:
            if sdata.ph.shape[0] == geom.surface_mask.sum() and baseline.ph.shape[0] == geom.surface_mask.sum():
                ph_gm = global_surface_mean(sdata.ph, geom)
                ph_base_gm = global_surface_mean(baseline.ph, geom)
                dph_global = ph_gm - ph_base_gm
            elif sdata.ph.shape[0] == 1 and baseline.ph.shape[0] == 1:
                dph_global = (sdata.ph - baseline.ph).reshape(-1)
                pH_method = "global pH series anomaly (no spatial weighting available)"
            else:
                dph_global = np.full(nt, np.nan)
                pH_method = "incompatible pH diagnostics shape across scenarios"
        else:
            dph_global = np.full(nt, np.nan)
            pH_method = "missing pH diagnostics in one or more files"

        tag = scen.replace("dredge_", "")
        ts[f"delta_pco2_ppm_{tag}"] = dpco2
        ts[f"delta_flux_global_GtCyr_{tag}"] = dflux_gtc
        ts[f"delta_cum_uptake_PgCO2_{tag}"] = cum_uptake_pgco2
        ts[f"delta_pH_global_{tag}"] = dph_global

    ts_parquet = os.path.join(postproc_dir, "timeseries_dredging.parquet")
    ts_csv = os.path.join(postproc_dir, "timeseries_dredging.csv")
    try:
        ts.to_parquet(ts_parquet, index=False)
        print(f"[INFO] wrote {ts_parquet}")
    except Exception as exc:
        print(f"[WARN] parquet write failed ({exc}); CSV only")
    ts.to_csv(ts_csv, index=False)
    print(f"[INFO] wrote {ts_csv}")

    maps_scen = args.maps_scenario
    s_map = scenario_data[maps_scen]
    dflux_mol_m2 = jco2_to_flux_mol_m2(s_map.jco2, geom) - jco2_to_flux_mol_m2(baseline.jco2, geom)

    dt_2020 = dt[mask_2020].reshape(1, -1)
    dt_2050 = dt[mask_2050].reshape(1, -1)

    cum_df_2020_ns = np.sum(dflux_mol_m2[:, mask_2020] * dt_2020, axis=1)
    cum_em_2020_ns = np.sum(np.maximum(dflux_mol_m2[:, mask_2020], 0.0) * dt_2020, axis=1)
    cum_df_2050_ns = np.sum(dflux_mol_m2[:, mask_2050] * dt_2050, axis=1)
    cum_em_2050_ns = np.sum(np.maximum(dflux_mol_m2[:, mask_2050], 0.0) * dt_2050, axis=1)

    cum_df_2020 = _surface_vec_to_map(cum_df_2020_ns, geom)
    cum_em_2020 = _surface_vec_to_map(cum_em_2020_ns, geom)
    cum_df_2050 = _surface_vec_to_map(cum_df_2050_ns, geom)
    cum_em_2050 = _surface_vec_to_map(cum_em_2050_ns, geom)

    nc_2020 = os.path.join(postproc_dir, "maps_cumulative_1996_2020.nc")
    nc_2050 = os.path.join(postproc_dir, "maps_cumulative_1996_2050.nc")
    write_maps_netcdf(
        nc_2020,
        geom.lon2d,
        geom.lat2d,
        cum_em_2020,
        cum_df_2020,
        attrs={"window": "1996-2020", "scenario": maps_scen},
    )
    write_maps_netcdf(
        nc_2050,
        geom.lon2d,
        geom.lat2d,
        cum_em_2050,
        cum_df_2050,
        attrs={"window": "1996-2050", "scenario": maps_scen},
    )
    print(f"[INFO] wrote {nc_2020}")
    print(f"[INFO] wrote {nc_2050}")

    idx_2020 = np.where(np.isclose(year, 2020))[0]
    if idx_2020.size == 0:
        raise ValueError("Year 2020 not found in time axis for Fig4")
    t2020 = int(idx_2020[0])

    # OCIM does not output upper-1000m integrated pH; pH_out is the surface layer (k=1).
    # The map variable retains the Atwood-style name for figure consistency.
    if (
        baseline.ph is None
        or s_map.ph is None
        or baseline.ph.shape[0] != geom.surface_mask.sum()
        or s_map.ph.shape[0] != geom.surface_mask.sum()
    ):
        dph_map = np.full((geom.ni, geom.nj), np.nan)
        dph_method_map = "surface pH anomaly: pH_out not available or wrong shape — map is all NaN"
    else:
        dph_surface_ns = s_map.ph[:, t2020] - baseline.ph[:, t2020]
        dph_map = _surface_vec_to_map(dph_surface_ns, geom)
        dph_method_map = "surface pH anomaly from pH_out at year 2020 (OCIM surface layer k=1)"

    if "surface pH anomaly" in dph_method_map:
        fig4_stem = "Fig4_dpH_surface_2020_atwood_style"
        fig4_title = "Delta surface ocean pH due to dredging (2020) - OCIM (surface fallback)"
    else:
        fig4_stem = "Fig4_dpH_2020_atwood_style"
        fig4_title = "Delta upper-1000m pH due to dredging (2020) - OCIM"

    nc_dph = os.path.join(postproc_dir, "map_dpH_upper1000m_2020.nc")
    write_dph_netcdf(nc_dph, geom.lon2d, geom.lat2d, dph_map, dph_method_map, maps_scen)
    print(f"[INFO] wrote {nc_dph}")

    make_fig2(year, ts, postproc_dir, args.dpi)
    make_fig3(
        geom.lon2d,
        geom.lat2d,
        {"cum_emissions": cum_em_2020, "cum_df": cum_df_2020},
        {"cum_emissions": cum_em_2050, "cum_df": cum_df_2050},
        postproc_dir,
        args.dpi,
        scenario_label=maps_scen,
    )
    make_fig4(
        geom.lon2d,
        geom.lat2d,
        dph_map,
        postproc_dir,
        args.dpi,
        scenario_label=maps_scen,
        fig4_title=fig4_title,
        fig4_stem=fig4_stem,
    )
    make_fig_hotspot_zoom(
        geom.lon2d,
        geom.lat2d,
        {"cum_emissions": cum_em_2020, "cum_df": cum_df_2020},
        {"cum_emissions": cum_em_2050, "cum_df": cum_df_2050},
        dph_map,
        postproc_dir,
        args.dpi,
        scenario_label=maps_scen,
    )
    make_fig_named_hotspots(
        geom.lon2d,
        geom.lat2d,
        cum_em_2020,
        cum_df_2020,
        dph_map,
        postproc_dir,
        args.dpi,
        scenario_label=maps_scen,
    )

    forcing_vec, _ = _load_forcing_vector(files["forcing"])
    if forcing_vec.size == geom.iocn_1based.size:
        make_fig_forcing_diagnostic(
            forcing_vec, geom, postproc_dir, args.dpi,
        )
    else:
        print(f"[WARN] forcing vector size mismatch — skipping forcing diagnostic figure")

    print("[INFO] figures saved (PDF/SVG/PNG)")

    if forcing_vec.size != geom.iocn_1based.size:
        print(
            f"[WARN] forcing length ({forcing_vec.size}) != ocean points ({geom.iocn_1based.size}); "
            "forcing check may be invalid."
        )
        forcing_total_tgc_yr = float("nan")
    else:
        v_ocean = geom.ocean_cell_volume_vec
        forcing_total_g_yr = np.sum(forcing_vec * v_ocean * RHO_SW * MOLAR_MASS_C * 1e-6)
        forcing_total_tgc_yr = forcing_total_g_yr / 1e12

    def _map_stats(a: np.ndarray, is_em: bool = False) -> Dict[str, float]:
        s = {
            "min": float(np.nanmin(a)),
            "max": float(np.nanmax(a)),
        }
        if is_em:
            s["emissions_nonneg"] = float(np.nanmin(a) >= -1e-15)
        else:
            s["has_pos"] = float(np.nanmax(a) > 0)
            s["has_neg"] = float(np.nanmin(a) < 0)
        return s

    ocean_surface = geom.surface_mask
    ocean_cells = int(np.sum(ocean_surface))
    map_stats = {
        "cum_emissions_1996_2020": _map_stats(cum_em_2020, is_em=True),
        "cum_df_1996_2020": _map_stats(cum_df_2020, is_em=False),
        "cum_emissions_1996_2050": _map_stats(cum_em_2050, is_em=True),
        "cum_df_1996_2050": _map_stats(cum_df_2050, is_em=False),
    }
    map_ocean_active = {
        "cum_emissions_1996_2020": float(np.sum((np.abs(cum_em_2020) > 0) & ocean_surface) / max(ocean_cells, 1)),
        "cum_df_1996_2020": float(np.sum((np.abs(cum_df_2020) > 0) & ocean_surface) / max(ocean_cells, 1)),
        "cum_emissions_1996_2050": float(np.sum((np.abs(cum_em_2050) > 0) & ocean_surface) / max(ocean_cells, 1)),
        "cum_df_1996_2050": float(np.sum((np.abs(cum_df_2050) > 0) & ocean_surface) / max(ocean_cells, 1)),
    }

    qc_path = os.path.join(postproc_dir, "qc_postproc.md")
    write_qc_report(
        qc_path=qc_path,
        files=files,
        forcing_total_tgc_yr=forcing_total_tgc_yr,
        forcing_active=int(np.count_nonzero(forcing_vec)),
        year=year,
        ts=ts,
        map_stats=map_stats,
        ph_method=pH_method,
        ocean_cells=ocean_cells,
        map_ocean_active=map_ocean_active,
    )
    print(f"[INFO] wrote {qc_path}")

    print("\n=== OUTPUT PATHS ===")
    print(f"- {ts_parquet} (or CSV fallback)")
    print(f"- {ts_csv}")
    print(f"- {nc_2020}")
    print(f"- {nc_2050}")
    print(f"- {nc_dph}")
    print(f"- {os.path.join(postproc_dir, 'Fig2_timeseries_atwood_style.pdf')}")
    print(f"- {os.path.join(postproc_dir, 'Fig3_maps_atwood_style.pdf')}")
    print(f"- {os.path.join(postproc_dir, fig4_stem + '.pdf')}")
    print(f"- {qc_path}")
    print("=== DONE ===")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
