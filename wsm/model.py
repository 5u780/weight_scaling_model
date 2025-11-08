"""Core modeling helpers for polars and single-speed performance.

This module centralizes computations that were previously duplicated across
`run_model.py`, `plots.py`, and calibration scripts. Keeping them here reduces
drift if aerodynamic formulas or coefficient handling change.

Functions
---------
compute_total_drag(person, wing_max, speed_mps, cd_wing_ref,
                   include_equalizer_drag=True) -> dict
    Compute pilot aerodynamic drag components (with or without equalizer drag),
    add suspension line drag, compute wing drag and glide ratio.

compute_polar(person, wing_max, speeds_mps, cd_wing_ref,
              include_equalizer_drag=True) -> pandas.DataFrame
    Vectorized helper building a polar (speed, sink, glide, drag breakdown).

summarize_polar(df) -> dict
    Extract min sink and best glide points.

get_cd_wing_ref() -> float
    Memoized retrieval of calibrated wing Cd (projected), excluding equalizers.

Notes
-----
Wing area modeling uses flat area = wing_max_kg / 4.5 and projected area via
`cfg.PROJECTED_TO_FLAT`. Pilot drag uses the existing iso_model + drag_all_iso
functions for consistency.
"""

from __future__ import annotations

from typing import Iterable, Dict
import pandas as pd

from . import config as cfg
from .iso_model import frontal_areas_iso, drag_all_iso
from .aero import lines_area_from_wing_size, calculate_wing_drag, calibrate_cd_wing_ref

_CD_WING_REF_CACHE: float | None = None

def get_cd_wing_ref() -> float:
    """Return calibrated wing Cd (projected area) with memoization.

    Calibration excludes equalizer drag to avoid biasing the baseline.
    """
    global _CD_WING_REF_CACHE
    if _CD_WING_REF_CACHE is None:
        _CD_WING_REF_CACHE = calibrate_cd_wing_ref(CD_EQUALIZERS_FOR_CAL=0.0)
    return _CD_WING_REF_CACHE

def compute_total_drag(person, wing_max_kg: float, speed_mps: float, cd_wing_ref: float,
                       include_equalizer_drag: bool = True) -> Dict[str, float]:
    """Compute single-speed aerodynamic breakdown.

    Returns a dictionary containing pilot component drags, line drag, wing drag,
    total drag, glide ratio, sink rate and key areas.
    """
    areas = frontal_areas_iso(person)
    Cd_eq = cfg.CD_EQUALIZERS if include_equalizer_drag else 0.0
    D_pilot = drag_all_iso(
        areas,
        speed_mps=speed_mps,
        rho_air=cfg.RHO_AIR,
        Cd_streamlined=cfg.CD_STREAMLINED,
        Cd_head=cfg.CD_HEAD,
        head_exposed_fraction=cfg.HEAD_EXPOSED_FRACTION,
        Cd_arms=cfg.CD_ARMS,
        Cd_equalizers=Cd_eq,
    )
    q = 0.5 * cfg.RHO_AIR * speed_mps**2
    A_lines = lines_area_from_wing_size(wing_max_kg)
    Drag_lines_N = q * cfg.CD_LINES * A_lines
    pilot_drag_N = D_pilot["Drag_total_N"] + Drag_lines_N
    W = calculate_wing_drag(wing_max_kg, pilot_drag_N, cd_wing_ref, speed_mps, cfg.RHO_AIR, e_oswald=cfg.E_OSWALD)
    glide = W["glide_ratio"]
    sink = speed_mps / glide if glide > 0 else float("nan")
    lift_force_N = wing_max_kg * cfg.GRAVITY
    return {
        "speed_mps": speed_mps,
        "speed_kmh": speed_mps * 3.6,
        "pilot_drag_N": pilot_drag_N,
        "wing_drag_N": W["drag_wing_N"],
        "drag_total_N": W["drag_total_N"],
        "glide_ratio": glide,
        "sink_mps": sink,
        "CL": W["CL"],
        "Cd0_flat": W["Cd0_flat"],
        "CDi": W["CDi"],
        "CD_total": W["CD_total"],
        # Wing geometry/aero
        "wing_area_m2": W["wing_area_m2_proj"],
        "wing_flat_area_m2": W["wing_area_m2_flat"],
        "k_induced": W["k_induced"],
        "Reynolds": W["Reynolds"],
        "Reynolds_ref": W["Reynolds_ref"],
        "re_factor": W["re_factor"],
        # Lines
        "A_lines_m2": A_lines,
        "Drag_lines_N": Drag_lines_N,
        # Areas (pilot)
        "A_total_m2": areas.get("A_total_m2", 0.0),
        "A_streamlined_m2": areas.get("A_streamlined_m2", 0.0),
        "A_head_m2": areas.get("A_head_m2", 0.0),
        "A_upper_arms_m2": areas.get("A_upper_arms_m2", 0.0),
        "A_forearms_m2": areas.get("A_forearms_m2", 0.0),
        "A_arms_m2": areas.get("A_arms_total_m2", 0.0),
        "A_equalizers_m2": areas.get("A_equalizers_m2", 0.0),
        # Pilot drag components
        "Drag_streamlined_N": D_pilot.get("Drag_streamlined_N", 0.0),
        "Drag_head_N": D_pilot.get("Drag_head_N", 0.0),
        "Drag_upper_arms_N": D_pilot.get("Drag_upper_arms_N", 0.0),
        "Drag_forearms_N": D_pilot.get("Drag_forearms_N", 0.0),
        "Drag_arms_N": D_pilot.get("Drag_arms_N", 0.0),
        "Drag_equalizers_N": D_pilot.get("Drag_equalizers_N", 0.0),
        # Misc
        "equalizer_length_mm": areas.get("equalizer_length_mm", 0.0),
        "lift_force_N": lift_force_N,
    }

def compute_polar(person, wing_max_kg: float, speeds_mps: Iterable[float], cd_wing_ref: float,
                  include_equalizer_drag: bool = True) -> pd.DataFrame:
    rows = [
        compute_total_drag(person, wing_max_kg, v, cd_wing_ref, include_equalizer_drag=include_equalizer_drag)
        for v in speeds_mps
    ]
    return pd.DataFrame(rows)

def summarize_polar(df: pd.DataFrame) -> Dict[str, float]:
    idx_min_sink = df["sink_mps"].idxmin()
    idx_best_ld = df["glide_ratio"].idxmax()
    # Convert via item() for type-checker friendliness
    v_min_sink = float(pd.Series(df.loc[idx_min_sink, "speed_kmh"]).item())
    sink_min = float(pd.Series(df.loc[idx_min_sink, "sink_mps"]).item())
    v_best_ld = float(pd.Series(df.loc[idx_best_ld, "speed_kmh"]).item())
    best_ld = float(pd.Series(df.loc[idx_best_ld, "glide_ratio"]).item())
    return {
        "v_min_sink_kmh": v_min_sink,
        "sink_min_mps": sink_min,
        "v_best_ld_kmh": v_best_ld,
        "best_ld": best_ld,
    }
