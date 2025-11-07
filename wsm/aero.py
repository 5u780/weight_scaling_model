"""
Aerodynamics: wing polar model, lines area model, calibration of reference Cd.
"""
import math
from typing import Dict

from . import config as cfg
from .iso_model import create_person_from_height_weight, frontal_areas_iso, drag_all_iso


def lines_area_from_wing_size(wing_max_kg: float) -> float:
    """Return total line frontal area (m^2) for a given wing size rating (kg).
    The formula is a linearisation from Luc Armants spreadsheet: A_lines = C0 + m*(size - 85).
    """
    
    line_area = cfg.LINES_AREA_C0 + cfg.LINES_AREA_M * (wing_max_kg - cfg.LINES_AREA_REF_KG)
    return max(0.0, line_area)


def calculate_wing_drag(
    wing_max_kg: float,
    pilot_drag_N: float,
    Cd_wing_ref: float,               # total Cd at reference (projected area)
    speed_mps: float = None,
    rho_air: float = None,
    e_oswald: float = cfg.E_OSWALD,
) -> Dict[str, float]:
    if speed_mps is None:
        speed_mps = cfg.SPEED_MPS
    if rho_air is None:
        rho_air = cfg.RHO_AIR

    q = 0.5 * rho_air * speed_mps**2

    wing_flat_area_m2 = wing_max_kg / 4.5
    wing_projected_area_m2 = cfg.PROJECTED_TO_FLAT * wing_flat_area_m2

    chord_m = (wing_flat_area_m2 / cfg.WING_ASPECT_RATIO) ** 0.5
    reynolds = speed_mps * chord_m / cfg.KINEMATIC_VISCOSITY

    ref_flat_area = cfg.REFERENCE_WING_MAX_KG / 4.5
    ref_chord = (ref_flat_area / cfg.WING_ASPECT_RATIO) ** 0.5
    reynolds_ref = cfg.REFERENCE_SPEED_MPS * ref_chord / cfg.KINEMATIC_VISCOSITY

    re_factor = (reynolds_ref / reynolds) ** cfg.REYNOLDS_EXPONENT

    ref_projected_area = cfg.PROJECTED_TO_FLAT * ref_flat_area
    q_ref = 0.5 * cfg.RHO_AIR * cfg.REFERENCE_SPEED_MPS**2
    W_ref_N = cfg.REFERENCE_WING_MAX_KG * cfg.GRAVITY

    CL_ref_proj = W_ref_N / (q_ref * ref_projected_area)
    AR_projected = cfg.WING_ASPECT_RATIO / cfg.PROJECTED_TO_FLAT
    k_induced_proj = 1.0 / (math.pi * e_oswald * AR_projected)
    CDi_ref_proj = k_induced_proj * (CL_ref_proj ** 2)
    Cd0_ref_proj = max(Cd_wing_ref - CDi_ref_proj, 0.0)

    Cd0_flat_ref = Cd0_ref_proj * (ref_projected_area / ref_flat_area)
    Cd0_corrected_flat = Cd0_flat_ref * re_factor

    total_weight_N = wing_max_kg * cfg.GRAVITY
    CL = total_weight_N / (q * wing_flat_area_m2)

    k_induced = 1.0 / (math.pi * e_oswald * cfg.WING_ASPECT_RATIO)

    CDi = k_induced * CL**2
    CD_total = Cd0_corrected_flat + CDi

    wing_drag_N = q * wing_flat_area_m2 * CD_total

    total_drag_N = pilot_drag_N + wing_drag_N
    glide_ratio = total_weight_N / total_drag_N if total_drag_N > 0 else 0.0

    return {
        "wing_area_m2_proj": wing_projected_area_m2,
        "wing_area_m2_flat": wing_flat_area_m2,
        "CL": CL,
        "Cd0_flat": Cd0_corrected_flat,
        "CDi": CDi,
        "CD_total": CD_total,
        "k_induced": k_induced,
        "Reynolds": reynolds,
        "Reynolds_ref": reynolds_ref,
        "re_factor": re_factor,
        "drag_wing_N": wing_drag_N,
        "drag_total_N": total_drag_N,
        "glide_ratio": glide_ratio,
    }


def calibrate_cd_wing_ref(CD_EQUALIZERS_FOR_CAL: float = 0.0) -> float:
    """Compute Cd_wing_ref on projected area so that the reference condition
    (CAL_REF_WING_KG at CAL_REF_SPEED_MPS) yields CAL_REF_SINK_MPS for CAL_REF_PILOT.

    Note: Calibration is performed with an explicit equalizer Cd passed via
    CD_EQUALIZERS_FOR_CAL (default 0.0) to avoid baking current cfg.CD_EQUALIZERS
    into the baseline wing Cd.
    """
    height_cm, weight_kg, sex = cfg.CAL_REF_PILOT

    person = create_person_from_height_weight(
        height_cm, weight_kg, sex,
        pilot_max_weight_kg=cfg.CAL_REF_WING_KG,
        back_angle_deg=cfg.BACK_ANGLE_DEG,
        upper_arm_angle_deg=cfg.UPPER_ARM_ANGLE_DEG,
        protector_height_mm=cfg.PROTECTOR_HEIGHT_MM,
        equalizer_diameter_mm=cfg.EQUALIZER_DIAMETER_MM,
    )
    areas = frontal_areas_iso(person)

    v = cfg.CAL_REF_SPEED_MPS
    q = 0.5 * cfg.RHO_AIR * v**2

    D_body = drag_all_iso(
        areas,
        speed_mps=v,
        rho_air=cfg.RHO_AIR,
        Cd_streamlined=cfg.CD_STREAMLINED,
        Cd_head=cfg.CD_HEAD,
        head_exposed_fraction=cfg.HEAD_EXPOSED_FRACTION,
        Cd_arms=cfg.CD_ARMS,
        Cd_equalizers=CD_EQUALIZERS_FOR_CAL,
    )
    A_lines = lines_area_from_wing_size(cfg.CAL_REF_WING_KG)
    D_lines = q * cfg.CD_LINES * A_lines
    D_pilot_total = D_body["Drag_total_N"] + D_lines

    S_flat = cfg.CAL_REF_WING_KG / 4.5
    L_total = cfg.CAL_REF_WING_KG * cfg.GRAVITY

    glide_target = v / cfg.CAL_REF_SINK_MPS
    D_total_target = L_total / glide_target
    D_wing_target = max(D_total_target - D_pilot_total, 1e-9)

    CD_total_req = D_wing_target / (q * S_flat)

    k_flat = 1.0 / (math.pi * cfg.E_OSWALD * cfg.WING_ASPECT_RATIO)
    CL_flat = L_total / (q * S_flat)
    CDi_flat = k_flat * CL_flat**2

    Cd0_flat_ref = max(CD_total_req - CDi_flat, 0.0)

    S_proj_ref = cfg.PROJECTED_TO_FLAT * S_flat
    Cd0_ref_proj = Cd0_flat_ref / cfg.PROJECTED_TO_FLAT

    AR_projected = cfg.WING_ASPECT_RATIO / cfg.PROJECTED_TO_FLAT
    k_proj = 1.0 / (math.pi * cfg.E_OSWALD * AR_projected)
    CL_ref_proj = L_total / (q * S_proj_ref)
    CDi_ref_proj = k_proj * CL_ref_proj**2

    Cd_total_ref_proj = Cd0_ref_proj + CDi_ref_proj
    return Cd_total_ref_proj
