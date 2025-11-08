"""
ISO/TR 7250-2 anthropometric scaling and drag breakdown for pilot/body.
"""
from dataclasses import dataclass
import math
from typing import Dict

from . import config as cfg

@dataclass
class ISOFields:
    sex: str  # "male" / "female"

    # Core ISO/TR 7250-2 inputs
    body_weight_kg: float
    cervicale_height_sitting_mm: float
    shoulder_breadth_biacromial_mm: float
    hip_breadth_sitting_mm: float

    head_circumference_mm: float
    face_length_mm: float

    shoulder_elbow_length_mm: float
    elbow_grip_length_mm: float

    arm_circumference_flexed_mm: float
    forearm_circumference_flexed_mm: float

    thigh_circumference_mm: float

    # New: Thigh clearance (sitting) used for lower-body frontal area
    thigh_clearance_sitting_mm: float

    # Options
    use_shoulder_width: bool = True
    back_angle_deg: float = cfg.BACK_ANGLE_DEG
    upper_arm_angle_deg: float = cfg.UPPER_ARM_ANGLE_DEG
    protector_heigth_mm: float = cfg.PROTECTOR_HEIGHT_MM

    # Head breadth from circumference: breadth ≈ k * (C/π)
    k_head_breadth_from_circ: float = 0.85

    # Equalizer parameters (risers)
    pilot_max_weight_kg: float = 90.0
    equalizer_diameter_mm: float = cfg.EQUALIZER_DIAMETER_MM

# DIN ISO/TR 7250-2 German data - P5, P50, P95 values (in mm)
DIN_DATA = {
    "body_height": {"male": (1650, 1750, 1855), "female": (1535, 1625, 1720)},
    "body_weight_kg": {"male": (64, 79, 100), "female": (52, 66, 87)},
    "cervicale_height_sitting": {"male": (810, 860, 910), "female": (810, 860, 910)},
    "shoulder_breadth_biacromial": {"male": (370, 405, 435), "female": (345, 370, 400)},
    "hip_breadth_sitting": {"male": (350, 375, 420), "female": (360, 390, 460)},
    "thigh_clearance_sitting": {"male": (130, 150, 180), "female": (125, 145, 175)},
    "head_circumference": {"male": (545, 570, 600), "female": (520, 545, 570)},
    "face_length": {"male": (105, 115, 130), "female": (95, 110, 125)},
    "shoulder_elbow_length": {"male": (330, 365, 400), "female": (290, 320, 350)},
    "elbow_grip_length": {"male": (325, 350, 390), "female": (295, 315, 350)},
    "arm_circumference_flexed": {"male": (280, 320, 380), "female": (230, 270, 320)},
    "forearm_circumference": {"male": (245, 265, 285), "female": (225, 245, 260)},
    "thigh_circumference": {"male": (490, 570, 640), "female": (485, 565, 670)},
}


def calculate_sd_from_percentiles(p5: float, p50: float, p95: float):
    """
    Calculate the mean and standard deviation assuming a normal distribution,
    given the 5th, 50th, and 95th percentiles.

    The standard deviation is estimated as (p95 - p5) / (2 * 1.645),
    where 1.645 is the z-score for the 5th and 95th percentiles.
    """
    mean = p50
    sd = (p95 - p5) / (2 * 1.645)
    return mean, sd


def estimate_z_score(height_cm: float, sex: str) -> float:
    """
    Estimate the z-score (number of standard deviations from the mean) for a given body height (in cm)
    and sex, based on DIN ISO/TR 7250-2 population data.

    The returned z-score represents how many standard deviations the given height is above (positive)
    or below (negative) the population mean for the specified sex.
    """
    h_data = DIN_DATA["body_height"][sex]
    h_p5, h_p50, h_p95 = [h/10 for h in h_data]

    if height_cm <= h_p5:
        h_sd = (h_p95 - h_p5) / (2 * 1.645)
        z_height = (height_cm - h_p50) / h_sd
    elif height_cm <= h_p50:
        pct_range = (height_cm - h_p5) / (h_p50 - h_p5)
        z_height = -1.645 + pct_range * 1.645
    elif height_cm <= h_p95:
        pct_range = (height_cm - h_p50) / (h_p95 - h_p50)
        z_height = pct_range * 1.645
    else:
        h_sd = (h_p95 - h_p5) / (2 * 1.645)
        z_height = (height_cm - h_p50) / h_sd

    return z_height

def calculate_dimension_from_z(dimension_name: str, sex: str, z_score: float) -> float:
    """
    Calculate the value of a specified anthropometric dimension for a given sex and z-score,
    using DIN ISO/TR 7250-2 data.

    Parameters:
        dimension_name (str): The name of the anthropometric dimension (e.g., 'body_weight_kg').
        sex (str): 'male' or 'female'.
        z_score (float): The z-score representing the number of standard deviations from the mean.

    Returns:
        float: The estimated value of the specified dimension.
    """
    data = DIN_DATA[dimension_name][sex]
    mean, sd = calculate_sd_from_percentiles(*data)
    return mean + z_score * sd


from typing import Optional

def create_person_from_height_weight(height_cm: float, weight_kg: float, sex: str,
                                     pilot_max_weight_kg: Optional[float] = None,
                                     back_angle_deg: Optional[float] = None,
                                     upper_arm_angle_deg: Optional[float] = None,
                                     protector_height_mm: Optional[float] = None,
                                     equalizer_diameter_mm: Optional[float] = None) -> ISOFields:
    z = estimate_z_score(height_cm, sex)

    body_weight = calculate_dimension_from_z("body_weight_kg", sex, z)
    cervicale = calculate_dimension_from_z("cervicale_height_sitting", sex, z)
    shoulder_breadth = calculate_dimension_from_z("shoulder_breadth_biacromial", sex, z)
    hip_breadth = calculate_dimension_from_z("hip_breadth_sitting", sex, z)
    thigh_clearance = calculate_dimension_from_z("thigh_clearance_sitting", sex, z)
    head_circ = calculate_dimension_from_z("head_circumference", sex, z)
    face_len = calculate_dimension_from_z("face_length", sex, z)
    shoulder_elbow = calculate_dimension_from_z("shoulder_elbow_length", sex, z)
    elbow_grip = calculate_dimension_from_z("elbow_grip_length", sex, z)
    arm_circ = calculate_dimension_from_z("arm_circumference_flexed", sex, z)
    forearm_circ = calculate_dimension_from_z("forearm_circumference", sex, z)
    thigh_circ = calculate_dimension_from_z("thigh_circumference", sex, z)

    return ISOFields(
        sex=sex,
        body_weight_kg=body_weight,
        cervicale_height_sitting_mm=cervicale,
        shoulder_breadth_biacromial_mm=shoulder_breadth,
        hip_breadth_sitting_mm=hip_breadth,
        thigh_clearance_sitting_mm=thigh_clearance,
        head_circumference_mm=head_circ,
        face_length_mm=face_len,
        shoulder_elbow_length_mm=shoulder_elbow,
        elbow_grip_length_mm=elbow_grip,
        arm_circumference_flexed_mm=arm_circ,
        forearm_circumference_flexed_mm=forearm_circ,
        thigh_circumference_mm=thigh_circ,
        back_angle_deg=back_angle_deg or cfg.BACK_ANGLE_DEG,
        upper_arm_angle_deg=upper_arm_angle_deg or cfg.UPPER_ARM_ANGLE_DEG,
        protector_heigth_mm=protector_height_mm or cfg.PROTECTOR_HEIGHT_MM,
        pilot_max_weight_kg=pilot_max_weight_kg or cfg.CAL_REF_WING_KG,
        equalizer_diameter_mm=equalizer_diameter_mm or cfg.EQUALIZER_DIAMETER_MM,
    )


def frontal_areas_iso(fields: ISOFields) -> Dict[str, float]:
    θ = math.radians(fields.back_angle_deg)

    torso_width_mm = fields.shoulder_breadth_biacromial_mm if fields.use_shoulder_width else fields.hip_breadth_sitting_mm

    A_upper_body = (torso_width_mm/1000) * ((fields.cervicale_height_sitting_mm/1000) * math.sin(θ))
    A_protector = (fields.protector_heigth_mm/1000) * (torso_width_mm/1000)

    head_breadth_mm = fields.k_head_breadth_from_circ * (fields.head_circumference_mm / math.pi)
    head_height_mm = 2.0 * fields.face_length_mm
    A_head = (head_breadth_mm/1000) * ((head_height_mm/1000) * math.sin(θ))

    upper_arm_thickness_mm = fields.arm_circumference_flexed_mm / math.pi
    forearm_thickness_mm   = fields.forearm_circumference_flexed_mm / math.pi

    L_shoulder_elbow_mm = fields.shoulder_elbow_length_mm
    A_upper_arm_one = (upper_arm_thickness_mm/1000) * ((L_shoulder_elbow_mm/1000) * math.cos(math.radians(fields.upper_arm_angle_deg)))

    A_forearm_one = (forearm_thickness_mm/1000) * (fields.elbow_grip_length_mm/1000)

    A_upper_arms_total = 2.0 * A_upper_arm_one
    A_forearms_total = 2.0 * A_forearm_one
    A_arms_total = A_upper_arms_total + A_forearms_total

    # Lower body frontal area
    A_lower_body = (fields.thigh_clearance_sitting_mm/1000) * (fields.hip_breadth_sitting_mm/1000)

    # Base equalizer length scales with wing size above ~90 kg baseline
    base_equalizer_length_mm = max(0.0, (fields.pilot_max_weight_kg - 90.0) * 22.0)
    equalizer_length_mm = base_equalizer_length_mm
    # Gear-weight based increase via configured mode
    try:
        gear_weight_kg = max(0.0, fields.pilot_max_weight_kg - fields.body_weight_kg)
        mode = getattr(cfg, 'EQUALIZER_GEAR_MODE', 'LEGACY_PER_KG')
        stage1_thresh = getattr(cfg, 'EQUALIZER_GEAR_STAGE1_THRESHOLD_KG', 1e9)
        stage2_thresh = getattr(cfg, 'EQUALIZER_GEAR_STAGE2_THRESHOLD_KG', 1e9)

        extra_len_mm = 0.0
        if mode == 'FIXED_STEPS':
            # Binary steps: add fixed amount based on threshold crossed
            step1 = getattr(cfg, 'FIXED_STEP_STAGE1_MM', 0.0)
            step2 = getattr(cfg, 'FIXED_STEP_STAGE2_MM', 0.0)
            if gear_weight_kg > stage2_thresh:
                extra_len_mm = step2
            elif gear_weight_kg > stage1_thresh:
                extra_len_mm = step1
        else:
            # LEGACY_PER_KG: per-kg piecewise linear increments
            stage1_mm_per_kg = getattr(cfg, 'EQUALIZER_GEAR_STAGE1_MM_PER_KG', 0.0)
            stage2_mm_per_kg = getattr(cfg, 'EQUALIZER_GEAR_STAGE2_MM_PER_KG', 0.0)
            if gear_weight_kg > stage2_thresh:
                extra_len_mm = stage2_mm_per_kg * (gear_weight_kg - stage2_thresh)
            elif gear_weight_kg > stage1_thresh:
                extra_len_mm = stage1_mm_per_kg * (gear_weight_kg - stage1_thresh)

        equalizer_length_mm = base_equalizer_length_mm + max(0.0, extra_len_mm)

        # Optional cap
        max_len = getattr(cfg, 'EQUALIZER_MAX_LENGTH_MM', None)
        if isinstance(max_len, (int, float)) and max_len is not None:
            equalizer_length_mm = min(equalizer_length_mm, float(max_len))
    except Exception:
        # Fail silently, keep base length
        equalizer_length_mm = base_equalizer_length_mm
    A_equalizer_one = (fields.equalizer_diameter_mm/1000) * (equalizer_length_mm/1000)
    A_equalizers = 2.0 * A_equalizer_one

    # A_streamlined = A_upper_body + A_protector + A_lower_body (per requirement)
    A_streamlined = A_upper_body + A_protector + A_lower_body

    return {
        "A_upper_body_m2": A_upper_body,
        "A_protector_m2": A_protector,
        "A_head_m2": A_head,
        "A_upper_arms_m2": A_upper_arms_total,
        "A_forearms_m2": A_forearms_total,
        "A_arms_total_m2": A_arms_total,
        "A_lower_body_m2": A_lower_body,
        "A_equalizers_m2": A_equalizers,
        "A_streamlined_m2": A_streamlined,
        "A_total_m2": A_upper_body + A_protector + A_head + A_arms_total + A_lower_body + A_equalizers,
        # echoes
        "torso_width_used_mm": torso_width_mm,
        "head_breadth_mm_est": head_breadth_mm,
        "head_height_mm_est": head_height_mm,
        "upper_arm_thickness_mm": upper_arm_thickness_mm,
        "forearm_thickness_mm": forearm_thickness_mm,
        "equalizer_length_mm": equalizer_length_mm,
        "thigh_clearance_sitting_mm": fields.thigh_clearance_sitting_mm,
    }


def drag_all_iso(areas: dict,
                 speed_mps: Optional[float] = None,
                 rho_air: Optional[float] = None,
                 Cd_streamlined: Optional[float] = None,
                 Cd_head: Optional[float] = None,
                 head_exposed_fraction: Optional[float] = None,
                 Cd_arms: Optional[float] = None,
                 Cd_equalizers: Optional[float] = None) -> dict:
    if speed_mps is None:
        speed_mps = cfg.SPEED_MPS
    if rho_air is None:
        rho_air = cfg.RHO_AIR
    if Cd_streamlined is None:
        Cd_streamlined = cfg.CD_STREAMLINED
    if Cd_head is None:
        Cd_head = cfg.CD_HEAD
    if head_exposed_fraction is None:
        head_exposed_fraction = cfg.HEAD_EXPOSED_FRACTION
    if Cd_arms is None:
        Cd_arms = cfg.CD_ARMS
    if Cd_equalizers is None:
        Cd_equalizers = cfg.CD_EQUALIZERS

    q = 0.5 * rho_air * speed_mps**2
    A_stream = areas["A_streamlined_m2"]
    A_head   = areas["A_head_m2"]
    A_upper_arms = areas["A_upper_arms_m2"]
    A_forearms = areas["A_forearms_m2"]
    A_arms   = areas["A_arms_total_m2"]
    A_equalizers = areas["A_equalizers_m2"]

    D_stream = q * Cd_streamlined * A_stream
    D_head   = q * (head_exposed_fraction * Cd_head) * A_head
    D_upper_arms = q * Cd_arms * A_upper_arms
    D_forearms = q * Cd_arms * A_forearms
    D_arms   = D_upper_arms + D_forearms
    D_equalizers = q * Cd_equalizers * A_equalizers
    D_total  = D_stream + D_head + D_arms + D_equalizers

    return {
        "q_Pa": q,
        "Cd_streamlined": Cd_streamlined,
        "Cd_head": Cd_head,
        "head_exposed_fraction": head_exposed_fraction,
        "Cd_arms": Cd_arms,
        "Cd_equalizers": Cd_equalizers,
        "Drag_streamlined_N": D_stream,
        "Drag_head_N": D_head,
        "Drag_upper_arms_N": D_upper_arms,
        "Drag_forearms_N": D_forearms,
        "Drag_arms_N": D_arms,
        "Drag_equalizers_N": D_equalizers,
        "Drag_total_N": D_total,
    }
