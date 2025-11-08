"""Standalone ISO/TR 7250-2 anthropometry + frontal area estimation extracted from `wsm.iso_model`.

This module focuses ONLY on generating body dimensions and frontal areas for a pilot given
height, weight, and sex. Aerodynamic drag coefficients or wing logic are intentionally omitted
so external models can plug in their own physics.

Main entry points:
    create_person_from_height_weight(...): percentile-aware dimension synthesis (DIN ISO/TR 7250-2)
    frontal_areas_iso(fields): frontal area components (m^2) + equalizer length
    compute_equalizer_length_mm(fields): equalizer strap length with mode-based extras

Differences vs original:
    - Drag calculation function removed (keep your own Cd assignments externally).
    - Configuration references are redirected to local `anthro.config`.

DISCLAIMER: Population scaling assumes normal distributions derived from P5/P50/P95 DIN values.
Replace with full percentile tables for production-grade ergonomic simulations.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Optional
import math

import config as cfg

# DIN ISO/TR 7250-2 German data - P5, P50, P95 values (in mm; weight in kg)
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


def _percentile_mean_sd(p5: float, p50: float, p95: float):
    mean = p50
    sd = (p95 - p5) / (2 * 1.645)
    return mean, sd


def estimate_height_z(height_cm: float, sex: str) -> float:
    data = DIN_DATA["body_height"][sex]
    p5, p50, p95 = [h / 10 for h in data]
    if height_cm <= p5:
        sd = (p95 - p5) / (2 * 1.645)
        return (height_cm - p50) / sd
    if height_cm <= p50:
        frac = (height_cm - p5) / (p50 - p5)
        return -1.645 + frac * 1.645
    if height_cm <= p95:
        frac = (height_cm - p50) / (p95 - p50)
        return frac * 1.645
    sd = (p95 - p5) / (2 * 1.645)
    return (height_cm - p50) / sd


def dimension_from_z(name: str, sex: str, z: float) -> float:
    p5, p50, p95 = DIN_DATA[name][sex]
    mean, sd = _percentile_mean_sd(p5, p50, p95)
    return mean + z * sd


@dataclass
class ISOFields:
    sex: str
    body_weight_kg: float
    cervicale_height_sitting_mm: float
    shoulder_breadth_biacromial_mm: float
    hip_breadth_sitting_mm: float
    thigh_clearance_sitting_mm: float
    head_circumference_mm: float
    face_length_mm: float
    shoulder_elbow_length_mm: float
    elbow_grip_length_mm: float
    arm_circumference_flexed_mm: float
    forearm_circumference_flexed_mm: float
    thigh_circumference_mm: float
    back_angle_deg: float
    upper_arm_angle_deg: float
    protector_height_mm: float
    pilot_max_weight_kg: float | None
    equalizer_diameter_mm: float
    k_head_breadth_from_circ: float = 0.85
    use_shoulder_width: bool = True


def create_person_from_height_weight(
    height_cm: float,
    weight_kg: float,
    sex: str,
    pilot_max_weight_kg: Optional[float] = None,
    back_angle_deg: Optional[float] = None,
    upper_arm_angle_deg: Optional[float] = None,
    protector_height_mm: Optional[float] = None,
    equalizer_diameter_mm: Optional[float] = None,
) -> ISOFields:
    sex = sex.lower()
    if sex not in ("male", "female"):
        raise ValueError("sex must be 'male' or 'female'")
    z = estimate_height_z(height_cm, sex)
    fields = ISOFields(
        sex=sex,
        body_weight_kg=dimension_from_z("body_weight_kg", sex, z),
        cervicale_height_sitting_mm=dimension_from_z("cervicale_height_sitting", sex, z),
        shoulder_breadth_biacromial_mm=dimension_from_z("shoulder_breadth_biacromial", sex, z),
        hip_breadth_sitting_mm=dimension_from_z("hip_breadth_sitting", sex, z),
        thigh_clearance_sitting_mm=dimension_from_z("thigh_clearance_sitting", sex, z),
        head_circumference_mm=dimension_from_z("head_circumference", sex, z),
        face_length_mm=dimension_from_z("face_length", sex, z),
        shoulder_elbow_length_mm=dimension_from_z("shoulder_elbow_length", sex, z),
        elbow_grip_length_mm=dimension_from_z("elbow_grip_length", sex, z),
        arm_circumference_flexed_mm=dimension_from_z("arm_circumference_flexed", sex, z),
        forearm_circumference_flexed_mm=dimension_from_z("forearm_circumference", sex, z),
        thigh_circumference_mm=dimension_from_z("thigh_circumference", sex, z),
        back_angle_deg=back_angle_deg if back_angle_deg is not None else cfg.BACK_ANGLE_DEG,
        upper_arm_angle_deg=upper_arm_angle_deg if upper_arm_angle_deg is not None else cfg.UPPER_ARM_ANGLE_DEG,
        protector_height_mm=protector_height_mm if protector_height_mm is not None else cfg.PROTECTOR_HEIGHT_MM,
        pilot_max_weight_kg=pilot_max_weight_kg,
        equalizer_diameter_mm=equalizer_diameter_mm if equalizer_diameter_mm is not None else cfg.EQUALIZER_DIAMETER_MM,
    )
    return fields


def compute_equalizer_length_mm(fields: ISOFields) -> float:
    """Compute equalizer strap length using original variable names (mirrors wsm.iso_model logic).
    Base: (pilot_max_weight_kg - 90) * 22  clipped at 0.
    Mode LEGACY_PER_KG: piecewise per-kg slopes above stage thresholds.
    Mode FIXED_STEPS: binary fixed extra lengths when thresholds crossed.
    """
    if fields.pilot_max_weight_kg is None:
        return 0.0
    base_equalizer_length_mm = max(0.0, (fields.pilot_max_weight_kg - 90.0) * 22.0)
    gear_weight_kg = max(0.0, fields.pilot_max_weight_kg - fields.body_weight_kg)
    mode = getattr(cfg, "EQUALIZER_GEAR_MODE", "LEGACY_PER_KG")
    stage1 = getattr(cfg, "EQUALIZER_GEAR_STAGE1_THRESHOLD_KG", 25.0)
    stage2 = getattr(cfg, "EQUALIZER_GEAR_STAGE2_THRESHOLD_KG", 35.0)
    extra_len_mm = 0.0
    if mode == "FIXED_STEPS":
        if gear_weight_kg > stage2:
            extra_len_mm = getattr(cfg, "FIXED_STEP_STAGE2_MM", 0.0)
        elif gear_weight_kg > stage1:
            extra_len_mm = getattr(cfg, "FIXED_STEP_STAGE1_MM", 0.0)
    else:  # LEGACY_PER_KG
        s1 = getattr(cfg, "EQUALIZER_GEAR_STAGE1_MM_PER_KG", 0.0)
        s2 = getattr(cfg, "EQUALIZER_GEAR_STAGE2_MM_PER_KG", 0.0)
        if gear_weight_kg > stage2:
            extra_len_mm = s2 * (gear_weight_kg - stage2)
        elif gear_weight_kg > stage1:
            extra_len_mm = s1 * (gear_weight_kg - stage1)
    length = base_equalizer_length_mm + max(0.0, extra_len_mm)
    cap = getattr(cfg, "EQUALIZER_MAX_LENGTH_MM", None)
    if isinstance(cap, (int, float)) and cap is not None:
        length = min(length, float(cap))
    return length


def frontal_areas_iso(fields: ISOFields) -> Dict[str, float]:
    theta = math.radians(fields.back_angle_deg)
    torso_width_mm = fields.shoulder_breadth_biacromial_mm if fields.use_shoulder_width else fields.hip_breadth_sitting_mm

    A_upper_body = (torso_width_mm/1000) * ((fields.cervicale_height_sitting_mm/1000) * math.sin(theta))
    A_protector = (fields.protector_height_mm/1000) * (torso_width_mm/1000)

    head_breadth_mm = fields.k_head_breadth_from_circ * (fields.head_circumference_mm / math.pi)
    head_height_mm = 2.0 * fields.face_length_mm
    A_head = (head_breadth_mm/1000) * ((head_height_mm/1000) * math.sin(theta))

    upper_arm_thickness_mm = fields.arm_circumference_flexed_mm / math.pi
    forearm_thickness_mm = fields.forearm_circumference_flexed_mm / math.pi
    A_upper_arm_one = (upper_arm_thickness_mm/1000) * ((fields.shoulder_elbow_length_mm/1000) * math.cos(math.radians(fields.upper_arm_angle_deg)))
    A_forearm_one = (forearm_thickness_mm/1000) * (fields.elbow_grip_length_mm/1000)
    A_upper_arms_total = 2.0 * A_upper_arm_one
    A_forearms_total = 2.0 * A_forearm_one
    A_arms_total = A_upper_arms_total + A_forearms_total

    A_lower_body = (fields.thigh_clearance_sitting_mm/1000) * (fields.hip_breadth_sitting_mm/1000)

    equalizer_length_mm = compute_equalizer_length_mm(fields)
    A_equalizer_one = (fields.equalizer_diameter_mm/1000) * (equalizer_length_mm/1000)
    A_equalizers = 2.0 * A_equalizer_one

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
        "equalizer_length_mm": equalizer_length_mm,
        # echo useful derived dimensions
        "torso_width_used_mm": torso_width_mm,
        "head_breadth_mm_est": head_breadth_mm,
        "head_height_mm_est": head_height_mm,
        "upper_arm_thickness_mm": upper_arm_thickness_mm,
        "forearm_thickness_mm": forearm_thickness_mm,
        "thigh_clearance_sitting_mm": fields.thigh_clearance_sitting_mm,
    }

__all__ = [
    "ISOFields",
    "create_person_from_height_weight",
    "frontal_areas_iso",
    "compute_equalizer_length_mm",
]


if __name__ == "__main__":
    print("=" * 70)
    print("ISO/TR 7250-2 Anthropometry Demo")
    print("=" * 70)
    
    # Example 1: Male pilot, P50 height/weight
    print("\nExample 1: Male pilot (175cm, 79kg)")
    person_male = create_person_from_height_weight(height_cm=175, weight_kg=79, sex="male")
    areas_male = frontal_areas_iso(person_male)
    
    print(f"  Body dimensions generated: {len(person_male.__dict__)} fields")
    print(f"  Total frontal area: {areas_male['A_total_m2']:.4f} m²")
    print(f"  Equalizer length: {areas_male['equalizer_length_mm']:.0f} mm")
    
    # Example 2: Female pilot, custom dimensions
    print("\nExample 2: Female pilot (162.5cm, 66kg)")
    person_female = create_person_from_height_weight(height_cm=162.5, weight_kg=66, sex="female")
    areas_female = frontal_areas_iso(person_female)
    
    print(f"  Total frontal area: {areas_female['A_total_m2']:.4f} m²")
    print(f"  Upper body area: {areas_female['A_upper_body_m2']:.4f} m²")
    print(f"  Head area: {areas_female['A_head_m2']:.4f} m²")
    print(f"  Arms total area: {areas_female['A_arms_total_m2']:.4f} m²")
    print(f"  Lower body area: {areas_female['A_lower_body_m2']:.4f} m²")
    
    print("\n" + "=" * 70)
    print("Demo complete. Import this module to use in your own code.")
    print("=" * 70)
