"""
Centralized configuration for the Weight Scaling Model.
Edit values here to change inputs without touching the computation code.
"""
from typing import Tuple, List
import os

# Aerodynamic parameters
SPEED_MPS: float = 10.0                    # Speed in m/s (36 km/h)
RHO_AIR: float = 1.225                     # Air density in kg/m^3
BACK_ANGLE_DEG: float = 30.0               # Back angle from vertical (aero position)
UPPER_ARM_ANGLE_DEG: float = 60.0          # Upper arm rotation forward from horizontal
HEAD_EXPOSED_FRACTION: float = 0.30        # Fraction of head exposed in aero position

# Plot/export options
EXPORT_CSV: bool = False                   # Save polar data to CSV files

# Figure saving
# When True, figures will be written to disk (PNG by default) in PLOTS_DIR.
# Files will be overwritten on subsequent runs using the same names.
SAVE_PLOTS: bool = True
PLOTS_DIR: str = os.path.join("plots")    # Relative to current working directory
SAVE_FORMAT: str = "png"                 # e.g., "png", "pdf", "svg"
SAVE_DPI: int = 800                        # Image DPI for raster formats

# When False, plots won't block during generation; a final block can be enabled separately.
SHOW_BLOCKING: bool = False
# If True, block once at the very end so all plot windows stay open. Set to False for non-blocking CI/VS Code runs.
BLOCK_AT_END: bool = False

# Geometry
PROTECTOR_HEIGHT_MM: float = 200.0         # Protector height under hips
EQUALIZER_DIAMETER_MM: float = 25.0        # Equalizer tube diameter

###############################################
# Equalizer gear-weight scaling modes
#
# Base formula shared by both modes (Luc Armant's version):
#   base_length_mm = max(0, (wing_max_kg - 90) * 22)
#   gear_weight_kg = max(0, wing_max_kg - body_weight_kg)
#
# Mode 1 (LEGACY_PER_KG): piecewise per-kg slopes (old behavior)
#   extra_length_mm =
#       if gear_weight_kg > STAGE2_THRESHOLD:  STAGE2_MM_PER_KG * (gear_weight_kg - STAGE2_THRESHOLD)
#       elif gear_weight_kg > STAGE1_THRESHOLD: STAGE1_MM_PER_KG * (gear_weight_kg - STAGE1_THRESHOLD)
#       else: 0
#
# Mode 2 (FIXED_STEPS): two fixed additive jumps based on thresholds (new request)
#   extra_length_mm =
#       if gear_weight_kg > STAGE2_THRESHOLD:  FIXED_STEP_STAGE2_MM
#       elif gear_weight_kg > STAGE1_THRESHOLD: FIXED_STEP_STAGE1_MM
#       else: 0
#   (No scaling with how far above threshold; binary application.)
#
# Select mode via EQUALIZER_GEAR_MODE. Accepted values: "LEGACY_PER_KG", "FIXED_STEPS"
###############################################
EQUALIZER_GEAR_MODE: str = "FIXED_STEPS"  # change to "LEGACY_PER_KG" to restore per-kg slopes

# Thresholds shared by both modes
EQUALIZER_GEAR_STAGE1_THRESHOLD_KG: float = 25.0
EQUALIZER_GEAR_STAGE2_THRESHOLD_KG: float = 35.0

# Legacy per-kg increments (kept for backward compatibility)
EQUALIZER_GEAR_STAGE1_MM_PER_KG: float = 72
EQUALIZER_GEAR_STAGE2_MM_PER_KG: float = 131  # rounded integer slope (previously 130)

# Fixed step lengths (requested): average 5 kg span * slope1, slope2
# stage1: 5 * 72 = 360 mm ; stage2: 5 * 141 = 705 mm (using provided 141 slope for second stage)
FIXED_STEP_STAGE1_MM: float = 360.0
FIXED_STEP_STAGE2_MM: float = 705.0

# Optional cap on total equalizer length (None disables)
EQUALIZER_MAX_LENGTH_MM: float | None = None
# Drag coefficients
CD_STREAMLINED: float = 0.23               # Torso + protector + legs
CD_HEAD: float = 0.40                      # Head with helmet
CD_ARMS: float = 1.00                      # Arms (cylinders perpendicular to flow)
# CD_EQUALIZERS: float = 1.00                # Equalizers (cylinders perpendicular to flow)
CD_EQUALIZERS: float = 0.00                # Equalizers (cylinde8rs perpendicular to flow)
CD_LINES: float = 1.00                     # Lines drag coefficient (cylinders approx perpendicular)

GRAVITY: float = 9.81                      # Gravitational acceleration (m/s^2)

# Wing parameters for drag calculation
REFERENCE_WING_MAX_KG: float = 115.0       # Reference wing size for Cd calculation (calibrated)
REFERENCE_GLIDE_RATIO: float = 11.0        # Legacy field (not used when calibrating by sink)
WING_ASPECT_RATIO: float = 5.5             # Typical CCC wing aspect ratio
KINEMATIC_VISCOSITY: float = 1.5e-5        # Air kinematic viscosity (m^2/s) at sea level
REYNOLDS_EXPONENT: float = 0.12            # Reynolds correction exponent (0.1-0.15 typical for paragliders)
REFERENCE_SPEED_MPS: float = 63.0/3.6      # Calibration speed used to derive CD_WING_REF (63 km/h)
E_OSWALD: float = 0.75                     # Oswald efficiency factor (typ. 0.7–0.85)
PROJECTED_TO_FLAT: float = 0.60            # Projected-to-flat area fraction (S_proj = 0.60 * S_flat)

# Lines equivalent frontal area model for suspension lines
# Source: derived from a Luc Armants spreadsheet that aggregates manufacturer weight optimized 
#         line plans and measured line sets across multiple wing sizes. The coefficients below
#         come from an ordinary least-squares fit of total line frontal area vs.
#         rated wing size (kg), referenced to REF = 85 kg.
# Model:  A_lines [m^2] = C0 + m * (size_kg - REF)
# Units:  C0 [m^2], m [m^2/kg], REF [kg]
LINES_AREA_C0: float = 0.149
LINES_AREA_M: float = 0.00104
LINES_AREA_REF_KG: float = 85.0

# Calibration target: make 115 kg wing the reference, with sink rate 1.93 m/s at 63 km/h
CAL_REF_WING_KG: float = 115.0
CAL_REF_SPEED_MPS: float = 63.0 / 3.6
CAL_REF_SINK_MPS: float = 1.93
CAL_REF_PILOT: Tuple[int, int, str] = (180, 90, "male")  # (height_cm, weight_kg, sex)
CAL_REF_LABEL: str = ("Man" if CAL_REF_PILOT[2] == "male" else "Woman") + f" {int(CAL_REF_PILOT[0])}cm {int(CAL_REF_PILOT[1])}kg"

# Pilot examples with wing sizes (men only), using default gear weight = 20kg
# For each wing size X, pilot_body_weight ≈ X - 20 kg and height derived from ISO mapping
# Approximate male height mapping used here:
#   70 kg -> ~170 cm, 75 kg -> ~173 cm, 85 kg -> ~178 cm,
#   95 kg -> ~184 cm, 105 kg -> ~190 cm, 115 kg -> ~196 cm
###############################
# Pilot populations (men only)
# Updated: capped at max wing size 135 kg.
# Base body weights: 65, 80, 95, 110 kg
# Gear weight scenarios implemented: +20 kg, +30 kg, +40 kg (dropping any >135 total)
# Tuple: (label, height_cm, body_weight_kg, sex, wing_max_weight_kg)
###############################

# Default group (20 kg gear)
PILOT_EXAMPLES: List[Tuple[str, int, int, str, int]] = [
    ("Man 167cm 65kg", 167, 65, "male", 85),
    ("Man 176cm 80kg", 176, 80, "male", 100),
    ("Man 184cm 95kg", 184, 95, "male", 115),
    ("Man 193cm 110kg", 193, 110, "male", 130),
]

# 30 kg gear group (exclude >135)
PILOT_EXAMPLES_30: List[Tuple[str, int, int, str, int]] = [
    ("Man 167cm 65kg", 167, 65, "male", 95),
    ("Man 176cm 80kg", 176, 80, "male", 110),
    ("Man 184cm 95kg", 184, 95, "male", 125),
    # 110kg +30 = 140 (>135) excluded
]

# 40 kg gear group (exclude >135)
PILOT_EXAMPLES_40: List[Tuple[str, int, int, str, int]] = [
    ("Man 167cm 65kg", 167, 65, "male", 105),
    ("Man 176cm 80kg", 176, 80, "male", 120),
    # 110kg +40 = 150 (>135) excluded
]
