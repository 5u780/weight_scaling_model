"""
Anthropometry configuration (standalone).
Only non-aero parameters live here.
"""
from __future__ import annotations

"""
Mirror the original variable names from wsm/config.py to avoid divergence.
Only non-aero parameters are duplicated here for standalone use.
"""

# Posture angles (degrees) — align with wsm/config.py defaults
BACK_ANGLE_DEG: float = 30.0
UPPER_ARM_ANGLE_DEG: float = 60.0

# Head exposure fraction for area (0..1)
HEAD_EXPOSED_FRACTION: float = 0.75

# Protector: additional effective height in mm
PROTECTOR_HEIGHT_MM: float = 200.0

# Equalizer geometry
EQUALIZER_DIAMETER_MM: float = 25.0

# Equalizer gear-weight scaling modes (names match wsm/config.py)
EQUALIZER_GEAR_MODE: str = "FIXED_STEPS"  # or "LEGACY_PER_KG"
EQUALIZER_GEAR_STAGE1_THRESHOLD_KG: float = 25.0
EQUALIZER_GEAR_STAGE2_THRESHOLD_KG: float = 35.0
# Fixed-step extra lengths
FIXED_STEP_STAGE1_MM: float = 360.0
FIXED_STEP_STAGE2_MM: float = 705.0
# Optional cap
EQUALIZER_MAX_LENGTH_MM: float | None = None

# Legacy per-kg slopes (kept for compatibility)
EQUALIZER_GEAR_STAGE1_MM_PER_KG: float = 72
EQUALIZER_GEAR_STAGE2_MM_PER_KG: float = 131  # rounded integer per previous discussion
