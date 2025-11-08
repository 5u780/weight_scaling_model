# Anthropometry Extraction Module

This folder contains a minimal, self-contained anthropometry and frontal area estimation layer extracted from `wsm/iso_model.py` so that other projects can reuse realistic pilot (human) dimensions without coupling to the aerodynamic or glider-specific logic.

## Goals
- Provide functions to create a simplified person profile from height, weight, and sex.
- Estimate key body segment dimensions (shoulders, chest, hips, limbs) scaled from ISO/TR 7250-2 / DIN ergonomic references.
- Compute approximate frontal areas for torso, head, arms, and optional protective equipment.
- Offer equalizer (harness side strap) length estimation logic with configurable modes.

## Non-Goals
- No aerodynamic drag coefficients or performance calculations.
- No plotting or polar generation.
- No dependency on glider wing size or calibration logic.

## Overview
```
create_person(height_cm, weight_kg, sex, **options) -> Person dict
compute_frontal_areas(person) -> Areas dict (m^2)
get_equalizer_length_mm(person) -> float
```

## Configuration Concepts
- Back angle and arm angle affect projected frontal area (reductions based on posture).
- Protector (back protector) can increase effective torso height and area.
- Equalizer gear mode:
  - FIXED_STEPS: adds one of two fixed extra lengths when gear weight exceeds thresholds.
  - LEGACY_PER_KG: piecewise per-kg slope based additive length.

All numeric parameters are centralized in `anthro/config.py`.

## Usage Example
```python
from anthro import person, areas
p = person.create_person(height_cm=180, weight_kg=78, sex="M")
a = areas.compute_frontal_areas(p)
print(a["torso_area_m2"], a["head_area_m2"], a["equalizer_length_mm"])
```

## Extension Ideas
- Replace linear scaling with allometric or percentile-based interpolation.
- Add limb circumference estimation for more detailed drag modeling.
- Introduce dataclasses for stronger typing.

