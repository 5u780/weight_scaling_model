# Weight Scaling Model (structured)

This repo now includes a small package (`wsm/`) to separate configuration, calculations, and plotting. 
The original scripts (e.g., `weight_scaling_model_v5.py`) still work; the new entry point is `run_model.py`.

## Structure

- `wsm/config.py` — all user-editable inputs in one place (CDs, geometry, speed, plotting flags, pilot list, calibration target).
- `wsm/iso_model.py` — ISO/TR 7250-2 anthropometric scaling and body drag breakdown.
- `wsm/aero.py` — wing polar model, lines area model, and reference calibration.
- `wsm/plots.py` — plotting utilities for reference polar and multi-pilot polars.
- `run_model.py` — orchestrates calibration, tabular outputs, and plots using the modules above.

## How to run

From the project folder:

```powershell
# Run the structured entry point
python .\run_model.py
```

If matplotlib is missing, install it (along with pandas/numpy):

```powershell
pip install matplotlib pandas numpy
```

## Change inputs

Edit `wsm/config.py`:

- Aerodynamics: `SPEED_MPS`, `RHO_AIR`, `HEAD_EXPOSED_FRACTION`, etc.
- Geometry: `PROTECTOR_HEIGHT_MM`, `EQUALIZER_DIAMETER_MM`.
- Drag coefficients: `CD_STREAMLINED`, `CD_HEAD`, `CD_ARMS`, `CD_EQUALIZERS`, `CD_LINES`.
- Wing/aero model: `REFERENCE_WING_MAX_KG`, `WING_ASPECT_RATIO`, `E_OSWALD`, `PROJECTED_TO_FLAT`, etc.
- Calibration target: `CAL_REF_PILOT`, `CAL_REF_WING_KG`, `CAL_REF_SPEED_MPS`, `CAL_REF_SINK_MPS`.
- Pilot list: `PILOT_EXAMPLES`.
- Plot/export flags: `EXPORT_CSV`.

## Notes

- `CD_WING_REF` is computed at runtime via `calibrate_cd_wing_ref()`.
- The plots include line drag and equalizer status. Speed sweep is capped at 65 km/h and sink axis is inverted.
- The multi-pilot deltas are referenced to `(CAL_REF_PILOT, REFERENCE_WING_MAX_KG)`.
