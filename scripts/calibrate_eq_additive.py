"""calibrate_eq_additive.py
=================================

Purpose
-------
Derive additive equalizer length slopes (mm/kg) for staged gear-weight thresholds
so that different gear loads (20 / 30 / 40 kg) on the *same* reference wing size
(115 kg) produce approximately equal TOTAL aerodynamic drag at the nominal
modeling speed (`cfg.SPEED_MPS`). Total drag here means: pilot aerodynamic drag
(streamlined + head + arms) + suspension line drag + wing drag, with equalizer
drag set to zero for the baseline comparison. By adding equalizer length we
increase equalizer drag, raising lower-drag gear cases toward the highest-drag
reference case, seeking parity in total system drag.

Conceptual Model
----------------
Base equalizer length (already in the core model):
    base_len_mm = max(0, (wing_max_kg - 90) * 22)

Additive staged increment (implemented elsewhere using the outputs of this script):
    extra_len_mm =
        stage2_mm_per_kg * (gear - T2)   if gear > T2
        stage1_mm_per_kg * (gear - T1)   elif gear > T1
        0                                otherwise

Calibration Goal
----------------
Choose stage1_mm_per_kg and stage2_mm_per_kg such that the total drag
    D_total(gear=20) ≈ D_total(gear=30) ≈ D_total(gear=40)
after adding the drag from the *extra* equalizer length. Instead of iterating
inside the aerodynamic core, we compute the total drag DIFFERENCES without any
equalizer drag and analytically back-solve the required mm/kg slopes.

Derivation
----------
Drag contribution from a pair of equalizers (cylinders, normal flow) at speed V:
    ΔD = q * Cd_eq * A_eq_pair
        = q * Cd_eq * 2 * d * L_add / 1000^2   (d, L_add in mm)
        = q * Cd_eq * 2 * d_mm/1000 * (L_add_mm/1000)

If per-kg slope is k (mm/kg) and we apply it over Δgear kg above a threshold:
    L_add_mm = k * Δgear
    ΔD = q * Cd_eq * 2 * d_m * (k * Δgear / 1000)

Solve for k given desired ΔD (drag difference we want to cancel):
    k = ΔD * 1000 / (q * Cd_eq * 2 * d_m * Δgear)

We use Δgear = 5 kg for each stage (30-25, 40-35) with thresholds T1=25, T2=35.
For code brevity the constant 2*Δgear is baked into a simplified rearranged
form (see inline comment near k1/k2 computation).

Usage
-----
Run directly:
    python scripts/calibrate_eq_additive.py
Then copy the printed k1 and k2 into `wsm/config.py`:
    EQUALIZER_GEAR_STAGE1_MM_PER_KG = k1
    EQUALIZER_GEAR_STAGE2_MM_PER_KG = k2

Limitations / Assumptions
-------------------------
1. Linear per-kg slopes; ignores Reynolds variation with added drag (small effect).
2. Uses body drag only (excludes existing equalizer drag) to isolate compensation.
3. Assumes thresholds of 25 kg and 35 kg and 5 kg span for calibration increments.
4. One-speed calibration; multi-speed neutrality would require further fitting.

Outputs
-------
Prints the raw total aerodynamic drag (pilot + lines + wing; equalizers disabled)
for 20/30/40 gear cases and suggested stage1 & stage2 mm/kg slopes. After slopes
are applied in the main model, extra equalizer length adds drag to lower-total-drag
gear cases to align them with the 20 kg reference.
"""

import sys, pathlib

# Ensure project root on path for "wsm" imports when run directly
ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from wsm import config as cfg
from wsm.iso_model import create_person_from_height_weight, frontal_areas_iso, drag_all_iso

# Reference calibration cases (same wing size 115 kg, different gear loads)
cases = [
    {"wing_max": 115, "gear": 20, "pilot_weight": 95, "height_cm": 184},
    {"wing_max": 115, "gear": 30, "pilot_weight": 85, "height_cm": 178},
    {"wing_max": 115, "gear": 40, "pilot_weight": 75, "height_cm": 173},
]

records = []
for c in cases:
    person = create_person_from_height_weight(
        c["height_cm"], c["pilot_weight"], "male",
        pilot_max_weight_kg=c["wing_max"],
        back_angle_deg=cfg.BACK_ANGLE_DEG,
        upper_arm_angle_deg=cfg.UPPER_ARM_ANGLE_DEG,
        protector_height_mm=cfg.PROTECTOR_HEIGHT_MM,
        equalizer_diameter_mm=cfg.EQUALIZER_DIAMETER_MM,
    )
    areas = frontal_areas_iso(person)
    D = drag_all_iso(
        areas,
        speed_mps=cfg.SPEED_MPS,
        rho_air=cfg.RHO_AIR,
        Cd_streamlined=cfg.CD_STREAMLINED,
        Cd_head=cfg.CD_HEAD,
        head_exposed_fraction=cfg.HEAD_EXPOSED_FRACTION,
        Cd_arms=cfg.CD_ARMS,
        Cd_equalizers=0.0,  # zero equalizer drag for clean baseline differences
    )
    # Pilot aerodynamic drag (excluding equalizers) + line drag + wing drag
    pilot_clean_drag = D["Drag_streamlined_N"] + D["Drag_head_N"] + D["Drag_arms_N"]
    q_case = 0.5 * cfg.RHO_AIR * (cfg.SPEED_MPS ** 2)
    # Lines frontal area from wing size
    from wsm.aero import lines_area_from_wing_size, calculate_wing_drag
    A_lines_m2 = lines_area_from_wing_size(c["wing_max"])
    Drag_lines_N = q_case * cfg.CD_LINES * A_lines_m2
    # Wing drag using current reference Cd (calibration had been performed elsewhere)
    # We supply pilot_drag_N without equalizers so baseline isn't biased.
    pilot_drag_no_eq = pilot_clean_drag + Drag_lines_N
    # Use existing wing Cd from config reference (approx); for simplicity reuse REFERENCE_WING_MAX_KG if needed.
    # We call calculate_wing_drag with wing_max from case and pilot drag.
    from wsm.aero import calibrate_cd_wing_ref
    # Calibrate wing Cd once (equalizers off) then reuse.
    if 'CD_WING_REF_CACHE' not in globals():
        globals()['CD_WING_REF_CACHE'] = calibrate_cd_wing_ref(CD_EQUALIZERS_FOR_CAL=0.0)
    W = calculate_wing_drag(c["wing_max"], pilot_drag_no_eq, globals()['CD_WING_REF_CACHE'], cfg.SPEED_MPS, cfg.RHO_AIR, e_oswald=cfg.E_OSWALD)
    total_drag_no_eq = pilot_drag_no_eq + W["drag_wing_N"]
    records.append({**c, "total_drag_N": total_drag_no_eq})

ref = next(r for r in records if r["gear"] == 20)  # 20 kg gear reference (total drag)

q = 0.5 * cfg.RHO_AIR * (cfg.SPEED_MPS ** 2)  # dynamic pressure at calibration speed
Cd_eq = cfg.CD_EQUALIZERS
d_m = cfg.EQUALIZER_DIAMETER_MM / 1000.0  # equalizer diameter in meters

# Stage spans (assumed 5 kg each over their thresholds)
DELTA_GEAR_STAGE = 5.0

r30 = next(r for r in records if r["gear"] == 30)
DeltaD30 = ref["total_drag_N"] - r30["total_drag_N"]  # positive: 30kg case total drag lower, needs added drag

r40 = next(r for r in records if r["gear"] == 40)
DeltaD40 = ref["total_drag_N"] - r40["total_drag_N"]

denom = q * Cd_eq * 2.0 * d_m * (DELTA_GEAR_STAGE / 1000.0)  # full denominator from derivation
if denom > 0:
    k1 = DeltaD30 / denom  # mm/kg (float)
    k2 = DeltaD40 / denom  # mm/kg (float)
else:
    k1 = k2 = 0.0

# Round slopes to whole millimeters per kg for use in config
k1_rounded = int(round(k1))
k2_rounded = int(round(k2))

print("Calibration (target: equal TOTAL drag across gear within 115 kg wing)\n")
for r in records:
    print(f"  {r['gear']}kg gear: total_drag_no_eq = {r['total_drag_N']:.4f} N")
print("")
print(f"Suggested EQUALIZER_GEAR_STAGE1_MM_PER_KG = {k1:.2f} mm/kg")
print(f"Suggested EQUALIZER_GEAR_STAGE2_MM_PER_KG = {k2:.2f} mm/kg")

print("\nCopy these integer values into wsm/config.py:")
print(f"  EQUALIZER_GEAR_STAGE1_MM_PER_KG = {k1_rounded}  # mm/kg")
print(f"  EQUALIZER_GEAR_STAGE2_MM_PER_KG = {k2_rounded}  # mm/kg")

# Optional visualization: show body drag before/after adding equalizer increments
try:
    import os
    import matplotlib.pyplot as plt

    # Sort by gear for stable plotting
    records_sorted = sorted(records, key=lambda r: r["gear"])
    gears = [r["gear"] for r in records_sorted]
    total_vals = [r["total_drag_N"] for r in records_sorted]

    # Extra equalizer length applied per gear (mm); 20 has none
    extra_len_mm = {
        20: 0.0,
        30: max(0.0, float(k1_rounded)) * 5.0,
        40: max(0.0, float(k2_rounded)) * 5.0,
    }
    # Drag added by the extra equalizer length (pair of equalizers)
    D_added = [q * Cd_eq * 2.0 * d_m * (extra_len_mm[g] / 1000.0) for g in gears]
    total_corrected = [total_vals[i] + D_added[i] for i in range(len(gears))]

    plt.figure(figsize=(9, 4.8))
    ax1 = plt.subplot(1, 2, 1)
    ax1.set_title("Total drag (pilot+lines+wing, eq off) before vs corrected")
    ax1.plot(gears, total_vals, "o-", label="baseline (no eq add)")
    ax1.plot(gears, total_corrected, "s-", label="after eq add")
    ax1.axhline(records_sorted[0]["total_drag_N"], color="gray", linestyle=":", label="20kg reference")
    for g, d0, d1 in zip(gears, total_vals, total_corrected):
        ax1.annotate(f"{d0:.2f}", (g, d0), textcoords="offset points", xytext=(0, 6), ha="center", fontsize=8)
        ax1.annotate(f"{d1:.2f}", (g, d1), textcoords="offset points", xytext=(0, -12), ha="center", fontsize=8, color="tab:orange")
    ax1.set_xlabel("gear weight (kg)")
    ax1.set_ylabel("Drag (N)")
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=8)

    ax2 = plt.subplot(1, 2, 2)
    ax2.set_title("Extra equalizer length applied (mm)")
    bars = ax2.bar([str(g) for g in gears], [extra_len_mm[g] for g in gears], color=["tab:blue", "tab:orange", "tab:green"]) 
    for b in bars:
        ax2.annotate(f"{b.get_height():.0f} mm", (b.get_x() + b.get_width()/2, b.get_height()),
                     ha="center", va="bottom", fontsize=8, xytext=(0, 3), textcoords="offset points")
    ax2.set_xlabel("gear weight (kg)")
    ax2.set_ylabel("Added length (mm)")
    ax2.grid(True, axis="y", alpha=0.3)
    ax2.text(0.5, -0.18,
             f"k1 = {k1_rounded} mm/kg over 25→30 kg;  k2 = {k2_rounded} mm/kg over 35→40 kg\n"
             "Added equalizer drag raises lower total-drag cases to match 20 kg baseline.",
             transform=ax2.transAxes, ha="center", va="top", fontsize=8)

    plt.tight_layout()
    # Save alongside other plots
    try:
        os.makedirs(cfg.PLOTS_DIR, exist_ok=True)
        out = os.path.join(cfg.PLOTS_DIR, f"calibrate_eq_additive.{cfg.SAVE_FORMAT}")
        plt.savefig(out, dpi=cfg.SAVE_DPI, format=cfg.SAVE_FORMAT, bbox_inches="tight")
        print(f"[saved] {out}")
    except Exception:
        pass

    # Show non-blocking consistent with project defaults
    try:
        plt.show(block=getattr(cfg, "SHOW_BLOCKING", False))
        if not getattr(cfg, "SHOW_BLOCKING", False):
            plt.pause(0.001)
    except Exception:
        try:
            plt.show()
        except Exception:
            pass
except Exception as exc:
    print(f"[note] Skipping visualization (reason: {exc})")
