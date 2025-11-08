"""
Structured runner for the Weight Scaling Model.
- Configuration: wsm/config.py
- Anthropometrics and body drag: wsm/iso_model.py
- Aerodynamics and calibration: wsm/aero.py
- Plotting: wsm/plots.py

This script mirrors the outputs from weight_scaling_model_v5.py but with clearer
separation of concerns and easier configuration.
"""
import pandas as pd

from wsm import config as cfg
from wsm.iso_model import create_person_from_height_weight, frontal_areas_iso, drag_all_iso, estimate_z_score
from wsm.aero import calculate_wing_drag, lines_area_from_wing_size
from wsm.model import get_cd_wing_ref, compute_total_drag
from wsm.plots import (
    plot_polar_for_reference,
    plot_polars_for_all_pilots,
    plot_115kg_ballast_comparison,
)


def main():
    # Calibrated wing Cd (cached, excludes equalizer drag for baseline neutrality)
    CD_WING_REF = get_cd_wing_ref()

    print("\n" + "="*80)
    print("SCALING MODEL: Dimensions calculated from height & weight using SD")
    print("="*80)

    print("\n" + "="*80)
    print("USER INPUT VARIABLES")
    print("="*80)

    print("\nAerodynamic Configuration:")
    print(f"  Speed: {cfg.SPEED_MPS} m/s ({cfg.SPEED_MPS * 3.6:.0f} km/h)")
    print(f"  Air density: {cfg.RHO_AIR} kg/m^3")
    print(f"  Back angle from vertical: {cfg.BACK_ANGLE_DEG} degrees")
    print(f"  Upper arm angle forward: {cfg.UPPER_ARM_ANGLE_DEG} degrees")
    print(f"  Head exposed fraction: {cfg.HEAD_EXPOSED_FRACTION * 100:.0f}%")

    print("\nGeometry:")
    print(f"  Protector height: {cfg.PROTECTOR_HEIGHT_MM} mm")
    print(f"  Equalizer tube diameter: {cfg.EQUALIZER_DIAMETER_MM} mm")

    print("\nDrag Coefficients:")
    print(f"  Streamlined body (torso + protector + legs): {cfg.CD_STREAMLINED}")
    print(f"  Head (with helmet): {cfg.CD_HEAD}")
    print(f"  Arms: {cfg.CD_ARMS}")
    print(f"  Equalizers: {cfg.CD_EQUALIZERS}")
    print(f"  Lines: {cfg.CD_LINES}")

    print("\nWing Parameters:")
    print(f"  Reference wing size: {cfg.REFERENCE_WING_MAX_KG} kg")
    print(f"  Reference glide ratio: {cfg.REFERENCE_GLIDE_RATIO}:1")
    print(f"  Wing aspect ratio: {cfg.WING_ASPECT_RATIO}")
    print(f"  Reynolds exponent: {cfg.REYNOLDS_EXPONENT}")
    print(f"  Reference speed (calibration): {cfg.REFERENCE_SPEED_MPS} m/s")
    print(f"  Oswald efficiency (e): {cfg.E_OSWALD}")
    print(f"  Projected/flat area fraction: {cfg.PROJECTED_TO_FLAT:.2f}")
    print(f"  Gravity: {cfg.GRAVITY} m/s^2")

    print("\nPilot Examples:")
    for label, height_cm, weight_kg, sex, wing_max in cfg.PILOT_EXAMPLES:
        person = create_person_from_height_weight(
            height_cm, weight_kg, sex,
            pilot_max_weight_kg=wing_max,
            back_angle_deg=cfg.BACK_ANGLE_DEG,
            upper_arm_angle_deg=cfg.UPPER_ARM_ANGLE_DEG,
            protector_height_mm=cfg.PROTECTOR_HEIGHT_MM,
            equalizer_diameter_mm=cfg.EQUALIZER_DIAMETER_MM,
        )
        areas = frontal_areas_iso(person)
        eq_len = areas.get("equalizer_length_mm", 0.0)
        base_len = max(0.0, (wing_max - 90.0) * 22.0)
        gear_weight = max(0.0, wing_max - weight_kg)
        extra = eq_len - base_len
        print(f"  {label}: Wing max {wing_max}kg, Gear {gear_weight:.0f}kg, Equalizer {eq_len:.0f}mm (base {base_len:.0f} + extra {extra:.0f})")

    print("\n" + "="*80)
    print("CALCULATED RESULTS")
    print("="*80)
    print(f"\nReference Cd_wing (fixed baseline value): {CD_WING_REF:.6f}")
    print(
        "Current drag coefficients used: "
        f"CD_STREAMLINED={cfg.CD_STREAMLINED}, CD_HEAD={cfg.CD_HEAD}, "
        f"CD_ARMS={cfg.CD_ARMS}, CD_EQUALIZERS={cfg.CD_EQUALIZERS}"
    )
    print("Note: Wing drag is fixed based on reference calibration, but pilot drag changes with CD values.")

    # Calculate all results
    scaled_rows = []
    for label, height_cm, weight_kg, sex, wing_max_weight in cfg.PILOT_EXAMPLES:
        person = create_person_from_height_weight(
            height_cm, weight_kg, sex,
            pilot_max_weight_kg=wing_max_weight,
            back_angle_deg=cfg.BACK_ANGLE_DEG,
            upper_arm_angle_deg=cfg.UPPER_ARM_ANGLE_DEG,
            protector_height_mm=cfg.PROTECTOR_HEIGHT_MM,
            equalizer_diameter_mm=cfg.EQUALIZER_DIAMETER_MM,
        )
        z_score = estimate_z_score(height_cm, sex)
        result = compute_total_drag(person, wing_max_weight, cfg.SPEED_MPS, CD_WING_REF, include_equalizer_drag=True)
        scaled_rows.append({
            "label": label,
            "height_cm": height_cm,
            "pilot_weight_kg": weight_kg,
            "body_weight_kg": person.body_weight_kg,
            "total_weight_kg": wing_max_weight,
            "wing_max_kg": wing_max_weight,
            "sex": sex,
            "z_score": z_score,
            "cervicale_mm": person.cervicale_height_sitting_mm,
            "shoulder_breadth_mm": person.shoulder_breadth_biacromial_mm,
            "hip_breadth_mm": person.hip_breadth_sitting_mm,
            "head_circ_mm": person.head_circumference_mm,
            "face_length_mm": person.face_length_mm,
            "shoulder_elbow_mm": person.shoulder_elbow_length_mm,
            "elbow_grip_mm": person.elbow_grip_length_mm,
            "arm_circ_mm": person.arm_circumference_flexed_mm,
            "forearm_circ_mm": person.forearm_circumference_flexed_mm,
            "thigh_circ_mm": person.thigh_circumference_mm,
            "equalizer_length_mm": result["equalizer_length_mm"],
            # Areas
            "A_total_m2": result["A_total_m2"],
            "A_streamlined_m2": result["A_streamlined_m2"],
            "A_head_m2": result["A_head_m2"],
            "A_upper_arms_m2": result["A_upper_arms_m2"],
            "A_forearms_m2": result["A_forearms_m2"],
            "A_arms_m2": result["A_arms_m2"],
            "A_equalizers_m2": result["A_equalizers_m2"],
            "A_lines_m2": result["A_lines_m2"],
            # Pilot drag components
            "Drag_streamlined_N": result["Drag_streamlined_N"],
            "Drag_head_N": result["Drag_head_N"],
            "Drag_upper_arms_N": result["Drag_upper_arms_N"],
            "Drag_forearms_N": result["Drag_forearms_N"],
            "Drag_arms_N": result["Drag_arms_N"],
            "Drag_equalizers_N": result["Drag_equalizers_N"],
            "Drag_lines_N": result["Drag_lines_N"],
            "Drag_pilot_N": result["pilot_drag_N"],
            "Drag_wing_N": result["wing_drag_N"],
            "Drag_total_N": result["drag_total_N"],
            # Wing aero/geom
            "wing_area_m2": result["wing_area_m2"],
            "wing_flat_area_m2": result["wing_flat_area_m2"],
            "CL": result["CL"],
            "Cd0_flat": result["Cd0_flat"],
            "CDi": result["CDi"],
            "CD_total": result["CD_total"],
            "Reynolds": result["Reynolds"],
            "re_factor": result["re_factor"],
            # Performance
            "lift_force_N": result["lift_force_N"],
            "glide_ratio": result["glide_ratio"],
        })

    df_scaled = pd.DataFrame(scaled_rows)

    ref_mask = (df_scaled["wing_max_kg"] == cfg.REFERENCE_WING_MAX_KG) & (df_scaled["label"] == cfg.CAL_REF_LABEL)
    if ref_mask.any():
        ref_row = df_scaled.loc[ref_mask].iloc[0]
        df_scaled["glide_delta_vs_ref"] = df_scaled["glide_ratio"] - ref_row["glide_ratio"]
        df_scaled["drag_delta_vs_ref"] = df_scaled["Drag_total_N"] - ref_row["Drag_total_N"]
    else:
        ref_mask_wing = df_scaled["wing_max_kg"] == cfg.REFERENCE_WING_MAX_KG
        if ref_mask_wing.any():
            ref_row = df_scaled.loc[ref_mask_wing].iloc[0]
            df_scaled["glide_delta_vs_ref"] = df_scaled["glide_ratio"] - ref_row["glide_ratio"]
            df_scaled["drag_delta_vs_ref"] = df_scaled["Drag_total_N"] - ref_row["Drag_total_N"]
        else:
            df_scaled["glide_delta_vs_ref"] = float("nan")
            df_scaled["drag_delta_vs_ref"] = float("nan")

    print("\nBody dimensions (all in mm except where noted):")
    print(df_scaled[["label", "height_cm", "pilot_weight_kg", "body_weight_kg", "wing_max_kg", "z_score", "cervicale_mm",
                     "shoulder_breadth_mm", "hip_breadth_mm"]].round(1).to_string(index=False))
    print("\n" + df_scaled[["label", "head_circ_mm", "face_length_mm", "shoulder_elbow_mm",
                            "elbow_grip_mm"]].round(1).to_string(index=False))
    print("\n" + df_scaled[["label", "arm_circ_mm", "forearm_circ_mm",
                            "thigh_circ_mm"]].round(1).to_string(index=False))
    print("\nEqualizers:")
    print(df_scaled[["label", "wing_max_kg", "equalizer_length_mm"]].round(1).to_string(index=False))
    print("\nAreas (m^2):")
    print(df_scaled[["label", "A_streamlined_m2", "A_head_m2", "A_upper_arms_m2",
                     "A_forearms_m2", "A_arms_m2", "A_equalizers_m2", "A_lines_m2", "A_total_m2"]].round(4).to_string(index=False))
    print("\nArm dimensions and drag (at current speed):")
    print(df_scaled[["label", "arm_circ_mm", "A_upper_arms_m2", "Drag_upper_arms_N",
                     "forearm_circ_mm", "A_forearms_m2", "Drag_forearms_N"]].round(4).to_string(index=False))

    print("\nDrag breakdown (at current speed):")
    print(df_scaled[["label", "Drag_streamlined_N", "Drag_head_N", "Drag_arms_N",
                     "Drag_equalizers_N", "Drag_lines_N", "Drag_pilot_N", "Drag_wing_N", "Drag_total_N"]].round(4).to_string(index=False))

    print("\nWing Aerodynamics:")
    print(df_scaled[["label", "wing_max_kg", "wing_flat_area_m2", "wing_area_m2"]].round(4).to_string(index=False))

    print("\nWing Drag Coefficients & Reynolds:")
    print(df_scaled[["label", "CL", "Cd0_flat", "CDi", "CD_total", "Reynolds", "re_factor"]].round(6).to_string(index=False))

    print("\nGlide Performance:")
    print(df_scaled[["label", "pilot_weight_kg", "total_weight_kg", "lift_force_N",
                     "Drag_total_N", "drag_delta_vs_ref", "glide_ratio", "glide_delta_vs_ref"]].round(2).to_string(index=False))

    # Calibration check
    try:
        cal_h, cal_w, cal_sex = cfg.CAL_REF_PILOT
        cal_person = create_person_from_height_weight(
            cal_h, cal_w, cal_sex,
            pilot_max_weight_kg=cfg.CAL_REF_WING_KG,
            back_angle_deg=cfg.BACK_ANGLE_DEG,
            upper_arm_angle_deg=cfg.UPPER_ARM_ANGLE_DEG,
            protector_height_mm=cfg.PROTECTOR_HEIGHT_MM,
            equalizer_diameter_mm=cfg.EQUALIZER_DIAMETER_MM,
        )
        cal_areas = frontal_areas_iso(cal_person)
        v_cal = cfg.CAL_REF_SPEED_MPS
        # Re-run reference condition WITHOUT equalizers to check calibration target
        D_body = drag_all_iso(
            cal_areas,
            speed_mps=v_cal,
            rho_air=cfg.RHO_AIR,
            Cd_streamlined=cfg.CD_STREAMLINED,
            Cd_head=cfg.CD_HEAD,
            head_exposed_fraction=cfg.HEAD_EXPOSED_FRACTION,
            Cd_arms=cfg.CD_ARMS,
            Cd_equalizers=0.0,
        )
        q_cal = 0.5 * cfg.RHO_AIR * v_cal**2
        A_lines_cal = lines_area_from_wing_size(cfg.CAL_REF_WING_KG)
        D_lines_cal = q_cal * cfg.CD_LINES * A_lines_cal
        pilot_drag_cal = D_body["Drag_total_N"] + D_lines_cal
        W_cal = calculate_wing_drag(cfg.CAL_REF_WING_KG, pilot_drag_cal, CD_WING_REF, v_cal, cfg.RHO_AIR, e_oswald=cfg.E_OSWALD)
        glide_cal = W_cal["glide_ratio"]
        sink_cal = v_cal / glide_cal if glide_cal > 0 else float("nan")
        print(f"\nCalibration check (115 kg @ {cfg.CAL_REF_SPEED_MPS*3.6:.0f} km/h, no equalizers): sink ≈ {sink_cal:.2f} m/s (target {cfg.CAL_REF_SINK_MPS:.2f})")
        if cfg.CD_EQUALIZERS > 0:
            # Show effect when applying current equalizer Cd at calibration speed
            D_body_eq = drag_all_iso(
                cal_areas,
                speed_mps=v_cal,
                rho_air=cfg.RHO_AIR,
                Cd_streamlined=cfg.CD_STREAMLINED,
                Cd_head=cfg.CD_HEAD,
                head_exposed_fraction=cfg.HEAD_EXPOSED_FRACTION,
                Cd_arms=cfg.CD_ARMS,
                Cd_equalizers=cfg.CD_EQUALIZERS,
            )
            q_cal_eq = 0.5 * cfg.RHO_AIR * v_cal**2
            A_lines_cal_eq = lines_area_from_wing_size(cfg.CAL_REF_WING_KG)
            D_lines_cal_eq = q_cal_eq * cfg.CD_LINES * A_lines_cal_eq
            pilot_drag_cal_eq = D_body_eq["Drag_total_N"] + D_lines_cal_eq
            W_cal_eq = calculate_wing_drag(cfg.CAL_REF_WING_KG, pilot_drag_cal_eq, CD_WING_REF, v_cal, cfg.RHO_AIR, e_oswald=cfg.E_OSWALD)
            glide_cal_eq = W_cal_eq["glide_ratio"]
            sink_cal_eq = v_cal / glide_cal_eq if glide_cal_eq > 0 else float("nan")
            print(f"  With equalizers Cd={cfg.CD_EQUALIZERS}: sink ≈ {sink_cal_eq:.2f} m/s (higher sink expected)")
    except Exception:
        pass

    # Plots
    plot_polar_for_reference(CD_WING_REF)
    plot_polars_for_all_pilots(CD_WING_REF)
    plot_115kg_ballast_comparison(CD_WING_REF)

    # Keep all figures open at the end of the run only if configured
    try:
        import matplotlib.pyplot as plt
        if getattr(cfg, "BLOCK_AT_END", False) or getattr(cfg, "SHOW_BLOCKING", False):
            plt.show()
    except Exception:
        pass


if __name__ == "__main__":
    main()
