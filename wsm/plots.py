"""
Plotting utilities for polars. These functions depend on the model but keep
matplotlib-specific code out of the core computations.
"""
import os
import numpy as np
import pandas as pd

from . import config as cfg
from .iso_model import create_person_from_height_weight, frontal_areas_iso, drag_all_iso
from .aero import lines_area_from_wing_size, calculate_wing_drag


def plot_polar_for_reference(CD_WING_REF: float):
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        print("\n[polar] matplotlib is not available. Install it to see the polar plot:")
        print("  python -m pip install matplotlib")
        print(f"  (reason: {exc})")
        return

    ref_label = cfg.CAL_REF_LABEL
    ref_h, ref_w, ref_sex = cfg.CAL_REF_PILOT
    ref_wing = cfg.REFERENCE_WING_MAX_KG

    ref_person = create_person_from_height_weight(
        ref_h, ref_w, ref_sex,
        pilot_max_weight_kg=ref_wing,
        back_angle_deg=cfg.BACK_ANGLE_DEG,
        upper_arm_angle_deg=cfg.UPPER_ARM_ANGLE_DEG,
        protector_height_mm=cfg.PROTECTOR_HEIGHT_MM,
        equalizer_diameter_mm=cfg.EQUALIZER_DIAMETER_MM,
    )

    areas = frontal_areas_iso(ref_person)

    speeds = np.linspace(8.0, 65.0/3.6, 57)

    rows = []
    A_lines_m2 = lines_area_from_wing_size(ref_wing)
    for v in speeds:
        D = drag_all_iso(
            areas,
            speed_mps=v,
            rho_air=cfg.RHO_AIR,
            Cd_streamlined=cfg.CD_STREAMLINED,
            Cd_head=cfg.CD_HEAD,
            head_exposed_fraction=cfg.HEAD_EXPOSED_FRACTION,
            Cd_arms=cfg.CD_ARMS,
            Cd_equalizers=cfg.CD_EQUALIZERS,
        )
        q_v = 0.5 * cfg.RHO_AIR * v**2
        Drag_lines_N = q_v * cfg.CD_LINES * A_lines_m2
        pilot_drag_N = D["Drag_total_N"] + Drag_lines_N

        W = calculate_wing_drag(ref_wing, pilot_drag_N, CD_WING_REF, v, cfg.RHO_AIR, e_oswald=cfg.E_OSWALD)
        glide = W["glide_ratio"]
        sink = v / glide if glide > 0 else float("nan")

        rows.append({
            "speed_mps": v,
            "speed_kmh": v * 3.6,
            "pilot_drag_N": pilot_drag_N,
            "wing_drag_N": W["drag_wing_N"],
            "drag_total_N": W["drag_total_N"],
            "glide_ratio": glide,
            "sink_mps": sink,
            "CL": W["CL"],
            "Cd0_flat": W["Cd0_flat"],
            "CDi": W["CDi"],
            "CD_total": W["CD_total"],
            "A_lines_m2": A_lines_m2,
            "Drag_lines_N": Drag_lines_N,
        })

    df_polar = pd.DataFrame(rows)

    idx_min_sink = df_polar["sink_mps"].idxmin()
    idx_best_ld = df_polar["glide_ratio"].idxmax()

    # Cast to float for type-checker friendliness
    v_min_sink = float(df_polar.loc[idx_min_sink, "speed_kmh"])
    sink_min = float(df_polar.loc[idx_min_sink, "sink_mps"])
    v_best_ld = float(df_polar.loc[idx_best_ld, "speed_kmh"])
    best_ld = float(df_polar.loc[idx_best_ld, "glide_ratio"])

    print(f"\nPolar summary (reference {ref_label} on {int(ref_wing)} kg wing):")
    print(f"  Min sink  ~ {sink_min:.2f} m/s at ~ {v_min_sink:.0f} km/h")
    print(f"  Best L/D  ~ {best_ld:.2f}:1 at ~ {v_best_ld:.0f} km/h")

    if cfg.EXPORT_CSV:
        try:
            df_polar.to_csv("polar_reference.csv", index=False)
            print("  [saved] polar_reference.csv")
        except Exception as exc:
            print(f"  [warn] failed to save polar_reference.csv: {exc}")

    plt.figure(figsize=(8, 5))
    plt.plot(df_polar["speed_kmh"], df_polar["sink_mps"], label="sink (m/s)", linewidth=2)
    # Use list wrappers for compatibility with various backends/type checkers
    plt.scatter([v_min_sink], [sink_min], color="red", zorder=5, label="min sink")
    plt.annotate(f"min sink\n{sink_min:.2f} m/s @ {v_min_sink:.0f} km/h",
                 xy=(v_min_sink, sink_min), xytext=(v_min_sink+3, sink_min+0.15),
                 arrowprops=dict(arrowstyle="->", color="red"), fontsize=8, color="red")
    sink_best_ld = df_polar.loc[idx_best_ld, "sink_mps"]
    plt.scatter(v_best_ld, sink_best_ld, color="green", zorder=5, label="best L/D")
    plt.annotate(f"best L/D {best_ld:.2f}:1\n@ {v_best_ld:.0f} km/h",
                 xy=(v_best_ld, sink_best_ld), xytext=(v_best_ld+3, sink_best_ld+0.15),
                 arrowprops=dict(arrowstyle="->", color="green"), fontsize=8, color="green")
    ax = plt.gca()
    ax.invert_yaxis()
    eq_status = "True" if cfg.CD_EQUALIZERS >= 0.5 else "False"
    plt.title("Polar: sink rate vs speed (reference pilot/wing)")
    plt.xlabel("Speed (km/h)")
    plt.ylabel("Sink (m/s)")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.text(0.98, 0.02, f"Equalizers: {eq_status}", transform=ax.transAxes,
             ha="right", va="bottom", fontsize=9,
             bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.7))
    plt.tight_layout()

    # Optional: save figure
    if getattr(cfg, "SAVE_PLOTS", False):
        try:
            os.makedirs(cfg.PLOTS_DIR, exist_ok=True)
            equ_flag = "True" if cfg.CD_EQUALIZERS >= 0.5 else "False"
            out_path = os.path.join(cfg.PLOTS_DIR, f"polar_reference_equ={equ_flag}.{cfg.SAVE_FORMAT}")
            plt.savefig(out_path, dpi=cfg.SAVE_DPI, format=cfg.SAVE_FORMAT, bbox_inches="tight")
            print(f"[saved] {out_path}")
        except Exception as exc:
            print(f"[warn] failed to save figure: {exc}")
    try:
        plt.show(block=cfg.SHOW_BLOCKING)
        if not cfg.SHOW_BLOCKING:
            plt.pause(0.001)
    except Exception:
        # Be robust to backend quirks on non-blocking shows
        try:
            plt.show()
        except Exception:
            pass


def plot_polars_for_all_pilots(CD_WING_REF: float):
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        print("\n[polar] matplotlib is not available. Install it to see the polar plot:")
        print("  python -m pip install matplotlib")
        print(f"  (reason: {exc})")
        return

    speeds = np.linspace(8.0, 65.0/3.6, 57)
    plt.figure(figsize=(9, 6))
    all_rows = []

    # Group A: 20kg gear weight (existing list) -> cool spectrum
    groupA = getattr(cfg, 'PILOT_EXAMPLES', [])
    nA = len(groupA)
    # Wider color sweep for better separation
    colorsA = [plt.get_cmap("Blues")(x) for x in np.linspace(0.20, 0.85, max(nA, 1))]
    linestyles = ['-', '--', '-.', ':']
    markers = ['o', 's', 'D', '^', 'v', 'P', 'X']

    for idx, (label, height_cm, weight_kg, sex, wing_max) in enumerate(groupA):
        person = create_person_from_height_weight(
            height_cm, weight_kg, sex,
            pilot_max_weight_kg=wing_max,
            back_angle_deg=cfg.BACK_ANGLE_DEG,
            upper_arm_angle_deg=cfg.UPPER_ARM_ANGLE_DEG,
            protector_height_mm=cfg.PROTECTOR_HEIGHT_MM,
            equalizer_diameter_mm=cfg.EQUALIZER_DIAMETER_MM,
        )
        areas = frontal_areas_iso(person)
        eq_len_mm = float(areas.get("equalizer_length_mm", 0.0))

        rows = []
        A_lines_m2 = lines_area_from_wing_size(wing_max)
        for v in speeds:
            D = drag_all_iso(
                areas,
                speed_mps=v,
                rho_air=cfg.RHO_AIR,
                Cd_streamlined=cfg.CD_STREAMLINED,
                Cd_head=cfg.CD_HEAD,
                head_exposed_fraction=cfg.HEAD_EXPOSED_FRACTION,
                Cd_arms=cfg.CD_ARMS,
                Cd_equalizers=cfg.CD_EQUALIZERS,
            )
            q_v = 0.5 * cfg.RHO_AIR * v**2
            Drag_lines_N = q_v * cfg.CD_LINES * A_lines_m2
            pilot_drag_N = D["Drag_total_N"] + Drag_lines_N
            W = calculate_wing_drag(wing_max, pilot_drag_N, CD_WING_REF, v, cfg.RHO_AIR, e_oswald=cfg.E_OSWALD)
            glide = W["glide_ratio"]
            sink = v / glide if glide > 0 else float("nan")
            rows.append({"speed_kmh": v * 3.6, "sink_mps": sink})

        df = pd.DataFrame(rows)
        color = colorsA[idx % len(colorsA)]
        ls = linestyles[idx % len(linestyles)]
        mk = markers[idx % len(markers)]
        lw = max(0.4, 1.2 - 0.1 * idx)
        ms = max(2.0, 4.0 - 0.3 * idx)
        plt.plot(
            df["speed_kmh"], df["sink_mps"],
            label=f"{label} ({wing_max}kg, 20kg gear, eq {eq_len_mm:.0f}mm)",
            color=color, linestyle=ls, linewidth=lw,
            marker=mk, markevery=6, markersize=ms
        )
        for r in rows:
            all_rows.append({
                "label": label,
                "wing_max_kg": wing_max,
                "gear_weight": 20,
                **r
            })

    # Additional gear weight groups: 30kg and 40kg if defined in config (new scheme)
    cmap_oranges = plt.get_cmap("Oranges")
    cmap_greens = plt.get_cmap("Greens")
    extra_groups = [
        (getattr(cfg, 'PILOT_EXAMPLES_30', []), 30, cmap_oranges, '30kg'),
        (getattr(cfg, 'PILOT_EXAMPLES_40', []), 40, cmap_greens, '40kg'),
    ]

    # Flag to indicate if additive scaling is active; used for legend annotation
    additive_active = getattr(cfg, 'EQUALIZER_GEAR_ADDITIVE_ENABLED', False)

    for group_list, ballast_value, cmap, ballast_label in extra_groups:
        if not group_list:
            continue
        nG = len(group_list)
        colorsG = cmap(np.linspace(0.20, 0.85, max(nG, 1)))
        for idx, (label, height_cm, weight_kg, sex, wing_max) in enumerate(group_list):
            person = create_person_from_height_weight(
                height_cm, weight_kg, sex,
                pilot_max_weight_kg=wing_max,
                back_angle_deg=cfg.BACK_ANGLE_DEG,
                upper_arm_angle_deg=cfg.UPPER_ARM_ANGLE_DEG,
                protector_height_mm=cfg.PROTECTOR_HEIGHT_MM,
                equalizer_diameter_mm=cfg.EQUALIZER_DIAMETER_MM,
            )
            areas = frontal_areas_iso(person)
            eq_len_mm = float(areas.get("equalizer_length_mm", 0.0))

            rows = []
            A_lines_m2 = lines_area_from_wing_size(wing_max)
            for v in speeds:
                D = drag_all_iso(
                    areas,
                    speed_mps=v,
                    rho_air=cfg.RHO_AIR,
                    Cd_streamlined=cfg.CD_STREAMLINED,
                    Cd_head=cfg.CD_HEAD,
                    head_exposed_fraction=cfg.HEAD_EXPOSED_FRACTION,
                    Cd_arms=cfg.CD_ARMS,
                    Cd_equalizers=cfg.CD_EQUALIZERS,
                )
                q_v = 0.5 * cfg.RHO_AIR * v**2
                Drag_lines_N = q_v * cfg.CD_LINES * A_lines_m2
                pilot_drag_N = D["Drag_total_N"] + Drag_lines_N
                W = calculate_wing_drag(wing_max, pilot_drag_N, CD_WING_REF, v, cfg.RHO_AIR, e_oswald=cfg.E_OSWALD)
                glide = W["glide_ratio"]
                sink = v / glide if glide > 0 else float("nan")
                rows.append({"speed_kmh": v * 3.6, "sink_mps": sink})

            df = pd.DataFrame(rows)
            color = colorsG[idx % len(colorsG)]
            ls = linestyles[idx % len(linestyles)]
            mk = markers[idx % len(markers)]
            # Later groups are plotted after Group A; make them thinner/smaller so they don't hide earlier lines
            lw_ex = max(0.3, 0.8 - 0.1 * idx - (0.1 if ballast_value >= 40 else 0.0))
            ms_ex = max(1.8, 3.0 - 0.2 * idx - (0.2 if ballast_value >= 40 else 0.0))
            plt.plot(
                df["speed_kmh"], df["sink_mps"],
                label=f"{label} ({wing_max}kg, {ballast_label} gear, eq {eq_len_mm:.0f}mm)",
                color=color, linestyle=ls, linewidth=lw_ex,
                marker=mk, markevery=6, markersize=ms_ex
            )
            for r in rows:
                all_rows.append({
                    "label": label,
                    "wing_max_kg": wing_max,
                    "gear_weight": ballast_value,
                    **r
                })

    ax = plt.gca()
    ax.invert_yaxis()
    eq_status = "True" if cfg.CD_EQUALIZERS >= 0.5 else "False"
    plt.title("Polars: sink vs speed (all pilots, gear weight groups)")
    plt.xlabel("Speed (km/h)")
    plt.ylabel("Sink (m/s)")
    plt.grid(True, alpha=0.3)
    legend_title = "Pilot groups" + (" (additive eq len)" if additive_active else "")
    plt.legend(loc="best", fontsize=8, title=legend_title)
    plt.text(0.98, 0.02, f"Equalizers: {eq_status}", transform=ax.transAxes,
             ha="right", va="bottom", fontsize=9,
             bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.7))
    plt.tight_layout()

    # Optional: save figure
    if getattr(cfg, "SAVE_PLOTS", False):
        try:
            os.makedirs(cfg.PLOTS_DIR, exist_ok=True)
            equ_flag = "True" if cfg.CD_EQUALIZERS >= 0.5 else "False"
            out_path = os.path.join(cfg.PLOTS_DIR, f"polars_all_pilots_equ={equ_flag}.{cfg.SAVE_FORMAT}")
            plt.savefig(out_path, dpi=cfg.SAVE_DPI, format=cfg.SAVE_FORMAT, bbox_inches="tight")
            print(f"[saved] {out_path}")
        except Exception as exc:
            print(f"[warn] failed to save figure: {exc}")
    try:
        plt.show(block=cfg.SHOW_BLOCKING)
        if not cfg.SHOW_BLOCKING:
            plt.pause(0.001)
    except Exception:
        try:
            plt.show()
        except Exception:
            pass

    if cfg.EXPORT_CSV and all_rows:
        try:
            pd.DataFrame(all_rows).to_csv("polars_all_pilots.csv", index=False)
            print("[saved] polars_all_pilots.csv")
        except Exception as exc:
            print(f"[warn] failed to save polars_all_pilots.csv: {exc}")


def plot_115kg_ballast_comparison(CD_WING_REF: float):
    """Plot a focused gear weight comparison around the M-size wing.

        Curves included:
            - 115 kg wing at 20/30/40 kg gear weight using:
          • 95 kg @ 184 cm (20 kg)
          • 85 kg @ 178 cm (30 kg)
          • 75 kg @ 173 cm (40 kg)
            - 125 kg wing at 20 kg gear using 105 kg @ 190 cm
            - For the 40 kg gear case pilot (75 kg @ 173 cm), also plot his 20 kg gear size:
                95 kg wing at 20 kg gear

    Notes:
    - Legend matches the multi-pilot plot: "Man {height}cm {weight}kg ({wing}kg, {gear}kg gear)".
    - Colors: Blues for 115 kg cases (shade by gear weight), Oranges for 125@20, Greens for 95@20.
    - Writes CSV when cfg.EXPORT_CSV is True (file: polar_115kg_gear_comparison.csv).
    """
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        print("\n[polar] matplotlib is not available. Install it to see the gear weight comparison plot:")
        print("  python -m pip install matplotlib")
        print(f"  (reason: {exc})")
        return

    cases = [
        {"wing_max": 115, "gear": 20, "pilot_weight": 95, "height_cm": 184},
        {"wing_max": 115, "gear": 30, "pilot_weight": 85, "height_cm": 178},
        {"wing_max": 115, "gear": 40, "pilot_weight": 75, "height_cm": 173},
        {"wing_max": 125, "gear": 20, "pilot_weight": 105, "height_cm": 190},
        # Same 40kg gear weight pilot (75kg/173cm) on his 20kg gear size (95kg wing)
        {"wing_max": 95,  "gear": 20, "pilot_weight": 75, "height_cm": 173},
    ]
    speeds = np.linspace(8.0, 65.0/3.6, 57)

    plt.figure(figsize=(8.5, 5.5))
    markers = ['o', 's', 'D', '^', 'v']
    all_rows = []

    for idx, c in enumerate(cases):
        wing_max = int(c["wing_max"])
        person = create_person_from_height_weight(
            c["height_cm"], c["pilot_weight"], "male",
            pilot_max_weight_kg=wing_max,
            back_angle_deg=cfg.BACK_ANGLE_DEG,
            upper_arm_angle_deg=cfg.UPPER_ARM_ANGLE_DEG,
            protector_height_mm=cfg.PROTECTOR_HEIGHT_MM,
            equalizer_diameter_mm=cfg.EQUALIZER_DIAMETER_MM,
        )
        areas = frontal_areas_iso(person)
        eq_len_mm = float(areas.get("equalizer_length_mm", 0.0))
        # Legend naming aligned with the multi-pilot plot, include eq len
        label = (
            f"Man {c['height_cm']}cm {c['pilot_weight']}kg "
            f"({wing_max}kg, {c['gear']}kg gear, eq {eq_len_mm:.0f}mm)"
        )
        rows = []
        A_lines_m2 = lines_area_from_wing_size(wing_max)
        for v in speeds:
            D = drag_all_iso(
                areas,
                speed_mps=v,
                rho_air=cfg.RHO_AIR,
                Cd_streamlined=cfg.CD_STREAMLINED,
                Cd_head=cfg.CD_HEAD,
                head_exposed_fraction=cfg.HEAD_EXPOSED_FRACTION,
                Cd_arms=cfg.CD_ARMS,
                Cd_equalizers=cfg.CD_EQUALIZERS,
            )
            q_v = 0.5 * cfg.RHO_AIR * v**2
            Drag_lines_N = q_v * cfg.CD_LINES * A_lines_m2
            pilot_drag_N = D["Drag_total_N"] + Drag_lines_N
            W = calculate_wing_drag(wing_max, pilot_drag_N, CD_WING_REF, v, cfg.RHO_AIR, e_oswald=cfg.E_OSWALD)
            glide = W["glide_ratio"]
            sink = v / glide if glide > 0 else float("nan")
            row = {
                "wing_max_kg": wing_max,
                "gear_weight_kg": c["gear"],
                "pilot_weight_kg": c["pilot_weight"],
                "height_cm": c["height_cm"],
                "speed_kmh": v * 3.6,
                "sink_mps": sink,
            }
            rows.append(row)
            all_rows.append({"label": label, **row})

        # Distinct colors per wing size; vary shade by gear weight for 115 kg wing
        if wing_max == 115:
            shade = {20: 0.25, 30: 0.55, 40: 0.85}.get(c["gear"], 0.6)
            color = plt.cm.get_cmap("Blues")(shade)
        elif wing_max == 125:
            color = plt.cm.get_cmap("Oranges")(0.7)
        elif wing_max == 95:
            color = plt.cm.get_cmap("Greens")(0.65)
        else:
            color = plt.cm.get_cmap("Greys")(0.5)

        df = pd.DataFrame(rows)
        plt.plot(
            df["speed_kmh"], df["sink_mps"],
            label=label,
            color=color, linestyle='-', linewidth=0.2,
            marker=markers[idx % len(markers)], markevery=6, markersize=3.5
        )

    ax = plt.gca()
    ax.invert_yaxis()
    eq_status = "True" if cfg.CD_EQUALIZERS >= 0.5 else "False"
    plt.title("M-size (115 kg) gear weight comparison + 125@20 + 95@20 cases")
    plt.xlabel("Speed (km/h)")
    plt.ylabel("Sink (m/s)")
    plt.grid(True, alpha=0.3)
    plt.legend(loc="best", fontsize=8)
    plt.text(0.98, 0.02, f"Equalizers: {eq_status}", transform=ax.transAxes,
             ha="right", va="bottom", fontsize=9,
             bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.7))
    plt.tight_layout()

    # Optional: save figure
    if getattr(cfg, "SAVE_PLOTS", False):
        try:
            os.makedirs(cfg.PLOTS_DIR, exist_ok=True)
            equ_flag = "True" if cfg.CD_EQUALIZERS >= 0.5 else "False"
            out_path = os.path.join(cfg.PLOTS_DIR, f"polar_115kg_gear_comparison_equ={equ_flag}.{cfg.SAVE_FORMAT}")
            plt.savefig(out_path, dpi=cfg.SAVE_DPI, format=cfg.SAVE_FORMAT, bbox_inches="tight")
            print(f"[saved] {out_path}")
        except Exception as exc:
            print(f"[warn] failed to save figure: {exc}")
    try:
        plt.show(block=cfg.SHOW_BLOCKING)
        if not cfg.SHOW_BLOCKING:
            plt.pause(0.001)
    except Exception:
        try:
            plt.show()
        except Exception:
            pass

    if cfg.EXPORT_CSV and all_rows:
        try:
            pd.DataFrame(all_rows).to_csv("polar_115kg_gear_comparison.csv", index=False)
            print("[saved] polar_115kg_gear_comparison.csv")
        except Exception as exc:
            print(f"[warn] failed to save polar_115kg_gear_comparison.csv: {exc}")
