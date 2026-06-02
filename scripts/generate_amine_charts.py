from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = ROOT / "outputs" / "evaluation_exploration_csv"
OUTPUT_DIR = ROOT / "outputs" / "evaluation_exploration_charts"

CLASS_ORDER = [
    "primary_amine",
    "secondary_amine",
    "tertiary_amine",
    "cyclic_amine",
    "diamine",
    "polyamine",
    "amino_silane",
    "silylamine",
    "boron_amide",
    "metal_amide",
    "inorganic_nitrogen_source",
]

PFAS_ORDER = ["no", "possible", "yes"]
PFAS_COLORS = {"no": "#2E7D32", "possible": "#F9A825", "yes": "#C62828"}
FLUOR_COLORS = {"no": "#4E79A7", "yes": "#59A14F"}
PHASE_COLORS = {
    "gas": "#4E79A7",
    "gas_or_supercritical": "#76B7B2",
    "liquid": "#59A14F",
    "solid": "#E15759",
    "unknown": "#9D9D9D",
}

CLASS_COLORS = {
    "primary_amine": "#4E79A7",
    "secondary_amine": "#F28E2B",
    "tertiary_amine": "#59A14F",
    "cyclic_amine": "#E15759",
    "diamine": "#76B7B2",
    "polyamine": "#EDC948",
    "amino_silane": "#B07AA1",
    "silylamine": "#FF9DA7",
    "boron_amide": "#9C755F",
    "metal_amide": "#BAB0AC",
    "inorganic_nitrogen_source": "#8CD17D",
}


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    summary = pd.read_csv(INPUT_DIR / "summary.csv")
    amines = pd.read_csv(INPUT_DIR / "amine_summary.csv")
    class_summary = pd.read_csv(INPUT_DIR / "amine_class_summary.csv")
    breakdown = pd.read_csv(INPUT_DIR / "candidate_breakdown.csv")

    summary = _prepare_summary(summary)
    amines = _prepare(amines)
    class_summary = _sort_class_summary(class_summary)

    index_rows = []
    _save(index_rows, "01_amine_class_fluorination.png", "Amine class counts split by fluorinated candidate status.", plot_class_fluorination, class_summary)
    _save(index_rows, "02_amine_class_pfas.png", "PFAS screening result by amine class.", plot_class_pfas, class_summary)
    _save(index_rows, "03_property_coverage.png", "Availability of key physical properties in the amine summary.", plot_property_coverage, amines)
    _save(index_rows, "04_boiling_point_by_class.png", "Boiling point distribution by amine class.", plot_boiling_point_by_class, amines)
    _save(index_rows, "05_vapor_pressure_25C_by_class.png", "Estimated or curated vapor pressure at 25C by amine class, log scale.", plot_vapor_pressure_by_class, amines)
    _save(index_rows, "06_phase_supply_matrix.png", "Phase versus supply class distribution for amines.", plot_phase_supply_matrix, amines)
    _save(index_rows, "07_mw_vs_boiling_pfas.png", "Molecular weight versus boiling point, colored by PFAS flag.", plot_mw_vs_boiling, amines)
    _save(index_rows, "08_precursor_family_counts.png", "Curated semiconductor amine precursor families excluding template-only unknown.", plot_precursor_family_counts, amines)
    _save(index_rows, "09_fluorination_pfas_matrix.png", "Fluorination level versus PFAS screening result.", plot_fluorination_pfas_matrix, amines)
    _save(index_rows, "10_process_role_counts.png", "Semiconductor process role tags on curated amine candidates.", plot_process_role_counts, amines)
    _save(index_rows, "11_property_scatter_matrix_by_class.png", "Scatter-matrix view of key amine properties colored by amine class.", plot_property_scatter_matrix, amines)
    _save(index_rows, "12_tb_vs_pvap25_by_class.png", "Boiling point versus 25C vapor pressure colored by amine class.", plot_tb_vs_pvap25, amines)
    _save(index_rows, "13_mw_vs_pvap25_by_class.png", "Molecular weight versus 25C vapor pressure colored by amine class.", plot_mw_vs_pvap25, amines)
    _save(index_rows, "14_tc_vs_pc_by_class.png", "Critical temperature versus critical pressure colored by amine class.", plot_tc_vs_pc, amines)
    _save(index_rows, "15_tm_vs_tb_by_class.png", "Melting point versus boiling point colored by amine class.", plot_tm_vs_tb, amines)
    _save(index_rows, "16_pvap25_vs_pvap60_by_class.png", "Vapor pressure at 25C versus 60C colored by amine class.", plot_pvap25_vs_pvap60, amines)
    _save(index_rows, "17_gwp_availability_by_family.png", "GWP100 AR6 availability rate by candidate family.", plot_gwp_availability_by_family, summary)
    _save(index_rows, "18_gwp100_ar6_known_values.png", "Known local GWP100 AR6 values merged into the screening output.", plot_gwp_known_values, summary)
    _save(index_rows, "19_gwp_by_pfas_flag.png", "Known GWP100 AR6 values grouped by PFAS screening flag.", plot_gwp_by_pfas_flag, summary)
    _save(index_rows, "20_gwp_vs_molecular_weight.png", "Known GWP100 AR6 values versus molecular weight, colored by PFAS screen.", plot_gwp_vs_molecular_weight, summary)
    _save(index_rows, "21_mw_vs_melting_pfas.png", "Molecular weight versus melting point, colored by PFAS flag.", plot_mw_vs_melting_pfas, amines)
    _save(index_rows, "22_mw_vs_critical_temperature_pfas.png", "Molecular weight versus critical temperature, colored by PFAS flag.", plot_mw_vs_critical_temperature_pfas, amines)
    _save(index_rows, "23_mw_vs_critical_pressure_pfas.png", "Molecular weight versus critical pressure, colored by PFAS flag.", plot_mw_vs_critical_pressure_pfas, amines)
    _save(index_rows, "24_mw_vs_pvap25_pfas.png", "Molecular weight versus 25C vapor pressure, colored by PFAS flag.", plot_mw_vs_pvap25_pfas, amines)
    _save(index_rows, "25_mw_vs_pvap40_pfas.png", "Molecular weight versus 40C vapor pressure, colored by PFAS flag.", plot_mw_vs_pvap40_pfas, amines)
    _save(index_rows, "26_mw_vs_pvap60_pfas.png", "Molecular weight versus 60C vapor pressure, colored by PFAS flag.", plot_mw_vs_pvap60_pfas, amines)

    with (OUTPUT_DIR / "plots_index.csv").open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["file", "purpose"])
        writer.writeheader()
        writer.writerows(index_rows)

    return 0


def _prepare(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in [
        "molecular_weight",
        "tm_C",
        "tb_C",
        "tc_C",
        "pc_MPa",
        "pvap_25C_kPa",
        "pvap_40C_kPa",
        "pvap_60C_kPa",
        "gwp100_ar6",
    ]:
        if col in out:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    out["amine_class_label"] = out["amine_class_label"].fillna("unknown_amine")
    out["pfas_flag"] = out["pfas_flag"].fillna("unknown")
    out["fluorinated_amine"] = out["fluorinated_amine"].fillna("unknown")
    out["amine_fluorination_level"] = out["amine_fluorination_level"].fillna("unknown")
    return out


def _prepare_summary(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in ["gwp100_ar6", "molecular_weight", "tb_C", "pvap_25C_kPa"]:
        if col in out:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    out["candidate_family"] = out["candidate_family"].fillna("unknown")
    out["preferred_name"] = out["preferred_name"].fillna("unknown")
    out["gwp100_ar6_status"] = out.get("gwp100_ar6_status", pd.Series(["missing"] * len(out))).fillna("missing")
    return out


def _sort_class_summary(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["amine_class_label"] = out["amine_class_label"].fillna("unknown_amine")
    order = {label: idx for idx, label in enumerate(CLASS_ORDER)}
    out["_order"] = out["amine_class_label"].map(order).fillna(999)
    return out.sort_values(["_order", "amine_class_label"]).drop(columns=["_order"])


def _save(index_rows: list[dict], filename: str, purpose: str, fn, *args) -> None:
    fig = fn(*args)
    fig.savefig(OUTPUT_DIR / filename, dpi=180, bbox_inches="tight")
    plt.close(fig)
    index_rows.append({"file": filename, "purpose": purpose})


def _setup(title: str, width: float = 11.0, height: float = 6.4):
    fig, ax = plt.subplots(figsize=(width, height))
    ax.set_title(title, fontsize=14, pad=12)
    ax.grid(axis="x", alpha=0.22)
    return fig, ax


def _footnote(fig):
    return None


def _ordered_classes(df: pd.DataFrame) -> list[str]:
    present = set(df["amine_class_label"].dropna())
    ordered = [label for label in CLASS_ORDER if label in present]
    ordered.extend(sorted(present - set(ordered)))
    return ordered


def _scatter_by_class(
    ax,
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    x_label: str,
    y_label: str,
    xlog: bool = False,
    ylog: bool = False,
    legend: bool = True,
) -> None:
    data = df[[x_col, y_col, "amine_class_label"]].dropna().copy()
    if xlog:
        data = data[data[x_col] > 0]
    if ylog:
        data = data[data[y_col] > 0]
    for label in _ordered_classes(data):
        subset = data[data["amine_class_label"] == label]
        if subset.empty:
            continue
        ax.scatter(
            subset[x_col],
            subset[y_col],
            s=24,
            alpha=0.62,
            color=CLASS_COLORS.get(label, "#666666"),
            label=label,
            edgecolors="none",
        )
    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.grid(alpha=0.22)
    if legend:
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8, frameon=False)


def plot_class_fluorination(class_summary: pd.DataFrame):
    fig, ax = _setup("Amine candidate count by class and fluorination")
    labels = class_summary["amine_class_label"].tolist()
    no = class_summary["non_fluorinated_count"].to_numpy()
    yes = class_summary["fluorinated_count"].to_numpy()
    y = np.arange(len(labels))
    ax.barh(y, no, color=FLUOR_COLORS["no"], label="non_fluorinated")
    ax.barh(y, yes, left=no, color=FLUOR_COLORS["yes"], label="fluorinated")
    for idx, total in enumerate(class_summary["count"].to_numpy()):
        ax.text(total + max(class_summary["count"]) * 0.01, idx, str(int(total)), va="center", fontsize=9)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("candidate count")
    ax.legend(loc="lower right")
    _footnote(fig)
    return fig


def plot_class_pfas(class_summary: pd.DataFrame):
    fig, ax = _setup("PFAS screen by amine class")
    labels = class_summary["amine_class_label"].tolist()
    y = np.arange(len(labels))
    left = np.zeros(len(labels))
    columns = {"no": "pfas_no_count", "possible": "pfas_possible_count", "yes": "pfas_yes_count"}
    for flag in PFAS_ORDER:
        values = class_summary[columns[flag]].to_numpy()
        ax.barh(y, values, left=left, color=PFAS_COLORS[flag], label=f"pfas_{flag}")
        left += values
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("candidate count")
    ax.legend(loc="lower right")
    _footnote(fig)
    return fig


def plot_property_coverage(amines: pd.DataFrame):
    properties = ["tm_C", "tb_C", "tc_C", "pc_MPa", "pvap_25C_kPa", "pvap_40C_kPa", "pvap_60C_kPa", "gwp100_ar6"]
    counts = [int(amines[prop].notna().sum()) for prop in properties]
    total = len(amines)
    fig, ax = _setup("Physical property coverage for amine candidates", width=10.5, height=5.8)
    colors = ["#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#76B7B2", "#EDC948", "#B07AA1"]
    ax.bar(properties, counts, color=colors)
    for idx, count in enumerate(counts):
        ax.text(idx, count + total * 0.01, f"{count}/{total}", ha="center", fontsize=9)
    ax.set_ylabel("present value count")
    ax.set_ylim(0, total * 1.10)
    ax.tick_params(axis="x", rotation=25)
    _footnote(fig)
    return fig


def plot_boiling_point_by_class(amines: pd.DataFrame):
    return _boxplot_by_class(
        amines,
        "tb_C",
        "Boiling point distribution by amine class",
        "normal boiling point (C)",
        linear=True,
    )


def plot_vapor_pressure_by_class(amines: pd.DataFrame):
    data = amines[amines["pvap_25C_kPa"] > 0].copy()
    fig = _boxplot_by_class(
        data,
        "pvap_25C_kPa",
        "Vapor pressure at 25C by amine class",
        "vapor pressure at 25C (kPa, log scale)",
        linear=False,
    )
    return fig


def _boxplot_by_class(amines: pd.DataFrame, column: str, title: str, ylabel: str, linear: bool):
    labels = [label for label in CLASS_ORDER if label in set(amines["amine_class_label"])]
    labels += sorted(set(amines["amine_class_label"]) - set(labels))
    groups = [amines.loc[amines["amine_class_label"] == label, column].dropna().to_numpy() for label in labels]
    labels = [label for label, values in zip(labels, groups) if len(values)]
    groups = [values for values in groups if len(values)]
    fig, ax = _setup(title, width=12, height=6.5)
    box = ax.boxplot(groups, tick_labels=labels, patch_artist=True, showfliers=False)
    palette = ["#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#76B7B2", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC", "#8CD17D"]
    for patch, color in zip(box["boxes"], palette):
        patch.set_facecolor(color)
        patch.set_alpha(0.72)
    rng = np.random.default_rng(7)
    for idx, values in enumerate(groups, start=1):
        if len(values) > 120:
            values = rng.choice(values, size=120, replace=False)
        jitter = rng.normal(idx, 0.045, size=len(values))
        ax.scatter(jitter, values, s=12, alpha=0.30, color="#333333", linewidths=0)
    if not linear:
        ax.set_yscale("log")
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", rotation=35)
    _footnote(fig)
    return fig


def plot_phase_supply_matrix(amines: pd.DataFrame):
    phase_order = ["gas", "gas_or_supercritical", "liquid", "solid", "unknown"]
    supply_order = [
        "compressed_or_liquefied_gas",
        "bubbler_or_direct_liquid_source",
        "heated_source_likely",
        "heated_source_required",
        "solid_or_sublimation_source_review",
        "unknown_review_required",
    ]
    table = pd.crosstab(amines["phase_25C_1atm"], amines["supply_class"]).reindex(index=phase_order, columns=supply_order, fill_value=0)
    fig, ax = plt.subplots(figsize=(12.5, 5.8))
    ax.set_title("Phase at 25C versus supply class", fontsize=14, pad=12)
    image = ax.imshow(table.to_numpy(), cmap="YlGnBu", aspect="auto")
    ax.set_xticks(np.arange(len(supply_order)), supply_order, rotation=35, ha="right")
    ax.set_yticks(np.arange(len(phase_order)), phase_order)
    for i in range(table.shape[0]):
        for j in range(table.shape[1]):
            value = int(table.iat[i, j])
            if value:
                ax.text(j, i, str(value), ha="center", va="center", fontsize=9, color="#111111")
    fig.colorbar(image, ax=ax, label="candidate count")
    _footnote(fig)
    return fig


def plot_mw_vs_boiling(amines: pd.DataFrame):
    fig, ax = _setup("Molecular weight versus boiling point", width=10.5, height=6.5)
    for flag in PFAS_ORDER:
        subset = amines[(amines["pfas_flag"] == flag) & amines["molecular_weight"].notna() & amines["tb_C"].notna()]
        ax.scatter(
            subset["molecular_weight"],
            subset["tb_C"],
            s=28,
            alpha=0.65,
            label=f"pfas_{flag}",
            color=PFAS_COLORS[flag],
            edgecolors="none",
        )
    ax.set_xlabel("molecular weight")
    ax.set_ylabel("normal boiling point (C)")
    ax.legend(loc="upper left")
    _footnote(fig)
    return fig


def _scatter_mw_vs_property_pfas(
    amines: pd.DataFrame,
    y_col: str,
    title: str,
    y_label: str,
    ylog: bool = False,
):
    data = amines[["molecular_weight", y_col, "pfas_flag"]].dropna().copy()
    if ylog:
        data = data[data[y_col] > 0]
    fig, ax = _setup(title, width=10.5, height=6.5)
    for flag in PFAS_ORDER:
        subset = data[data["pfas_flag"] == flag]
        if subset.empty:
            continue
        ax.scatter(
            subset["molecular_weight"],
            subset[y_col],
            s=28,
            alpha=0.65,
            label=f"pfas_{flag}",
            color=PFAS_COLORS[flag],
            edgecolors="none",
        )
    if ylog:
        ax.set_yscale("log")
    ax.set_xlabel("molecular weight")
    ax.set_ylabel(y_label)
    ax.legend(loc="upper left")
    _footnote(fig)
    return fig


def plot_mw_vs_melting_pfas(amines: pd.DataFrame):
    return _scatter_mw_vs_property_pfas(
        amines,
        "tm_C",
        "Molecular weight versus melting point",
        "normal melting point (C)",
    )


def plot_mw_vs_critical_temperature_pfas(amines: pd.DataFrame):
    return _scatter_mw_vs_property_pfas(
        amines,
        "tc_C",
        "Molecular weight versus critical temperature",
        "critical temperature (C)",
    )


def plot_mw_vs_critical_pressure_pfas(amines: pd.DataFrame):
    return _scatter_mw_vs_property_pfas(
        amines,
        "pc_MPa",
        "Molecular weight versus critical pressure",
        "critical pressure (MPa)",
    )


def plot_mw_vs_pvap25_pfas(amines: pd.DataFrame):
    return _scatter_mw_vs_property_pfas(
        amines,
        "pvap_25C_kPa",
        "Molecular weight versus vapor pressure at 25C",
        "vapor pressure at 25C (kPa, log scale)",
        ylog=True,
    )


def plot_mw_vs_pvap40_pfas(amines: pd.DataFrame):
    return _scatter_mw_vs_property_pfas(
        amines,
        "pvap_40C_kPa",
        "Molecular weight versus vapor pressure at 40C",
        "vapor pressure at 40C (kPa, log scale)",
        ylog=True,
    )


def plot_mw_vs_pvap60_pfas(amines: pd.DataFrame):
    return _scatter_mw_vs_property_pfas(
        amines,
        "pvap_60C_kPa",
        "Molecular weight versus vapor pressure at 60C",
        "vapor pressure at 60C (kPa, log scale)",
        ylog=True,
    )


def plot_precursor_family_counts(amines: pd.DataFrame):
    series = amines["precursor_family"].replace("", np.nan).dropna()
    series = series[series != "unknown"]
    counts = series.value_counts().sort_values()
    fig, ax = _setup("Curated semiconductor amine precursor families", width=11.5, height=5.8)
    ax.barh(counts.index, counts.values, color="#4E79A7")
    for idx, value in enumerate(counts.values):
        ax.text(value + max(counts.values) * 0.03, idx, str(int(value)), va="center", fontsize=9)
    ax.set_xlabel("candidate count")
    _footnote(fig)
    return fig


def plot_fluorination_pfas_matrix(amines: pd.DataFrame):
    fluor_order = ["none", "single_fluorinated_substituent", "multiple_fluorinated_substituents"]
    table = pd.crosstab(amines["amine_fluorination_level"], amines["pfas_flag"]).reindex(index=fluor_order, columns=PFAS_ORDER, fill_value=0)
    fig, ax = _setup("Fluorination level versus PFAS screen", width=10.5, height=5.8)
    y = np.arange(len(fluor_order))
    left = np.zeros(len(fluor_order))
    for flag in PFAS_ORDER:
        values = table[flag].to_numpy()
        ax.barh(y, values, left=left, color=PFAS_COLORS[flag], label=f"pfas_{flag}")
        left += values
    ax.set_yticks(y, fluor_order)
    ax.invert_yaxis()
    ax.set_xlabel("candidate count")
    ax.legend(loc="lower right")
    _footnote(fig)
    return fig


def plot_process_role_counts(amines: pd.DataFrame):
    roles: Counter[str] = Counter()
    for value in amines["semiconductor_process_roles"].fillna(""):
        for role in [item.strip() for item in str(value).split(";") if item.strip()]:
            roles[role] += 1
    counts = pd.Series(roles).sort_values()
    fig, ax = _setup("Semiconductor process role tags", width=11.5, height=6.2)
    ax.barh(counts.index, counts.values, color="#F28E2B")
    if not counts.empty:
        for idx, value in enumerate(counts.values):
            ax.text(value + max(counts.values) * 0.03, idx, str(int(value)), va="center", fontsize=9)
    ax.set_xlabel("tag count")
    _footnote(fig)
    return fig


def plot_property_scatter_matrix(amines: pd.DataFrame):
    data = amines.copy()
    data["log10_pvap_25C_kPa"] = np.where(data["pvap_25C_kPa"] > 0, np.log10(data["pvap_25C_kPa"]), np.nan)
    variables = [
        ("molecular_weight", "MW"),
        ("tm_C", "Tm_C"),
        ("tb_C", "Tb_C"),
        ("tc_C", "Tc_C"),
        ("pc_MPa", "Pc_MPa"),
        ("log10_pvap_25C_kPa", "log10 Pvap25"),
    ]
    n = len(variables)
    fig, axes = plt.subplots(n, n, figsize=(15, 15), sharex="col", sharey="row")
    fig.suptitle("Amine property scatter matrix by class", fontsize=16, y=0.995)
    classes = _ordered_classes(data)
    for row, (y_col, y_label) in enumerate(variables):
        for col, (x_col, x_label) in enumerate(variables):
            ax = axes[row, col]
            if row == col:
                ax.text(0.5, 0.5, y_label, ha="center", va="center", fontsize=10, transform=ax.transAxes)
                ax.set_xticks([])
                ax.set_yticks([])
            elif row > col:
                subset = data[[x_col, y_col, "amine_class_label"]].dropna()
                for label in classes:
                    points = subset[subset["amine_class_label"] == label]
                    if points.empty:
                        continue
                    ax.scatter(
                        points[x_col],
                        points[y_col],
                        s=8,
                        alpha=0.42,
                        color=CLASS_COLORS.get(label, "#666666"),
                        edgecolors="none",
                    )
                ax.grid(alpha=0.16)
            else:
                ax.axis("off")
            if row == n - 1:
                ax.set_xlabel(x_label, fontsize=8)
            if col == 0:
                ax.set_ylabel(y_label, fontsize=8)
            ax.tick_params(labelsize=7)
    handles = [
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=CLASS_COLORS.get(label, "#666666"), label=label, markersize=6)
        for label in classes
    ]
    fig.legend(handles=handles, loc="center right", bbox_to_anchor=(1.09, 0.5), frameon=False, fontsize=8)
    fig.tight_layout(rect=[0, 0, 0.91, 0.98])
    return fig


def plot_tb_vs_pvap25(amines: pd.DataFrame):
    fig, ax = _setup("Boiling point vs vapor pressure at 25C by amine class", width=10.8, height=6.4)
    _scatter_by_class(
        ax,
        amines,
        "tb_C",
        "pvap_25C_kPa",
        "normal boiling point (C)",
        "vapor pressure at 25C (kPa, log scale)",
        ylog=True,
    )
    return fig


def plot_mw_vs_pvap25(amines: pd.DataFrame):
    fig, ax = _setup("Molecular weight vs vapor pressure at 25C by amine class", width=10.8, height=6.4)
    _scatter_by_class(
        ax,
        amines,
        "molecular_weight",
        "pvap_25C_kPa",
        "molecular weight",
        "vapor pressure at 25C (kPa, log scale)",
        ylog=True,
    )
    return fig


def plot_tc_vs_pc(amines: pd.DataFrame):
    fig, ax = _setup("Critical temperature vs critical pressure by amine class", width=10.8, height=6.4)
    _scatter_by_class(
        ax,
        amines,
        "tc_C",
        "pc_MPa",
        "critical temperature (C)",
        "critical pressure (MPa)",
    )
    return fig


def plot_tm_vs_tb(amines: pd.DataFrame):
    fig, ax = _setup("Melting point vs boiling point by amine class", width=10.8, height=6.4)
    _scatter_by_class(
        ax,
        amines,
        "tm_C",
        "tb_C",
        "normal melting point (C)",
        "normal boiling point (C)",
    )
    return fig


def plot_pvap25_vs_pvap60(amines: pd.DataFrame):
    fig, ax = _setup("Vapor pressure at 25C vs 60C by amine class", width=10.8, height=6.4)
    _scatter_by_class(
        ax,
        amines,
        "pvap_25C_kPa",
        "pvap_60C_kPa",
        "vapor pressure at 25C (kPa, log scale)",
        "vapor pressure at 60C (kPa, log scale)",
        xlog=True,
        ylog=True,
    )
    return fig


def plot_gwp_availability_by_family(summary: pd.DataFrame):
    data = summary.copy()
    data["gwp_available"] = np.where(data["gwp100_ar6"].notna(), "available", "missing")
    table = pd.crosstab(data["candidate_family"], data["gwp_available"])
    for col in ["available", "missing"]:
        if col not in table:
            table[col] = 0
    table = table[["available", "missing"]]
    table["total"] = table["available"] + table["missing"]
    table["available_pct"] = np.where(table["total"] > 0, table["available"] / table["total"] * 100.0, 0.0)
    table["missing_pct"] = 100.0 - table["available_pct"]
    table = table.sort_values(["available_pct", "total"], ascending=[True, True])

    fig, ax = _setup("GWP100 AR6 availability rate by candidate family", width=10.8, height=5.8)
    y = np.arange(len(table.index))
    ax.barh(y, table["available_pct"], color="#4E79A7", label="available")
    ax.barh(y, table["missing_pct"], left=table["available_pct"], color="#BAB0AC", label="missing")
    for idx, row in enumerate(table.itertuples()):
        ax.text(
            min(99, row.available_pct + 1.2),
            idx,
            f"{int(row.available)}/{int(row.total)} ({row.available_pct:.1f}%)",
            va="center",
            fontsize=9,
        )
    ax.set_yticks(y, table.index)
    ax.set_xlim(0, 100)
    ax.set_xlabel("share of candidates (%)")
    ax.legend(loc="lower right")
    return fig


def plot_gwp_known_values(summary: pd.DataFrame):
    data = summary[summary["gwp100_ar6"].notna()].copy()
    data = data.sort_values("gwp100_ar6")
    fig, ax = _setup("Known local GWP100 AR6 values", width=11.2, height=6.2)
    labels = data["preferred_name"].astype(str).tolist()
    values = data["gwp100_ar6"].to_numpy()
    colors = ["#C62828" if value >= 10000 else "#F28E2B" if value >= 1000 else "#4E79A7" for value in values]
    y = np.arange(len(labels))
    ax.barh(y, values, color=colors)
    for idx, value in enumerate(values):
        ax.text(value * 1.02 if value > 0 else 1, idx, f"{value:g}", va="center", fontsize=9)
    ax.set_yticks(y, labels)
    ax.set_xlabel("GWP100 AR6 (kg CO2e / kg)")
    if len(values) and values.max() / max(values.min(), 1) > 100:
        ax.set_xscale("log")
    return fig


def plot_gwp_by_pfas_flag(summary: pd.DataFrame):
    data = summary[summary["gwp100_ar6"].notna()].copy()
    present_flags = [flag for flag in PFAS_ORDER if flag in set(data["pfas_flag"])]
    positions = {flag: idx for idx, flag in enumerate(present_flags)}
    fig, ax = _setup("Known GWP100 AR6 by PFAS screen", width=10.8, height=6.4)
    rng = np.random.default_rng(11)
    for flag in present_flags:
        subset = data[data["pfas_flag"] == flag].sort_values("gwp100_ar6")
        x = np.full(len(subset), positions[flag], dtype=float) + rng.normal(0, 0.035, len(subset))
        ax.scatter(
            x,
            subset["gwp100_ar6"],
            s=64,
            alpha=0.80,
            color=PFAS_COLORS[flag],
            label=f"pfas_{flag}",
            edgecolors="white",
            linewidths=0.7,
        )
        for x_value, row in zip(x, subset.itertuples()):
            ax.annotate(str(row.preferred_name), (x_value, row.gwp100_ar6), xytext=(5, 2), textcoords="offset points", fontsize=8)
    ax.set_xticks(list(positions.values()), [f"pfas_{flag}" for flag in present_flags])
    ax.set_yscale("log")
    ax.set_xlabel("PFAS screening flag")
    ax.set_ylabel("GWP100 AR6 (kg CO2e / kg, log scale)")
    ax.legend(loc="lower right")
    return fig


def plot_gwp_vs_molecular_weight(summary: pd.DataFrame):
    data = summary[summary["gwp100_ar6"].notna() & summary["molecular_weight"].notna()].copy()
    fig, ax = _setup("GWP100 AR6 versus molecular weight", width=10.8, height=6.4)
    for flag in PFAS_ORDER:
        subset = data[data["pfas_flag"] == flag]
        if subset.empty:
            continue
        ax.scatter(
            subset["molecular_weight"],
            subset["gwp100_ar6"],
            s=58,
            alpha=0.78,
            color=PFAS_COLORS[flag],
            label=f"pfas_{flag}",
            edgecolors="white",
            linewidths=0.7,
        )
    for row in data.itertuples():
        ax.annotate(str(row.preferred_name), (row.molecular_weight, row.gwp100_ar6), xytext=(5, 2), textcoords="offset points", fontsize=8)
    ax.set_yscale("log")
    ax.set_xlabel("molecular weight")
    ax.set_ylabel("GWP100 AR6 (kg CO2e / kg, log scale)")
    ax.legend(loc="lower right")
    return fig


if __name__ == "__main__":
    raise SystemExit(main())
