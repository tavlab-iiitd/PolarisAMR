"""Draw the TrACSS audit and investment-versus-burden figures from calculated results."""
# Copyright 2026 Tavpritesh Sethi and Jasmine Kaur (Tavlab, IIITD)
# Licensed under the Apache License, Version 2.0 (see LICENSE)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
import pandas as pd
from common import ROOT


def plot_results():
    """Recreate Figures 1 and 2 using the matching analysis columns."""
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.edgecolor": "#888", "axes.linewidth": 0.8, "figure.dpi": 140,
                         "savefig.dpi": 160, "font.family": "DejaVu Sans"})

    OI = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73", "red": "#D55E00",
          "yellow": "#F0E442", "black": "#222222"}

    os.makedirs("figures", exist_ok=True)

    # figure 1 - TrACSS audit
    lab = pd.read_csv("outputs/tracss/audit_lab_items_by_wave_sector.csv", index_col=0)
    # The 1-versus-4 comparison holds in 2020-21, not in every survey wave.
    n_human = int(lab.loc["2020-21", "human"])
    n_animal = int(lab.loc["2020-21", "animal_food"])
    lad = pd.read_csv("outputs/tracss/audit_ladder_comparisons.csv").set_index("ladder_stability")["n"]
    total = int(lad["total"])

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [1, 1.25]})
    a1.bar(["Human\nclinical lab", "Animal / food\nlab battery"], [n_human, n_animal],
           color=[OI["blue"], OI["red"]], width=0.6, zorder=3)
    for x, v in zip([0, 1], [n_human, n_animal]):
        a1.text(x, v + 0.08, str(v), ha="center", fontweight="bold", fontsize=12)
    a1.set_ylabel("Gradeable A–E lab items (2020–21)")
    a1.set_ylim(0, n_animal + 0.6)
    a1.set_title("A. Sector coverage of laboratory items", loc="left", fontweight="bold", fontsize=10)

    left = 0
    for cat, label, col in [("stable", "stable", OI["green"]),
                            ("reworded_same_meaning", "reworded", OI["yellow"]),
                            ("redefined", "redefined", OI["red"])]:
        v = int(lad.get(cat, 0))
        a2.barh(0, v, left=left, color=col, zorder=3, height=0.5)
        a2.text(left + v / 2, 0, f"{label}\n{v}", ha="center", va="center", fontsize=9,
                fontweight="bold", color="#222")
        left += v
    a2.set_xlim(0, total)
    a2.set_ylim(-0.5, 0.9)
    a2.set_yticks([])
    a2.set_xlabel(f"Matched item comparisons across waves (n={total})")
    a2.set_title(f"B. A–E ladders redefined between waves: {int(lad.get('redefined', 0))}/{total}",
                 loc="left", fontweight="bold", fontsize=10)
    a2.text(0, 0.55, 'National action plan, top grade "E": implemented (2019) → '
            "funded in budget (2022) → evaluated (2025)", fontsize=8, color="#555", style="italic")
    plt.tight_layout()
    plt.savefig("figures/figure1_tracss_audit.png", bbox_inches="tight")
    plt.close()

    # figure 2 - R&D vs burden
    h = pd.read_csv("outputs/rd_investment/investment_vs_burden_by_country.csv")
    h = h[(h.hub_record_missing == False) & (h.hosted_usd_per_capita > 0)].copy()
    rho = pd.read_csv("outputs/rd_investment/spearman_correlations.csv").set_index("label")
    rho = rho.loc["burden vs hosted USD per capita", "rho"]
    ter = pd.read_csv("outputs/rd_investment/burden_terciles.csv").set_index("burden_tercile")

    fig, ax = plt.subplots(figsize=(8.2, 6))
    sc = ax.scatter(h.burden_mean_baseline_pctR, h.hosted_usd_per_capita, s=90,
                    c=h.share_not_improving, cmap="YlOrRd", vmin=0.5, vmax=1.0,
                    edgecolor="#333", linewidth=0.7, zorder=3)
    ax.set_yscale("log")
    for _, r in h.iterrows():
        if r.iso3 in ["TUR", "DNK", "GRC", "THA", "USA", "GBR", "MEX", "PHL", "IND", "ZAF", "BRA", "CHN"]:
            ax.annotate(r.iso3, (r.burden_mean_baseline_pctR, r.hosted_usd_per_capita), fontsize=8,
                        xytext=(4, 3), textcoords="offset points", fontweight="bold")
    tur = h[h.iso3 == "TUR"].iloc[0]
    dnk = h[h.iso3 == "DNK"].iloc[0]
    ax.annotate(f"Türkiye\n{tur.burden_mean_baseline_pctR:.0f}% R, "
                f"{100 * tur.share_not_improving:.0f}% not improving\n${tur.hosted_usd_per_capita:.3f}/capita",
                (tur.burden_mean_baseline_pctR, tur.hosted_usd_per_capita), fontsize=8,
                xytext=(-20, -45), textcoords="offset points", ha="right", color=OI["red"],
                fontweight="bold",
                arrowprops=dict(arrowstyle="->", color=OI["red"]))
    ax.annotate(f"Denmark\n{dnk.burden_mean_baseline_pctR:.0f}% R, ${dnk.hosted_usd_per_capita:.1f}/capita",
                (dnk.burden_mean_baseline_pctR, dnk.hosted_usd_per_capita), fontsize=8,
                xytext=(10, -40), textcoords="offset points", color=OI["blue"], fontweight="bold",
                arrowprops=dict(arrowstyle="->", color=OI["blue"]))
    cb = plt.colorbar(sc, ax=ax, pad=0.02)
    cb.set_label("Share of trajectories not improving")
    ax.set_xlabel("Measured resistance burden (mean baseline % resistant)")
    ax.set_ylabel("AMR R&D investment hosted, US$ per capita (log)")
    ax.set_title(f"Spearman ρ = {rho:.2f}; highest-burden tercile: "
                 f"{ter.loc['high', 'pct_of_population']:.1f}% of population, "
                 f"{ter.loc['high', 'pct_of_hosted_usd']:.1f}% of hosted R&D",
                 fontsize=10, loc="left")
    plt.tight_layout()
    plt.savefig("figures/figure2_investment_vs_burden.png", bbox_inches="tight")
    plt.close()


def main():
    os.chdir(ROOT)
    plot_results()


if __name__ == "__main__":
    main()
