"""Compare hosted R&D investment with observed resistance burden."""
# Copyright 2026 Tavpritesh Sethi and Jasmine Kaur (Tavlab, IIITD)
# Licensed under the Apache License, Version 2.0 (see LICENSE)

import numpy as np
import os
import pandas as pd
from common import ROOT
from scipy.stats import pearsonr, rankdata, spearmanr


def analyse_investment():
    """Calculate correlations, burden terciles and the GDP-adjusted sensitivity."""
    # hosted R&D vs resistance burden: spearman, terciles, sensitivity, partial corr | GDP

    OUT = "outputs/rd_investment"
    os.makedirs("outputs/rd_investment", exist_ok=True)

    RNG = np.random.default_rng(20260723)
    N_BOOT_SPEARMAN = 5000
    N_BOOT_PARTIAL = 10000

    def spearman_ci(x, y, nboot=N_BOOT_SPEARMAN):
        m = (~pd.isna(x)) & (~pd.isna(y))
        x, y = np.asarray(x)[m], np.asarray(y)[m]
        n = len(x)
        rho, p = spearmanr(x, y)
        boots = []
        for _ in range(nboot):
            idx = RNG.integers(0, n, n)
            if len(np.unique(x[idx])) < 3 or len(np.unique(y[idx])) < 3:
                continue
            boots.append(spearmanr(x[idx], y[idx])[0])
        lo, hi = np.nanpercentile(boots, [2.5, 97.5])
        return dict(n=int(n), rho=float(rho), ci_lo=float(lo), ci_hi=float(hi), p=float(p))

    def partial_spearman(x, y, z):
        a, b, c = rankdata(x), rankdata(y), rankdata(z)
        Z = np.column_stack([np.ones(len(c)), c])
        ra = a - Z @ np.linalg.lstsq(Z, a, rcond=None)[0]
        rb = b - Z @ np.linalg.lstsq(Z, b, rcond=None)[0]
        return pearsonr(ra, rb)

    def ols_r2(y, X):
        X = np.column_stack([np.ones(len(y))] + list(X))
        beta = np.linalg.lstsq(X, y, rcond=None)[0]
        resid = y - X @ beta
        return 1 - np.sum(resid ** 2) / np.sum((y - y.mean()) ** 2), beta

    df = pd.read_csv("outputs/dataset/rd_country_dataset.csv")
    df.to_csv(os.path.join(OUT, "investment_vs_burden_by_country.csv"), index=False)
    present = df[~df["hub_record_missing"]].copy()
    print(f"countries: {len(df)}; with a Hub record: {len(present)}; missing: "
          f"{df.loc[df.hub_record_missing, 'iso3'].tolist()}")

    # 1. Spearman correlations (order fixed: bootstrap CIs share one RNG stream)
    corr = []
    for label, x, y in [
        ("burden vs hosted USD (total)", "burden_mean_baseline_pctR", "usd_nominal_hosted"),
        ("burden vs hosted USD per capita", "burden_mean_baseline_pctR", "hosted_usd_per_capita"),
        ("burden vs share of hosted USD", "burden_mean_baseline_pctR", "share_of_global_hosted"),
        ("burden vs funded USD (funder country)", "burden_mean_baseline_pctR", "usd_nominal_funded"),
        ("share improving vs hosted USD per capita", "share_improving", "hosted_usd_per_capita"),
        ("share not improving vs hosted USD per capita", "share_not_improving", "hosted_usd_per_capita"),
        ("GDP pc vs hosted USD per capita", "gdp_pc_ppp", "hosted_usd_per_capita"),
        ("GDP pc vs burden", "gdp_pc_ppp", "burden_mean_baseline_pctR"),
    ]:
        corr.append(dict(label=label, x=x, y=y, **spearman_ci(present[x], present[y])))
    corr = pd.DataFrame(corr)
    corr.to_csv(os.path.join(OUT, "spearman_correlations.csv"), index=False)
    print(corr[["label", "n", "rho", "ci_lo", "ci_hi", "p"]].round(3).to_string(index=False))

    # 2. burden terciles
    present["burden_tercile"] = pd.qcut(present["burden_mean_baseline_pctR"], 3,
                                        labels=["low", "mid", "high"])
    ter = (present.groupby("burden_tercile", observed=True)
           .agg(n_countries=("iso3", "size"),
                mean_burden_pctR=("burden_mean_baseline_pctR", "mean"),
                hosted_usd=("usd_nominal_hosted", "sum"),
                mean_share_improving=("share_improving", "mean"),
                total_pop=("population_2020", "sum"))
           .reset_index())
    ter["pct_of_hosted_usd"] = 100 * ter["hosted_usd"] / ter["hosted_usd"].sum()
    ter["pct_of_population"] = 100 * ter["total_pop"] / ter["total_pop"].sum()
    ter["hosted_usd_per_capita"] = ter["hosted_usd"] / ter["total_pop"]
    ter.to_csv(os.path.join(OUT, "burden_terciles.csv"), index=False)
    print(ter.round(2).to_string(index=False))

    # 3. sensitivity analyses and 4. partial correlation
    sens = []

    def add(analysis, n, estimate, ci_lo=np.nan, ci_hi=np.nan, p=np.nan):
        sens.append(dict(analysis=analysis, n_countries=n, estimate=estimate,
                         ci_lo=ci_lo, ci_hi=ci_hi, p=p))

    r, p = spearmanr(present["burden_standardised"], present["hosted_usd_per_capita"])
    add("Spearman: standardised burden vs hosted USD per capita", len(present), r, p=p)
    r, p = spearmanr(present["burden_standardised"], present["burden_mean_baseline_pctR"])
    add("Spearman: standardised vs unstandardised burden (country ranking)", len(present), r, p=p)

    ex = present[~present.iso3.isin(["USA", "GBR"])]
    r, p = spearmanr(ex["burden_mean_baseline_pctR"], ex["hosted_usd_per_capita"])
    add("Spearman: burden vs hosted USD per capita, excluding USA and GBR", len(ex), r, p=p)

    z = df.copy()
    z["usd0"] = z["usd_nominal_hosted"].fillna(0.0)
    r, p = spearmanr(z["burden_mean_baseline_pctR"], z["usd0"] / z["population_2020"])
    add("Spearman: burden vs hosted USD per capita, no-record countries as zero", len(z), r, p=p)

    x = present["burden_mean_baseline_pctR"].values
    y = np.log(present["hosted_usd_per_capita"].values)
    g = np.log(present["gdp_pc_ppp"].values)
    pr, pp = partial_spearman(x, y, g)
    rng = np.random.default_rng(20260723)
    boots = []
    for _ in range(N_BOOT_PARTIAL):
        i = rng.integers(0, len(x), len(x))
        boots.append(partial_spearman(x[i], y[i], g[i])[0])
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    add("Partial Spearman: burden vs hosted USD per capita | GDP per capita",
        len(x), pr, lo, hi, pp)

    r2_gdp, _ = ols_r2(y, [g])
    r2_both, beta = ols_r2(y, [x, g])
    add("OLS R2: log hosted USD per capita ~ log GDP per capita", len(x), r2_gdp)
    add("OLS R2: log hosted USD per capita ~ burden + log GDP per capita", len(x), r2_both)

    sens = pd.DataFrame(sens)
    sens.to_csv(os.path.join(OUT, "sensitivity_and_partial.csv"), index=False)
    print(sens.round(3).to_string(index=False))


def main():
    os.chdir(ROOT)
    analyse_investment()


if __name__ == "__main__":
    main()
