"""Screen readiness measures, compare organism groups and check testing-volume trends."""
# Copyright 2026 Tavpritesh Sethi and Jasmine Kaur (Tavlab, IIITD)
# Licensed under the Apache License, Version 2.0 (see LICENSE)

import numpy as np
import os
import pandas as pd
import time
import warnings
from common import ROOT
from scipy.stats import chi2, norm, pearsonr, spearmanr


ABX_CLASS = {
    "Amikacin": "Aminoglycoside",
    "Amoxycillin.clavulanate": "Penicillin_BLI",
    "Ampicillin": "Penicillin",
    "Aztreonam": "Monobactam",
    "Cefepime": "Cephalosporin",
    "Ceftaroline": "Cephalosporin",
    "Ceftazidime": "Cephalosporin",
    "Ceftazidime.avibactam": "Cephalosporin_BLI",
    "Imipenem": "Carbapenem",
    "Meropenem": "Carbapenem",
    "Levofloxacin": "Fluoroquinolone",
    "Piperacillin.tazobactam": "Penicillin_BLI",
}

COVARIATES = ["cov_log_gdp_pc_ppp", "cov_uhc_index", "cov_log_atlas_testing_volume"]

def firth_fit(X, y, max_iter=500, tol=1e-7):
    n, p = X.shape
    beta = np.zeros(p)
    for _ in range(max_iter):
        eta = X @ beta
        pr = 1.0 / (1.0 + np.exp(-eta))
        pr = np.clip(pr, 1e-12, 1 - 1e-12)
        W = pr * (1 - pr)
        info = X.T @ (X * W[:, None])
        try:
            info_inv = np.linalg.inv(info)
        except np.linalg.LinAlgError:
            return None, None, False
        h = W * np.einsum("ij,ij->i", X @ info_inv, X)       # hat-matrix diagonal
        U = X.T @ (y - pr + h * (0.5 - pr))                   # Firth-modified score
        try:
            step = np.linalg.solve(info, U)
        except np.linalg.LinAlgError:
            return None, None, False
        beta = beta + step
        if not np.all(np.isfinite(beta)):
            return None, None, False
        if np.max(np.abs(step)) < tol:
            return beta, info_inv, True
    return beta, info_inv, True

def firth_fit_constrained(X, y, fix_idx, max_iter=500, tol=1e-7):
    n, p = X.shape
    beta = np.zeros(p)
    free = [j for j in range(p) if j != fix_idx]
    for _ in range(max_iter):
        eta = X @ beta
        pr = 1.0 / (1.0 + np.exp(-eta))
        pr = np.clip(pr, 1e-12, 1 - 1e-12)
        W = pr * (1 - pr)
        info = X.T @ (X * W[:, None])
        try:
            info_inv = np.linalg.inv(info)
        except np.linalg.LinAlgError:
            return None, False
        h = W * np.einsum("ij,ij->i", X @ info_inv, X)
        U = X.T @ (y - pr + h * (0.5 - pr))
        try:
            step = np.linalg.solve(info[np.ix_(free, free)], U[free])
        except np.linalg.LinAlgError:
            return None, False
        beta[free] += step
        if not np.all(np.isfinite(beta)):
            return None, False
        if np.max(np.abs(step)) < tol:
            return beta, True
    return beta, True

def pen_loglik(X, y, beta):
    """Penalised log-likelihood: log L + 0.5 log|I(beta)|."""
    eta = X @ beta
    pr = 1.0 / (1.0 + np.exp(-eta))
    pr = np.clip(pr, 1e-12, 1 - 1e-12)
    ll = np.sum(y * np.log(pr) + (1 - y) * np.log(1 - pr))
    W = pr * (1 - pr)
    _, logdet = np.linalg.slogdet(X.T @ (X * W[:, None]))
    return ll + 0.5 * logdet

def country_standardise(sub, expo):
    cl = sub.groupby("iso3")[expo].first()
    return cl.mean(), cl.std(ddof=0)

def design_matrix(sub, expo, mu, sd, covs, extra=None):
    if sd == 0 or not np.isfinite(sd):
        return None
    z = (sub[expo] - mu) / sd
    parts = [pd.Series(1.0, index=sub.index, name="const"), z.rename("z")]
    for c in covs:
        parts.append(sub[c].rename(c))
    if extra:
        for nm, s in extra.items():
            parts.append(s.rename(nm))
    parts.append(pd.get_dummies(sub["abx_class"], prefix="abx", drop_first=True).astype(float))
    Xdf = pd.concat(parts, axis=1)
    return Xdf, sub["improving"].values.astype(float), sub["iso3"].values

def cluster_bootstrap(X, y, iso, target_idx, always_keep, n_boot, seed):
    countries = np.unique(iso)
    k = len(countries)
    idx_by_c = {c: np.where(iso == c)[0] for c in countries}
    rng = np.random.default_rng(seed)
    coefs, n_fail = [], 0
    ak = set(always_keep)
    for _ in range(n_boot):
        pick = rng.choice(countries, size=k, replace=True)
        rows = np.concatenate([idx_by_c[c] for c in pick])
        Xb, yb = X[rows], y[rows]
        keep = sorted(ak | {j for j in range(X.shape[1]) if j not in ak and np.ptp(Xb[:, j]) > 0})
        ti = keep.index(target_idx)
        beta, _, ok = firth_fit(Xb[:, keep], yb)
        if (not ok) or beta is None or not np.isfinite(beta[ti]):
            n_fail += 1
            continue
        coefs.append(beta[ti])
    return np.array(coefs), n_fail

def percentile_ci(coefs):
    if len(coefs) == 0:
        return np.nan, np.nan
    return float(np.exp(np.percentile(coefs, 2.5))), float(np.exp(np.percentile(coefs, 97.5)))

def bootstrap_p(coefs):
    b = len(coefs)
    if b == 0:
        return np.nan
    p = 2 * min(np.mean(coefs <= 0), np.mean(coefs >= 0))
    return float(max(min(p, 1.0), 1.0 / (b + 1)))

def bh_fdr(pvals):
    """Benjamini-Hochberg adjusted q-values (same order as input)."""
    pv = np.asarray(pvals, dtype=float)
    m = len(pv)
    order = np.argsort(pv)
    ranked = pv[order]
    q = np.empty(m)
    prev = 1.0
    for rank in range(m - 1, -1, -1):
        prev = min(prev, ranked[rank] * m / (rank + 1))
        q[rank] = prev
    out = np.empty(m)
    out[order] = q
    return out


def screen_measures():
    """Fit 46 adjusted Firth models with 2,000 country bootstrap replicates."""
    # screen of the 46 readiness measures
    # improving ~ z(measure) + log GDP pc + UHC + log testing volume + abx class, Firth logistic
    # CIs from 2000 country bootstrap resamples, BH across the 46

    # Antibiotic -> mechanism class, used as a fixed-effect covariate.

    warnings.filterwarnings("ignore")

    OUT = "outputs/readiness_screen"
    os.makedirs("outputs/readiness_screen", exist_ok=True)
    N_BOOT = 2000
    SEED = 20260722
    UNDETERMINED_SECTOR_ITEMS = {"ITEM_009", "ITEM_054_s2"}

    def build_exposure_set(reg):
        sel = reg[(reg.role == "primary") & (reg.n_countries >= 30)].copy()
        item_of = sel.name.str.replace("exp_item_", "", regex=False)
        sel = sel[~(sel.name.str.startswith("exp_item_") & item_of.isin(UNDETERMINED_SECTOR_ITEMS))]
        # exp_lab_diagnostic_2020_21 duplicates exp_item_ITEM_063; keep one copy only
        sel = sel[sel.name != "exp_lab_diagnostic_2020_21"]
        return sel.reset_index(drop=True)

    def load_data():
        df = pd.read_csv("outputs/dataset/analysis_dataset.csv")
        df["abx_class"] = df["antibiotic"].map(ABX_CLASS)
        assert df["abx_class"].isna().sum() == 0, "antibiotic without a class"
        reg = pd.read_csv("outputs/dataset/exposure_registry.csv")
        return df, reg

    def fit_measure(df, expo, seed):
        sub = df.dropna(subset=[expo] + COVARIATES).copy()
        mu, sd = country_standardise(sub, expo)
        built = design_matrix(sub, expo, mu, sd, COVARIATES)
        if built is None:
            return dict(status="FAILED_design_zero_variance", n_rows=0), None
        Xdf, y, iso = built
        X = Xdf.values.astype(float)
        zi = list(Xdf.columns).index("z")

        beta, info_inv, ok = firth_fit(X, y)
        if (not ok) or beta is None:
            return dict(status="FAILED_nonconvergence", n_rows=len(sub)), None
        logOR, se = beta[zi], np.sqrt(info_inv[zi, zi])

        bc, okc = firth_fit_constrained(X, y, zi)
        if okc and bc is not None:
            stat = 2 * (pen_loglik(X, y, beta) - pen_loglik(X, y, bc))
            plr_p = float(chi2.sf(max(stat, 0.0), 1))
        else:
            plr_p = np.nan

        coefs, n_fail = cluster_bootstrap(X, y, iso, zi, [0, zi], N_BOOT, seed)
        ci_low, ci_high = percentile_ci(coefs)
        res = dict(N=sub["iso3"].nunique(), status="ok", OR_per_SD=float(np.exp(logOR)),
                   CI_low=ci_low, CI_high=ci_high, p=plr_p, logOR=logOR, se=se,
                   wald_p=float(2 * norm.sf(abs(logOR / se))), n_rows=len(sub),
                   n_boot_success=int(len(coefs)), n_boot_fail=int(n_fail),
                   direction="positive" if logOR > 0 else "negative",
                   boot_p=bootstrap_p(coefs))
        return res, coefs

    t0 = time.time()
    df, reg = load_data()
    exposures = build_exposure_set(reg)
    print(f"screened measures: {len(exposures)}")
    assert len(exposures) == 46

    rows, draws = [], {}
    for i, srow in exposures.iterrows():
        expo = srow["name"]
        res, coefs = fit_measure(df, expo, SEED + i)
        rows.append(dict(exposure=expo, sector=srow.sector, level=srow.level, **res))
        if coefs is not None:
            draws[expo] = coefs
        if res["status"] == "ok":
            print(f"[{i + 1}/{len(exposures)}] {expo}: OR={res['OR_per_SD']:.3f} "
                  f"CI=[{res['CI_low']:.3f},{res['CI_high']:.3f}] p={res['p']:.3g} "
                  f"({time.time() - t0:.0f}s)", flush=True)

    res = pd.DataFrame(rows)
    ok = res["p"].notna()
    res.loc[ok, "q"] = bh_fdr(res.loc[ok, "p"].values)
    res["survives_q05"] = res["q"] < 0.05
    res["survives_q10"] = res["q"] < 0.10
    okb = res["boot_p"].notna()
    res.loc[okb, "boot_q"] = bh_fdr(res.loc[okb, "boot_p"].values)
    res["boot_CI_excludes_1"] = (res["CI_low"] > 1) | (res["CI_high"] < 1)
    res["boot_survives_q05"] = res["boot_q"] < 0.05
    res["boot_survives_q10"] = res["boot_q"] < 0.10
    res = res.sort_values("q", na_position="last").reset_index(drop=True)

    cols = ["exposure", "sector", "level", "N", "OR_per_SD", "CI_low", "CI_high", "p", "q",
            "survives_q05", "survives_q10", "direction", "logOR", "se", "wald_p",
            "n_rows", "n_boot_success", "n_boot_fail", "status",
            "boot_p", "boot_q", "boot_CI_excludes_1", "boot_survives_q05", "boot_survives_q10"]
    res[cols].to_csv(os.path.join(OUT, "screen_results.csv"), index=False)
    np.savez_compressed(os.path.join(OUT, "bootstrap_logOR_draws.npz"), **draws)

    surv = res[res.boot_survives_q05]
    print(f"\nmeasures with boot_q < 0.05: {len(surv)}")
    print(surv[["exposure", "sector", "OR_per_SD", "CI_low", "CI_high", "boot_q"]].to_string(index=False))
    print(f"total {time.time() - t0:.0f}s")


def compare_organism_groups():
    """Fit ITEM_058 by organism group and estimate its interaction."""

    # Antibiotic -> mechanism class, used as a fixed-effect covariate.

    warnings.filterwarnings("ignore")

    OUT = "outputs/organism_stratified"
    os.makedirs("outputs/organism_stratified", exist_ok=True)
    N_BOOT = 2000
    SEED = 20260722
    EXPOSURE = "exp_item_ITEM_058"
    ZOONOTIC = {"Escherichia coli", "Klebsiella pneumoniae"}

    def fit_point(X, y, idx):
        beta, info_inv, ok = firth_fit(X, y)
        if (not ok) or beta is None or not np.isfinite(beta[idx]):
            return dict(status="FAILED_nonconvergence", OR=np.nan, logOR=np.nan, se=np.nan, p=np.nan)
        logOR, se = float(beta[idx]), float(np.sqrt(info_inv[idx, idx]))
        bc, okc = firth_fit_constrained(X, y, idx)
        p = (float(chi2.sf(max(2 * (pen_loglik(X, y, beta) - pen_loglik(X, y, bc)), 0.0), 1))
             if (okc and bc is not None) else np.nan)
        return dict(status="ok", OR=float(np.exp(logOR)), logOR=logOR, se=se, p=p,
                    wald_p=float(2 * norm.sf(abs(logOR / se))))

    df = pd.read_csv("outputs/dataset/analysis_dataset.csv")
    df["abx_class"] = df["antibiotic"].map(ABX_CLASS)
    df["grp"] = df["organism"].apply(lambda o: "zoonotic" if o in ZOONOTIC else "hospital")
    assert set(df["organism"]) <= (ZOONOTIC | {"Pseudomonas aeruginosa"})

    cc = df.dropna(subset=[EXPOSURE] + COVARIATES).copy()
    mu, sd = country_standardise(cc, EXPOSURE)
    rows = []

    # (a), (b) stratum-specific models
    for gi, grp in enumerate(["zoonotic", "hospital"]):
        sub = cc[cc["grp"] == grp].copy()
        Xdf, y, iso = design_matrix(sub, EXPOSURE, mu, sd, COVARIATES)
        X = Xdf.values.astype(float)
        zi = list(Xdf.columns).index("z")
        pt = fit_point(X, y, zi)
        coefs, nf = cluster_bootstrap(X, y, iso, zi, [0, zi], N_BOOT, SEED + gi)
        lo, hi = percentile_ci(coefs)
        rows.append(dict(model=grp, organisms=("E. coli + K. pneumoniae" if grp == "zoonotic"
                                               else "P. aeruginosa"),
                         n_countries=sub["iso3"].nunique(), n_rows=len(sub),
                         OR_per_SD=pt["OR"], CI_low=lo, CI_high=hi, p=pt["p"],
                         n_boot_success=len(coefs), n_boot_fail=nf))

    # (c) pooled model with interaction (reference = E. coli + K. pneumoniae)
    g_hosp = (cc["grp"] == "hospital").astype(float)
    z = (cc[EXPOSURE] - mu) / sd
    Xdf, y, iso = design_matrix(cc, EXPOSURE, mu, sd, COVARIATES,
                                extra={"grp_hospital": g_hosp, "z_x_hospital": z * g_hosp})
    names = list(Xdf.columns)
    X = Xdf.values.astype(float)
    ti, zi = names.index("z_x_hospital"), names.index("z")
    pt = fit_point(X, y, ti)
    coefs, nf = cluster_bootstrap(X, y, iso, ti, [0, zi, names.index("grp_hospital"), ti],
                                  N_BOOT, SEED + 5)
    lo, hi = percentile_ci(coefs)
    rows.append(dict(model="interaction", organisms="ratio of ORs: P. aeruginosa vs E. coli + K. pneumoniae",
                     n_countries=cc["iso3"].nunique(), n_rows=len(cc),
                     OR_per_SD=pt["OR"], CI_low=lo, CI_high=hi, p=pt["p"],
                     n_boot_success=len(coefs), n_boot_fail=nf))

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(OUT, "item058_by_organism_group.csv"), index=False)
    print(out.round(3).to_string(index=False))


def check_testing_trend():
    """Relate the share improving to changes in ATLAS testing volume."""
    # share improving vs trend in testing volume (country level)

    OUT = "outputs/testing_volume"
    os.makedirs("outputs/testing_volume", exist_ok=True)

    def testing_trend():
        cy = pd.read_csv("outputs/atlas/atlas_testing_volume_country_year.csv")
        w = cy[(cy.year >= 2014) & (cy.year <= 2022)]
        rows = []
        for iso, g in w.groupby("iso3"):
            g = g.sort_values("year")
            if len(g) < 3:
                rows.append(dict(iso3=iso, n_years_trend=len(g), slope_raw=np.nan, slope_log=np.nan))
                continue
            x, y = g.year.values.astype(float), g.isolates.values.astype(float)
            rows.append(dict(iso3=iso, n_years_trend=len(g),
                             slope_raw=float(np.polyfit(x, y, 1)[0]),
                             slope_log=float(np.polyfit(x, np.log(y), 1)[0]) if (y > 0).all() else np.nan))
        return pd.DataFrame(rows)

    df = pd.read_csv("outputs/dataset/analysis_dataset.csv")
    trend = testing_trend()
    cty = (df.groupby("iso3")
           .agg(share_improving=("improving", "mean"), n_trajectories=("improving", "size"),
                log_vol=("cov_log_atlas_testing_volume", "first"))
           .reset_index()
           .merge(trend, on="iso3", how="left"))
    cty.to_csv(os.path.join(OUT, "testing_volume_trend_by_country.csv"), index=False)

    rows = []
    for label, col in [("testing_volume_trend_raw_slope", "slope_raw"),
                       ("testing_volume_trend_log_slope", "slope_log"),
                       ("testing_volume_level_log", "log_vol")]:
        sub = cty.dropna(subset=[col, "share_improving"])
        pr, pp = pearsonr(sub["share_improving"], sub[col])
        sr, sp = spearmanr(sub["share_improving"], sub[col])
        rows.append(dict(y="share_improving", x=label, n_countries=len(sub),
                         pearson_r=pr, pearson_p=pp, spearman_rho=sr, spearman_p=sp))
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(OUT, "improving_vs_testing_volume.csv"), index=False)
    print(out.round(3).to_string(index=False))


def main():
    os.chdir(ROOT)
    screen_measures()
    compare_organism_groups()
    check_testing_trend()


if __name__ == "__main__":
    main()
