"""Select the analysis cohort and join readiness, surveillance, investment and covariates."""
# Copyright 2026 Tavpritesh Sethi and Jasmine Kaur (Tavlab, IIITD)
# Licensed under the Apache License, Version 2.0 (see LICENSE)

import csv
import datetime as dt
import hashlib
import json
import math
import numpy as np
import os
import pandas as pd
import statistics
from collections import Counter, OrderedDict, defaultdict
from common import ATLAS_ALIASES, DROP_BY_RULE, NAME_TO_ISO3, ROOT, WINDOW_LABELS, normalise_country
from datetime import datetime, timezone


def select_combinations():
    """Retain combinations meeting the country and isolate-count thresholds."""

    OUT = "outputs/trajectories"
    os.makedirs("outputs/trajectories", exist_ok=True)
    MIN_ISOLATES = 30
    MIN_WINDOWS = 2
    N_GRID = [5, 10, 15, 20, 25, 30, 35]
    N_CHOSEN = 20

    traj = defaultdict(set)
    with open("outputs/trajectories/trajectory_countries.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            traj[(r["organism"], r["source"], r["antibiotic"])].add(r["country"])

    with open("outputs/atlas/density_table.csv", newline="") as fh:
        dens = {(r["organism"], r["source"], r["antibiotic"]): r for r in csv.DictReader(fh)}
    assert len(dens) == 362, f"density table has {len(dens)} rows, expected 362"

    good = defaultdict(set)
    cwr = defaultdict(list)
    with open("outputs/atlas/density_country_window.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            raw = r["isolates"]
            if not raw.strip().isdigit():
                raise AssertionError(f"unparseable isolate count {raw!r}")
            cwr[(r["organism"], r["source"], r["antibiotic"], r["country"])].append(int(raw))
    for (o, s, a, c), counts in cwr.items():
        if sum(1 for n in counts if n >= MIN_ISOLATES) >= MIN_WINDOWS:
            good[(o, s, a)].add(c)

    qual_all = {k: len(good.get(k, set())) for k in dens}
    qual_traj = {k: len(good.get(k, set()) & traj[k]) for k in dens}

    sweep = []
    for N in N_GRID:
        ka = [k for k in dens if qual_all[k] >= N]
        kt = [k for k in dens if qual_traj[k] >= N]
        sweep.append(dict(
            N=N,
            retained_ALL=len(ka), dropped_ALL=362 - len(ka),
            organisms_ALL=len(set(k[0] for k in ka)), antibiotics_ALL=len(set(k[2] for k in ka)),
            retained_TRAJ=len(kt), dropped_TRAJ=362 - len(kt),
            organisms_TRAJ=len(set(k[0] for k in kt)), antibiotics_TRAJ=len(set(k[2] for k in kt)),
        ))

    rows = []
    for k in sorted(dens):
        d = dens[k]
        qa, qt = qual_all[k], qual_traj[k]
        retained = qt >= N_CHOSEN
        if retained:
            reason = "retained"
        elif int(d["n_countries_any_data"]) == 0:
            reason = "dropped: no ATLAS isolates for this organism-source-antibiotic"
        elif qa == 0:
            reason = f"dropped: no country has >={MIN_WINDOWS} windows with >={MIN_ISOLATES} isolates"
        elif qa >= N_CHOSEN and qt < N_CHOSEN:
            reason = (f"dropped: {qa} countries meet density but only {qt} of them are in "
                      f"the trajectory country set (need >={N_CHOSEN})")
        else:
            reason = (f"dropped: only {qt} trajectory countries have >={MIN_WINDOWS} windows "
                      f"with >={MIN_ISOLATES} isolates (need >={N_CHOSEN})")
        rows.append(dict(
            organism=k[0], source=k[1], antibiotic=k[2],
            retained=int(retained), drop_reason=reason,
            n_qualifying_countries_TRAJ=qt, n_qualifying_countries_ALL=qa,
            n_countries_any_data=d["n_countries_any_data"],
            n_countries_ge2_valid_windows=d["n_countries_ge2_valid_windows"],
            n_trajectory_countries=d["n_trajectory_countries"],
            isolates_per_cw_median=d["isolates_per_cw_median"],
            isolates_per_cw_q1=d["isolates_per_cw_q1"],
            isolates_per_cw_q3=d["isolates_per_cw_q3"],
            pct_cw_below_10=d["pct_cw_below_10"],
            pct_cw_below_30=d["pct_cw_below_30"],
            year_min=d["year_min"], year_max=d["year_max"],
            n_years_with_data=d["n_years_with_data"],
            rule=(f"retain iff >= {N_CHOSEN} trajectory countries have >= {MIN_WINDOWS} "
                  f"rolling windows with >= {MIN_ISOLATES} valid-AST isolates"),
        ))

    path = os.path.join(OUT, "basket_frozen.csv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(OUT, "basket_threshold_sweep.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(sweep[0].keys()))
        w.writeheader()
        w.writerows(sweep)

    kept = [r for r in rows if r["retained"]]
    manifest = dict(
        artefact="basket_frozen.csv",
        sha256=hashlib.sha256(open(path, "rb").read()).hexdigest(),
        created_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
        rule=rows[0]["rule"], N=N_CHOSEN, min_isolates=MIN_ISOLATES, min_windows=MIN_WINDOWS,
        total_combinations=len(rows), retained=len(kept), dropped=len(rows) - len(kept),
        retained_by_organism={o: sum(1 for r in kept if r["organism"] == o)
                              for o in sorted(set(r["organism"] for r in kept))},
        sweep=sweep,
    )
    with open(os.path.join(OUT, "basket_frozen_manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"retained {len(kept)} of {len(rows)} combinations (N={N_CHOSEN})")


def classify_trajectories():
    """Compare the first and last available AMROrbit windows."""
    # Improving = level down AND velocity down; Turning = velocity down only

    OUT = "outputs/trajectories"
    os.makedirs("outputs/trajectories", exist_ok=True)

    CANON_FIRST, CANON_LAST = WINDOW_LABELS[0], WINDOW_LABELS[-1]

    # this one only matches global.csv up to floating point noise (not in final set anyway)
    NOT_EXACTLY_VERIFIED = {("Escherichia coli", "Blood", "Tigecycline")}

    def wkey(w):
        a, b = w.split("-")
        return (int(a), int(b))

    def parse_float(raw, ctx):
        if raw is None or raw.strip() == "":
            raise AssertionError(f"empty numeric value: {ctx}")
        try:
            return float(raw)
        except ValueError:
            raise AssertionError(f"unparseable float {raw!r}: {ctx}")

    retained = set()
    with open(os.path.join(OUT, "basket_frozen.csv"), newline="") as fh:
        for r in csv.DictReader(fh):
            if r["retained"] == "1":
                retained.add((r["organism"], r["source"], r["antibiotic"]))

    series = defaultdict(dict)   # (org, src, ab, country) -> {window: values}
    verifiable = set()
    with open(os.path.join(OUT, "reconstructed_windows.csv"), newline="") as fh:
        for r in csv.DictReader(fh):
            org = r["organism_global"]
            if org == "UNMAPPED":
                continue
            key = (org, r["source"], r["antibiotic"])
            verifiable.add(key)
            series[key + (r["country"],)][r["window"]] = dict(
                amp=parse_float(r["amplitude"], f"{key} {r['country']} {r['window']} amp"),
                vel=parse_float(r["velocity"], f"{key} {r['country']} {r['window']} vel"),
                prov_file=r["prov_file"], prov_row_index=r["prov_row_index"],
                prov_raw_intercept=r["prov_raw_intercept"], prov_raw_slope=r["prov_raw_slope"])

    analysis_combos = retained & verifiable
    rows = []
    excluded_lt2, partial_span = defaultdict(int), defaultdict(int)
    class_counts = defaultdict(int)

    for (org, src, ab, country), wd in series.items():
        key = (org, src, ab)
        if key not in analysis_combos:
            continue
        real_ws = sorted(wd.keys(), key=wkey)
        if len(real_ws) < 2:
            excluded_lt2[key] += 1
            continue
        fw, lw = real_ws[0], real_ws[-1]
        first, last = wd[fw], wd[lw]
        d_level = last["amp"] - first["amp"]
        d_velocity = last["vel"] - first["vel"]
        if d_level < 0 and d_velocity < 0:
            cls = "Improving"
        elif d_velocity < 0 and d_level >= 0:
            cls = "Turning"
        else:
            cls = "Not improving"
        class_counts[cls] += 1
        spans_grid = (fw == CANON_FIRST and lw == CANON_LAST)
        if not spans_grid:
            partial_span[key] += 1
        rows.append(dict(
            organism=org, source=src, antibiotic=ab, country=country,
            n_real_windows=len(real_ws), first_window=fw, last_window=lw,
            spans_canonical_grid=int(spans_grid),
            amplitude_first=repr(first["amp"]), amplitude_last=repr(last["amp"]),
            velocity_first=repr(first["vel"]), velocity_last=repr(last["vel"]),
            d_level=repr(d_level), d_velocity=repr(d_velocity),
            outcome_class=cls, improving=int(cls == "Improving"),
            combo_not_exactly_verified=int(key in NOT_EXACTLY_VERIFIED),
            prov_first_file=first["prov_file"], prov_first_row=first["prov_row_index"],
            prov_first_raw_intercept=first["prov_raw_intercept"],
            prov_first_raw_slope=first["prov_raw_slope"],
            prov_last_file=last["prov_file"], prov_last_row=last["prov_row_index"],
            prov_last_raw_intercept=last["prov_raw_intercept"],
            prov_last_raw_slope=last["prov_raw_slope"],
        ))

    rows.sort(key=lambda r: (r["organism"], r["source"], r["antibiotic"], r["country"]))
    with open(os.path.join(OUT, "outcome_within_country.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    excl_rows = []
    for key in sorted(analysis_combos):
        excl_rows.append(dict(
            organism=key[0], source=key[1], antibiotic=key[2],
            country_trajectories_with_outcome=sum(
                1 for r in rows if (r["organism"], r["source"], r["antibiotic"]) == key),
            excluded_lt2_real_windows=excluded_lt2.get(key, 0),
            included_but_partial_span=partial_span.get(key, 0),
            not_exactly_verified=int(key in NOT_EXACTLY_VERIFIED),
        ))
    with open(os.path.join(OUT, "outcome_exclusions_by_combo.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(excl_rows[0].keys()))
        w.writeheader()
        w.writerows(excl_rows)

    summary = dict(
        retained_combinations=len(retained),
        verifiable_reconstructed_combos=len(verifiable),
        analysis_combos_retained_and_verifiable=len(analysis_combos),
        country_combination_trajectories=len(rows),
        excluded_country_trajectories_lt2_windows=sum(excluded_lt2.values()),
        included_trajectories_partial_grid_span=sum(partial_span.values()),
        outcome_distribution=dict(class_counts),
        outcome_distribution_pct={k: round(100.0 * v / len(rows), 1) for k, v in class_counts.items()},
    )
    with open(os.path.join(OUT, "outcome_within_country_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))


def exclude_near_zero_resistance():
    """Apply the original 2% resistance and near-zero amplitude rules."""

    csv.field_size_limit(10 ** 7)

    OUT = "outputs/trajectories"
    os.makedirs("outputs/trajectories", exist_ok=True)
    DUST = 1e-6
    FLOORS = [1.0, 2.0, 5.0]
    CHOSEN_FLOOR = 2.0

    def analysis_combos():
        verifiable = set()
        with open(os.path.join(OUT, "reconstructed_windows.csv"), newline="") as fh:
            for r in csv.DictReader(fh):
                if r["organism_global"] != "UNMAPPED":
                    verifiable.add((r["organism_global"], r["source"], r["antibiotic"]))
        retained = set()
        with open(os.path.join(OUT, "basket_frozen.csv"), newline="") as fh:
            for r in csv.DictReader(fh):
                if r["retained"] == "1":
                    retained.add((r["organism"], r["source"], r["antibiotic"]))
        return retained & verifiable

    def amplitude_metrics(combos):
        amps = defaultdict(list)
        with open(os.path.join(OUT, "reconstructed_windows.csv"), newline="") as fh:
            for r in csv.DictReader(fh):
                key = (r["organism_global"], r["source"], r["antibiotic"])
                if key in combos:
                    amps[key].append(float(r["amplitude"]))
        out = {}
        for k, vs in amps.items():
            av = [abs(x) for x in vs]
            out[k] = dict(
                n_country_windows=len(vs),
                median_amplitude=statistics.median(vs),
                median_abs_amplitude=statistics.median(av),
                max_abs_amplitude=max(av),
                frac_amplitude_below_dust=round(sum(1 for x in av if x < DUST) / len(av), 4),
            )
        return out

    combos = analysis_combos()
    amp = amplitude_metrics(combos)
    res = {}
    with open("outputs/atlas/resistance_by_combo.csv", newline="") as fh:
        for row in csv.DictReader(fh):
            key = (row["organism"], row["source"], row["antibiotic"])
            res[key] = {c: (float(row[c]) if row[c] else None)
                        for c in ["median_pctR", "mean_pctR", "max_pctR"]}
            res[key]["n_country_windows_atlas"] = int(row["n_country_windows_atlas"])

    rows = []
    for k in sorted(combos):
        a, r = amp[k], res[k]
        rows.append(dict(
            organism=k[0], source=k[1], antibiotic=k[2],
            median_pctR_atlas=("" if r["median_pctR"] is None else round(r["median_pctR"], 3)),
            mean_pctR_atlas=r["mean_pctR"], max_pctR_atlas=r["max_pctR"],
            n_country_windows_atlas=r["n_country_windows_atlas"],
            median_amplitude=round(a["median_amplitude"], 6),
            median_abs_amplitude=round(a["median_abs_amplitude"], 6),
            max_abs_amplitude=round(a["max_abs_amplitude"], 6),
            frac_amplitude_below_dust=a["frac_amplitude_below_dust"],
            n_country_windows_model=a["n_country_windows"],
        ))

    def is_dust(row):
        return row["frac_amplitude_below_dust"] >= 0.80 or row["max_abs_amplitude"] < DUST

    def is_low(row, floor):
        mr = row["median_pctR_atlas"]
        return mr != "" and mr < floor

    sweep = []
    for floor in FLOORS:
        dropped = [(r["organism"], r["source"], r["antibiotic"]) for r in rows
                   if is_low(r, floor) or is_dust(r)]
        sweep.append(dict(floor_pctR=floor, dropped=len(dropped), retained=len(rows) - len(dropped)))

    dropped_final = []
    for row in rows:
        low, dust = is_low(row, CHOSEN_FLOOR), is_dust(row)
        row["degenerate"] = int(low or dust)
        reasons = []
        if low:
            reasons.append(f"median %R {row['median_pctR_atlas']} < {CHOSEN_FLOOR}% floor")
        if dust:
            reasons.append(f"amplitude ~0 (frac<1e-6={row['frac_amplitude_below_dust']}, "
                           f"max|amp|={row['max_abs_amplitude']})")
        row["degeneracy_reason"] = "; ".join(reasons)
        if low or dust:
            dropped_final.append((row["organism"], row["source"], row["antibiotic"]))

    with open(os.path.join(OUT, "degeneracy_screen.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    summary = dict(
        analysis_combos_in=len(rows), chosen_floor_pctR=CHOSEN_FLOOR,
        dust_rule="share(|amp| < 1e-6) >= 0.80 OR max|amp| < 1e-6",
        dropped_degenerate=len(dropped_final),
        surviving_non_degenerate=len(rows) - len(dropped_final),
        dropped_list=[list(d) for d in sorted(dropped_final)],
        floor_sweep=sweep,
    )
    with open(os.path.join(OUT, "degeneracy_screen_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"{len(rows)} combinations screened; {len(dropped_final)} excluded; "
          f"{len(rows) - len(dropped_final)} retained")


def save_final_combinations():
    """Save the 85 combinations used in the submission."""

    OUT = "outputs/trajectories"
    os.makedirs("outputs/trajectories", exist_ok=True)

    KEY = ["organism", "source", "antibiotic"]

    deg = pd.read_csv(os.path.join(OUT, "degeneracy_screen.csv"))
    oc = pd.read_csv(os.path.join(OUT, "outcome_within_country.csv"))

    keep = deg[deg.degenerate == 0][KEY + ["median_pctR_atlas"]]
    oc = oc.merge(keep[KEY], on=KEY)

    per = (oc.groupby(KEY)
           .agg(n_country_trajectories=("country", "size"),
                n_improving=("outcome_class", lambda s: int((s == "Improving").sum())),
                n_turning=("outcome_class", lambda s: int((s == "Turning").sum())),
                n_not_improving=("outcome_class", lambda s: int((s == "Not improving").sum())),
                combo_not_exactly_verified=("combo_not_exactly_verified", "max"))
           .reset_index())
    per["pct_improving"] = (100 * per.n_improving / per.n_country_trajectories).round(1)
    per = per.merge(keep, on=KEY, how="left")
    per = per[KEY + ["n_country_trajectories", "n_improving", "n_turning", "n_not_improving",
                     "pct_improving", "median_pctR_atlas", "combo_not_exactly_verified"]]
    per = per.sort_values(KEY).reset_index(drop=True)

    path = os.path.join(OUT, "analysis_set_final.csv")
    per.to_csv(path, index=False)

    n_retained = int(pd.read_csv(os.path.join(OUT, "basket_frozen.csv")).retained.sum())
    summary = {
        "analysis_set_final_combos": int(len(per)),
        "funnel": {"retained_on_density": n_retained,
                   "retained_and_verifiable": int(len(deg)),
                   "minus_degenerate_2pct": int(len(per))},
        "country_combination_trajectories_final": int(len(oc)),
        "outcome_distribution_final": oc.outcome_class.value_counts().to_dict(),
        "surviving_by_organism": per.organism.value_counts().to_dict(),
        "surviving_by_source": per.source.value_counts().to_dict(),
        "sha256": hashlib.sha256(open(path, "rb").read()).hexdigest(),
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    with open(os.path.join(OUT, "analysis_set_final_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))


def merge_readiness():
    """Join the trajectory rows to TrACSS and World Bank covariates."""
    # TrACSS reference wave = 2020-21 (or latest earlier wave)

    OUT = "outputs/dataset"
    os.makedirs("outputs/dataset", exist_ok=True)

    def check(cond, msg):
        if not cond:
            raise AssertionError(msg)

    NAME_MAP = dict(NAME_TO_ISO3)

    NAME_MAP.update(ATLAS_ALIASES)

    def to_iso3(raw):
        n = normalise_country(raw)
        if n == "" or n.lower() == "nan":
            return None, n, "blank"
        if n in DROP_BY_RULE:
            return None, n, "drop_by_rule"
        i = NAME_MAP.get(n)
        return (i, n, "mapped") if i else (None, n, "UNMAPPED")

    PRE_WAVES = ["2016-17", "2017-18", "2018-19", "2019-20", "2020-21"]
    BOUNDARY_WAVE = "2020-21"

    POST_WAVES = ["2022", "2023", "2024", "2025"]   # post-date the trajectories: descriptive only

    # 1. outcome rows for the final combinations
    final_combos = set()

    for r in csv.DictReader(open("outputs/trajectories/analysis_set_final.csv")):
        final_combos.add((r["organism"], r["source"], r["antibiotic"]))

    check(len(final_combos) == 85, f"expected 85 final combinations, got {len(final_combos)}")

    outcome_rows, outcome_iso = [], set()

    n_dropped_non_member = 0

    for r in csv.DictReader(open("outputs/trajectories/outcome_within_country.csv")):
        if (r["organism"], r["source"], r["antibiotic"]) not in final_combos:
            continue
        iso3, name, status = to_iso3(r["country"])
        if status == "drop_by_rule":          # Taiwan / Hong Kong: not in TrACSS
            n_dropped_non_member += 1
            continue
        check(status == "mapped", f"outcome country not in lookup: {r['country']!r}")
        outcome_iso.add(iso3)
        outcome_rows.append({
            "organism": r["organism"], "source": r["source"], "antibiotic": r["antibiotic"],
            "country": name, "iso3": iso3,
            "n_real_windows": r["n_real_windows"], "first_window": r["first_window"],
            "last_window": r["last_window"], "d_level": r["d_level"], "d_velocity": r["d_velocity"],
            "outcome_class": r["outcome_class"], "improving": r["improving"],
            "combo_not_exactly_verified": r["combo_not_exactly_verified"],
        })

    analysis_iso = sorted(outcome_iso)

    print(f"outcome rows={len(outcome_rows)}  countries={len(analysis_iso)}  "
          f"non-member rows dropped={n_dropped_non_member}")

    # 2. item-level exposures
    grade, item_waves, item_sector = {}, defaultdict(set), {}

    for r in csv.DictReader(open("outputs/tracss/tracss_analysis_long.csv")):
        if r["parse_status"] != "parsed":
            continue
        g = r["grade_ordinal"]
        if g in ("", "nan"):
            continue
        iso, item, wave = r["iso3"], r["item_id"], r["wave"]
        grade[(iso, item, wave)] = int(g)
        item_waves[item].add(wave)
        item_sector[item] = r["sector"]

    # item-level ladder verdict = worst verdict across the item's waves
    _LADDER_RANK = {"stable": 0, "reworded_same_meaning": 1, "redefined": 2}

    _cw_rows = defaultdict(list)

    for r in csv.DictReader(open("outputs/tracss/tracss_item_crosswalk.csv")):
        _cw_rows[r["item_id"]].append(r)

    cw_meta = {}

    for iid, rs in _cw_rows.items():
        cw_meta[iid] = {
            "construct_label": rs[0]["construct_label"],
            "ladder_stability": max(rs, key=lambda r: _LADDER_RANK.get(r["ladder_stability"], 0))["ladder_stability"],
            "longitudinal_eligible_primary_0.70": rs[0]["longitudinal_eligible_primary_0.70"],
        }

    def reference_wave(item):
        w = item_waves[item]
        pre = [x for x in PRE_WAVES if x in w]
        if pre:
            return pre[-1], "predictive"
        post = [x for x in POST_WAVES if x in w]
        if post:
            return post[0], "descriptive"
        return None, "none"

    item_specs, item_value = [], {}

    for item in sorted(item_waves):
        ref_wave, kind = reference_wave(item)
        if ref_wave is None:
            continue
        meta = cw_meta.get(item, {})
        ladder = meta.get("ladder_stability", "")
        if kind == "descriptive":
            role, prefix, xsec = "descriptive", "desc_item_", False
        elif ladder in ("stable", "reworded_same_meaning"):
            role, prefix, xsec = "primary", "exp_item_", False
        else:
            role, prefix, xsec = "cross_sectional", "exp_item_", True
        col = prefix + item
        item_specs.append(OrderedDict(
            name=col, sector=item_sector.get(item, "undetermined"),
            construct_label=meta.get("construct_label", ""), reference_wave=ref_wave,
            boundary_wave=(ref_wave == BOUNDARY_WAVE), ladder_stability=ladder,
            longitudinal_eligible=meta.get("longitudinal_eligible_primary_0.70", "") == "True",
            cross_sectional_only=xsec, role=role))
        for iso in analysis_iso:
            v = grade.get((iso, item, ref_wave))
            if v is not None:
                item_value[(col, iso)] = v

    # 3. domain exposures at 2020-21, total readiness and lean contrasts
    DOMAIN_ITEMS = {
        "one_health_coordination": ["ITEM_036", "ITEM_007"],
        "nap_maturity": ["ITEM_005"],
        "human_training": ["ITEM_004", "ITEM_065"],
        "human_am_use_monitoring": ["ITEM_003", "ITEM_066", "ITEM_074"],
        "human_amr_surveillance": ["ITEM_000"],
        "human_ipc": ["ITEM_001"],
    }

    DOMAINS = list(DOMAIN_ITEMS)

    domain_item_2021 = {}

    for dom, items in DOMAIN_ITEMS.items():
        at = [i for i in items if BOUNDARY_WAVE in item_waves.get(i, set())]
        check(len(at) == 1, f"domain {dom}: expected one item at {BOUNDARY_WAVE}, got {at}")
        domain_item_2021[dom] = at[0]

    domain_value, domain_sector = defaultdict(dict), {}

    for dom in DOMAINS:
        it = domain_item_2021[dom]
        domain_sector[dom] = item_sector.get(it, "undetermined")
        for iso in analysis_iso:
            v = grade.get((iso, it, BOUNDARY_WAVE))
            if v is not None:
                domain_value[dom][iso] = v

    total_readiness, lean_value = {}, defaultdict(dict)

    for iso in analysis_iso:
        vals = {d: domain_value[d].get(iso) for d in DOMAINS}
        if all(v is not None for v in vals.values()):
            total_readiness[iso] = sum(vals.values()) / len(DOMAINS)
            for d in DOMAINS:
                others = [vals[o] for o in DOMAINS if o != d]
                lean_value[d][iso] = vals[d] - sum(others) / len(others)

    # 4. human lab composite (2020-21) and single diagnostic item
    lab_composite, lab_diag = {}, {}

    for r in csv.DictReader(open("outputs/tracss/tracss_human_lab_composite_2020_21.csv")):
        if r["lab_composite_equalweight"] not in ("", "nan"):
            lab_composite[r["iso3"]] = float(r["lab_composite_equalweight"])
        if r["lab_diagnostic_techniques"] not in ("", "nan"):
            lab_diag[r["iso3"]] = float(r["lab_diagnostic_techniques"])

    # 5. covariates (latest available year 2014-2020)
    COV_YEARS = list(range(2014, 2021))

    covariates = list(csv.DictReader(open("outputs/world_bank/covariates.csv")))
    def cov_values(name):
        return {r["iso3"]: (float(r[name]), int(float(r[name + "_year"])))
                for r in covariates if r[name]}
    gdp_ppp = cov_values("gdp_pc_ppp")
    che_pc = cov_values("che_pc_usd")
    uhc = cov_values("uhc_index")

    atlas_vol = {r["iso3"]: int(r["isolates_2014_2022"])
                 for r in csv.DictReader(open("outputs/atlas/atlas_testing_volume_country.csv"))}

    # 6. exposure registry
    _ANALYSIS_SET = set(analysis_iso)

    registry = []

    def reg(name, sector, level, waves_used, iso_present, ladder, role, long_elig="", note=""):
        n = len(set(iso_present) & _ANALYSIS_SET)
        miss = len(analysis_iso) - n
        registry.append(OrderedDict(
            name=name, sector=sector, level=level, waves_used=waves_used,
            n_countries=n, n_missing=miss, missing_frac=round(miss / len(analysis_iso), 3),
            ladder_stability=ladder, longitudinal_eligible=long_elig, role=role, note=note))

    item_iso_present = defaultdict(set)

    for (col, iso) in item_value:
        item_iso_present[col].add(iso)

    for spec in item_specs:
        reg(spec["name"], spec["sector"], "item", spec["reference_wave"],
            item_iso_present.get(spec["name"], set()), spec["ladder_stability"], spec["role"],
            long_elig=spec["longitudinal_eligible"],
            note=("boundary_wave" if spec["boundary_wave"] else "") +
                 ("; cross_sectional_only" if spec["cross_sectional_only"] else "") +
                 ("; post-2021 descriptive only" if spec["role"] == "descriptive" else ""))

    for dom in DOMAINS:
        reg("exp_domain_" + dom, domain_sector[dom], "domain", BOUNDARY_WAVE, set(domain_value[dom]),
            "domain_composite", "primary", note=f"{domain_item_2021[dom]} at {BOUNDARY_WAVE}")

    reg("exp_total_readiness", "cross", "derived", BOUNDARY_WAVE, set(total_readiness),
        "domain_mean", "primary", note="mean of 6 domains at 2020-21 (complete case)")

    for dom in DOMAINS:
        reg("exp_lean_" + dom, domain_sector[dom], "lean", BOUNDARY_WAVE, set(lean_value[dom]),
            "domain_lean", "primary", note=f"{dom} minus mean of the other 5 domains")

    reg("exp_lab_composite_2020_21", "human", "composite", BOUNDARY_WAVE, set(lab_composite),
        "human_lab_composite", "primary", note="2020-21 human clinical-lab composite (complete case)")

    reg("exp_lab_diagnostic_2020_21", "human", "item", BOUNDARY_WAVE, set(lab_diag),
        "human_lab_single_item", "sensitivity", note="diagnostic techniques item alone (same as ITEM_063)")

    # 7. wide analysis dataset
    outcome_cols = ["organism", "source", "antibiotic", "country", "iso3", "n_real_windows",
                    "first_window", "last_window", "d_level", "d_velocity", "outcome_class",
                    "improving", "combo_not_exactly_verified"]

    exp_item_cols = [s["name"] for s in item_specs]

    exp_domain_cols = ["exp_domain_" + d for d in DOMAINS]

    exp_lean_cols = ["exp_lean_" + d for d in DOMAINS]

    exp_derived_cols = ["exp_total_readiness"] + exp_lean_cols + \
                       ["exp_lab_composite_2020_21", "exp_lab_diagnostic_2020_21"]

    cov_cols = ["cov_gdp_pc_ppp", "cov_log_gdp_pc_ppp", "cov_gdp_pc_ppp_year",
                "cov_che_pc_usd", "cov_log_che_pc_usd", "cov_che_pc_usd_year",
                "cov_uhc_index", "cov_uhc_index_year", "cov_income_group",
                "cov_atlas_testing_volume_2014_2022", "cov_log_atlas_testing_volume"]

    header = outcome_cols + exp_domain_cols + exp_derived_cols + exp_item_cols + cov_cols

    BLANK = ""

    def logv(x):
        return round(math.log(x), 6) if (x is not None and x > 0) else BLANK

    def cov_for(iso):
        g, h, u, a = gdp_ppp.get(iso), che_pc.get(iso), uhc.get(iso), atlas_vol.get(iso)
        return {
            "cov_gdp_pc_ppp": g[0] if g else BLANK,
            "cov_log_gdp_pc_ppp": logv(g[0]) if g else BLANK,
            "cov_gdp_pc_ppp_year": g[1] if g else BLANK,
            "cov_che_pc_usd": h[0] if h else BLANK,
            "cov_log_che_pc_usd": logv(h[0]) if h else BLANK,
            "cov_che_pc_usd_year": h[1] if h else BLANK,
            "cov_uhc_index": u[0] if u else BLANK,
            "cov_uhc_index_year": u[1] if u else BLANK,
            "cov_income_group": "not_available",
            "cov_atlas_testing_volume_2014_2022": a if a is not None else BLANK,
            "cov_log_atlas_testing_volume": logv(a) if a is not None else BLANK,
        }

    rows_out = []

    for r in outcome_rows:
        iso = r["iso3"]
        row = {c: r.get(c, BLANK) for c in outcome_cols}
        for dom in DOMAINS:
            row["exp_domain_" + dom] = domain_value[dom].get(iso, BLANK)
            row["exp_lean_" + dom] = round(lean_value[dom][iso], 6) if iso in lean_value[dom] else BLANK
        row["exp_total_readiness"] = round(total_readiness[iso], 6) if iso in total_readiness else BLANK
        row["exp_lab_composite_2020_21"] = round(lab_composite[iso], 6) if iso in lab_composite else BLANK
        row["exp_lab_diagnostic_2020_21"] = lab_diag.get(iso, BLANK)
        for col in exp_item_cols:
            row[col] = item_value.get((col, iso), BLANK)
        row.update(cov_for(iso))
        rows_out.append(row)

    seen = set()

    for r in rows_out:
        k = (r["organism"], r["source"], r["antibiotic"], r["iso3"])
        check(k not in seen, f"duplicate country-combination key {k}")
        seen.add(k)

    ds_path = os.path.join(OUT, "analysis_dataset.csv")

    with open(ds_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header)
        w.writeheader()
        w.writerows(rows_out)

    reg_path = os.path.join(OUT, "exposure_registry.csv")

    with open(reg_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(registry[0].keys()))
        w.writeheader()
        w.writerows(registry)

    # 8. missingness
    miss_path = os.path.join(OUT, "missingness_by_exposure_wave.csv")

    with open(miss_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["exposure", "level", "sector", "reference_wave", "role",
                    "n_countries_present", "n_missing", "missing_frac"])

        def line(name, level, sector, wave, role, present):
            n = len(present)
            w.writerow([name, level, sector, wave, role, n, len(analysis_iso) - n,
                        round((len(analysis_iso) - n) / len(analysis_iso), 3)])

        for spec in item_specs:
            line(spec["name"], "item", spec["sector"], spec["reference_wave"], spec["role"],
                 item_iso_present.get(spec["name"], set()))
        for dom in DOMAINS:
            line("exp_domain_" + dom, "domain", domain_sector[dom], BOUNDARY_WAVE, "primary",
                 set(domain_value[dom]))
        for nm, present in [("exp_total_readiness", set(total_readiness) & _ANALYSIS_SET),
                            ("exp_lab_composite_2020_21", set(lab_composite) & _ANALYSIS_SET),
                            ("exp_lab_diagnostic_2020_21", set(lab_diag) & _ANALYSIS_SET)]:
            line(nm, "derived", "-", BOUNDARY_WAVE, "-", present)
        for nm, d in [("cov_gdp_pc_ppp", gdp_ppp), ("cov_che_pc_usd", che_pc), ("cov_uhc_index", uhc),
                      ("cov_atlas_testing_volume_2014_2022", atlas_vol)]:
            line(nm, "covariate", "-", "-", "-", {iso for iso in analysis_iso if iso in d})

    def sha256(path):
        return hashlib.sha256(open(path, "rb").read()).hexdigest()

    manifest = OrderedDict(
        created_utc=datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        n_country_combination_rows=len(rows_out),
        n_countries=len(analysis_iso),
        countries=analysis_iso,
        n_final_combinations=len(final_combos),
        non_member_rows_dropped=n_dropped_non_member,
        n_exposures=len(registry),
        exposures_by_sector=dict(Counter(r["sector"] for r in registry)),
        exposures_by_level=dict(Counter(r["level"] for r in registry)),
        exposures_by_role=dict(Counter(r["role"] for r in registry)),
        reference_wave=BOUNDARY_WAVE,
        covariate_years=f"latest available in {COV_YEARS[0]}-{COV_YEARS[-1]}",
        artefacts={os.path.basename(p): {"sha256": sha256(p), "bytes": os.path.getsize(p)}
                   for p in (ds_path, reg_path, miss_path)},
    )

    with open(os.path.join(OUT, "assembly_manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)

    print(f"dataset rows={len(rows_out)} countries={len(analysis_iso)} exposures={len(registry)}")


def merge_investment():
    """Prepare country-level burden and investment for the R&D analysis."""
    KEY = ["organism", "source", "antibiotic", "country"]
    def load_country_table():
        ad = pd.read_csv("outputs/dataset/analysis_dataset.csv")
        rw = pd.read_csv("outputs/trajectories/reconstructed_windows.csv")
        cov = pd.read_csv("outputs/hub/hub_by_country.csv")

        pop = pd.read_csv("outputs/world_bank/population_2020.csv")
        gdp = (ad.groupby("iso3")["cov_gdp_pc_ppp"].first().reset_index()
               .rename(columns={"cov_gdp_pc_ppp": "gdp_pc_ppp"}))

        # baseline amplitude = amplitude in the earliest reconstructed window
        rw = rw.rename(columns={"organism_global": "organism"})
        rw["win_start"] = rw["window"].str.slice(0, 4).astype(int)
        baseline = (rw.sort_values("win_start")
                    .groupby(KEY, as_index=False).first()[KEY + ["amplitude"]]
                    .rename(columns={"amplitude": "baseline_amplitude"}))
        adm = ad.merge(baseline, on=KEY, how="left")
        assert adm["baseline_amplitude"].notna().all()

        g = adm.groupby(["organism", "source", "antibiotic"])["baseline_amplitude"]
        adm["baseline_z"] = (adm["baseline_amplitude"] - g.transform("mean")) / g.transform("std")

        burden = (adm.groupby(["iso3", "country"])
                  .agg(burden_mean_baseline_pctR=("baseline_amplitude", "mean"),
                       burden_standardised=("baseline_z", "mean"),
                       n_combos=("baseline_amplitude", "size"),
                       share_improving=("improving", "mean"))
                  .reset_index())
        burden["share_not_improving"] = 1 - burden["share_improving"]

        hub = cov[["iso3", "in_hub_as_institution", "n_projects_hosted", "usd_nominal_hosted",
                   "in_hub_as_funder", "usd_nominal_funded"]]
        df = (burden.merge(hub, on="iso3", how="left").merge(pop, on="iso3", how="left")
              .merge(gdp, on="iso3", how="left"))

        # An absent country has no Hub record; its investment remains missing below.
        for col in ["in_hub_as_institution", "in_hub_as_funder"]:
            df[col] = df[col].fillna(False).astype(bool)

        # no Hub record = missing, not a recorded zero
        df["hub_record_missing"] = df["in_hub_as_institution"].isna() | (df["in_hub_as_institution"] == False)
        for c in ["n_projects_hosted", "usd_nominal_hosted", "usd_nominal_funded"]:
            df.loc[df["hub_record_missing"], c] = np.nan
        total = df.loc[~df["hub_record_missing"], "usd_nominal_hosted"].sum()
        df["share_of_global_hosted"] = df["usd_nominal_hosted"] / total
        df["hosted_usd_per_capita"] = df["usd_nominal_hosted"] / df["population_2020"]
        return df.sort_values("burden_mean_baseline_pctR", ascending=False).reset_index(drop=True)
    load_country_table().to_csv("outputs/dataset/rd_country_dataset.csv", index=False)


def main():
    os.chdir(ROOT)
    select_combinations()
    classify_trajectories()
    exclude_near_zero_resistance()
    save_final_combinations()
    merge_readiness()
    merge_investment()


if __name__ == "__main__":
    main()
