"""Count ATLAS isolates, summarise resistance and calculate testing volume."""
# Copyright 2026 Tavpritesh Sethi and Jasmine Kaur (Tavlab, IIITD)
# Licensed under the Apache License, Version 2.0 (see LICENSE)

import csv
import json
import os
import pandas as pd
import statistics
from collections import Counter, defaultdict
from common import (
    ATLAS_ALIASES, DROP_BY_RULE, NAME_TO_ISO3, ROOT, WINDOWS, WINDOW_LABELS, atlas_column, normalise_country, organism_groups, parse_country_list
)


def count_isolates():
    """Count valid S/I/R results by combination, country and rolling window."""

    SOURCES = {"Blood", "Sputum", "Urine", "Wound"}

    VALID_AST = {"Susceptible", "Intermediate", "Resistant"}

    csv.field_size_limit(10 ** 7)

    OUT = "outputs/atlas"
    os.makedirs("outputs/atlas", exist_ok=True)
    ATLAS = "data/atlas/atlas_vivli_2004_2024.csv"

    UNPARSED = []

    def note_unparsed(row_index, column, raw, reason):
        UNPARSED.append(dict(file=os.path.basename(ATLAS), row_index=row_index, column=column,
                             raw_cell_value="" if raw is None else raw, reason=reason))

    def quantile(vals, p):
        vs = sorted(vals)
        if not vs:
            return ""
        k = (len(vs) - 1) * p
        lo, hi = int(k), min(int(k) + 1, len(vs) - 1)
        return vs[lo] + (vs[hi] - vs[lo]) * (k - lo)

    with open("data/amrorbit/global.csv", newline="") as fh:
        grows = list(csv.DictReader(fh))
    assert len(grows) == 362, f"expected 362 global.csv rows, got {len(grows)}"

    targets = {}  # (org, src, ab) -> set of trajectory countries
    for r in grows:
        cs = set()
        for col in ("Spiral In", "Spiral Out", "Constant", "Other"):
            cs |= set(parse_country_list(r[col]))
        targets[(r["Organism"], r["Source"], r["Antibiotic"])] = cs

    with open(ATLAS, newline="") as fh:
        header = next(csv.reader(fh))
    ab_col = {}
    for (_, _, ab) in targets:
        col = atlas_column(ab)
        if col not in header:
            raise AssertionError(f"no ATLAS column for antibiotic {ab!r} (tried {col!r})")
        ab_col[ab] = col

    need = defaultdict(set)
    for (org, src, ab) in targets:
        need[(org, src)].add(ab)

    # single streaming pass over ATLAS
    counts = defaultdict(int)       # (org, src, ab, country, year) -> isolates
    years_seen = defaultdict(set)
    years_in_win = defaultdict(set)
    staph_no_oxacillin = bad_year = n_rows = 0
    bad_ast = defaultdict(int)

    with open(ATLAS, newline="") as fh:
        for i, row in enumerate(csv.DictReader(fh)):
            n_rows += 1
            src = (row.get("Source") or "").strip()
            if src not in SOURCES:
                continue
            species = (row.get("Species") or "").strip()
            ox = (row.get("Oxacillin_I") or "").strip()
            if species == "Staphylococcus aureus" and ox not in ("Resistant", "Susceptible"):
                staph_no_oxacillin += 1
                if ox not in ("", "Intermediate"):
                    note_unparsed(i, "Oxacillin_I", ox, "unexpected_oxacillin_value")
            orgs = organism_groups(species, ox)
            if not orgs:
                continue
            raw_year = (row.get("Year") or "").strip()
            if not raw_year.isdigit():
                note_unparsed(i, "Year", row.get("Year"), "not_an_integer_year")
                bad_year += 1
                continue
            year = int(raw_year)
            country = (row.get("Country") or "").strip()
            if country == "":
                note_unparsed(i, "Country", row.get("Country"), "empty_country")
                continue
            for org in orgs:
                for ab in need.get((org, src), ()):
                    v = (row.get(ab_col[ab]) or "").strip()
                    if v == "":
                        continue
                    if v not in VALID_AST:
                        note_unparsed(i, ab_col[ab], v, "unexpected_AST_interpretation")
                        bad_ast[v] += 1
                        continue
                    counts[(org, src, ab, country, year)] += 1
                    years_seen[(org, src, ab)].add(year)
                    if 2014 <= year <= 2022:
                        years_in_win[(org, src, ab)].add(year)

    # roll years into windows
    cw = defaultdict(int)
    for (org, src, ab, country, year), n in counts.items():
        for (a, b), wl in zip(WINDOWS, WINDOW_LABELS):
            if a <= year <= b:
                cw[(org, src, ab, country, wl)] += n

    density_rows, cwrows = [], []
    for key in sorted(targets):
        org, src, ab = key
        per_country = defaultdict(dict)
        for (o, s, a, c, w), n in cw.items():
            if (o, s, a) == key:
                per_country[c][w] = n
        for c, wd in per_country.items():
            for w, n in sorted(wd.items()):
                cwrows.append(dict(organism=org, source=src, antibiotic=ab, country=c, window=w, isolates=n))
        all_counts = [n for wd in per_country.values() for n in wd.values()]
        ys = sorted(years_seen.get(key, ()))
        yw = sorted(years_in_win.get(key, ()))
        wins_with_data = sorted({wl for (o, s_, a, c, wl) in cw if (o, s_, a) == key})
        traj = targets[key]
        density_rows.append(dict(
            organism=org, source=src, antibiotic=ab, atlas_column=ab_col[ab],
            n_countries_any_data=len(per_country),
            n_countries_ge2_valid_windows=sum(1 for wd in per_country.values() if len(wd) >= 2),
            n_countries_ge2_windows_ge10_isolates=sum(
                1 for wd in per_country.values() if sum(1 for n in wd.values() if n >= 10) >= 2),
            n_countries_ge2_windows_ge30_isolates=sum(
                1 for wd in per_country.values() if sum(1 for n in wd.values() if n >= 30) >= 2),
            n_country_windows=len(all_counts),
            isolates_total=sum(all_counts),
            isolates_per_cw_median=(statistics.median(all_counts) if all_counts else ""),
            isolates_per_cw_q1=quantile(all_counts, 0.25),
            isolates_per_cw_q3=quantile(all_counts, 0.75),
            isolates_per_cw_iqr=("" if not all_counts
                                 else quantile(all_counts, 0.75) - quantile(all_counts, 0.25)),
            pct_cw_below_10=("" if not all_counts else
                             round(100.0 * sum(1 for n in all_counts if n < 10) / len(all_counts), 2)),
            pct_cw_below_30=("" if not all_counts else
                             round(100.0 * sum(1 for n in all_counts if n < 30) / len(all_counts), 2)),
            year_min=(ys[0] if ys else ""), year_max=(ys[-1] if ys else ""),
            n_years_with_data=len(ys),
            year_coverage_all=",".join(str(y) for y in ys),
            year_min_in_windows=(yw[0] if yw else ""),
            year_max_in_windows=(yw[-1] if yw else ""),
            n_years_in_windows=len(yw),
            n_windows_with_data=len(wins_with_data),
            windows_with_data="|".join(wins_with_data),
            n_trajectory_countries=len(traj),
            n_countries_overlap_trajectory=len(traj & set(per_country)),
            n_trajectory_countries_not_in_atlas=len(traj - set(per_country)),
        ))

    def dump(name, rows, fields=None):
        with open(os.path.join(OUT, name), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields or list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)

    dump("density_table.csv", density_rows)
    dump("density_country_window.csv", sorted(
        cwrows, key=lambda r: (r["organism"], r["source"], r["antibiotic"], r["country"], r["window"])))
    dump("atlas_unparsed.csv", UNPARSED, ["file", "row_index", "column", "raw_cell_value", "reason"])

    prov = dict(
        atlas_file=os.path.basename(ATLAS), atlas_rows_read=n_rows,
        windows=WINDOW_LABELS, sources=sorted(SOURCES), valid_ast_values=sorted(VALID_AST),
        staph_aureus_isolates_without_oxacillin_result=staph_no_oxacillin,
        rows_with_unparseable_year=bad_year,
        unexpected_ast_values=dict(bad_ast),
        unparsed_records=len(UNPARSED),
    )
    with open(os.path.join(OUT, "density_provenance.json"), "w") as fh:
        json.dump(prov, fh, indent=2)
    print(json.dumps(prov, indent=2))


def summarise_resistance():
    """Save resistance summaries for the near-zero resistance exclusion."""
    VALID_AST = {"Susceptible", "Intermediate", "Resistant"}
    ATLAS = "data/atlas/atlas_vivli_2004_2024.csv"
    def atlas_resistance(combos):
        with open(ATLAS, newline="") as fh:
            header = next(csv.reader(fh))
        need, ab_col = defaultdict(set), {}
        for (org, src, ab) in combos:
            col = atlas_column(ab)
            assert col in header, f"no ATLAS column {col!r}"
            ab_col[ab] = col
            need[(org, src)].add(ab)

        agg = defaultdict(lambda: [0, 0])   # (org, src, ab, country, window) -> [n_valid, n_R]
        with open(ATLAS, newline="") as fh:
            for row in csv.DictReader(fh):
                src = (row.get("Source") or "").strip()
                species = (row.get("Species") or "").strip()
                ox = (row.get("Oxacillin_I") or "").strip()
                for org in organism_groups(species, ox):
                    if (org, src) not in need:
                        continue
                    raw_year = (row.get("Year") or "").strip()
                    if not raw_year.isdigit():
                        continue
                    year = int(raw_year)
                    win = next((wl for (a, b), wl in zip(WINDOWS, WINDOW_LABELS) if a <= year <= b), None)
                    if win is None:
                        continue
                    country = (row.get("Country") or "").strip()
                    if not country:
                        continue
                    for ab in need[(org, src)]:
                        v = (row.get(ab_col[ab]) or "").strip()
                        if v not in VALID_AST:
                            continue
                        cell = agg[(org, src, ab, country, win)]
                        cell[0] += 1
                        if v == "Resistant":
                            cell[1] += 1

        pct_r = defaultdict(list)
        for (org, src, ab, c, w), (nv, nr) in agg.items():
            if nv > 0:
                pct_r[(org, src, ab)].append(100.0 * nr / nv)
        out = {}
        for k in combos:
            vs = pct_r.get(k, [])
            out[k] = dict(
                n_country_windows_atlas=len(vs),
                median_pctR=(statistics.median(vs) if vs else None),
                mean_pctR=(round(statistics.fmean(vs), 3) if vs else None),
                max_pctR=(max(vs) if vs else None),
            )
        return out
    combos = set(pd.read_csv("outputs/atlas/density_table.csv")[["organism", "source", "antibiotic"]].itertuples(index=False, name=None))
    res = atlas_resistance(combos)
    rows = [dict(organism=k[0], source=k[1], antibiotic=k[2], **res[k]) for k in sorted(combos)]
    pd.DataFrame(rows).to_csv("outputs/atlas/resistance_by_combo.csv", index=False)


def count_testing_volume():
    """Count all ATLAS rows by country and year for the adjustment covariate."""

    ATLAS = "data/atlas/atlas_vivli_2004_2024.csv"
    OUT = "outputs/atlas"
    os.makedirs("outputs/atlas", exist_ok=True)

    NAME_MAP = dict(NAME_TO_ISO3)

    NAME_MAP.update(ATLAS_ALIASES)

    def resolve(raw):
        name = normalise_country(raw)
        if name == "" or name.lower() == "nan":
            return None, name, "blank"
        if name in DROP_BY_RULE:
            return None, name, "drop_by_rule"
        iso3 = NAME_MAP.get(name)
        return (iso3, name, "mapped") if iso3 else (None, name, "UNMAPPED")

    cy = Counter()
    name_status, unmapped_names, dropped_names, bad_year = Counter(), Counter(), Counter(), Counter()
    nrows = 0
    with open(ATLAS, newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            nrows += 1
            iso3, name, status = resolve(row.get("Country"))
            name_status[status] += 1
            if status == "UNMAPPED":
                unmapped_names[name] += 1
                continue
            if status == "drop_by_rule":
                dropped_names[name] += 1
                continue
            if status == "blank":
                continue
            yr_raw = (row.get("Year") or "").strip()
            try:
                yr = int(yr_raw)
            except (ValueError, TypeError):
                bad_year[yr_raw] += 1
                continue
            cy[(iso3, yr)] += 1

    if unmapped_names:
        print("Unmapped ATLAS country names (rows):")
        for n, c in unmapped_names.most_common():
            print(f"  {n!r}: {c}")

    with open(os.path.join(OUT, "atlas_testing_volume_country_year.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["iso3", "year", "isolates"])
        for (iso3, yr), n in sorted(cy.items()):
            w.writerow([iso3, yr, n])

    traj_years = range(2014, 2023)
    per_country = defaultdict(lambda: {"isolates_2014_2022": 0, "isolates_all_years": 0,
                                       "years_with_data_2014_2022": set()})
    for (iso3, yr), n in cy.items():
        per_country[iso3]["isolates_all_years"] += n
        if yr in traj_years:
            per_country[iso3]["isolates_2014_2022"] += n
            per_country[iso3]["years_with_data_2014_2022"].add(yr)
    with open(os.path.join(OUT, "atlas_testing_volume_country.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["iso3", "isolates_2014_2022", "n_years_with_data_2014_2022",
                    "mean_isolates_per_year_2014_2022", "isolates_all_years"])
        for iso3 in sorted(per_country):
            d = per_country[iso3]
            ny = len(d["years_with_data_2014_2022"])
            w.writerow([iso3, d["isolates_2014_2022"], ny,
                        round(d["isolates_2014_2022"] / ny if ny else 0, 3), d["isolates_all_years"]])

    prov = {
        "atlas_rows_read": nrows,
        "name_status_counts": dict(name_status),
        "unmapped_names": dict(unmapped_names),
        "dropped_by_rule_names": dict(dropped_names),
        "unparseable_year_values": dict(bad_year),
        "n_countries_mapped": len(per_country),
        "definition": "isolates per country-year = count of ATLAS rows (all species and sources)",
    }
    with open(os.path.join(OUT, "atlas_testing_volume_provenance.json"), "w") as fh:
        json.dump(prov, fh, indent=2)
    print(f"rows read: {nrows} | countries: {len(per_country)} | "
          f"unmapped rows: {sum(unmapped_names.values())}")


def main():
    os.chdir(ROOT)
    count_isolates()
    summarise_resistance()
    count_testing_volume()


if __name__ == "__main__":
    main()
