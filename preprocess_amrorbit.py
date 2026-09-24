"""Reconstruct the supplied AMROrbit estimates and check the scorecard."""
# Copyright 2026 Tavpritesh Sethi and Jasmine Kaur (Tavlab, IIITD)
# Licensed under the Apache License, Version 2.0 (see LICENSE)

import csv
import glob
import json
import os
import statistics
from collections import defaultdict
from common import ROOT, WINDOW_LABELS, parse_country_list


def reconstruct_estimates():
    """Reconstruct window levels and slopes from the saved model estimates."""
    # amplitude = Global intercept + country deviation, velocity = Global slope + country deviation

    OUT = "outputs/trajectories"
    os.makedirs("outputs/trajectories", exist_ok=True)
    ME_DIR = "data/amrorbit/Model_estimates"
    Q_DIR = "data/amrorbit/Quadrant"

    CANONICAL_WINDOWS = WINDOW_LABELS

    UNPARSED = []

    def rel(path):
        return os.path.relpath(path, "data")

    def parse_float(raw, file, row_index, column):
        if raw is None:
            UNPARSED.append(dict(file=rel(file), row_index=row_index, column=column,
                                 raw_cell_value="", reason="missing_column"))
            return None, False
        s = raw.strip()
        if s == "":
            UNPARSED.append(dict(file=rel(file), row_index=row_index, column=column,
                                 raw_cell_value=raw, reason="empty_string"))
            return None, False
        try:
            return float(s), True
        except ValueError:
            UNPARSED.append(dict(file=rel(file), row_index=row_index, column=column,
                                 raw_cell_value=raw, reason="not_a_float"))
            return None, False

    def parse_int(raw, file, row_index, column):
        if raw is None or raw.strip() == "":
            UNPARSED.append(dict(file=rel(file), row_index=row_index, column=column,
                                 raw_cell_value=raw or "", reason="empty_string"))
            return None, False
        try:
            return int(raw.strip()), True
        except ValueError:
            UNPARSED.append(dict(file=rel(file), row_index=row_index, column=column,
                                 raw_cell_value=raw, reason="not_an_int"))
            return None, False

    def window_key(fname):
        parts = os.path.basename(fname)[:-4].split("_")
        return (int(parts[-2]), int(parts[-1]))

    def load_window(path):
        with open(path, newline="") as fh:
            rows = list(csv.DictReader(fh))
        g, dev, prov = None, {}, {}
        for i, r in enumerate(rows):
            country = r["Country"].strip()
            icept, ok1 = parse_float(r.get("intercept"), path, i, "intercept")
            slope, ok2 = parse_float(r.get("slope"), path, i, "slope")
            prov[country] = dict(file=path, row_index=i, raw_intercept=r.get("intercept"),
                                 raw_slope=r.get("slope"),
                                 status="ok" if (ok1 and ok2) else "unparsed")
            if not (ok1 and ok2):
                continue
            if country == "Global":
                g = (icept, slope)
            else:
                dev[country] = (icept, slope)
        if g is None:
            raise AssertionError(f"no Global row in {path}")
        return g, dev, prov

    def quadrant(amp, vel, amp_med, vel_med):
        amp_high = amp >= amp_med
        vel_high = vel >= vel_med
        if not amp_high and not vel_high:
            return 1
        if amp_high and not vel_high:
            return 2
        if not amp_high and vel_high:
            return 3
        return 4

    def build_combo(ab_dir):
        files = sorted(glob.glob(os.path.join(ab_dir, "*.csv")), key=window_key)
        windows, prov_all = {}, {}
        for f in files:
            wk = window_key(f)
            wlabel = f"{wk[0]}-{wk[1]}"
            if wlabel not in CANONICAL_WINDOWS:
                raise AssertionError(f"non-canonical window {wlabel} in {f}")
            g, dev, prov = load_window(f)
            if not dev:
                continue
            absv = {c: (g[0] + d[0], g[1] + d[1]) for c, d in dev.items()}
            amp_med = statistics.median(v[0] for v in absv.values())
            vel_med = statistics.median(v[1] for v in absv.values())
            windows[wlabel] = dict(
                amp_med=amp_med, vel_med=vel_med, global_intercept=g[0], global_slope=g[1], file=f,
                countries={c: (a, v, quadrant(a, v, amp_med, vel_med)) for c, (a, v) in absv.items()})
            prov_all[wlabel] = prov
        return windows, prov_all

    def trajectory_label(first_q, final_q):
        if first_q == final_q:
            return "Constant"
        if final_q == 1:
            return "Spiral In"
        if final_q == 4:
            return "Spiral Out"
        return "Other"

    def classify_combo(windows):
        countries = set()
        for w in windows.values():
            countries |= set(w["countries"])
        labels, detail = {}, {}
        first_w, final_w = CANONICAL_WINDOWS[0], CANONICAL_WINDOWS[-1]
        for c in sorted(countries):
            grid = [windows.get(w, {}).get("countries", {}).get(c, (None, None, 0))[2]
                    for w in CANONICAL_WINDOWS]
            fq, lq = grid[0], grid[-1]
            detail[c] = dict(grid=grid, n_valid_windows=sum(1 for q in grid if q != 0),
                             first_window=first_w, final_window=final_w, first_q=fq, final_q=lq)
            labels[c] = None if (fq == 0 and lq == 0) else trajectory_label(fq, lq)
        return labels, detail

    def load_global():
        with open("data/amrorbit/global.csv", newline="") as fh:
            rows = list(csv.DictReader(fh))
        assert len(rows) == 362, f"expected 362 global.csv rows, got {len(rows)}"
        idx = {}
        for r in rows:
            d = {}
            for col in ("Spiral In", "Spiral Out", "Constant", "Other"):
                for c in parse_country_list(r[col]):
                    d[c] = col
            idx[(r["Organism"], r["Source"], r["Antibiotic"])] = d
        return rows, idx

    def quadrant_grid_labels(qpath):
        with open(qpath, newline="") as fh:
            rows = list(csv.DictReader(fh))
        cols = [c for c in rows[0].keys() if c != "Country"]
        out = {}
        for i, r in enumerate(rows):
            vals = []
            for c in cols:
                v, ok = parse_int(r.get(c), qpath, i, c)
                if not ok:
                    raise AssertionError(f"unparseable quadrant cell {qpath} row {i} col {c}")
                vals.append(v)
            if vals[0] == 0 and vals[-1] == 0:
                continue
            out[r["Country"].strip()] = trajectory_label(vals[0], vals[-1])
        return out

    def resolve_organism_folders(me_combos, gidx):
        g_orgs = sorted(set(k[0] for k in gidx))
        folders = sorted(set(k[0] for k in me_combos))
        mapping, log = {}, []
        for f in folders:
            cands = [g for g in g_orgs if g.split()[0] == f.split()[0]]
            scores = {}
            for cand in cands:
                exact = total = 0
                for (fo, src, ab) in me_combos:
                    if fo != f or (cand, src, ab) not in gidx:
                        continue
                    qp = os.path.join(Q_DIR, f, src, f"{f}_{ab}_{src}.csv")
                    if not os.path.exists(qp):
                        continue
                    total += 1
                    if quadrant_grid_labels(qp) == gidx[(cand, src, ab)]:
                        exact += 1
                scores[cand] = dict(exact_label_match=exact, comparable=total)
            best = max(scores, key=lambda c: (scores[c]["exact_label_match"],
                                              scores[c]["comparable"])) if scores else None
            confident = bool(best and scores[best]["comparable"] > 0 and
                             scores[best]["exact_label_match"] == scores[best]["comparable"])
            mapping[f] = best if confident else None
            log.append(dict(folder=f, candidate_scores=scores, chosen=mapping[f], confident=confident))
        return mapping, log

    me_combos = {}
    for p in sorted(glob.glob(os.path.join(ME_DIR, "*", "*", "*"))):
        if not os.path.isdir(p):
            continue
        parts = p.split(os.sep)
        # one antibiotic directory name carries a trailing space; strip it
        org, src, abd = parts[-3].strip(), parts[-2].strip(), parts[-1].strip()
        ab = abd[:-2] if abd.endswith("_I") else abd
        me_combos[(org, src, ab)] = p

    grows, gidx = load_global()
    orgmap, orglog = resolve_organism_folders(me_combos, gidx)

    results, window_rows, traj_rows = {}, [], []
    for (folder_org, src, ab), path in sorted(me_combos.items()):
        org = orgmap.get(folder_org)
        windows, prov = build_combo(path)
        if not windows:
            continue
        labels, detail = classify_combo(windows)
        results[(folder_org, src, ab)] = dict(windows=windows, labels=labels, detail=detail,
                                              global_key=(org, src, ab) if org else None, path=path)
        for w, wd in windows.items():
            for c, (a, v, q) in wd["countries"].items():
                window_rows.append(dict(
                    organism_folder=folder_org, organism_global=org or "UNMAPPED",
                    source=src, antibiotic=ab, window=w, country=c,
                    amplitude=repr(a), velocity=repr(v), quadrant=q,
                    amp_median=repr(wd["amp_med"]), vel_median=repr(wd["vel_med"]),
                    global_intercept=repr(wd["global_intercept"]),
                    global_slope=repr(wd["global_slope"]),
                    prov_file=rel(wd["file"]), prov_row_index=prov[w][c]["row_index"],
                    prov_raw_intercept=prov[w][c]["raw_intercept"],
                    prov_raw_slope=prov[w][c]["raw_slope"]))
        for c, lab in labels.items():
            d = detail[c]
            traj_rows.append(dict(
                organism_folder=folder_org, organism_global=org or "UNMAPPED",
                source=src, antibiotic=ab, country=c,
                label=lab if lab else "EXCLUDED_absent_at_both_endpoints",
                quadrant_grid="|".join(str(x) for x in d["grid"]),
                n_valid_windows=d["n_valid_windows"],
                first_window=d["first_window"], final_window=d["final_window"],
                first_quadrant=d["first_q"], final_quadrant=d["final_q"],
                prov_dir=rel(path)))

    # verification against global.csv
    verified, mismatches, missing_data = [], [], []
    for key in gidx:
        cand = [k for k, v in results.items() if v["global_key"] == key]
        if not cand:
            missing_data.append(key)
            continue
        r = results[cand[0]]
        got = {c: l for c, l in r["labels"].items() if l}
        if got == gidx[key]:
            verified.append(key)
        else:
            diffs = [dict(country=c, expected=gidx[key].get(c), got=got.get(c),
                          grid="|".join(str(x) for x in r["detail"].get(c, {}).get("grid", [])))
                     for c in sorted(set(gidx[key]) | set(got)) if gidx[key].get(c) != got.get(c)]
            mismatches.append(dict(key=list(key), n_diff=len(diffs), diffs=diffs))
    unmapped = [list(k) for k, v in results.items() if v["global_key"] is None]

    # per-window check against Quadrant/
    cells_cmp = cells_bad = 0
    bad_by_combo = defaultdict(int)
    for (folder_org, src, ab), r in results.items():
        qp = os.path.join(Q_DIR, folder_org, src, f"{folder_org}_{ab}_{src}.csv")
        if not os.path.exists(qp):
            continue
        with open(qp, newline="") as fh:
            qrows = list(csv.DictReader(fh))
        for i, row in enumerate(qrows):
            c = row["Country"].strip()
            for w in CANONICAL_WINDOWS:
                if w not in row:
                    continue
                exp, ok = parse_int(row.get(w), qp, i, w)
                if not ok or exp == 0:
                    continue
                mine = r["windows"].get(w, {}).get("countries", {}).get(c)
                if mine is None:
                    continue
                cells_cmp += 1
                if mine[2] != exp:
                    cells_bad += 1
                    bad_by_combo[(folder_org, src, ab)] += 1

    def dump(name, rows, fields=None):
        with open(os.path.join(OUT, name), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields or list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)

    dump("reconstructed_windows.csv", window_rows)
    dump("reconstructed_trajectories.csv", traj_rows)
    dump("unparsed_values.csv", UNPARSED, ["file", "row_index", "column", "raw_cell_value", "reason"])
    dump("verified_rows.csv", [dict(organism=k[0], source=k[1], antibiotic=k[2]) for k in sorted(verified)])

    report = dict(
        global_csv_rows=len(grows),
        model_estimates_combos=len(me_combos),
        organism_folder_resolution=orglog,
        rows_with_source_data=len(verified) + len(mismatches),
        rows_verified_exact=len(verified),
        rows_mismatched=len(mismatches),
        rows_no_source_data=len(missing_data),
        me_combos_unmapped_to_global=len(unmapped),
        quadrant_cells_compared=cells_cmp,
        quadrant_cells_mismatched=cells_bad,
        quadrant_combos_with_cell_mismatch={f"{k[0]}|{k[1]}|{k[2]}": v
                                            for k, v in sorted(bad_by_combo.items())},
        unparsed_value_count=len(UNPARSED),
        mismatch_detail=mismatches,
    )
    with open(os.path.join(OUT, "verification_report.json"), "w") as fh:
        json.dump(report, fh, indent=2, default=str)

    print(f"global.csv rows with source data: {len(verified) + len(mismatches)}; "
          f"reproduced exactly: {len(verified)}; mismatched: {len(mismatches)}")
    print(f"Quadrant/ cells compared: {cells_cmp}; mismatched: {cells_bad}")
    for e in orglog:
        print(f"  {e['folder']} -> {e['chosen']} (confident={e['confident']})")

    membership = []
    for r in grows:
        for label in ("Spiral In", "Spiral Out", "Constant", "Other"):
            for country in parse_country_list(r[label]):
                membership.append(dict(organism=r["Organism"], source=r["Source"],
                                       antibiotic=r["Antibiotic"], country=country, label=label))
    dump("trajectory_countries.csv", membership)


def main():
    os.chdir(ROOT)
    reconstruct_estimates()


if __name__ == "__main__":
    main()
