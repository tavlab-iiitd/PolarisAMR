"""Clean TrACSS 2016-2025 and build the question crosswalk."""
# Copyright 2026 Tavpritesh Sethi and Jasmine Kaur (Tavlab, IIITD)
# Licensed under the Apache License, Version 2.0 (see LICENSE)

import difflib
import hashlib
import itertools
import json
import os
import pandas as pd
import re
from collections import Counter, defaultdict
from common import DROP_BY_RULE, NAME_TO_ISO3, ROOT, norm_ws, normalise_country
from datetime import datetime, timezone


WAVES = [
    ("2016-17", "AMR-self-assessment-survey-country-responses-2016-17.xlsx", 0, 2),
    ("2017-18", "AMR-self-assessment-survey-country-responses-2017-18.xlsx", 0, 2),
    ("2018-19", "AMR-self-assessment-survey-country-responses-2018-19.xls", 0, 1),
    ("2019-20", "AMR-self-assessment-survey-country-responses-2019-20.xls", 0, 0),
    ("2020-21", "AMR-self-assessment-survey-country-responses-2020-21.xlsx", 0, 2),
    ("2022", "AMR-self-assessment-survey-responses-2022.xlsx", 0, 0),
    ("2023", "AMR-self-assessment-survey-responses-TrACSS-2023.xlsx", 0, 0),
    ("2024", "TrACSS_Survey_2024_Dataset_v3.xlsx", 0, 0),
    ("2025", "TrACSS-2025-Data-export-01102025.xlsx", 0, 0),
]
WAVE_ORDER = [w for w, _, _, _ in WAVES]


def extract_responses():
    """Read the nine survey workbooks and parse their responses."""
    # Read all 9 TrACSS waves into one table


    # answers look like "C - Formalized multisector ..." - only accept a dash here,
    # otherwise things like "E. coli" or "A) Pesticides" get read as grades
    _LABELLED = re.compile(r"^([A-Ea-e])\s*-\s*(.+)$")

    _BARE = re.compile(r"^([A-Ea-e])$")

    def parse_ordinal(value_norm: str):
        m = _BARE.match(value_norm)
        if m:
            return m.group(1).upper(), ""
        m = _LABELLED.match(value_norm)
        if m:
            return m.group(1).upper(), m.group(2).strip()
        return None, None

    _YESNO = {"YES", "NO"}

    def classify_value(value_norm: str) -> str:
        """Type of a single value: blank | ordinal_bare | ordinal_labelled | binary | numeric | text."""
        if value_norm == "":
            return "blank"
        g, opt = parse_ordinal(value_norm)
        if g is not None:
            return "ordinal_bare" if opt == "" else "ordinal_labelled"
        if value_norm.upper() in _YESNO:
            return "binary"
        try:
            float(value_norm.replace(",", ""))
            return "numeric"
        except ValueError:
            return "text"

    DOMINANCE = 0.90

    def column_response_type(value_types) -> str:
        vals = [t for t in value_types if t != "blank"]
        if not vals:
            return "blank"
        n = len(vals)
        c = Counter(vals)
        ordinal = c["ordinal_labelled"] + c["ordinal_bare"]
        if ordinal / n >= DOMINANCE:
            return "ordinal_labelled" if c["ordinal_labelled"] >= c["ordinal_bare"] else "ordinal_bare"
        if c["binary"] / n >= DOMINANCE:
            return "binary"
        if c["numeric"] / n >= DOMINANCE:
            return "numeric"
        if (c["binary"] + c["text"]) / n >= DOMINANCE and c["binary"]:
            return "binary_or_text"
        if (c["numeric"] + c["text"]) / n >= DOMINANCE and c["numeric"]:
            return "numeric_or_text"
        return "text"

    ORDINAL_TYPES = ("ordinal_bare", "ordinal_labelled")

    def to_iso3(raw):
        name = normalise_country(raw)
        if name == "" or name.lower() == "nan":
            return None, name, "blank"
        if name in DROP_BY_RULE:
            return None, name, "drop_by_rule"
        iso3 = NAME_TO_ISO3.get(name)
        if iso3:
            return iso3, name, "mapped"
        return None, name, "UNMAPPED"

    TRACSS_DIR = "data/tracss"
    OUT_DIR = "outputs/tracss"
    os.makedirs("outputs/tracss", exist_ok=True)
    HEADER_SCAN_ROWS = 15

    # header-row detection
    def detect_header_row(raw: pd.DataFrame) -> tuple[int, dict]:
        scores = {}
        for i in range(min(HEADER_SCAN_ROWS, len(raw))):
            cells = [norm_ws(v) for v in raw.iloc[i].tolist()]
            scores[i] = len({c for c in cells if len(c) >= 8})
        best = max(scores, key=lambda i: (scores[i], -i))
        return best, scores

    def find_country_column(wave: str, header: list[str]) -> int:
        lowered = [h.lower() for h in header]
        cand = [i for i, h in enumerate(lowered) if "countr" in h and len(h) < 60]
        if not cand:
            raise AssertionError(
                f"[{wave}] detected header row contains no country-like column: {header[:8]}"
            )
        # prefer the leftmost bare country identifier over a question mentioning countries
        for i in cand:
            h = lowered[i]
            if h.startswith("country") or h.endswith("country") or "country name" in h:
                if not any(ch.isdigit() for ch in h):
                    return i
        return cand[0]

    def find_iso3_column(header: list[str]):
        for i, h in enumerate(header):
            if h.strip().lower().replace(" ", "") in ("iso3", "iso3code", "iso3_code"):
                return i
        return None

    # extraction of one wave
    def extract_wave(wave, filename, sheet_idx, expected_header):
        path = os.path.join(TRACSS_DIR, filename)
        xl = pd.ExcelFile(path)
        sheet_name = xl.sheet_names[sheet_idx]
        raw = pd.read_excel(path, sheet_name=sheet_idx, header=None, dtype=object)

        header_row, scores = detect_header_row(raw)
        header = [norm_ws(v) for v in raw.iloc[header_row].tolist()]
        country_col = find_country_column(wave, header)
        iso3_col = find_iso3_column(header)

        prov = {
            "wave": wave,
            "file": filename,
            "sheet": sheet_name,
            "sheet_index": sheet_idx,
            "header_row_detected": header_row,
            "header_row_expected": expected_header,
            "header_row_matches_expected": header_row == expected_header,
            "header_detection_scores": scores,
            "n_rows_raw": int(raw.shape[0]),
            "n_data_rows": int(raw.shape[0] - header_row - 1),
            "n_columns": int(raw.shape[1]),
            "country_column_index": country_col,
            "country_column_header": header[country_col],
            "native_iso3_column_index": iso3_col,
            "native_iso3_column_header": header[iso3_col] if iso3_col is not None else None,
        }

        body = raw.iloc[header_row + 1:].reset_index(drop=True)

        # resolve countries; non-country rows (blank, footnotes) are dropped by rule
        row_iso3, row_name, row_status = [], [], []
        for _, r in body.iterrows():
            iso3, name, status = to_iso3(r.iloc[country_col])
            row_iso3.append(iso3)
            row_name.append(name)
            row_status.append(status)

        unmapped = sorted({n for n, s in zip(row_name, row_status) if s == "UNMAPPED"})
        if unmapped:
            raise AssertionError(
                f"[{wave}] {len(unmapped)} country name(s) missing from the ISO3 lookup: {unmapped}"
            )

        drop_counts = Counter(s for s in row_status if s != "mapped")

        # cross-check against a native ISO3 column where one exists
        iso3_disagreements = []
        if iso3_col is not None:
            for i, (iso3, name, status) in enumerate(zip(row_iso3, row_name, row_status)):
                if status != "mapped":
                    continue
                native = norm_ws(body.iloc[i, iso3_col]).upper()
                if native and native != iso3:
                    iso3_disagreements.append(
                        {"wave": wave, "country": name, "lookup_iso3": iso3, "native_iso3": native}
                    )

        keep = [i for i, s in enumerate(row_status) if s == "mapped"]
        iso_counts = Counter(row_iso3[i] for i in keep)
        dup = {k: v for k, v in iso_counts.items() if v > 1}

        # pass 1: type every column from its observed values
        col_value_types = defaultdict(list)
        for ci in range(raw.shape[1]):
            if ci == country_col or ci == iso3_col:
                continue
            for i in keep:
                v = norm_ws(body.iloc[i, ci])
                col_value_types[ci].append(classify_value(v))
        col_type = {ci: column_response_type(ts) for ci, ts in col_value_types.items()}

        # A-E columns: record any values outside A-E
        residue_rows = []
        for ci, ctype in col_type.items():
            if ctype not in ORDINAL_TYPES:
                continue
            bad = Counter()
            for i in keep:
                v = norm_ws(body.iloc[i, ci])
                if v == "":
                    continue
                g, _ = parse_ordinal(v)
                if g is None:
                    bad[v] += 1
            if bad:
                residue_rows.append(
                    {
                        "wave": wave,
                        "column_index": ci,
                        "raw_header": header[ci],
                        "response_type": ctype,
                        "n_residue": sum(bad.values()),
                        "residue_values": json.dumps(bad.most_common(20)),
                    }
                )

        # pass 2: emit long rows
        records, unparsed_rows = [], []
        for ci in range(raw.shape[1]):
            if ci == country_col or ci == iso3_col:
                continue
            ctype = col_type[ci]
            hdr = header[ci]
            for i in keep:
                v = norm_ws(body.iloc[i, ci])
                if v == "":
                    continue
                grade, option_text, parse_status = None, None, "n/a"
                if ctype in ORDINAL_TYPES:
                    grade, option_text = parse_ordinal(v)
                    if grade is None:
                        parse_status = "unparsed"
                        unparsed_rows.append(
                            {
                                "wave": wave,
                                "file": filename,
                                "sheet": sheet_name,
                                "header_row": header_row,
                                "column_index": ci,
                                "raw_header": hdr,
                                "iso3": row_iso3[i],
                                "expected_type": ctype,
                                "raw_cell_value": v,
                            }
                        )
                    else:
                        parse_status = "parsed"
                records.append(
                    {
                        "iso3": row_iso3[i],
                        "country": row_name[i],
                        "wave": wave,
                        "file": filename,
                        "sheet": sheet_name,
                        "header_row": header_row,
                        "column_index": ci,
                        "raw_header": hdr,
                        "response_type": ctype,
                        "parse_status": parse_status,
                        "grade": grade,
                        "option_text": option_text,
                        "raw_value": v,
                        "raw_cell_value": v,
                    }
                )

        summary = {
            **{k: v for k, v in prov.items() if k != "header_detection_scores"},
            "n_rows_kept": len(keep),
            "n_rows_dropped": int(len(row_status) - len(keep)),
            "rows_dropped_by_reason": dict(drop_counts),
            "n_unique_iso3": len(set(row_iso3[i] for i in keep)),
            "duplicate_iso3_within_wave": dup,
            "n_long_rows": len(records),
            "n_unparsed": len(unparsed_rows),
            "n_ordinal_columns": sum(1 for t in col_type.values() if t in ORDINAL_TYPES),
            "column_type_counts": dict(Counter(col_type.values())),
            "native_iso3_disagreements": len(iso3_disagreements),
        }
        return records, unparsed_rows, residue_rows, summary, prov, iso3_disagreements

    all_records, all_unparsed, all_residue, summaries, provs, all_disagree = [], [], [], [], [], []
    failures = []

    for wave, fn, si, exp in WAVES:
        try:
            rec, unp, res, summ, prov, dis = extract_wave(wave, fn, si, exp)
        except AssertionError as e:
            failures.append(str(e))
            print(f"ERROR: {e}")
            continue
        all_records += rec
        all_unparsed += unp
        all_residue += res
        all_disagree += dis
        summaries.append(summ)
        provs.append(prov)
        print(
            f"{wave}: header={prov['header_row_detected']} "
            f"({'OK' if prov['header_row_matches_expected'] else 'MISMATCH'}) "
            f"cols={prov['n_columns']} kept={summ['n_rows_kept']} "
            f"iso3={summ['n_unique_iso3']} long_rows={summ['n_long_rows']} "
            f"unparsed={summ['n_unparsed']}"
        )

    if failures:
        raise SystemExit("Stopped: errors above; no outputs written.")

    long_df = pd.DataFrame(all_records)
    long_df.to_csv(os.path.join(OUT_DIR, "tracss_raw_long.csv"), index=False)
    pd.DataFrame(provs).drop(columns=["header_detection_scores"]).to_csv(
        os.path.join(OUT_DIR, "tracss_wave_provenance.csv"), index=False
    )
    unparsed_cols = ["wave", "file", "sheet", "header_row", "column_index", "raw_header",
                     "iso3", "expected_type", "raw_cell_value"]
    pd.DataFrame(all_unparsed, columns=unparsed_cols).to_csv(
        os.path.join(OUT_DIR, "tracss_unparsed_values.csv"), index=False
    )
    pd.DataFrame(all_residue).to_csv(
        os.path.join(OUT_DIR, "tracss_ordinal_residue.csv"), index=False
    )
    if all_disagree:
        pd.DataFrame(all_disagree).to_csv(
            os.path.join(OUT_DIR, "tracss_native_iso3_disagreements.csv"), index=False
        )

    # global checks
    checks = {}
    ordm = long_df[long_df.response_type.isin(ORDINAL_TYPES) & (long_df.parse_status == "parsed")]
    vocab = sorted(ordm.grade.dropna().unique())
    checks["ordinal_grade_vocabulary"] = vocab
    checks["ordinal_vocab_subset_ABCDE"] = set(vocab) <= set("ABCDE")
    checks["n_waves"] = len(summaries)
    checks["all_header_rows_match_expected"] = all(p["header_row_matches_expected"] for p in provs)
    checks["total_long_rows"] = int(len(long_df))
    checks["total_unparsed"] = len(all_unparsed)
    checks["n_ordinal_columns_with_residue"] = len(all_residue)
    checks["drop_by_rule_observed_in_tracss"] = sorted(set(long_df.country) & set(DROP_BY_RULE))

    with open(os.path.join(OUT_DIR, "tracss_extraction_report.json"), "w") as fh:
        json.dump({"per_wave": summaries, "global_checks": checks}, fh, indent=2, default=str)

    for k, v in checks.items():
        print(f"  {k}: {v}")
    assert checks["ordinal_vocab_subset_ABCDE"], f"ordinal vocabulary outside A-E: {vocab}"
    assert not checks["drop_by_rule_observed_in_tracss"], "a non-member entity appeared in TrACSS"
    print(f"wrote {len(long_df):,} rows -> {OUT_DIR}/tracss_raw_long.csv")


def inventory_questions():
    """Count responses and record the original question columns."""


    def to_iso3(raw):
        name = normalise_country(raw)
        if name == "" or name.lower() == "nan":
            return None, name, "blank"
        if name in DROP_BY_RULE:
            return None, name, "drop_by_rule"
        iso3 = NAME_TO_ISO3.get(name)
        if iso3:
            return iso3, name, "mapped"
        return None, name, "UNMAPPED"

    OUT_DIR = "outputs/tracss"
    os.makedirs("outputs/tracss", exist_ok=True)
    TRACSS_DIR = "data/tracss"

    QNUM = re.compile(r"^(\d+(?:\.\d+)*)")

    def question_prefix(header: str) -> str:
        m = QNUM.match(header.strip())
        return m.group(1) if m else ""

    def build_inventory() -> pd.DataFrame:
        d = pd.read_csv(os.path.join(OUT_DIR, "tracss_raw_long.csv"), dtype=str, keep_default_na=False)
        prov = pd.read_csv(os.path.join(OUT_DIR, "tracss_wave_provenance.csv"))
        ncols = dict(zip(prov.wave.astype(str), prov.n_columns))

        rows = []
        for (wave, ci), g in d.groupby(["wave", "column_index"], sort=False):
            vc = Counter(g.raw_value)
            grades = sorted(set(g.grade) - {""})
            hdr = g.raw_header.iloc[0]
            rows.append(
                {
                    "wave": wave,
                    "column_index": int(ci),
                    "raw_header": hdr,
                    "question_number_prefix": question_prefix(hdr),
                    "response_type": g.response_type.iloc[0],
                    "n_nonnull": int(len(g)),
                    "n_unique": int(len(vc)),
                    "n_countries": int(g.iso3.nunique()),
                    "grade_vocabulary": ",".join(grades),
                    "value_vocabulary_top10": json.dumps(
                        [[v[:160], c] for v, c in vc.most_common(10)], ensure_ascii=False
                    ),
                }
            )
        inv = pd.DataFrame(rows).sort_values(
            ["wave", "column_index"], key=lambda s: s.map(lambda x: x) if s.name != "wave" else s
        )

        covered = inv.groupby("wave").column_index.nunique().to_dict()
        print("column coverage (inventoried / total in workbook):")
        for w, _, _, _ in WAVES:
            print(f"  {w}: {covered.get(w, 0)} / {ncols.get(w)}")
        return inv

    def whitespace_diagnostic() -> pd.DataFrame:
        rows = []
        for wave, fn, si, _ in WAVES:
            raw = pd.read_excel(os.path.join(TRACSS_DIR, fn), sheet_name=si, header=None, dtype=object)
            prov = pd.read_csv(os.path.join(OUT_DIR, "tracss_wave_provenance.csv"))
            hr = int(prov.loc[prov.wave.astype(str) == wave, "header_row_detected"].iloc[0])
            header = [norm_ws(v) for v in raw.iloc[hr].tolist()]
            body = raw.iloc[hr + 1:]
            keep = [i for i in range(len(body)) if to_iso3(body.iloc[i, 0])[2] == "mapped"]
            for ci in range(raw.shape[1]):
                vals = [body.iloc[i, ci] for i in keep if norm_ws(body.iloc[i, ci]) != ""]
                if not vals:
                    continue
                rawu = {str(v) for v in vals}
                normu = {norm_ws(v).upper() for v in vals}
                if len(rawu) > len(normu):
                    rows.append(
                        {
                            "wave": wave,
                            "column_index": ci,
                            "raw_header": header[ci],
                            "n_unique_raw": len(rawu),
                            "n_unique_normalised": len(normu),
                            "collapsed_by": len(rawu) - len(normu),
                            "variants_removed": json.dumps(sorted(rawu - normu), ensure_ascii=False),
                        }
                    )
        return pd.DataFrame(rows)

    inv = build_inventory()
    inv.to_csv(os.path.join(OUT_DIR, "tracss_column_inventory.csv"), index=False)
    print(f"wrote {len(inv)} inventory rows")
    print(inv.response_type.value_counts().to_string())

    ws = whitespace_diagnostic()
    ws.to_csv(os.path.join(OUT_DIR, "tracss_whitespace_variants.csv"), index=False)
    print(f"whitespace/case normalisation repaired {len(ws)} columns")


def match_questions():
    """Find laboratory items and compare question wording across waves."""
    # stable >= 0.99 similarity, reworded >= 0.70, redefined < 0.70


    OUT_DIR = "outputs/tracss"
    os.makedirs("outputs/tracss", exist_ok=True)
    TRACSS_DIR = "data/tracss"

    # keep this high - at lower values it chained "...AMR in humans" onto
    # "...AMR in live terrestrial animals"
    FUZZY_THRESHOLD = 0.93

    LAB_PAT = re.compile(r"laborator|bacteriolog|diagnostic|susceptibilit|\bAST\b|\bEQA\b|in-vitro", re.I)

    ANIMAL_PAT = re.compile(
        r"animal health and food safety|in animals|terrestrial|aquatic|veterinar|food safety|"
        r"food (?:and|&) agricultur|livestock|farm",
        re.I,
    )

    HUMAN_PAT = re.compile(r"human health|clinical bacteriolog|patient management|in humans", re.I)

    QNUM_RE = re.compile(r"^(\d+(?:\.\d+)*)")

    def question_prefix(h: str) -> str:
        m = QNUM_RE.match(h.strip())
        return m.group(1) if m else ""

    def parent_headers() -> dict:
        prov = pd.read_csv(os.path.join(OUT_DIR, "tracss_wave_provenance.csv"))
        inv = pd.read_csv(os.path.join(OUT_DIR, "tracss_column_inventory.csv"))
        answered = {(str(r.wave), int(r.column_index)) for r in inv.itertuples()}

        out = {}
        for wave, fn, si, _ in WAVES:
            raw = pd.read_excel(os.path.join(TRACSS_DIR, fn), sheet_name=si, header=None, dtype=object)
            hr = int(prov.loc[prov.wave.astype(str) == wave, "header_row_detected"].iloc[0])
            header = [norm_ws(v) for v in raw.iloc[hr].tolist()]
            own_q = [question_prefix(h) for h in header]

            def in_scope(parent_col: int, parent_q: str, c: int) -> bool:
                if parent_q and own_q[c]:
                    return own_q[c] == parent_q or own_q[c].startswith(parent_q + ".")
                if own_q[c]:
                    return False
                return all(not own_q[k] for k in range(parent_col + 1, c + 1))

            # (1) banner in a row above the header row
            banner = [""] * raw.shape[1]
            for r in range(0, hr):
                for c0 in range(raw.shape[1]):
                    v = norm_ws(raw.iloc[r, c0])
                    if not v:
                        continue
                    pq = question_prefix(v)
                    for c in range(c0, raw.shape[1]):
                        if c > c0 and not in_scope(c0, pq, c):
                            break
                        if not banner[c]:
                            banner[c] = v

            # (2) section title: a header-row cell without answers of its own
            section = [""] * raw.shape[1]
            for c0 in range(raw.shape[1]):
                h = header[c0]
                if not h or (wave, c0) in answered:
                    continue
                pq = question_prefix(h)
                for c in range(c0 + 1, raw.shape[1]):
                    if not in_scope(c0, pq, c):
                        break
                    section[c] = h

            for c in range(raw.shape[1]):
                out[(wave, c)] = {"banner": banner[c], "section": section[c]}
        return out

    def sector_evidence(header: str, parent: dict) -> tuple[str, str, str]:
        own_a, own_h = bool(ANIMAL_PAT.search(header)), bool(HUMAN_PAT.search(header))
        if own_a and not own_h:
            return "animal_food", "high", f"own header: {header[:200]}"
        if own_h and not own_a:
            return "human", "high", f"own header: {header[:200]}"
        if own_a and own_h:
            return "AMBIGUOUS", "low", f"own header mentions BOTH sectors: {header[:200]}"

        for src in ("banner", "section"):
            p = parent.get(src, "")
            if not p:
                continue
            pa, ph = bool(ANIMAL_PAT.search(p)), bool(HUMAN_PAT.search(p))
            if pa and not ph:
                return "animal_food", "medium", f"parent {src}: {p[:200]}"
            if ph and not pa:
                return "human", "medium", f"parent {src}: {p[:200]}"
        return "UNDETERMINED", "none", f"no sector token in header or parent: {header[:200]}"

    def strip_qnum(h: str) -> str:
        t = re.sub(r"^\(?\d+(?:\.\d+)*\s*", "", h.strip())
        t = re.sub(r"^\(?[a-eA-E]\)\s*", "", t)
        t = re.sub(r"^[-.:)]\s*", "", t)
        return t.strip().lower()

    inv = pd.read_csv(os.path.join(OUT_DIR, "tracss_column_inventory.csv"))
    long = pd.read_csv(os.path.join(OUT_DIR, "tracss_raw_long.csv"), dtype=str, keep_default_na=False)
    parents = parent_headers()

    # 1. lab construct census
    lab = inv[inv.raw_header.str.contains(LAB_PAT)].copy()
    rows = []
    for r in lab.itertuples():
        p = parents.get((str(r.wave), int(r.column_index)), {})
        sec, conf, ev = sector_evidence(r.raw_header, p)
        rows.append(
            {
                "wave": r.wave,
                "column_index": r.column_index,
                "question_number_prefix": r.question_number_prefix,
                "raw_header": r.raw_header,
                "response_type": r.response_type,
                "n_nonnull": r.n_nonnull,
                "n_unique": r.n_unique,
                "is_AE_ordinal": r.response_type in ("ordinal_bare", "ordinal_labelled"),
                "sector_candidate": sec,
                "sector_confidence": conf,
                "sector_evidence": ev,
                "parent_banner": p.get("banner", ""),
                "parent_section": p.get("section", ""),
            }
        )
    census = pd.DataFrame(rows).sort_values(
        ["wave", "column_index"], key=lambda s: s.map(WAVE_ORDER.index) if s.name == "wave" else s
    )
    census.to_csv(os.path.join(OUT_DIR, "tracss_lab_item_census.csv"), index=False)
    print("A-E ordinal lab items by wave x sector:")
    print(census[census.is_AE_ordinal].groupby(["wave", "sector_candidate"]).size()
          .unstack(fill_value=0).to_string())

    # 2. candidate crosswalk
    ordinal = inv[inv.response_type.isin(["ordinal_bare", "ordinal_labelled"])].copy()
    ordinal["prompt"] = ordinal.raw_header.map(strip_qnum)
    groups = defaultdict(list)
    for r in ordinal.itertuples():
        groups[r.prompt].append(r)

    # exact-prompt groups first (largest first), then fuzzy-link remaining prompts
    keys = list(groups)
    assigned, item_of = {}, {}
    for k in sorted(keys, key=lambda x: -len(groups[x])):
        if k in assigned:
            continue
        item_id = f"ITEM_{len(item_of):03d}"
        item_of[item_id] = [k]
        assigned[k] = item_id
        for k2 in keys:
            if k2 in assigned or not k2 or not k:
                continue
            sim = difflib.SequenceMatcher(None, k, k2).ratio()
            if sim >= FUZZY_THRESHOLD:
                assigned[k2] = item_id
                item_of[item_id].append(k2)

    cw_rows = []
    for prompt, members in groups.items():
        iid = assigned[prompt]
        for r in members:
            p = parents.get((str(r.wave), int(r.column_index)), {})
            sec, sconf, ev = sector_evidence(r.raw_header, p)
            cw_rows.append(
                {
                    "item_id": iid,
                    "construct_prompt_normalised": prompt,
                    "wave": r.wave,
                    "column_index": r.column_index,
                    "question_number_prefix": r.question_number_prefix,
                    "raw_header": r.raw_header,
                    "n_nonnull": r.n_nonnull,
                    "sector_candidate": sec,
                    "sector_confidence": sconf,
                    "sector_evidence": ev,
                    "match_confidence": "exact" if len(item_of[iid]) == 1 else "semantic",
                    "REVIEW_REQUIRED": True,
                }
            )
    cw = pd.DataFrame(cw_rows)
    nwaves = cw.groupby("item_id").wave.nunique().rename("n_waves_spanned")
    cw = cw.merge(nwaves, on="item_id")

    # flag groups that look over-merged: >1 column in one wave, or mixed sectors
    dupe_wave = cw.groupby(["item_id", "wave"]).size().rename("n_cols_in_wave").reset_index()
    over = set(dupe_wave.loc[dupe_wave.n_cols_in_wave > 1, "item_id"])
    mixed = {
        i
        for i, g in cw.groupby("item_id")
        if len({s for s in g.sector_candidate if s not in ("UNDETERMINED", "AMBIGUOUS")}) > 1
    }
    cw["multiple_columns_in_a_wave"] = cw.item_id.isin(over)
    cw["mixed_sector_candidates"] = cw.item_id.isin(mixed)
    cw.loc[cw.item_id.isin(over | mixed), "match_confidence"] = "uncertain"
    cw = cw.merge(dupe_wave, on=["item_id", "wave"], how="left")
    cw["wave_order"] = cw.wave.map(WAVE_ORDER.index)
    cw = cw.sort_values(["n_waves_spanned", "item_id", "wave_order"], ascending=[False, True, True])
    cw.drop(columns=["wave_order"]).to_csv(
        os.path.join(OUT_DIR, "tracss_candidate_crosswalk.csv"), index=False
    )
    print(f"candidate crosswalk: {cw.item_id.nunique()} items over {len(cw)} wave-columns")

    # 3. ladder drift
    ord_long = long[long.response_type.isin(["ordinal_bare", "ordinal_labelled"])]
    opts = (
        ord_long[ord_long.option_text != ""]
        .groupby(["wave", "column_index", "grade"])
        .option_text.agg(lambda s: s.mode().iloc[0])
        .reset_index()
    )
    key = {(str(r.wave), int(r.column_index)): r.item_id for r in cw.itertuples()}
    opts["item_id"] = [key.get((str(w), int(c))) for w, c in zip(opts.wave, opts.column_index)]
    opts = opts[opts.item_id.notna()]

    drift = []
    for iid, g in opts.groupby("item_id"):
        waves = sorted(g.wave.unique(), key=WAVE_ORDER.index)
        if len(waves) < 2:
            continue
        base = waves[0]
        bmap = dict(zip(g[g.wave == base].grade, g[g.wave == base].option_text))
        for w in waves[1:]:
            wmap = dict(zip(g[g.wave == w].grade, g[g.wave == w].option_text))
            sims, detail = [], {}
            for gr in "ABCDE":
                if gr in bmap and gr in wmap:
                    s = difflib.SequenceMatcher(None, bmap[gr].lower(), wmap[gr].lower()).ratio()
                    sims.append(s)
                    detail[gr] = round(s, 3)
            if not sims:
                continue
            mn, mean = min(sims), sum(sims) / len(sims)
            verdict = ("stable" if mn >= 0.99
                       else "reworded_same_meaning" if mn >= 0.70
                       else "redefined")
            hdr = cw[(cw.item_id == iid) & (cw.wave == w)].raw_header
            drift.append(
                {
                    "item_id": iid,
                    "base_wave": base,
                    "compare_wave": w,
                    "raw_header_compare_wave": hdr.iloc[0] if len(hdr) else "",
                    "n_grades_compared": len(sims),
                    "min_similarity": round(mn, 3),
                    "mean_similarity": round(mean, 3),
                    "per_grade_similarity": json.dumps(detail),
                    "ladder_stability": verdict,
                    "eligible_for_longitudinal": verdict in ("stable", "reworded_same_meaning"),
                }
            )
    dr = pd.DataFrame(drift)
    dr.to_csv(os.path.join(OUT_DIR, "tracss_ladder_drift.csv"), index=False)
    print(f"ladder drift: {len(dr)} matched item comparisons")
    print(dr.ladder_stability.value_counts().to_string())


def finish_crosswalk():
    """Split ambiguous matches, assign sectors and build the lab composite."""
    # 2020-21 human lab composite


    OUT_DIR = "outputs/tracss"
    os.makedirs("outputs/tracss", exist_ok=True)
    GRADE_ORD = {"A": 1, "B": 2, "C": 3, "D": 4, "E": 5}

    def check(cond: bool, msg: str):
        if not cond:
            raise AssertionError(msg)

    def load():
        long = pd.read_csv(os.path.join(OUT_DIR, "tracss_raw_long.csv"), dtype=str, keep_default_na=False)
        cand = pd.read_csv(os.path.join(OUT_DIR, "tracss_candidate_crosswalk.csv"),
                           dtype=str, keep_default_na=False)
        check(len(long) == 139975, f"tracss_raw_long: expected 139975 rows, got {len(long)}")
        long["column_index"] = long["column_index"].astype(int)
        cand["column_index"] = cand["column_index"].astype(int)
        return long, cand

    # split over-merged candidate items
    def strip_qnum(h: str) -> str:
        t = re.sub(r"^\(?\d+(?:\.\d+)*\s*", "", h.strip())
        t = re.sub(r"^\(?[a-eA-E]\)\s*", "", t)
        t = re.sub(r"^[-.:)]\s*", "", t)
        return t.strip().lower()

    def split_key(h: str) -> str:
        t = strip_qnum(h)
        t = re.sub(r"[\*\+†‡§]", "", t)
        t = re.sub(r"\{[^}]*\}", "", t)
        return re.sub(r"\s+", " ", t).strip()

    def apply_splits(cand: pd.DataFrame):
        cand = cand.copy()
        cand["prompt_key"] = cand["raw_header"].map(split_key)
        to_split = sorted(cand.loc[cand.match_confidence == "uncertain", "item_id"].unique())
        log = []
        new_item_id = cand["item_id"].copy()
        for iid in to_split:
            sub = cand[cand.item_id == iid]
            distinct = sorted(sub.prompt_key.unique())
            if len(distinct) <= 1:
                log.append({"orig_item_id": iid, "action": "kept_no_distinct_construct",
                            "n_new_items": 1, "detail": distinct[0][:120] if distinct else ""})
                continue
            for k, pk in enumerate(distinct, 1):
                sub_idx = sub.index[sub.prompt_key == pk]
                sid = f"{iid}_s{k}"
                new_item_id.loc[sub_idx] = sid
                waves = ",".join(sorted(cand.loc[sub_idx, "wave"], key=WAVE_ORDER.index))
                log.append({"orig_item_id": iid, "action": "split", "new_item_id": sid,
                            "n_new_items": len(distinct), "waves": waves, "construct": pk[:160]})
        cand["item_id"] = new_item_id
        return cand.drop(columns=["prompt_key"]), pd.DataFrame(log)

    # sector assignment
    ANIMAL_PAT = re.compile(r"animal|terrestrial|aquatic|veterinar|livestock|woah|\boie\b", re.I)

    FOOD_PAT = re.compile(r"food safety|food production|food processing|in food|foods|infarm|"
                          r"and food|manufactur", re.I)

    PLANT_ENV_PAT = re.compile(r"plant|pesticid|environment|wastewater|bactericid|fungicid|crop", re.I)

    HUMAN_PAT = re.compile(r"human health|clinical bacteriolog|patient management|in humans|"
                           r"hospital|primary care|human$", re.I)

    CROSS_STRONG_PAT = re.compile(r"multi-sector|one health|national action plan|behaviour change", re.I)

    CROSS_WEAK_PAT = re.compile(r"raising awareness|awareness and understanding", re.I)

    # One Health coordination and national action plan items are cross-sectoral by design
    # (ITEM_036 is the 2016-17 predecessor of the One Health coordination item).
    FORCE_CROSS = {"ITEM_007", "ITEM_005", "ITEM_036"}
    CROSS_NOTE = "Cross-sectoral item (applies to all sectors)."

    def semantic_sector(header: str) -> str:
        h = header
        if CROSS_STRONG_PAT.search(h):
            return "cross"
        if HUMAN_PAT.search(h) and not ANIMAL_PAT.search(h):
            return "human"
        a, f, p = bool(ANIMAL_PAT.search(h)), bool(FOOD_PAT.search(h)), bool(PLANT_ENV_PAT.search(h))
        if a and f:
            return "animal_food"
        if a:
            return "animal"
        if f:
            return "food"
        if p:
            return "environment"
        if HUMAN_PAT.search(h):
            return "human"
        if CROSS_WEAK_PAT.search(h):
            return "cross"
        return "undetermined"

    def assign_sector(row) -> dict:
        iid = row.item_id.split("_s")[0]        # split children inherit the parent rule
        conf = row.sector_confidence            # high | medium | none | low
        cand_sec = row.sector_candidate         # human | animal_food | UNDETERMINED | AMBIGUOUS

        if iid in FORCE_CROSS:
            return dict(sector="cross", sector_status="fixed_cross",
                        sector_review_required=False, sector_note=CROSS_NOTE,
                        sector_evidence=row.sector_evidence)

        if conf == "high":
            # sector stated in the item's own header
            if cand_sec == "human":
                sec = "human"
            elif cand_sec == "animal_food":
                sec = semantic_sector(row.raw_header)
                if sec in ("undetermined", "human", "cross"):
                    sec = "animal_food"
            else:
                sec = semantic_sector(row.raw_header)
            return dict(sector=sec, sector_status="auto_accepted_high",
                        sector_review_required=False, sector_note="",
                        sector_evidence=row.sector_evidence)

        if conf == "medium":
            # sector found only in a parent header: use it, but flag for review
            rec_sec = cand_sec if cand_sec not in ("UNDETERMINED", "AMBIGUOUS") else "animal_food"
            return dict(sector=rec_sec, sector_status="medium_FLAGGED_for_review",
                        sector_review_required=True,
                        sector_note=("Medium confidence: sector inferred from the parent section "
                                     f"header, not the item's own header. Evidence: {row.sector_evidence}"),
                        sector_evidence=row.sector_evidence)

        sec = semantic_sector(row.raw_header)
        return dict(sector=sec, sector_status="none_best_effort_FLAGGED",
                    sector_review_required=True,
                    sector_note=("No sector wording in header or parent; semantic reading = "
                                 f"'{sec}'; flagged for review."),
                    sector_evidence=row.sector_evidence)

    # ladder stability per item-wave at several thresholds
    THRESHOLDS = {  # name -> (reworded_min, stable_min); below reworded_min = redefined
        "primary_0.70": (0.70, 0.99),
        "lenient_0.55": (0.55, 0.99),
        "strict_0.80": (0.80, 0.99),
    }

    def classify(minsim: float, reworded_min: float, stable_min: float) -> str:
        if minsim >= stable_min:
            return "stable"
        if minsim >= reworded_min:
            return "reworded_same_meaning"
        return "redefined"

    def ladder_perwave(cand: pd.DataFrame, opts: dict):
        rows = []
        for iid, g in cand.groupby("item_id"):
            waves = sorted(g.wave.unique(), key=WAVE_ORDER.index)
            with_opts = [w for w in waves if opts.get((iid, w))]
            base = with_opts[0] if with_opts else None
            bmap = opts.get((iid, base), {}) if base else {}
            for w in waves:
                rec = {"item_id": iid, "wave": w, "ladder_min_similarity": ""}
                if base is None or w == base or not opts.get((iid, w)):
                    for tname in THRESHOLDS:
                        rec[f"ladder_{tname}"] = "stable"
                else:
                    wmap = opts[(iid, w)]
                    sims = [difflib.SequenceMatcher(None, bmap[gr].lower(), wmap[gr].lower()).ratio()
                            for gr in "ABCDE" if gr in bmap and gr in wmap]
                    if not sims:
                        for tname in THRESHOLDS:
                            rec[f"ladder_{tname}"] = "stable"
                    else:
                        mn = min(sims)
                        rec["ladder_min_similarity"] = round(mn, 3)
                        for tname, (rmin, smin) in THRESHOLDS.items():
                            rec[f"ladder_{tname}"] = classify(mn, rmin, smin)
                rows.append(rec)
        return pd.DataFrame(rows)

    def longitudinal_eligible(perwave: pd.DataFrame, cand: pd.DataFrame, tname: str):
        """Eligible if the item spans >= 2 waves and no wave is 'redefined'."""
        nwaves = cand.groupby("item_id").wave.nunique()
        col = f"ladder_{tname}"
        return {iid: bool(nwaves[iid] >= 2 and (g[col] != "redefined").all())
                for iid, g in perwave.groupby("item_id")}

    def option_texts(long: pd.DataFrame, cand: pd.DataFrame):
        ord_long = long[long.response_type.isin(["ordinal_bare", "ordinal_labelled"])]
        ord_long = ord_long[ord_long.option_text != ""]
        key = {(r.wave, r.column_index): r.item_id for r in cand.itertuples()}
        opt = {}
        for r in ord_long.itertuples():
            iid = key.get((r.wave, r.column_index))
            if iid is None:
                continue
            opt.setdefault((iid, r.wave), {}).setdefault(r.grade, {})
            opt[(iid, r.wave)][r.grade][r.option_text] = opt[(iid, r.wave)][r.grade].get(r.option_text, 0) + 1
        return {k: {g: (max(v, key=v.get) if v else "") for g, v in grades.items()}
                for k, grades in opt.items()}

    def build_crosswalk(long, cand):
        cand, split_log = apply_splits(cand)

        sec_recs = cand.apply(assign_sector, axis=1, result_type="expand")
        cand = pd.concat([cand.reset_index(drop=True), sec_recs.reset_index(drop=True)], axis=1)

        labels = {}
        for iid, g in cand.groupby("item_id"):
            g2 = g.sort_values("wave", key=lambda s: s.map(WAVE_ORDER.index))
            labels[iid] = strip_qnum(g2.raw_header.iloc[0])[:180]
        cand["construct_label"] = cand.item_id.map(labels)

        # split items are single-question by construction
        conf = cand.match_confidence.copy()
        conf[cand.item_id.str.contains("_s")] = "semantic"
        cand["confidence"] = conf

        opts = option_texts(long, cand)
        for L in "ABCDE":
            cand[f"option_text_{L}"] = [opts.get((r.item_id, r.wave), {}).get(L, "")
                                        for r in cand.itertuples()]

        pw = ladder_perwave(cand, opts)
        cand = cand.merge(pw, on=["item_id", "wave"], how="left")

        elig_by_thr = {t: longitudinal_eligible(pw, cand, t) for t in THRESHOLDS}
        for t in THRESHOLDS:
            cand[f"longitudinal_eligible_{t}"] = cand.item_id.map(elig_by_thr[t])

        def make_note(r):
            parts = [r.sector_note] if r.sector_note else []
            if "_s" in r.item_id:
                parts.append(f"Split from {r.item_id.split('_s')[0]} (different question text).")
            if r.item_id.split("_s")[0] in FORCE_CROSS and CROSS_NOTE not in parts:
                parts.append(CROSS_NOTE)
            return " | ".join(parts)
        cand["notes"] = cand.apply(make_note, axis=1)
        cand["ladder_stability"] = cand["ladder_primary_0.70"]

        cand["wave_order"] = cand.wave.map(WAVE_ORDER.index)
        cand = cand.sort_values(["item_id", "wave_order"]).reset_index(drop=True)

        mandatory = ["item_id", "construct_label", "sector", "wave", "raw_header", "column_index",
                     "option_text_A", "option_text_B", "option_text_C", "option_text_D", "option_text_E",
                     "confidence", "ladder_stability", "notes"]
        extra = ["question_number_prefix", "n_nonnull", "sector_status", "sector_confidence",
                 "sector_review_required", "sector_evidence", "ladder_min_similarity",
                 "ladder_primary_0.70", "ladder_lenient_0.55", "ladder_strict_0.80",
                 "longitudinal_eligible_primary_0.70", "longitudinal_eligible_lenient_0.55",
                 "longitudinal_eligible_strict_0.80", "n_waves_spanned"]
        crosswalk = cand[mandatory + [c for c in extra if c in cand.columns]].copy()
        return crosswalk, cand, split_log, elig_by_thr

    def build_analysis_long(long, cand):
        key = {(r.wave, r.column_index): r for r in cand.itertuples()}
        ord_long = long[long.response_type.isin(["ordinal_bare", "ordinal_labelled"])].copy()
        rows, unparsed = [], 0
        for r in ord_long.itertuples():
            meta = key.get((r.wave, r.column_index))
            if meta is None:
                continue
            base = {"iso3": r.iso3, "wave": r.wave, "item_id": meta.item_id, "sector": meta.sector,
                    "ladder_stability": meta.ladder_stability, "confidence": meta.confidence,
                    "raw_value": r.raw_value}
            if r.grade not in GRADE_ORD:
                if r.grade == "" and r.raw_value.strip() == "":
                    continue
                unparsed += 1
                rows.append({**base, "grade_ordinal": "", "parse_status": "unparsed"})
                continue
            rows.append({**base, "grade_ordinal": GRADE_ORD[r.grade], "parse_status": "parsed"})
        cols = ["iso3", "wave", "item_id", "sector", "grade_ordinal", "ladder_stability",
                "confidence", "raw_value", "parse_status"]
        return pd.DataFrame(rows)[cols], unparsed

    # 2020-21 human lab composite, from these columns:
    #   column 49: national AMR laboratory network membership (Yes/No)
    #   column 50: national accreditation authority (Yes/No)
    #   column 51: accreditation to international standards (Yes/No; asked if col 50 = Yes)
    #   column 56: diagnostic (bacteriology) techniques, A-E
    def build_lab_composite(long):
        w = long[long.wave == "2020-21"].copy()

        def col(ci):
            return w[w.column_index == ci].set_index("iso3")["raw_value"]

        c49, c50, c51, c56 = col(49), col(50), col(51), col(56)
        isos = sorted(set(c49.index) | set(c50.index) | set(c51.index) | set(c56.index))
        DK = {"Don't know", "Don’t know", "Dont know"}

        def yn(v):
            v = (v or "").strip()
            if v in DK:
                return ("missing_dk", None)
            if v == "":
                return ("blank", None)
            if v.lower() == "yes":
                return ("Yes", 1)
            if v.lower() == "no":
                return ("No", 0)
            return ("unparsed:" + v, None)

        rows = []
        for iso in isos:
            v49, v50, v51, v56 = c49.get(iso, ""), c50.get(iso, ""), c51.get(iso, ""), c56.get(iso, "")
            s49, net = yn(v49)
            s50, _ = yn(v50)
            s51, _ = yn(v51)

            # accreditation: 0 = no authority, 1 = authority, 2 = authority + international
            if s50 == "No":
                acc, acc_status = 0, "observed_no_authority"
            elif s50 == "missing_dk":
                acc, acc_status = None, "genuine_missing_dk_col50"
            elif s50 == "blank":
                acc, acc_status = None, "genuine_missing_blank_col50"
            elif s50 == "Yes":
                if s51 == "Yes":
                    acc, acc_status = 2, "observed_authority_intl"
                elif s51 == "No":
                    acc, acc_status = 1, "observed_authority_not_intl"
                elif s51 == "missing_dk":
                    acc, acc_status = None, "genuine_missing_dk_col51"
                elif s51 == "blank":
                    acc, acc_status = None, "genuine_missing_blank_col51_where_col50_yes"
                else:
                    acc, acc_status = None, "unparsed_col51"
            else:
                acc, acc_status = None, "unparsed_col50"

            m = re.match(r"^([A-E])\b", v56.strip())
            diag = GRADE_ORD[m.group(1)] if m else None
            diag_status = "observed" if diag is not None else ("blank" if v56.strip() == "" else "unparsed")

            n_acc = (acc / 2) if acc is not None else None
            n_diag = ((diag - 1) / 4) if diag is not None else None
            comps = [net, n_acc, n_diag]
            composite = (sum(comps) / 3) if all(x is not None for x in comps) else None

            rows.append({
                "iso3": iso,
                "lab_network_membership": net, "network_status": s49, "col49_raw": v49,
                "lab_accreditation": acc, "accreditation_status": acc_status,
                "col50_raw": v50, "col51_raw": v51,
                "col51_structural_absent": s50 in ("No", "missing_dk", "blank"),
                "lab_diagnostic_techniques": diag, "diagnostic_status": diag_status,
                "col56_grade": (m.group(1) if m else ""),
                "n_network_norm": net, "acc_norm": n_acc, "diag_norm": n_diag,
                "lab_composite_equalweight": composite,
            })
        return pd.DataFrame(rows)

    def lab_intercorrelations(comp: pd.DataFrame):
        cols = {"network": "lab_network_membership",
                "accreditation": "lab_accreditation",
                "diagnostic": "lab_diagnostic_techniques"}
        rows = []
        for a, b in itertools.combinations(cols, 2):
            sub = comp[[cols[a], cols[b]]].dropna()
            n = len(sub)
            if n >= 3 and sub[cols[a]].nunique() > 1 and sub[cols[b]].nunique() > 1:
                sp = sub[cols[a]].rank().corr(sub[cols[b]].rank())
            else:
                sp = None
            rows.append({"pair": f"{a}__{b}", "n_pairwise_complete": n,
                         "spearman_rho": round(sp, 3) if sp is not None else None})
        return pd.DataFrame(rows)

    long, cand = load()
    crosswalk, full, split_log, elig = build_crosswalk(long, cand)

    check(crosswalk.sector.notna().all() and (crosswalk.sector != "").all(), "sector missing")
    check(crosswalk.confidence.isin(["exact", "semantic", "uncertain"]).all(), "bad confidence value")
    check(crosswalk.ladder_stability.isin(["stable", "reworded_same_meaning", "redefined"]).all(),
          "bad ladder_stability value")

    cw_path = os.path.join(OUT_DIR, "tracss_item_crosswalk.csv")
    crosswalk.to_csv(cw_path, index=False)
    sha = hashlib.sha256(open(cw_path, "rb").read()).hexdigest()
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    with open(os.path.join(OUT_DIR, "tracss_item_crosswalk_manifest.json"), "w") as fh:
        json.dump({"artifact": "tracss_item_crosswalk.csv", "sha256": sha, "created_utc": ts,
                   "n_rows": int(len(crosswalk)), "n_items": int(crosswalk.item_id.nunique())},
                  fh, indent=2)

    al, unparsed = build_analysis_long(long, full)
    al.to_csv(os.path.join(OUT_DIR, "tracss_analysis_long.csv"), index=False)

    split_log.to_csv(os.path.join(OUT_DIR, "tracss_crosswalk_splits.csv"), index=False)
    crosswalk[crosswalk.sector_status == "medium_FLAGGED_for_review"][
        ["item_id", "wave", "sector", "raw_header", "sector_evidence"]
    ].to_csv(os.path.join(OUT_DIR, "tracss_sector_medium_confidence_review.csv"), index=False)

    comp = build_lab_composite(long)
    comp.to_csv(os.path.join(OUT_DIR, "tracss_human_lab_composite_2020_21.csv"), index=False)
    lab_intercorrelations(comp).to_csv(
        os.path.join(OUT_DIR, "tracss_lab_component_intercorrelations.csv"), index=False)

    nwaves = full.groupby("item_id").wave.nunique()
    sens_rows = []
    for t in THRESHOLDS:
        e = elig[t]
        multi = [i for i in e if nwaves[i] >= 2]
        n_long = sum(e[i] for i in multi)
        sens_rows.append({"threshold": t, "reworded_min": THRESHOLDS[t][0],
                          "stable_min": THRESHOLDS[t][1], "n_multiwave_items": len(multi),
                          "n_longitudinal_eligible": n_long,
                          "n_cross_sectional_only": len(multi) - n_long})
    sens = pd.DataFrame(sens_rows)
    sens.to_csv(os.path.join(OUT_DIR, "tracss_ladder_threshold_sensitivity.csv"), index=False)

    print(f"crosswalk: {len(crosswalk)} item-wave rows, {crosswalk.item_id.nunique()} items")
    print(f"analysis_long: {len(al)} rows ({unparsed} unparsed)")
    print(sens.to_string(index=False))
    print(f"lab composite (complete case): {comp.lab_composite_equalweight.notna().sum()} countries")


def main():
    os.chdir(ROOT)
    extract_responses()
    inventory_questions()
    match_questions()
    finish_crosswalk()


if __name__ == "__main__":
    main()
