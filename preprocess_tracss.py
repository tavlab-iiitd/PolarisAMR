# preprocess tracss lab items by sector, match questions across waves, compare A-E wording
# stable >= 0.99 similarity, reworded >= 0.70, redefined < 0.70
# Copyright 2026 Tavpritesh Sethi and Jasmine Kaur (Tavlab, IIITD)
# Licensed under the Apache License, Version 2.0 (see LICENSE)

from collections import defaultdict
import difflib
import json
import os
import re
import unicodedata

import pandas as pd


# wave, file, sheet, expected header row (only used as a check - header is detected)
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


# whitespace / text normalisation
_NBSP = "\u00a0\u2007\u202f\u200b\ufeff"


def norm_ws(s) -> str:
    if s is None:
        return ""
    if isinstance(s, float) and s != s:  # NaN
        return ""
    if s is getattr(__import__("pandas"), "NaT", object()):
        return ""
    t = str(s)
    t = unicodedata.normalize("NFKC", t)
    for ch in _NBSP:
        t = t.replace(ch, " ")
    t = t.replace("\u2013", "-").replace("\u2014", "-").replace("\u2212", "-")
    t = t.replace("\r", " ").replace("\n", " ").replace("\t", " ")
    t = re.sub(r"\s+", " ", t)
    return t.strip()


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
