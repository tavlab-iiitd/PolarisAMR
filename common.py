"""Country names and ATLAS definitions used by more than one script."""
# Copyright 2026 Tavpritesh Sethi and Jasmine Kaur (Tavlab, IIITD)
# Licensed under the Apache License, Version 2.0 (see LICENSE)

import csv
import html
import re
import unicodedata
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

NAME_TO_ISO3, ATLAS_ALIASES, DROP_BY_RULE = {}, {}, set()
with open(ROOT / "lookups/country_iso3.csv", encoding="utf-8") as fh:
    for r in csv.DictReader(fh):
        if r["table"] == "tracss":
            NAME_TO_ISO3[r["name"]] = r["iso3"]
        elif r["table"] == "atlas":
            ATLAS_ALIASES[r["name"]] = r["iso3"]
        else:
            DROP_BY_RULE.add(r["name"])


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

def normalise_country(raw):
    return norm_ws(html.unescape(norm_ws(raw)))

WINDOWS = [(2014, 2017), (2015, 2018), (2016, 2019),
           (2017, 2020), (2018, 2021), (2019, 2022)]

WINDOW_LABELS = [f"{a}-{b}" for a, b in WINDOWS]

def parse_country_list(s):
    if s is None or s.strip() == "":
        return []
    parts = [p.strip() for p in s.split(",")]
    out, i = [], 0
    while i < len(parts):
        if parts[i] == "Korea" and i + 1 < len(parts) and parts[i + 1] in ("South", "North"):
            out.append(f"Korea, {parts[i + 1]}")
            i += 2
        else:
            out.append(parts[i])
            i += 1
    return [p for p in out if p]

def organism_groups(species: str, oxacillin_interp: str) -> list[str]:
    out = []
    simple = {
        "Acinetobacter baumannii", "Escherichia coli", "Enterococcus faecium",
        "Klebsiella pneumoniae", "Klebsiella aerogenes", "Pseudomonas aeruginosa",
    }
    if species in simple:
        out.append(species)
    elif species == "Staphylococcus aureus":
        if oxacillin_interp == "Resistant":
            out.append("Staphylococcus aureus (MRSA)")
        elif oxacillin_interp == "Susceptible":
            out.append("Staphylococcus aureus (MSSA)")
    elif species.startswith("Enterobacter"):
        out.append("Enterobacter species with cloacae")
        if species == "Enterobacter cloacae":
            out.append("Enterobacter cloacae")
        else:
            out.append("Enterobacter species without cloacae")
    return out

def atlas_column(antibiotic: str) -> str:
    return antibiotic.replace(".", " ") + "_I"
