"""Prepare GDP, health expenditure, UHC and population from the World Bank extracts."""
# Copyright 2026 Tavpritesh Sethi and Jasmine Kaur (Tavlab, IIITD)
# Licensed under the Apache License, Version 2.0 (see LICENSE)

import csv
import os
import pandas as pd
from common import ROOT
from pathlib import Path


YEARS = range(2014, 2021)
WDI = "data/world_bank/wdi/9fb3b7d2-e3cd-4fe3-bc94-30b4b56aa084_Data.csv"
HNP = "data/world_bank/hnp/75c52d32-65e4-482c-a9d6-4d0b8985ac57_Data.csv"

def latest_values(path, series_code):
    values = {}
    with open(path, encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            if row.get("Series Code") != series_code:
                continue
            for year in YEARS:
                raw = (row.get(f"{year} [YR{year}]") or "").strip()
                if raw not in ("", ".."):
                    values[row["Country Code"]] = (float(raw), year)
    return values

def main():
    os.chdir(ROOT)
    out = Path("outputs/world_bank")
    out.mkdir(parents=True, exist_ok=True)
    series = {
        "gdp_pc_ppp": latest_values(WDI, "NY.GDP.PCAP.PP.CD"),
        "che_pc_usd": latest_values(WDI, "SH.XPD.CHEX.PC.CD"),
        "uhc_index": latest_values(HNP, "SH.UHC.SRVS.CV.XD"),
    }
    countries = sorted(set().union(*(set(v) for v in series.values())))
    rows = []
    for iso in countries:
        row = {"iso3": iso}
        for name, values in series.items():
            row[name], row[name + "_year"] = values.get(iso, (None, None))
        rows.append(row)
    pd.DataFrame(rows).to_csv(out / "covariates.csv", index=False)
    raw = pd.read_csv(HNP)
    pop = raw.loc[raw["Series Code"] == "SP.POP.TOTL", ["Country Code", "2020 [YR2020]"]].copy()
    pop.columns = ["iso3", "population_2020"]
    pop["population_2020"] = pd.to_numeric(pop["population_2020"], errors="coerce")
    assert not pop.iso3.duplicated().any(), "duplicate population country"
    pop.to_csv(out / "population_2020.csv", index=False)
    print(f"Covariates: {len(rows)} countries; population: {pop.population_2020.notna().sum()} countries")

if __name__ == "__main__":
    main()
