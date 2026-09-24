"""Summarise laboratory coverage and changing response ladders for Figure 1."""
# Copyright 2026 Tavpritesh Sethi and Jasmine Kaur (Tavlab, IIITD)
# Licensed under the Apache License, Version 2.0 (see LICENSE)

import os
import pandas as pd
from common import ROOT
from preprocess_tracss import WAVE_ORDER


def summarise_audit():
    """Tabulate the laboratory item counts, ladder comparisons and NAP wording."""


    OUT = "outputs/tracss"
    os.makedirs("outputs/tracss", exist_ok=True)

    census = pd.read_csv(os.path.join(OUT, "tracss_lab_item_census.csv"))
    lab = (census[census.is_AE_ordinal]
           .groupby(["wave", "sector_candidate"]).size().unstack(fill_value=0)
           .reindex(WAVE_ORDER).fillna(0).astype(int))
    lab.to_csv(os.path.join(OUT, "audit_lab_items_by_wave_sector.csv"))
    print("Gradeable A-E laboratory items per wave:")
    print(lab.to_string())

    drift = pd.read_csv(os.path.join(OUT, "tracss_ladder_drift.csv"))
    ladder = drift.ladder_stability.value_counts().rename_axis("ladder_stability").reset_index(name="n")
    ladder.loc[len(ladder)] = ["total", int(ladder.n.sum())]
    ladder.to_csv(os.path.join(OUT, "audit_ladder_comparisons.csv"), index=False)
    print("\nLadder comparisons:")
    print(ladder.to_string(index=False))

    cw = pd.read_csv(os.path.join(OUT, "tracss_item_crosswalk.csv"), keep_default_na=False)
    nap = cw[cw.item_id == "ITEM_005"][["wave", "raw_header", "option_text_E", "ladder_stability"]]
    nap.to_csv(os.path.join(OUT, "audit_nap_top_grade_by_wave.csv"), index=False)
    print("\nNational action plan item, grade E wording by wave:")
    for r in nap.itertuples():
        print(f"  {r.wave}: {r.option_text_E[:110]}")


def main():
    os.chdir(ROOT)
    summarise_audit()


if __name__ == "__main__":
    main()
