# PolarisAMR

Python code for the Vivli AMR Data Challenge 2026 submission on resistance trajectories, policy readiness and hosted AMR R&D investment.

The analysis combines Pfizer ATLAS surveillance data and supplied AMROrbit estimates with TrACSS survey responses, Global AMR R&D Hub investment data and World Bank indicators. All eleven Python files are stored directly in the repository root

| File | Purpose |
| --- | --- |
| [preprocess_tracss.py](preprocess_tracss.py) | Cleans nine TrACSS survey waves, matches questions across years and builds the laboratory composite. |
| [preprocess_amrorbit.py](preprocess_amrorbit.py) | Reconstructs resistance levels and slopes from the supplied AMROrbit model estimates. |
| [preprocess_atlas.py](preprocess_atlas.py) | Calculates isolate counts, resistance summaries and country-level testing volume. |
| [preprocess_hub.py](preprocess_hub.py) | Cleans R&D project records and totals investment by host and funder country. |
| [preprocess_world_bank.py](preprocess_world_bank.py) | Prepares GDP, health expenditure, universal health coverage and population measures. |
| [merge_data.py](merge_data.py) | Selects the analysis combinations and joins the prepared datasets. |
| [analysis_readiness.py](analysis_readiness.py) | Screens readiness measures using Firth logistic regression and country bootstrap inference; compares organism groups and examines testing-volume trends. |
| [analysis_tracss.py](analysis_tracss.py) | Summarises laboratory-question coverage and changes in response ladders. |
| [analysis_rd.py](analysis_rd.py) | Examines investment–burden correlations, burden terciles and adjustment for GDP per capita. |
| [make_figures.py](make_figures.py) | Generates the TrACSS audit and investment–burden figures. |
| [common.py](common.py) | Provides shared country-name handling, organism mappings and rolling-window definitions. |

The workflow proceeds from preprocessing to merging, then the three analyses and figure generation. `common.py` supports the other scripts.

This repository contains the analysis code. Raw datasets, supplied model estimates, the country lookup table and generated outputs are not included. File paths refer to the original local analysis inputs. The AMROrbit preprocessing script uses saved model estimates rather than fitting the original models from scratch.
