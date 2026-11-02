# Your simple association project

Start with `correlation.py`. It retrieves one gene from one cBioPortal study, matches expression and copy-number samples, selects one sample per patient, adds patient clinical fields, and computes two introductory associations.

## Run

```powershell
python -m pip install requests pandas scipy
python correlation.py
```

CSV files open in Excel without openpyxl. The default is PTDSS1 in the pancreatic cancer PanCancer Atlas study, with up to 30 eligible patients. This is a convenience subset for learning, not a representative research sample. Requests may take time; progress appears in the terminal.

## Read the outputs

All files are in `correlation_output` beside the script:

| File | Purpose |
|---|---|
| `profiles.csv` | Available molecular datasets; inspect descriptions and units |
| `clinical_attributes.csv` | Available clinical field definitions |
| `matched_patients.csv` | One selected tumor specimen per patient, with gene and clinical data |
| `raw_data.json` | Original retrieved records for checking transformations |
| `missing_values.csv` | Missing entries per column in the matched table |
| `settings.json` | Study, profiles, gene, retrieval time, and selection limit |
| `analysis.txt` | Expression-age and expression-CNA Spearman results |

The script deliberately fails if it cannot select an unambiguous profile. Inspect `profiles.csv`, then set `EXPRESSION_PROFILE` to an exact expression profile ID. For ambiguous CNA/mutation profiles, explicitly select the ID in `main()` after reviewing the listed descriptions.

## Understand the code in six steps

1. `api()` sends one request. GET reads metadata; POST submits gene/sample filters.
2. `pages()` repeats requests until it has all returned records.
3. `main()` discovers IDs and intersects expression/CNA sample eligibility.
4. It picks the alphabetically first eligible sample for each patient. Review sample types before research use; this rule does not choose baseline or pretreatment tissue.
5. It joins measurements by sample ID and clinical data by study/patient IDs. No response-order matching occurs.
6. `run_analysis()` drops missing pairs separately for each question and reports the actual analysis size.

## Interpret the first analysis

Spearman rho describes rank association: positive means larger values tend to occur together; negative means they tend to move in opposite directions. A p-value is not an effect size or evidence of causation. Small groups, ties, selection, and confounding can make exploratory results unstable. Inspect distributions and sample sizes before conclusions. The two p-values are unadjusted exploratory tests.

Mutation types are preserved as returned annotations. A sample with no returned variant is not labeled wild type because gene-panel coverage has not been verified. The script does not yet test mutation groups statistically. Shallow/deep CNA deletions are putative categories, not confirmed zygosity.

## Optional survival lesson

First inspect `OS_STATUS`, `OS_MONTHS`, and their definitions. The supplied mapping recognizes only `0:LIVING` and `1:DECEASED`; other values remain missing. Change the mapping only after checking study documentation.

```powershell
python -m pip install lifelines
```

Set `RUN_SURVIVAL = True`. The script fits an unadjusted continuous-expression Cox model only when its basic learning safeguards pass. The minimum counts are guardrails, not a statistical power calculation. It saves `survival_cox.csv` and prints proportional-hazards checks. Expression is standardized, so the hazard ratio corresponds to one standard deviation. Review convergence warnings and assumptions; a saved model is not automatically reliable. See [lifelines Cox documentation](https://lifelines.readthedocs.io/en/latest/Survival%20Regression.html).

## Expand one thing at a time

1. Change `GENE` to another gene and compare output structures.
2. Set `MAX_PATIENTS = None` for the eligible cohort after validating the test table. Archive previous outputs first; each run overwrites files in the output folder.
3. Inspect the actual age/sex/race fields. Add sex-group expression summaries before introducing more tests.
4. Define baseline sample eligibility and verify gene coverage before making mutation-present/absent comparisons.
5. Add age and appropriately encoded sex to a survival model once event counts and assumptions support the larger model.
6. Add race or therapy only when field definitions, group sizes, missingness, and treatment timing support a meaningful question. Therapy response is not implemented here.

This first version is a teaching foundation, not the entire research analysis. It keeps the matched table visible so you can check each expansion. The API routes are based on the [public specification](https://www.cbioportal.org/api/v3/api-docs). Live network retrieval was unavailable during implementation.
