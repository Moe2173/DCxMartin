# GDC single-gene presentation demo

Install dependencies: `python -m pip install requests pandas numpy matplotlib`

Run all eligible primary-tumor RNA cases:

```powershell
python GDC_Test_API.py --gene PTDSS1 --project TCGA-SKCM
```

Small download check: add `--max-rna-cases 5`. This is a nonrandom subset and may contain no mutation-positive cases. Use `--max-rna-cases 0` for all cases. Use `--sample-type Metastatic` to analyze metastatic RNA separately. Change `--gene BRAF` to demonstrate a different gene.

Outputs are in `gdc_results/<project>_<gene>_<sample type>/`:

- `mutation_cohort.csv`: one row per case linked to masked somatic mutation files.
- `mutations.csv` and `mutation_occurrences.json`: reported gene-associated mutation occurrences, including consequences in the nested JSON.
- `rna_manifest.json`: RNA file metadata before selection.
- `matched_expression.csv`: mutation status, target expression, marker TPM values and scores for cases with selected RNA.
- `group_summary.csv`: descriptive counts, medians and means by mutation status.
- `comparisons.png`: expression and three illustrative immune-marker comparisons; dots are cases, black lines are medians.
- `run_metadata.json`: time, parameters, GDC status/release, counts and panel definitions.
- `rna_cache/`: reusable expression downloads.

## Interpretation

Mutation frequency counts unique cases with any reported gene-associated mutation among cases linked to masked mutation files. It does not restrict to protein-changing mutations or establish whether mutations are functional. File availability is a denominator proxy, not proof of adequate coverage at every base of the target gene. Cases without suitable files are excluded; a missing occurrence is labeled “no mutation reported,” not “wild type.”

Mutation status is case-level. The mutation and RNA specimens are not guaranteed to match, so plots cannot establish that a particular RNA specimen carries the mutation. RNA selection keeps one unambiguous file per case of the requested sample type, choosing the first file ID in sorted order when several qualify. Primary and metastatic specimens are not pooled. Cohort mutation frequency is reported across the project mutation cohort; RNA comparisons use its smaller overlapping subgroup.

Scores are the mean log2(TPM + 1) across each listed panel, requiring every marker:

- CD8 markers: CD8A, CD8B.
- Cytotoxic markers: PRF1, GZMB, NKG7.
- Inhibitory markers: PDCD1, LAG3, TIGIT, HAVCR2, TOX.

These are illustrative panels, not validated exhaustion classifiers or estimated cell fractions. Markers can be expressed by other cells. Scores are not comparable in magnitude across panels. The script does not adjust for purity, immune abundance, clinical factors, mutation burden or other confounders. It produces descriptive associations without causal claims or significance tests. It does not measure protein activity, phosphatidylserine exposure or T-cell function. Small mutation-positive groups may prevent useful comparisons.

## API references

- https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/
- https://docs.gdc.cancer.gov/API/Users_Guide/Data_Analysis/
- https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
