# Learn the cBioPortal API with Python

Written for this project on October 31, 2026.

This tutorial teaches a workflow you can adapt yourself: discover a study, inspect its available data, select samples and genes, retrieve results, and join those results correctly. Examples use your genes (`PTDSS1`, `BRAF`, `PISD`, `LPIN1`) and your existing pancreatic cancer study. Examples describe requests, not a guarantee that a study contains every requested measurement.

There **is official documentation**. The challenge is that it is split between an API reference and explanations of the underlying datasets:

- [API overview and clients](https://docs.cbioportal.org/web-api-and-clients/)
- [Interactive Swagger reference](https://www.cbioportal.org/api/swagger-ui/index.html)
- [Machine-readable OpenAPI specification](https://www.cbioportal.org/api/v3/api-docs)
- [Data format definitions](https://docs.cbioportal.org/file-formats/)
- [Interpretation FAQ](https://docs.cbioportal.org/user-guide/faq/)

Keep the reference open as you work. Endpoint availability and accepted parameters can change between portal instances. The teaching examples below were syntax checked, but live retrieval could not be tested in this environment.

## Contents

1. [Understand the data model](#1-understand-the-data-model)
2. [Set up Python](#2-set-up-python)
3. [Understand requests and responses](#3-understand-requests-and-responses)
4. [Build a reusable request helper](#4-build-a-reusable-request-helper)
5. [Discover studies](#5-discover-studies)
6. [Discover profiles and sample lists](#6-discover-profiles-and-sample-lists)
7. [Look up genes](#7-look-up-genes)
8. [Retrieve mutations](#8-retrieve-mutations)
9. [Retrieve copy-number calls](#9-retrieve-copy-number-calls)
10. [Retrieve expression and other numeric measurements](#10-retrieve-expression-and-other-numeric-measurements)
11. [Retrieve clinical data](#11-retrieve-clinical-data)
12. [Join datasets](#12-join-datasets)
13. [Save and analyze results](#13-save-and-analyze-results)
14. [Teach yourself any additional endpoint](#14-teach-yourself-any-additional-endpoint)
15. [Troubleshoot and scale up](#15-troubleshoot-and-scale-up)
16. [Practice exercises](#16-practice-exercises)

## 1. Understand the data model

Think of your analysis as a set of linked tables rather than one large response.

```text
Study
  |-- Patients
  |     |-- Clinical attributes and outcomes
  |     `-- Samples
  |           `-- Clinical sample attributes
  |-- Molecular profiles (measurement datasets)
  |     |-- Mutation records for samples and genes
  |     |-- Copy-number calls for samples and genes
  |     `-- Expression measurements for samples and genes
  `-- Sample lists (named subsets)
```

| Object | Question it answers | Identifier to retain |
|---|---|---|
| Study | Which cohort/publication am I using? | `studyId` |
| Patient | Which person does this record describe? | `patientId` plus study |
| Sample | Which specimen was measured? | `sampleId` plus study |
| Molecular profile | Which assay and processing produced these values? | `molecularProfileId` |
| Sample list | Which named subset do I want? | `sampleListId` |
| Gene | Which gene am I requesting? | `entrezGeneId`, with symbol |

An illustrative profile ID might end with `_mutations` or `_gistic`, but **discover the ID instead of constructing it from a suffix**. One study can have several expression profiles with different scales.

Important distinctions for your project:

- `BRAF` is the gene symbol you want, not `BRAF1`.
- cBioPortal molecular requests commonly use numeric Entrez gene IDs. Your GDC `ENSG...` IDs are a different identifier system.
- A sample is not necessarily a unique patient.
- An individual mutation record is not a patient count.
- GDC variant IDs cannot be used as cBioPortal sample IDs.

## 2. Set up Python

From the project terminal:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install requests pandas openpyxl
```

Select `.venv` as the IDE interpreter. You can run scripts without activating the environment:

```powershell
.\.venv\Scripts\python.exe your_script.py
```

Verify the environment:

```powershell
.\.venv\Scripts\python.exe -c "import sys, requests; print(sys.executable); print(requests.__version__)"
```

Use the same executable for installation and execution. If one Python installs a package and another runs the script, the package may appear to be missing.

## 3. Understand requests and responses

A request has four parts:

| Part | Example | Python |
|---|---|---|
| HTTP method | GET to read a resource | `session.get(...)` |
| URL path | `/studies/<study>/samples` | Part of the URL |
| Query parameters | Page size or projection | `params={...}` |
| JSON body | A list of sample and gene IDs | `json={...}` on POST |

Many `/fetch` endpoints use POST to submit a structured search. In those cases POST retrieves data; it does not upload a study or modify the portal.

```python
import requests

response = requests.get(
    "https://www.cbioportal.org/api/genes/BRAF",
    timeout=(10, 60),
)
response.raise_for_status()
gene = response.json()
print(gene)
```

`timeout=(10, 60)` sets connection and read timeouts. It is not a strict 70-second limit on an entire analysis. `raise_for_status()` prevents you from treating an HTTP error page as data. `.json()` converts JSON into Python dictionaries/lists.

To inspect an unfamiliar result:

```python
print(type(gene))
print(gene.keys())
print(gene.get("entrezGeneId"))
```

For a list response, inspect one record before building loops:

```python
if records:
    print(records[0])
else:
    print("No records returned")
```

Do not print thousands of records while debugging a schema.

## 4. Build a reusable request helper

Put the following in a scratch Python file. Later examples assume these definitions exist.

```python
import json
from pathlib import Path
import requests

BASE = "https://www.cbioportal.org/api"
STUDY = "paad_tcga_pan_can_atlas_2018"
GENES = ["PTDSS1", "BRAF", "PISD", "LPIN1"]
session = requests.Session()
session.headers.update({"Accept": "application/json"})


def request(method, path, *, params=None, body=None):
    response = session.request(
        method,
        BASE + path,
        params=params,
        json=body,
        timeout=(10, 60),
    )
    if not response.ok:
        raise RuntimeError(
            f"HTTP {response.status_code}: {response.url}\n"
            f"{response.text[:800]}"
        )
    return response.json()


def get(path, **params):
    return request("GET", path, params=params)


def post(path, body, **params):
    return request("POST", path, params=params, body=body)


def get_pages(path, page_size=100, **params):
    # Only use on endpoints whose reference lists these pagination parameters.
    page = 0
    while True:
        batch = get(
            path, pageSize=page_size, pageNumber=page, **params
        )
        if not isinstance(batch, list):
            raise TypeError("Expected a list response")
        yield from batch
        if len(batch) < page_size:
            break
        page += 1
```

`yield from` lets the caller process each record without first creating one giant list. `list(get_pages(...))` collects all pages when needed. If a result count is an exact multiple of the page size, this helper requests one final empty page before stopping.

Use a supported, stable sort when paging. A moving server dataset can still change while you retrieve it; record the retrieval date and save raw responses for reproducibility.

Some endpoints support `projection=DETAILED`, which requests additional fields. Projection controls response detail, not eligibility for an assay. Check the reference rather than attaching every parameter to every endpoint.

## 5. Discover studies

```python
studies = list(get_pages("/studies", projection="DETAILED"))

for study in studies:
    if "pancrea" in study.get("name", "").lower():
        print(study["studyId"], study["name"])

study = get(f"/studies/{STUDY}")
print(json.dumps(study, indent=2))
```

Before analyzing, write down why you chose the study: diagnosis, publication, assay coverage, cohort size, and desired clinical fields. Studies with similar names can contain overlapping people or differently processed measurements. Combining them requires checking overlap first.

**Checkpoint:** you can explain what your selected cohort represents without referring only to its ID.

## 6. Discover profiles and sample lists

```python
profiles = list(get_pages(f"/studies/{STUDY}/molecular-profiles"))
for profile in profiles:
    print(
        profile["molecularProfileId"],
        profile.get("molecularAlterationType"),
        profile.get("datatype"),
        profile.get("name"),
        profile.get("description"),
    )

sample_lists = list(get_pages(f"/studies/{STUDY}/sample-lists"))
for sample_list in sample_lists:
    print(sample_list["sampleListId"], sample_list.get("name"))
```

A profile contains measurements; a sample list contains IDs. A named sequenced subset does not imply every gene was measured in every sample. For targeted sequencing, consult gene panels and assay availability before assigning wild-type status.

Choose a profile explicitly after inspecting the list:

```python
mutation_candidates = [
    p for p in profiles
    if p.get("molecularAlterationType") == "MUTATION_EXTENDED"
]
if len(mutation_candidates) != 1:
    raise ValueError("Inspect mutation profiles and choose one explicitly")
mutation_profile = mutation_candidates[0]["molecularProfileId"]

mutation_samples = get(
    f"/molecular-profiles/{mutation_profile}/sample-ids"
)
test_samples = sorted(mutation_samples)[:5]
if not test_samples:
    raise ValueError("No mutation samples available")
```

Selecting the first five IDs is convenient for checking code. It is not random sampling and may produce no alterations in your chosen genes.

## 7. Look up genes

```python
gene_records = [get(f"/genes/{symbol}") for symbol in GENES]
gene_ids = [gene["entrezGeneId"] for gene in gene_records]
symbol_by_id = {
    gene["entrezGeneId"]: gene["hugoGeneSymbol"]
    for gene in gene_records
}
for gene in gene_records:
    print(gene["hugoGeneSymbol"], gene["entrezGeneId"])
```

Use IDs from the response rather than guessing them. An unexpected lookup is a reason to investigate a symbol or alias before proceeding. Keep the mapping in your exports so the numeric request IDs remain understandable.

## 8. Retrieve mutations

This is a bounded test request for the selected genes and five samples:

```python
mutations = post(
    f"/molecular-profiles/{mutation_profile}/mutations/fetch",
    {"sampleIds": test_samples, "entrezGeneIds": gene_ids},
    projection="DETAILED",
    pageSize=100,
    pageNumber=0,
)
for mutation in mutations[:5]:
    print(
        mutation.get("sampleId"),
        symbol_by_id.get(mutation.get("entrezGeneId")),
        mutation.get("mutationType"),
        mutation.get("variantType"),
        mutation.get("proteinChange"),
    )
```

For full retrieval, repeat supported POST pagination until a page contains fewer than `pageSize` records. Do not silently assume the first page is the complete result.

For your requested categories, retain **two separate classifications**:

```python
def classify_sequence_variant(record):
    consequence = str(record.get("mutationType", "")).lower()
    variant_type = str(record.get("variantType", "")).upper()
    labels = []
    if "missense" in consequence:
        labels.append("Missense")
    if "frame_shift" in consequence or "frameshift" in consequence:
        labels.append("Frameshift")
    if variant_type in {"SNP", "DNP", "TNP", "ONP"}:
        labels.append("Substitution")
    return labels
```

A variant can be both missense and substitution. Category counts therefore need not sum to the number of variants. Preserve the original fields alongside these derived labels and inspect the actual observed vocabulary before relying on a classifier.

An empty response is not automatically an error. Verify the profile, genes, selected samples, and assay coverage. A gene can have no returned variants in five samples while being altered elsewhere in the study.

## 9. Retrieve copy-number calls

Use samples belonging to the copy-number profile, rather than assuming the mutation sample list applies:

```python
cna_candidates = [
    p for p in profiles
    if p.get("molecularAlterationType") == "COPY_NUMBER_ALTERATION"
    and p.get("datatype") == "DISCRETE"
]
if len(cna_candidates) != 1:
    raise ValueError("Inspect discrete copy-number profiles and choose one")
cna_profile = cna_candidates[0]["molecularProfileId"]
cna_samples = sorted(get(
    f"/molecular-profiles/{cna_profile}/sample-ids"
))[:5]

cna = post(
    f"/molecular-profiles/{cna_profile}/discrete-copy-number/fetch",
    {"sampleIds": cna_samples, "entrezGeneIds": gene_ids},
    projection="DETAILED",
)
print(cna[:2])
```

Interpret discrete calls as follows:

| Value | Label |
|---|---|
| -2 | Deep deletion; possible homozygous deletion |
| -1 | Shallow deletion; possible heterozygous deletion |
| 0 | Neutral/diploid call |
| 1 | Low-level gain |
| 2 | Amplification |

These are inferred categories, not direct proof of absolute allele counts. Do not apply this mapping to continuous log-ratio data. Do not convert an absent record to zero. See the [copy-number FAQ](https://docs.cbioportal.org/user-guide/faq/).

If you want matched mutation/CNA comparisons, intersect the two profiles' eligible sample IDs **before** selecting your test cohort. Independently selecting five IDs from each profile need not select the same samples.

## 10. Retrieve expression and other numeric measurements

First inspect profiles for the measurement you need. Expression abundance and expression z-scores answer different questions. Do not call an arbitrary numeric expression value TPM or FPKM unless the profile documentation supports that interpretation.

```python
expression_candidates = [
    p for p in profiles
    if p.get("molecularAlterationType") == "MRNA_EXPRESSION"
]
for profile in expression_candidates:
    print(profile["molecularProfileId"], profile.get("description"))
```

After choosing an exact profile ID from that output:

```python
# Replace this with an ID you actually discovered.
expression_profile = "REPLACE_WITH_DISCOVERED_PROFILE_ID"
expression_samples = sorted(get(
    f"/molecular-profiles/{expression_profile}/sample-ids"
))[:5]
expression = post(
    f"/molecular-profiles/{expression_profile}/molecular-data/fetch",
    {"sampleIds": expression_samples, "entrezGeneIds": gene_ids},
    projection="DETAILED",
)
print(expression[:2])
```

For correlations in `correlation.py`, align measurements on the same sample IDs. Never assume two API arrays have matching order. Inspect the scale, missing values, and repeated patient samples before choosing Pearson or Spearman correlation.

Other numeric profiles can use similar retrieval patterns, but verify that the endpoint supports the profile type. Structural variants, gene sets, and generic assays have their own schemas; they should not be forced into a gene-expression request.

## 11. Retrieve clinical data

Your dedicated [clinical guide](CBIOPORTAL_CLINICAL_GUIDE.md) and [small exporter](APIs/cBioPortal_clinical_export.py) complement this section.

Discover metadata first:

```python
attributes = list(get_pages(
    f"/studies/{STUDY}/clinical-attributes", projection="DETAILED"
))
for attribute in attributes:
    print(
        attribute["clinicalAttributeId"],
        attribute.get("displayName"),
        attribute.get("patientAttribute"),
        attribute.get("datatype"),
        attribute.get("description"),
    )
```

Then request values at the correct level:

```python
patient_clinical = list(get_pages(
    f"/studies/{STUDY}/clinical-data",
    clinicalDataType="PATIENT",
))
sample_clinical = list(get_pages(
    f"/studies/{STUDY}/clinical-data",
    clinicalDataType="SAMPLE",
))
samples = list(get_pages(f"/studies/{STUDY}/samples"))
```

Those examples download the study's clinical records. For a quick test, use the existing five-patient exporter instead. Patient-level survival fields and sample-level pathology fields must not be joined on an arbitrary row number.

Clinical events are another structure: a patient can have multiple events. Discover a patient ID, then inspect their timeline:

```python
patient_id = samples[0]["patientId"]
events = list(get_pages(
    f"/studies/{STUDY}/patients/{patient_id}/clinical-events"
))
print(events[:2])
```

Interpret event-specific attributes according to the study. Missing treatment events do not establish that treatment never occurred. Survival needs a defined event indicator and time variable; preserve censoring rather than averaging follow-up as if every patient had an event.

## 12. Join datasets

Use pandas once you understand the raw records:

```python
import pandas as pd

sample_map = pd.DataFrame(samples)[
    ["studyId", "sampleId", "patientId"]
].drop_duplicates()

patient_long = pd.DataFrame(patient_clinical)
keys = ["studyId", "patientId", "clinicalAttributeId"]
if patient_long.duplicated(keys).any():
    raise ValueError("Resolve duplicate clinical values before pivoting")
patient_wide = patient_long.pivot(
    index=["studyId", "patientId"],
    columns="clinicalAttributeId",
    values="value",
).reset_index()

sample_with_clinical = sample_map.merge(
    patient_wide,
    on=["studyId", "patientId"],
    how="left",
    validate="many_to_one",
)
```

`how="left"` keeps samples even when clinical information is unavailable. `validate="many_to_one"` catches unexpected duplication on the clinical side. Handle an empty clinical response before this block: an empty DataFrame does not have the expected columns.

Attach molecular results using `(studyId, sampleId)`. If a molecular response does not include `studyId`, add it from the known selected study; do not infer it from a sample name. If both tables contain `patientId`, resolve duplicated columns explicitly and check that the identifiers agree.

Keep one row per sample-gene for expression/CNA when that is the intended table. Mutation tables can contain several variants per sample-gene; define an aggregation rule before pivoting. For a binary alteration indicator, count distinct sample-gene pairs, not raw mutation rows.

For patient-level analyses, select or aggregate repeated samples using a documented rule. The rule is an analysis decision, not something the API chooses for you.

## 13. Save and analyze results

Save raw data before transformations:

```python
Path("cbio_raw_mutations.json").write_text(
    json.dumps(mutations, indent=2), encoding="utf-8"
)
```

Save separate tables in Excel:

```python
with pd.ExcelWriter("cbio_learning_output.xlsx", engine="openpyxl") as writer:
    pd.DataFrame(gene_records).to_excel(writer, sheet_name="Genes", index=False)
    pd.json_normalize(mutations).to_excel(writer, sheet_name="Mutations", index=False)
    pd.json_normalize(cna).to_excel(writer, sheet_name="CNA", index=False)
    sample_map.to_excel(writer, sheet_name="Sample mapping", index=False)
    patient_wide.to_excel(writer, sheet_name="Patient clinical", index=False)
```

Close the workbook in Excel before overwriting it on Windows. Very large tables belong in CSV or a columnar format rather than a single worksheet. Keep raw text as text when preparing shared spreadsheets so values starting with `=` are not unintentionally interpreted as formulas.

Record these settings beside every export:

- Portal base URL and retrieval date.
- Study, profile, gene, and selected sample IDs.
- Request parameters and any page or sample limits.
- Profile descriptions and clinical field definitions.
- Missing-data rules, category mappings, and patient/sample counting rules.
- Errors and whether retrieval completed.

## 14. Teach yourself any additional endpoint

Use this repeatable method in Swagger:

1. Find a tag matching the resource: studies, genes, clinical data, mutations, molecular data, or another category.
2. Open a route and read the HTTP method and description.
3. Separate required path parameters from optional query parameters.
4. For POST, inspect the request-body schema. Expand referenced schemas rather than guessing field names.
5. Read the response schema to learn whether the endpoint returns a list, object, count, or nested records.
6. Try one valid study, one sample, and one gene where applicable.
7. Translate query inputs to `params=` and the body to `json=` in Python.
8. Inspect a raw response, then add transformations and pagination.

You can inspect the specification programmatically too:

```python
spec_response = session.get(
    "https://www.cbioportal.org/api/v3/api-docs", timeout=(10, 60)
)
spec_response.raise_for_status()
spec = spec_response.json()

for path, methods in spec["paths"].items():
    if "clinical" in path:
        print(path, list(methods))

# Example: inspect a route and its request body/parameters.
route = "/api/studies/{studyId}/clinical-data/fetch"
operation = spec["paths"][route]["post"]
print(json.dumps(operation.get("parameters", []), indent=2))
print(json.dumps(operation.get("requestBody", {}), indent=2))
```

A `$ref` such as `#/components/schemas/SomeFilter` points to `spec["components"]["schemas"]["SomeFilter"]`. Follow that reference and inspect required fields, allowed values, and nested references. This is how you learn an endpoint without relying on copied examples.

For advanced work, search the specification for structural variants, gene panels, generic assays, treatments, survival, and multi-study fetch routes. Their request bodies differ. Inspect the schema and start with a minimal request before scaling it up.

If a plot on the website is difficult to reproduce, your browser's Network panel can show the request used to populate it. Compare that request with the public reference before depending on it; a website can use routes that are not stable public contracts.

## 15. Troubleshoot and scale up

| Symptom | What to check |
|---|---|
| `ModuleNotFoundError` | Install into the interpreter running your script |
| HTTP 400 | Request schema, enum values, numeric gene IDs, parameter spelling |
| HTTP 401/403 | Authentication and access for that portal/study |
| HTTP 404 | Base URL, route, discovered IDs |
| HTTP 429 | Slow down; honor retry instructions |
| HTTP 5xx | Temporary server problem; retry with a bounded delay |
| Timeout | Reduce batch size; retry transient failures; inspect progress |
| JSON parsing fails | Inspect status, content type, and the beginning of response text |
| Empty list | Verify dataset availability, sample eligibility, genes, and filters |
| Missing dictionary key | Inspect the record; check projection and optional fields |
| Unexpectedly huge counts | Check duplicate joins and repeated samples/variants |
| Excel permission error | Close Excel; check destination permissions |

The public portal's openly accessible datasets can be queried without configuring a private-instance token. Other deployments may require authentication. Follow that deployment's instructions; do not put credentials into exported files or commit them to source control.

For larger jobs:

1. Cache discovery results instead of looking up the same genes and profiles repeatedly.
2. Select the cohort and required fields before requesting molecular records.
3. Batch modest numbers of samples/genes, then increase only after checking response size.
4. Use pagination only where the route supports it; batching and pagination solve different problems.
5. Retry connection failures and temporary server errors with bounded backoff. Do not retry malformed requests indefinitely.
6. Save completed batches and progress. A retry must not duplicate already saved records.
7. Limit concurrency and stop increasing load if the service throttles you.

No universal request-per-second limit is assumed here. Use the server's current guidance. Before starting a long run, verify a one-gene/one-sample request and record the expected keys.

## 16. Practice exercises

Work through these in order and inspect the raw results at every stage:

1. Find three pancreatic cancer studies. Explain why they might not be interchangeable.
2. List every molecular profile in your chosen study and describe what you would use it for.
3. Resolve your four genes to Entrez IDs and save the lookup table.
4. Retrieve five mutation-profile samples. Report the number of returned variants and distinct affected samples separately.
5. Retrieve five CNA-profile samples. Count amplification, shallow deletion, and deep deletion calls without treating missing records as neutral.
6. Select samples shared by mutation and CNA profiles, then repeat both requests on the same cohort.
7. Discover available clinical fields. Select two patient fields and two sample fields and explain the difference.
8. Join the molecular data to clinical information. Check the row count before and after each join.
9. Choose an expression profile, document its scale, and align two genes by sample ID before plotting them.
10. Change to another study using discovered IDs. List which assumptions or fields changed.

You have learned the core workflow when you can discover the required IDs, inspect a new request schema, make a small request, interpret its rows, and join those rows without silently changing the unit of analysis.

## How this relates to your existing scripts

| File | Role |
|---|---|
| `APIs/cBioPortal_API_Test.py` | Initial study and metadata exploration |
| `APIs/cBioPortal_clinical_export.py` | Bounded patient/sample clinical export |
| `mutations.py` | GDC sequence variants plus a separate cBioPortal CNA test |
| `correlation.py` | Review its inputs and sample alignment before using newly retrieved measurements |

The tutorial snippets build on shared variables and are learning examples, not independent scripts. In particular, choose a real expression profile before running that section. Keep the five-patient exporter as your first practical test.
