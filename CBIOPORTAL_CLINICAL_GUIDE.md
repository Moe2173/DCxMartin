# cBioPortal clinical outputs

Prepared October 31, 2026. Example study: `paad_tcga_pan_can_atlas_2018`, matching `APIs/cBioPortal_API_Test.py`.

## What you can retrieve

Clinical information describes patients and tumor samples. It can be joined to molecular results to compare alteration groups, investigate expression associations, or prepare outcome analyses. Availability depends on the study: there is no universal list of populated fields.

Start by retrieving the study's attribute dictionary. Treat its IDs, descriptions, data types, and patient/sample classification as the authority for that study.

| Output group | Example field IDs to look for, not guaranteed | Possible use |
|---|---|---|
| Demographics | `AGE`, `SEX` | Describe the cohort; consider confounding |
| Diagnosis | `CANCER_TYPE`, `CANCER_TYPE_DETAILED`, `ONCOTREE_CODE` | Select a diagnosis or subtype |
| Pathology | Study-specific stage, grade, histology fields | Compare disease characteristics |
| Overall survival | `OS_MONTHS`, `OS_STATUS` | Time-to-event analysis |
| Other outcomes | `DFS_MONTHS`, `DFS_STATUS`, study-specific progression fields | Analyze a defined endpoint |
| Sample information | `SAMPLE_TYPE`, study-specific collection or purity fields | Identify relevant tissue and samples |
| Treatment and timeline | Study-specific attributes or clinical events | Explore treatment history where recorded |

These are candidate fields to discover, not a claim that the example study provides all of them. Clinical file conventions are described in the [official format guide](https://docs.cbioportal.org/file-formats/).

## API routes

Base URL: `https://www.cbioportal.org/api`. Replace `{study}` with your study ID.

| Request | Purpose |
|---|---|
| `GET /studies/{study}/clinical-attributes` | Discover fields and their metadata |
| `GET /studies/{study}/clinical-data?clinicalDataType=PATIENT` | Patient attribute values |
| `GET /studies/{study}/clinical-data?clinicalDataType=SAMPLE` | Sample attribute values |
| `GET /studies/{study}/samples` | Sample-to-patient mapping |
| `GET /studies/{study}/patients/{patient}/clinical-events` | Patient timeline, when available |

List endpoints support pagination; page numbers start at zero. Use `projection=DETAILED` where supported to request additional metadata. Consult the [current API specification](https://www.cbioportal.org/api/v3/api-docs) for endpoint parameters and response schemas.

## How to read the outputs

The attribute dictionary explains what a field means. Its `patientAttribute` flag separates patient fields from sample fields. Clinical values arrive in long form: an entity, an attribute ID, and a value. Values may be numeric-looking strings; retain the original values before conversion.

An illustrative row, not downloaded patient data:

```json
{
  "studyId": "example_study",
  "patientId": "example_patient",
  "clinicalAttributeId": "AGE",
  "value": "65"
}
```

Missing fields, blank values, and study-specific missing-value labels should remain missing. They do not imply absence of disease or treatment. An empty events response means no events were returned, not that a patient received no treatment.

## Join clinical information to your gene results

1. Keep `studyId` and `sampleId` with each molecular result.
2. Join sample metadata using `(studyId, sampleId)`.
3. Use that mapping's `patientId` to join patient attributes using `(studyId, patientId)`.
4. Keep repeated samples visible. One patient can have multiple tumor samples; counting sample rows as patients would inflate the cohort size.

The current GDC variant list in `mutations.py` does not identify carriers. You cannot assign patient clinical information to those variants from gene name or mutation UUID alone. Retrieve case/sample occurrences and verify identifier mappings first. The cBioPortal copy-number output already retains sample and study IDs, making it a suitable starting point.

Do not interpret missing mutation data as wild type without confirming that the gene and sample were assayed. Do not combine survival records across repeated samples without choosing a patient-level analysis rule.

## Five-patient Excel test

The companion script selects up to five patients and one sample for each patient. It exports the field dictionary, original clinical records, sample mapping, and retrieval status. This is a convenience sample for testing retrieval, not a representative research cohort.

```powershell
python -m pip install requests openpyxl
python APIs/cBioPortal_clinical_export.py
```

Output: `cbioportal_clinical_sample.xlsx` in the project root. Change `STUDY` and `MAX_PATIENTS` near the top of the script. The script does not export timeline events; those require a separate request to the events route above.

For survival work, inspect both time and status definitions. Preserve censoring, verify time units and the event coding, and avoid correlating survival time alone with gene expression. Treatment availability does not establish response or causation.

API access overview: [cBioPortal API and clients](https://docs.cbioportal.org/web-api-and-clients/). The example exporter follows the public REST API; it has not been verified against live study data in this environment.
