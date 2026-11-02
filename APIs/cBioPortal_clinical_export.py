"""Export a small clinical dataset; see CBIOPORTAL_CLINICAL_GUIDE.md."""
from pathlib import Path
import json
import requests
from openpyxl import Workbook

BASE = "https://www.cbioportal.org/api"
STUDY = "paad_tcga_pan_can_atlas_2018"
MAX_PATIENTS = 5
OUTPUT = Path(__file__).resolve().parents[1] / "cbioportal_clinical_sample.xlsx"


def get(path, **params):
    response = requests.get(BASE + path, params=params, timeout=(10, 60))
    response.raise_for_status()
    return response.json()


def pages(path, **params):
    page = 0
    while True:
        batch = get(path, **params, pageSize=100, pageNumber=page)
        yield from batch
        if len(batch) < 100:
            break
        page += 1


def add_sheet(workbook, name, records):
    sheet = workbook.create_sheet(name)
    if not records:
        sheet.append(["No records returned"])
        return
    columns = list(dict.fromkeys(key for record in records for key in record))
    sheet.append(columns)
    for record in records:
        values = []
        for key in columns:
            value = record.get(key)
            if isinstance(value, (dict, list)):
                value = json.dumps(value)
            values.append(value)
        sheet.append(values)
        # Preserve API strings as text, including strings beginning with '='.
        for cell in sheet[sheet.max_row]:
            if isinstance(cell.value, str):
                cell.data_type = "s"
    sheet.freeze_panes = "A2"
    sheet.auto_filter.ref = sheet.dimensions


def main():
    attributes = list(pages(f"/studies/{STUDY}/clinical-attributes", projection="DETAILED"))
    patients = get(
        f"/studies/{STUDY}/patients", pageSize=MAX_PATIENTS,
        pageNumber=0, sortBy="patientId", direction="ASC",
    )
    patient_rows, sample_rows, mapping, status = [], [], [], []
    for patient in patients:
        patient_id = patient["patientId"]
        print(f"Retrieving clinical information: {patient_id}", flush=True)
        try:
            patient_rows.extend(pages(f"/studies/{STUDY}/patients/{patient_id}/clinical-data"))
            samples = get(
                f"/studies/{STUDY}/patients/{patient_id}/samples",
                pageSize=1, pageNumber=0, sortBy="sampleId", direction="ASC",
            )
            for sample in samples:
                mapping.append(sample)
                sample_rows.extend(pages(f"/studies/{STUDY}/samples/{sample['sampleId']}/clinical-data"))
            status.append({"patientId": patient_id, "status": "Complete", "selectedSamples": len(samples)})
        except (requests.RequestException, ValueError, KeyError) as error:
            status.append({"patientId": patient_id, "status": f"Incomplete: {error}"})

    workbook = Workbook()
    workbook.remove(workbook.active)
    for name, records in [
        ("Attributes", attributes), ("Patient clinical", patient_rows),
        ("Sample clinical", sample_rows), ("Sample mapping", mapping),
        ("Retrieval status", status),
        ("Run settings", [{"studyId": STUDY, "maximumPatients": MAX_PATIENTS,
                           "samplesPerPatient": 1, "selection": "First IDs alphabetically"}]),
    ]:
        add_sheet(workbook, name, records)
    workbook.save(OUTPUT)
    print(f"Saved: {OUTPUT}")


if __name__ == "__main__":
    main()
