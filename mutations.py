import json
import requests
import time
from pathlib import Path
from collections import defaultdict
from openpyxl import Workbook

BASE_URL = "https://api.gdc.cancer.gov"

GENES = ["PTDSS1", "BRAF", "PISD", "LPIN1"]
OUTPUT_FILE = Path(__file__).with_name("gene_mutations.xlsx")
MAX_VARIANTS_PER_GENE = 5
MAX_CNA_SAMPLES = 5
# Same study used in APIs/cBioPortal_API_Test.py.
CBIO_STUDY = "paad_tcga_pan_can_atlas_2018"


def gdc_get(endpoint, params):
    for attempt in range(1, 4):
        print(f"Requesting {endpoint} (attempt {attempt}/3)...", flush=True)
        try:
            response = requests.get(
                f"{BASE_URL}/{endpoint}",
                params=params,
                timeout=(10, 60),
            )
            response.raise_for_status()
            break
        except requests.RequestException as error:
            status = error.response.status_code if error.response is not None else None
            retryable = isinstance(error, (requests.Timeout, requests.ConnectionError)) or status in (429, 500, 502, 503, 504)
            if not retryable or attempt == 3:
                raise
            print(f"Temporary GDC request failure: {error}. Retrying...", flush=True)
            time.sleep(attempt * 2)

    payload = response.json()
    if payload.get("warnings"):
        print("GDC warnings:", payload["warnings"])

    return payload["data"]


def get_gene_id(symbol):
    """Find the Ensembl gene ID for a gene symbol."""
    filters = {
        "op": "in",
        "content": {
            "field": "symbol",
            "value": [symbol],
        },
    }

    data = gdc_get(
        "genes",
        {
            "filters": json.dumps(filters),
            "fields": "gene_id,symbol",
            "size": 100,
        },
    )

    for gene in data["hits"]:
        if gene.get("symbol") == symbol:
            return gene["gene_id"]

    raise ValueError(f"No gene found for symbol: {symbol}")


def get_gene_mutations(gene_id):
    """Retrieve all matching mutations using pagination."""
    filters = {
        "op": "in",
        "content": {
            "field": "consequence.transcript.gene.gene_id",
            "value": [gene_id],
        },
    }

    offset = 0
    page_size = MAX_VARIANTS_PER_GENE

    while True:
        data = gdc_get(
            "ssms",
            {
                "filters": json.dumps(filters),
                "fields": (
                    "ssm_id,genomic_dna_change,"
                    "mutation_subtype"
                ),
                "size": page_size,
                "from": offset,
            },
        )

        hits = data.get("hits", [])
        if not hits:
            break

        yield from hits

        offset += len(hits)
        if offset >= min(data["pagination"]["total"], MAX_VARIANTS_PER_GENE):
            break


def get_mutation_details(ssm_id):
    """Retrieve transcript consequences for one mutation."""
    return gdc_get(
        f"ssms/{ssm_id}",
        {"expand": "consequence.transcript.gene"},
    )


def get_copy_number_sample():
    """Fetch discrete CNA calls for a small, fixed sample cohort."""
    base = "https://www.cbioportal.org/api"

    def get(path, **params):
        response = requests.get(base + path, params=params, timeout=(10, 60))
        response.raise_for_status()
        return response.json()

    profiles = get(f"/studies/{CBIO_STUDY}/molecular-profiles")
    profile = next((p for p in profiles if p.get("molecularAlterationType") == "COPY_NUMBER_ALTERATION" and p.get("datatype") == "DISCRETE"), None)
    if profile is None:
        raise ValueError("No discrete copy-number profile available")
    profile_id = profile["molecularProfileId"]
    sample_ids = sorted(get(f"/molecular-profiles/{profile_id}/sample-ids"))[:MAX_CNA_SAMPLES]
    genes = [get(f"/genes/{symbol}") for symbol in GENES]
    response = requests.post(
        f"{base}/molecular-profiles/{profile_id}/discrete-copy-number/fetch",
        params={"projection": "DETAILED"},
        json={"sampleIds": sample_ids, "entrezGeneIds": [g["entrezGeneId"] for g in genes]},
        timeout=(10, 60),
    )
    response.raise_for_status()
    symbols = {g["entrezGeneId"]: g["hugoGeneSymbol"] for g in genes}
    labels = {
        -2: "Deep deletion (possible homozygous deletion)",
        -1: "Shallow deletion (possible heterozygous deletion)",
        0: "Neutral", 1: "Gain", 2: "Amplification",
    }
    calls = {(r["sampleId"], r["entrezGeneId"]): r["alteration"] for r in response.json()}
    return [
        [symbols[gene["entrezGeneId"]], sample, CBIO_STUDY, profile_id,
         calls.get((sample, gene["entrezGeneId"])),
         labels.get(calls.get((sample, gene["entrezGeneId"])), "Not returned / unavailable")]
        for sample in sample_ids for gene in genes
    ]


def export_excel(rows, statuses, cna_rows=None):
    """Save transcript annotations and counts of distinct variants by type."""
    workbook = Workbook()
    sheet = workbook.active
    sheet.title = "Mutation details"
    sheet.append([
        "Gene", "Gene ID", "Mutation ID", "DNA change", "DNA change type",
        "Mutation consequence", "Protein change", "Transcript", "Canonical",
    ])
    counts = defaultdict(set)
    for row in rows:
        sheet.append(row)
        counts[(row[0], row[5])].add(row[2])

    summary = workbook.create_sheet("Summary")
    summary.append(["Gene", "Mutation consequence", "Distinct mutation count"])
    for (symbol, consequence), mutation_ids in sorted(counts.items()):
        summary.append([symbol, consequence, len(mutation_ids)])

    status_sheet = workbook.create_sheet("Run status")
    status_sheet.append(["Gene", "Retrieved mutation count", "Status"])
    for status in statuses:
        status_sheet.append(status)
    status_sheet.append(["Counts describe variants, not patients."])
    status_sheet.append(["A variant can appear under multiple transcript consequences."])
    status_sheet.append(["Results cover GDC projects without a cancer-type filter."])
    status_sheet.append([f"TEST: at most {MAX_VARIANTS_PER_GENE} variants per gene; counts are not full totals."])
    status_sheet.append(["GDC variants and cBioPortal CNA samples are separate datasets, not matched samples."])
    cna_sheet = workbook.create_sheet("Copy number sample")
    cna_sheet.append(["Gene", "Sample", "Study", "Profile", "Discrete CNA value", "Copy-number category"])
    for row in cna_rows or []:
        cna_sheet.append(row)

    categories = workbook.create_sheet("Requested types")
    categories.append(["Gene", "Mutation ID", "Requested type", "Transcript"])
    seen = set()
    for row in rows:
        labels = []
        consequence = str(row[5]).lower()
        if "missense" in consequence:
            labels.append("Missense")
        if "frameshift" in consequence or "frame_shift" in consequence:
            labels.append("Frameshift")
        if "substitution" in str(row[4]).lower():
            labels.append("Substitution")
        for label in labels:
            key = (row[0], row[2], label, row[7])
            if key not in seen:
                categories.append(list(key))
                seen.add(key)

    for worksheet in workbook:
        worksheet.freeze_panes = "A2"
        worksheet.auto_filter.ref = worksheet.dimensions
        for column in worksheet.columns:
            width = max(len(str(cell.value or "")) for cell in column) + 2
            worksheet.column_dimensions[column[0].column_letter].width = min(width, 60)
    workbook.save(OUTPUT_FILE)
    print(f"\nExcel workbook saved to: {OUTPUT_FILE}")


def main():
    rows = []
    statuses = []
    for symbol in GENES:
        print(f"\n{'=' * 60}")
        print(f"Gene: {symbol}")

        mutation_count = 0
        try:
            gene_id = get_gene_id(symbol)
            print(f"Ensembl gene ID: {gene_id}")

            mutation_count = 0

            for mutation in get_gene_mutations(gene_id):
                mutation_count += 1
                ssm_id = mutation["ssm_id"]

                print(f"\nMutation ID: {ssm_id}")
                print(
                    "DNA change:",
                    mutation.get("genomic_dna_change", "Unknown"),
                )
                print(
                    "DNA change type:",
                    mutation.get("mutation_subtype", "Unknown"),
                )

                details = get_mutation_details(ssm_id)

                matched = False
                for consequence in details.get("consequence", []):
                    transcript = consequence.get("transcript", {})
                    # A variant may also overlap another gene.
                    if transcript.get("gene", {}).get("gene_id") != gene_id:
                        continue
                    matched = True
                    rows.append([
                        symbol, gene_id, ssm_id,
                        mutation.get("genomic_dna_change", ""),
                        mutation.get("mutation_subtype", ""),
                        transcript.get("consequence_type") or "Unknown",
                        transcript.get("aa_change", ""),
                        transcript.get("transcript_id", ""),
                        transcript.get("is_canonical"),
                    ])

                    print(
                        "  Transcript:",
                        transcript.get("transcript_id", "Unknown"),
                    )
                    print(
                        "  Canonical:",
                        transcript.get("is_canonical", "Unknown"),
                    )
                    print(
                        "  Consequence:",
                        transcript.get("consequence_type", "Unknown"),
                    )
                    print(
                        "  Protein change:",
                        transcript.get("aa_change") or "Not available",
                    )

                if not matched:
                    rows.append([
                        symbol, gene_id, ssm_id,
                        mutation.get("genomic_dna_change", ""),
                        mutation.get("mutation_subtype", ""),
                        "Unknown", "", "", None,
                    ])

            print(f"\nTotal mutations retrieved: {mutation_count}")
            statuses.append([symbol, mutation_count, "Complete"])

        except (requests.RequestException, ValueError, KeyError) as error:
            print(f"Failed to process {symbol}: {error}")
            statuses.append([symbol, mutation_count, f"Incomplete: {error}"])
        except KeyboardInterrupt:
            print("\nStopped by user. Saving results collected so far...")
            statuses.append([symbol, mutation_count, "Incomplete: interrupted by user"])
            break

    cna_rows = []
    try:
        print(f"Fetching copy-number calls for {MAX_CNA_SAMPLES} samples from {CBIO_STUDY}...", flush=True)
        cna_rows = get_copy_number_sample()
        statuses.append(["cBioPortal CNA", MAX_CNA_SAMPLES, "Complete"])
    except (requests.RequestException, ValueError, KeyError) as error:
        statuses.append(["cBioPortal CNA", 0, f"Incomplete: {error}"])
    except KeyboardInterrupt:
        statuses.append(["cBioPortal CNA", 0, "Incomplete: interrupted by user"])
    export_excel(rows, statuses, cna_rows)


if __name__ == "__main__":
    main()
