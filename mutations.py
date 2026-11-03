"""List every mutation and copy-number alteration in a set of genes (TCGA LUAD, cBioPortal).

Excel layout: one column per gene. Row 1 holds the gene names; every cell below it
is one altered sample: "sample | alteration type | protein change | DNA change".
"""
from collections import Counter
from pathlib import Path

import requests
from openpyxl import Workbook
from openpyxl.styles import Font

BASE = "https://www.cbioportal.org/api"
STUDY = "luad_tcga_pan_can_atlas_2018"
GENES = ["PTDSS1", "BRAF", "PISD", "LPIN1"]
OUTPUT_FILE = Path(__file__).with_name("gene_mutations.xlsx")

# Readable names for the MAF mutation types
MUTATION_LABELS = {
    "Missense_Mutation": "Missense",
    "Nonsense_Mutation": "Nonsense",
    "Frame_Shift_Del": "Frameshift deletion",
    "Frame_Shift_Ins": "Frameshift insertion",
    "In_Frame_Del": "In-frame deletion",
    "In_Frame_Ins": "In-frame insertion",
    "Splice_Site": "Splice site",
    "Splice_Region": "Splice region",
    "Nonstop_Mutation": "Nonstop",
    "Translation_Start_Site": "Start codon",
    "Silent": "Silent",
}

# How the DNA changed
VARIANT_LABELS = {
    "SNP": "single-base substitution",
    "DNP": "2-base substitution",
    "TNP": "3-base substitution",
    "ONP": "multi-base substitution",
    "DEL": "deletion",
    "INS": "insertion",
}

# GISTIC copy-number calls (putative, not measured copy counts); 0 = diploid is skipped
CNA_LABELS = {
    2: "Amplification",
    1: "Gain",
    -1: "Shallow deletion (likely heterozygous)",
    -2: "Deep deletion (likely homozygous)",
}


def get(path, **params):
    response = requests.get(BASE + path, params=params, timeout=(10, 120))
    response.raise_for_status()
    return response.json()


def post(path, body, **params):
    response = requests.post(BASE + path, json=body, params=params, timeout=(10, 120))
    response.raise_for_status()
    return response.json()


def get_gene_ids(symbols):
    """Map Entrez gene ID -> gene symbol."""
    genes = post("/genes/fetch", symbols, geneIdType="HUGO_GENE_SYMBOL")
    found = {gene["entrezGeneId"]: gene["hugoGeneSymbol"] for gene in genes}
    missing = set(symbols) - set(found.values())
    if missing:
        raise ValueError(f"Genes not found: {', '.join(sorted(missing))}")
    return found


def get_mutations(gene_ids):
    """Mutations in sequenced samples, as (gene, sample, type, protein change, DNA change)."""
    records = post(
        f"/molecular-profiles/{STUDY}_mutations/mutations/fetch",
        {"sampleListId": f"{STUDY}_sequenced", "entrezGeneIds": list(gene_ids)},
    )
    rows = []
    for record in records:
        mutation_type = MUTATION_LABELS.get(record["mutationType"], record["mutationType"])
        dna_change = VARIANT_LABELS.get(record.get("variantType"), record.get("variantType") or "")
        if record.get("referenceAllele") and record.get("variantAllele"):
            dna_change += f" {record['referenceAllele']}>{record['variantAllele']}"
        rows.append((gene_ids[record["entrezGeneId"]], record["sampleId"], mutation_type,
                     record.get("proteinChange") or "", dna_change))
    return rows


def get_copy_number(gene_ids):
    """Non-diploid copy-number calls, in the same format as get_mutations()."""
    records = post(
        f"/molecular-profiles/{STUDY}_gistic/discrete-copy-number/fetch",
        {"sampleListId": f"{STUDY}_cna", "entrezGeneIds": list(gene_ids)},
        discreteCopyNumberEventType="ALL",
    )
    return [(gene_ids[record["entrezGeneId"]], record["sampleId"], CNA_LABELS[record["alteration"]], "", "")
            for record in records if record["alteration"] in CNA_LABELS]


def export_excel(rows):
    workbook = Workbook()

    # Sheet 1: one column per gene, one altered sample per cell
    sheet = workbook.active
    sheet.title = "Alterations by gene"
    type_order = list(MUTATION_LABELS.values()) + list(CNA_LABELS.values())
    for column, gene in enumerate(GENES, start=1):
        sheet.cell(row=1, column=column, value=gene).font = Font(bold=True)
        gene_rows = sorted((r for r in rows if r[0] == gene),
                           key=lambda r: (type_order.index(r[2]) if r[2] in type_order else len(type_order), r[1]))
        for row_number, (_, sample, alteration, protein, dna) in enumerate(gene_rows, start=2):
            text = " | ".join(part for part in (sample, alteration, protein, dna) if part)
            sheet.cell(row=row_number, column=column, value=text)
        sheet.column_dimensions[sheet.cell(row=1, column=column).column_letter].width = 55
    sheet.freeze_panes = "A2"

    # Sheet 2: how many samples have each alteration type in each gene
    counts = Counter((gene, alteration) for gene, _, alteration, _, _ in rows)
    summary = workbook.create_sheet("Counts")
    summary.append(["Alteration type"] + GENES)
    for alteration in type_order:
        if any(counts[(gene, alteration)] for gene in GENES):
            summary.append([alteration] + [counts[(gene, alteration)] for gene in GENES])
    for cell in summary[1]:
        cell.font = Font(bold=True)
    summary.column_dimensions["A"].width = 40
    summary.append([])
    summary.append([f"Study: {STUDY} (cBioPortal)"])
    summary.append(["Mutations: sequenced samples. Copy number: GISTIC calls; diploid samples are not listed."])
    summary.append(["Heterozygous/homozygous labels for deletions are putative (GISTIC), not measured."])
    summary.append(["A sample can appear more than once in a gene (e.g. a mutation and a copy-number change)."])

    workbook.save(OUTPUT_FILE)
    print(f"Excel workbook saved to: {OUTPUT_FILE}")


def main():
    gene_ids = get_gene_ids(GENES)
    print(f"Fetching mutations for {', '.join(GENES)} in {STUDY}...", flush=True)
    rows = get_mutations(gene_ids)
    print("Fetching copy-number calls...", flush=True)
    rows += get_copy_number(gene_ids)

    for gene in GENES:
        counts = Counter(alteration for g, _, alteration, _, _ in rows if g == gene)
        print(f"{gene}: " + ", ".join(f"{alteration} {n}" for alteration, n in counts.most_common()))
    export_excel(rows)


if __name__ == "__main__":
    main()
