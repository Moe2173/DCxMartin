import requests
import json
import re


## Creating a gene dictionary to Get information on multiple genes at once
genes_dict = {
    'PTDSS1': 'ENSG00000156471',
    'BRAF1':  'ENSG00000157764',
    'PISD': 'ENSP00000391739',
    'LPIN1': 'ENSG00000134324'
}

## Pulls information about the gene
def gene_sequencing(gene_dictionary):
    for x in genes_dict:
        gene_value = gene_dictionary.get(x)
        genes_endpt = f"https://api.gdc.cancer.gov/genes/{gene_value}?pretty=true"
        response = requests.get(genes_endpt)
        print(response.text)

## Pulls information about the gene somatic mutation
def gene_somatic_mutations(gene_dictionary):
    for x in genes_dict:
        gene_value = gene_dictionary.get(x)
        ssms_endpt = f"https://api.gdc.cancer.gov/ssms/{gene_value}?pretty=true"
        response = requests.get(ssms_endpt)
        print(response.text)


gene_sequencing(genes_dict)
# gene_somatic_mutations(genes_dict)