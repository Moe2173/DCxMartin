import pandas as pd
import requests

GENE = 'PTDSS1'
GENE_ID = 'ENSG00000156471'
cancers = {
    'Pancreatic': 'TCGA-PAAD',
    'Melanoma': 'TCGA-SKCM',
    'Breast': 'TCGA-BRCA',
    'Lung adenocarcinoma': 'TCGA-LUAD',
    'Colon': 'TCGA-COAD',
}

BASE = 'https://api.gdc.cancer.gov'


# Search GDC using a set of matching fields.
def search(endpoint, conditions, fields):
    rules = [
        {'op': 'in', 'content': {'field': field, 'value': values}}
        for field, values in conditions.items()
    ]
    response = requests.post(f'{BASE}/{endpoint}', json={
        'filters': {'op': 'and', 'content': rules},
        'fields': fields,
        'size': 10000,
    }, timeout=120)
    response.raise_for_status()
    data = response.json()['data']
    if data['pagination']['total'] > len(data['hits']):
        raise ValueError('Too many results for this small demo.')
    return data['hits']

# Compare mutation frequency across the cancer projects.
results = []
for cancer, project in cancers.items():
    print(f'Checking {cancer}...', flush=True)

    cases = search('cases', {
        'project.project_id': [project],
        'files.data_type': ['Masked Somatic Mutation'],
    }, 'case_id')

    mutations = search('ssm_occurrences', {
        'case.project.project_id': [project],
        'ssm.consequence.transcript.gene.gene_id': [GENE_ID],
    }, 'case.case_id')

    # Sets count each case once, even if it has several mutations.
    tested_ids = {case['case_id'] for case in cases}
    mutated_ids = {item['case']['case_id'] for item in mutations} & tested_ids

    if not tested_ids:
        print('No mutation data available; skipping this project.')
        continue

    results.append({
        'Cancer': cancer,
        'Project': project,
        'Cases with mutation data': len(tested_ids),
        'Mutated cases': len(mutated_ids),
        'Mutation frequency (%)': 100 * len(mutated_ids) / len(tested_ids),
    })

## Followed a tutorial to set up this table 
if results:
    ranking = pd.DataFrame(results).sort_values('Mutation frequency (%)', ascending=False)
    print(f'\n{GENE} mutation frequency by cancer:')
    print(ranking.to_string(index=False, float_format='%.2f'))
    
else:
    print('No results found.')
