"""GDC presentation demo: case-level mutations and bulk RNA associations."""
import argparse
import io
import json
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

BASE = 'https://api.gdc.cancer.gov'
# Illustrative marker panels, not validated exhaustion or cell-fraction scores.
PANELS = {'CD8_markers': ['CD8A', 'CD8B'],
          'cytotoxic_markers': ['PRF1', 'GZMB', 'NKG7'],
          'inhibitory_markers': ['PDCD1', 'LAG3', 'TIGIT', 'HAVCR2', 'TOX']}


def rule(field, values):
    return {'op': 'in', 'content': {'field': field, 'value': values}}


def filters(*rules):
    return {'op': 'and', 'content': list(rules)}


class GDC:
    def __init__(self):
        self.session = requests.Session()
        retry = Retry(total=3, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504],
                      allowed_methods=['GET', 'POST'])
        self.session.mount('https://', HTTPAdapter(max_retries=retry))

    def request(self, method, endpoint, **kwargs):
        r = self.session.request(method, BASE + '/' + endpoint, timeout=(15, 120), **kwargs)
        if not r.ok:
            raise RuntimeError(f'GDC {endpoint}: HTTP {r.status_code}\n{r.text[:2000]}')
        return r

    def search(self, endpoint, query, fields):
        hits = []
        while True:
            data = self.request('POST', endpoint, json={'filters': query, 'fields': fields,
                                'size': 500, 'from': len(hits)}).json()['data']
            batch = data['hits']
            hits.extend(batch)
            if len(hits) >= data['pagination']['total']:
                return hits
            if not batch:
                raise RuntimeError('Incomplete pagination: ' + endpoint)


def extract(text, gene_id):
    table = pd.read_csv(io.StringIO(text), sep='\t', comment='#')
    target = table.loc[table.gene_id.str.split('.').str[0] == gene_id]
    if len(target) != 1:
        raise ValueError('Expected exactly one target gene row')
    tpm = float(target.iloc[0].tpm_unstranded)
    result = {'gene_TPM': tpm, 'gene_log2_TPM_plus1': float(np.log2(tpm + 1)),
              'gene_raw_count': float(target.iloc[0].unstranded)}
    for panel, markers in PANELS.items():
        levels = []
        for marker in markers:
            match = table.loc[table.gene_name == marker, 'tpm_unstranded']
            value = float(match.iloc[0]) if len(match) == 1 else np.nan
            result[marker + '_TPM'] = value
            levels.append(value)
        result[panel] = float(np.mean(np.log2(np.array(levels) + 1)))
    return result


def plot_results(table, output, gene):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    rng = np.random.default_rng(7)
    groups = ['no_mutation_reported', 'mutation_reported']
    for ax, metric in zip(axes.flat, ['gene_log2_TPM_plus1', *PANELS]):
        labels = []
        for x, group in enumerate(groups):
            values = table.loc[table.mutation_status == group, metric].dropna()
            ax.scatter(x + rng.uniform(-.12, .12, len(values)), values, alpha=.6, s=20)
            if len(values):
                ax.plot([x-.2, x+.2], [values.median()]*2, color='black')
            labels.append(group.replace('_', ' ') + f'\nn={len(values)}')
        ax.set_xticks([0, 1], labels)
        ax.set_title(metric.replace('gene_', gene + ' ').replace('_', ' '))
        ax.set_ylabel('log2(TPM + 1); panels = mean across markers')
    fig.suptitle('Exploratory bulk-tumor associations; case-level mutation status')
    fig.tight_layout()
    fig.savefig(output / 'comparisons.png', dpi=180)
    plt.close(fig)


def run(args):
    output = Path(args.output or f'gdc_results/{args.project}_{args.gene}_{args.sample_type.replace(" ", "_")}')
    output.mkdir(parents=True, exist_ok=True)
    cache = output / 'rna_cache'
    cache.mkdir(exist_ok=True)
    api = GDC()
    genes = api.search('genes', rule('symbol', [args.gene]), 'gene_id,symbol')
    if len(genes) != 1:
        raise ValueError('Use an official gene symbol that resolves to exactly one gene.')
    gene_id = genes[0]['gene_id']
    print('Finding cases with available masked somatic mutation files...', flush=True)
    cases = api.search('cases', filters(rule('project.project_id', [args.project]),
                       rule('files.data_type', ['Masked Somatic Mutation'])), 'case_id,submitter_id')
    eligible = {c['case_id']: c['submitter_id'] for c in cases}
    if not eligible:
        raise ValueError('No cases with available masked mutation files.')
    print('Finding gene mutation occurrences...', flush=True)
    occurrences = api.search('ssm_occurrences', filters(
        rule('case.project.project_id', [args.project]),
        rule('ssm.consequence.transcript.gene.gene_id', [gene_id])),
        'case.case_id,case.submitter_id,ssm.ssm_id,ssm.chromosome,ssm.start_position,'
        'ssm.reference_allele,ssm.tumor_allele,ssm.consequence.transcript.consequence_type,'
        'ssm.consequence.transcript.aa_change,ssm.consequence.transcript.gene.gene_id')
    (output / 'mutation_occurrences.json').write_text(json.dumps(occurrences, indent=2))
    mutated = {hit['case']['case_id'] for hit in occurrences}
    mutation_rows = []
    for hit in occurrences:
        ssm = hit.get('ssm', {})
        transcripts = [c.get('transcript', {}) for c in ssm.get('consequence', [])]
        transcripts = [t for t in transcripts if t.get('gene', {}).get('gene_id') == gene_id]
        mutation_rows.append({'case_id': hit['case']['case_id'],
            'ssm_id': ssm.get('ssm_id'), 'chromosome': ssm.get('chromosome'),
            'position': ssm.get('start_position'), 'reference': ssm.get('reference_allele'),
            'alternate': ssm.get('tumor_allele'),
            'consequence_types': ';'.join(sorted({t['consequence_type'] for t in transcripts if t.get('consequence_type')})),
            'protein_changes': ';'.join(sorted({t['aa_change'] for t in transcripts if t.get('aa_change')})),
            'mutation_json': json.dumps(ssm)})
    pd.DataFrame(mutation_rows, columns=['case_id', 'ssm_id', 'chromosome', 'position',
        'reference', 'alternate', 'consequence_types', 'protein_changes', 'mutation_json']).to_csv(
            output / 'mutations.csv', index=False)
    cohort = pd.DataFrame([{'case_id': cid, 'patient': patient,
        'mutation_status': 'mutation_reported' if cid in mutated else 'no_mutation_reported'}
        for cid, patient in sorted(eligible.items())])
    cohort.to_csv(output / 'mutation_cohort.csv', index=False)
    n_mutated = len(mutated & eligible.keys())
    print(f'{args.gene}: {n_mutated}/{len(eligible)} cases ({100*n_mutated/len(eligible):.2f}%)')
    files = api.search('files', filters(rule('cases.project.project_id', [args.project]),
        rule('data_type', ['Gene Expression Quantification']),
        rule('analysis.workflow_type', ['STAR - Counts']),
        rule('cases.samples.sample_type', [args.sample_type]), rule('access', ['open'])),
        'file_id,file_name,cases.case_id,cases.samples.sample_id,'
        'cases.samples.submitter_id,cases.samples.sample_type')
    (output / 'rna_manifest.json').write_text(json.dumps(files, indent=2))
    selected = {}
    for file in sorted(files, key=lambda f: f['file_id']):
        if len(file.get('cases', [])) != 1:
            continue
        case = file['cases'][0]
        samples = case.get('samples', [])
        if (case['case_id'] in eligible and len(samples) == 1
                and samples[0].get('sample_type') == args.sample_type):
            selected.setdefault(case['case_id'], file)
    items = sorted(selected.items())
    if args.max_rna_cases:
        items = items[:args.max_rna_cases]
    rows = []
    for i, (cid, file) in enumerate(items, 1):
        print(f'RNA {i}/{len(items)}: {eligible[cid]}', flush=True)
        cached = cache / (file['file_id'] + '.tsv')
        text = cached.read_text() if cached.exists() else api.request('GET', 'data/' + file['file_id']).text
        values = extract(text, gene_id)
        if not cached.exists():
            cached.write_text(text)
        sample = file['cases'][0]['samples'][0]
        rows.append({'case_id': cid, 'rna_file_id': file['file_id'],
                     'rna_sample_id': sample['sample_id'], 'rna_sample': sample.get('submitter_id'), **values})
    expression = pd.DataFrame(rows) if rows else pd.DataFrame(columns=['case_id'])
    matched = cohort.merge(expression, on='case_id', how='inner')
    matched.to_csv(output / 'matched_expression.csv', index=False)
    if rows:
        matched.groupby('mutation_status')[['gene_TPM', 'gene_log2_TPM_plus1', *PANELS]].agg(
            ['count', 'median', 'mean']).to_csv(output / 'group_summary.csv')
        plot_results(matched, output, args.gene)
    metadata = {'timestamp_utc': datetime.now(timezone.utc).isoformat(), 'arguments': vars(args),
        'gene_id': gene_id, 'gdc_status': api.request('GET', 'status').json(),
        'denominator_cases': len(eligible), 'mutated_cases': n_mutated,
        'mutation_percent': 100*n_mutated/len(eligible), 'rna_cases_available': len(selected),
        'rna_cases_analyzed': len(rows), 'marker_panels': PANELS}
    (output / 'run_metadata.json').write_text(json.dumps(metadata, indent=2))
    print(f'Results: {output.resolve()}')
    print('All reported gene-associated mutations included; case-level status, not RNA-specimen status.')
    print('Denominator uses masked mutation file availability, not gene-level callable coverage.')
    print('Marker scores are illustrative bulk RNA signals, not measured T-cell effects.')
    if args.max_rna_cases:
        print('RNA uses a limited nonrandom subset. Use --max-rna-cases 0 for all eligible RNA cases.')
    if not rows:
        print('No unambiguous RNA matches. Mutation results were saved.')
    elif (matched.mutation_status == 'mutation_reported').sum() < 5:
        print('Fewer than five mutated RNA cases: comparisons are especially limited.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gene', default='PTDSS1')
    parser.add_argument('--project', default='TCGA-SKCM')
    parser.add_argument('--sample-type', default='Primary Tumor')
    parser.add_argument('--max-rna-cases', type=int, default=0)
    parser.add_argument('--output')
    args = parser.parse_args()
    if args.max_rna_cases < 0:
        parser.error('--max-rna-cases must be zero or positive')
    try:
        run(args)
    except (requests.RequestException, RuntimeError, ValueError) as exc:
        raise SystemExit(f'Analysis stopped: {exc}') from exc
