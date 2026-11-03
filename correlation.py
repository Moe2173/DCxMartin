from datetime import datetime
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
import requests
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from scipy import stats

BASE = 'https://www.cbioportal.org/api'
studyId = 'luad_tcga_pan_can_atlas_2018'
geneSymbol = 'PTDSS1'
minGroupSize = 10
outputRoot = Path(__file__).with_name('correlation_output')

## Every variable the engine can correlate, and its type.
## To add a variable: add a column in buildPatientTable() and one line here.
variables = {
    'copy_number': 'ordinal',          # -2..2 GISTIC call
    'cna_group': 'categorical',        # Deletion / Diploid / Gain/Amp
    'expression_z': 'continuous',      # mRNA z-score
    'mutated': 'categorical',          # Mutated / Wild-type
    'mutation_type': 'categorical',    # Missense_Mutation, ... / Wild-type
    'age': 'continuous',
    'sex': 'categorical',
    'race': 'categorical',
    'stage': 'ordinal',                # 1..4
    'overall_survival': 'survival',
    'pfs': 'survival',
}

## Survival variables are two columns: time (months) and event (1 = death/progression, 0 = censored)
survivalColumns = {
    'overall_survival': ('os_months', 'os_event'),
    'pfs': ('pfs_months', 'pfs_event'),
}

## Variables measured from the same data; correlating them with each other is trivial
sameSource = [{'copy_number', 'cna_group'}, {'mutated', 'mutation_type'}]

## Pairs that always get a plot, even when not significant
keyPlots = [('cna_group', 'overall_survival'), ('expression_z', 'overall_survival'),
            ('mutated', 'overall_survival'), ('cna_group', 'expression_z')]

clinicalFields = ['OS_MONTHS', 'OS_STATUS', 'PFS_MONTHS', 'PFS_STATUS',
                  'AGE', 'SEX', 'RACE', 'AJCC_PATHOLOGIC_TUMOR_STAGE']

## Categorical colors in fixed order (validated colorblind-safe palette)
palette = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']
inkColor = '#0b0b0b'
mutedInk = '#52514e'
gridColor = '#e4e3df'


## Small helpers so every request has a timeout and raises on errors
def apiGet(path, **params):
    response = requests.get(f'{BASE}{path}', params=params, timeout=120)
    response.raise_for_status()
    return response.json()

def apiPost(path, body, **params):
    response = requests.post(f'{BASE}{path}', json=body, params=params, timeout=120)
    response.raise_for_status()
    return response.json()


# Search cancer database for eligble studies
def getStudies():
    response = requests.get(f'{BASE}/studies', timeout=30)
    response.raise_for_status()
    studies = response.json()

    for study in studies:
        if "lung adenocarcinoma" in study["name"].lower():
            print(study["studyId"], study["name"])

## Pull a molecular profile of the study (default: discrete copy number)
## idEndsWith picks between profiles of the same type, e.g. the different z-score profiles
def getStudyMolecularProfile(alterationType='COPY_NUMBER_ALTERATION', datatype='discrete', idEndsWith=None):
    response = requests.get(f'{BASE}/studies/{studyId}/molecular-profiles', timeout=30)
    response.raise_for_status()
    profiles = response.json()

    for profile in profiles:
        if(profile['molecularAlterationType'] == alterationType and profile['datatype'].lower() == datatype.lower()):
            if idEndsWith is None or profile['molecularProfileId'].endswith(idEndsWith):
                return(profile["molecularProfileId"])

    raise ValueError(f'No {alterationType}/{datatype} profile in {studyId}')

## identifies if the study has the clinical data I want to compare
def getClinicalAttributes():
    response = requests.get(f'{BASE}/studies/{studyId}/clinical-attributes', timeout=30)
    response.raise_for_status()
    clinicalData = response.json()

    available = {data["clinicalAttributeId"] for data in clinicalData}
    for field in clinicalFields:
        print(field, 'found' if field in available else 'MISSING')

## gets gene ID
def getGeneInfo():
    response = requests.get(f'{BASE}/genes/{geneSymbol}', timeout=30)
    response.raise_for_status()
    geneInfo = response.json()

    return geneInfo["entrezGeneId"]

##gets sample list (cna, sequenced, rna_seq_v2_mrna, ...)
def getSampleList(listSuffix='cna'):
    response = requests.get(f'{BASE}/studies/{studyId}/sample-lists', timeout=30)
    response.raise_for_status()
    sampleList = response.json()

    for x in sampleList:
        if x['sampleListId'] == f'{studyId}_{listSuffix}':
            return x["sampleListId"]

    raise ValueError(f'No {listSuffix} sample list in {studyId}')

##gets Id of samples from list
def getSampleId(listSuffix='cna'):
    sampleListId = getSampleList(listSuffix)
    response = requests.get(f'{BASE}/sample-lists/{sampleListId}/sample-ids', timeout=30)
    response.raise_for_status()
    sampleIds = response.json()

    return(sampleIds[0:])

## Get the copy number of the sample
def getCopyNumber():
    sampleIds = getSampleId()
    gene_id = getGeneInfo()
    molecularProfileId = getStudyMolecularProfile()

    body = {
        "sampleIds": sampleIds,
        "entrezGeneIds": [gene_id]
    }

    response = requests.post(f'{BASE}/molecular-profiles/{molecularProfileId}/discrete-copy-number/fetch', json=body,timeout=30, params={"discreteCopyNumberEventType": "ALL"})
    response.raise_for_status()
    copyNumber = response.json()

    return copyNumber

## Get mRNA expression of the gene for every sample with RNA-seq: {sampleId: value}
def getExpression(datatype, idEndsWith):
    molecularProfileId = getStudyMolecularProfile('MRNA_EXPRESSION', datatype, idEndsWith)
    body = {
        "sampleListId": getSampleList('rna_seq_v2_mrna'),
        "entrezGeneIds": [getGeneInfo()]
    }
    expression = apiPost(f'/molecular-profiles/{molecularProfileId}/molecular-data/fetch', body)

    return {record["sampleId"]: record["value"] for record in expression}

## Get mutations of the gene. Also returns which samples were sequenced,
## so "no mutation record" can be told apart from "never sequenced"
def getMutations():
    molecularProfileId = getStudyMolecularProfile('MUTATION_EXTENDED', 'maf')
    sampleListId = getSampleList('sequenced')
    body = {
        "sampleListId": sampleListId,
        "entrezGeneIds": [getGeneInfo()]
    }
    mutations = apiPost(f'/molecular-profiles/{molecularProfileId}/mutations/fetch', body)
    sequencedIds = apiGet(f'/sample-lists/{sampleListId}/sample-ids')

    return sequencedIds, mutations

## Get every sample (with its patient and tissue type)
def getSamples():
    return apiGet(f'/studies/{studyId}/samples', pageSize=100000)

## Get the clinical data for all patients in one request: {patientId: {attribute: value}}
def getAllClinicalData():
    records = apiGet(f'/studies/{studyId}/clinical-data', clinicalDataType='PATIENT', pageSize=100000)
    clinical = {}
    for x in records:
        if x["clinicalAttributeId"] in clinicalFields:
            clinical.setdefault(x["patientId"], {})[x["clinicalAttributeId"]] = x["value"]

    return clinical

## Fetch everything buildPatientTable() needs
def getAllData():
    print('Fetching samples...', flush=True)
    samples = getSamples()
    print('Fetching copy number...', flush=True)
    copyNumber = getCopyNumber()
    print('Fetching expression...', flush=True)
    expressionZ = getExpression('z-score', '_rna_seq_v2_mrna_median_all_sample_Zscores')
    expressionRaw = getExpression('continuous', '_rna_seq_v2_mrna')
    print('Fetching mutations...', flush=True)
    sequencedIds, mutations = getMutations()
    print('Fetching clinical data...', flush=True)
    clinical = getAllClinicalData()

    return samples, copyNumber, expressionZ, expressionRaw, sequencedIds, mutations, clinical


## Turn API text into numbers; anything unreadable, negative or infinite becomes NaN
def parseNumber(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return np.nan
    return number if np.isfinite(number) and number >= 0 else np.nan

## "1:DECEASED" / "1:PROGRESSION" -> 1, "0:LIVING" / "0:CENSORED" -> 0, anything else -> NaN
def parseEvent(value):
    if isinstance(value, str) and value.startswith('1:'):
        return 1
    if isinstance(value, str) and value.startswith('0:'):
        return 0
    return np.nan

## "STAGE IIIA" -> 3
def parseStage(value):
    if not isinstance(value, str) or not value.upper().startswith('STAGE '):
        return np.nan
    roman = value.upper().replace('STAGE ', '').rstrip('ABC')
    return {'I': 1, 'II': 2, 'III': 3, 'IV': 4}.get(roman, np.nan)

## One row per patient. Everything is joined by sample/patient ID, never by list order.
## Missing data stays NaN (a missing CNA call is NOT diploid).
def buildPatientTable(samples, copyNumber, expressionZ, expressionRaw, sequencedIds, mutations, clinical):
    ## one primary tumor sample per patient: first sample ID alphabetically
    sampleOf = {}
    for sample in sorted(samples, key=lambda s: s['sampleId']):
        if sample['sampleType'] == 'Primary Solid Tumor':
            sampleOf.setdefault(sample['patientId'], sample['sampleId'])

    cna = {record['sampleId']: record['alteration'] for record in copyNumber}
    mutationTypes = {}
    for record in mutations:
        mutationTypes.setdefault(record['sampleId'], set()).add(record['mutationType'])
    sequenced = set(sequencedIds)

    rows = []
    for patientId, sampleId in sorted(sampleOf.items()):
        info = clinical.get(patientId, {})
        if sampleId not in sequenced:
            mutated, mutationType = np.nan, np.nan
        elif sampleId in mutationTypes:
            mutated, mutationType = 'Mutated', ';'.join(sorted(mutationTypes[sampleId]))
        else:
            mutated, mutationType = 'Wild-type', 'Wild-type'
        raw = expressionRaw.get(sampleId)

        rows.append({
            'patientId': patientId,
            'sampleId': sampleId,
            'copy_number': cna.get(sampleId, np.nan),
            'expression_z': expressionZ.get(sampleId, np.nan),
            'expression_log2': np.log2(raw + 1) if raw is not None and raw >= 0 else np.nan,
            'mutated': mutated,
            'mutation_type': mutationType,
            'os_months': parseNumber(info.get('OS_MONTHS')),
            'os_event': parseEvent(info.get('OS_STATUS')),
            'pfs_months': parseNumber(info.get('PFS_MONTHS')),
            'pfs_event': parseEvent(info.get('PFS_STATUS')),
            'age': parseNumber(info.get('AGE')),
            'sex': info.get('SEX', np.nan),
            'race': info.get('RACE', np.nan),
            'stage': parseStage(info.get('AJCC_PATHOLOGIC_TUMOR_STAGE')),
        })

    table = pd.DataFrame(rows)
    table['cna_group'] = table['copy_number'].map({-2: 'Deletion', -1: 'Deletion', 0: 'Diploid', 1: 'Gain/Amp', 2: 'Gain/Amp'})
    median = table['expression_z'].median()
    table['expression_group'] = pd.Series(np.where(table['expression_z'] >= median, 'High', 'Low'),
                                          dtype=object).where(table['expression_z'].notna())

    return table


## Columns that hold a variable (survival variables use two)
def columnsOf(var):
    return list(survivalColumns.get(var, (var,)))

## Drop categories with too few patients; returns the data and the dropped categories
def dropSmallGroups(data, column):
    counts = data[column].value_counts()
    small = sorted(counts[counts < minGroupSize].index)
    return data[~data[column].isin(small)], small

## Median survival per group, as text
def medianSurvival(data, timeCol, eventCol, groupCol):
    parts = []
    for group, groupData in data.groupby(groupCol):
        kmf = KaplanMeierFitter().fit(groupData[timeCol], groupData[eventCol])
        median = kmf.median_survival_time_
        parts.append(f'{group}: {median:.1f}' if np.isfinite(median) else f'{group}: not reached')
    return '; '.join(parts)


def correlate(table, varA, varB, types=variables):
    row = {'varA': varA, 'varB': varB, 'test': '', 'n': 0, 'statistic': np.nan,
           'effect_measure': '', 'effect_size': np.nan, 'effect_detail': '', 'p_value': np.nan, 'note': ''}

    ## put survival first, then categorical, so each case is handled once
    order = {'survival': 0, 'categorical': 1, 'continuous': 2, 'ordinal': 2}
    if order[types[varA]] > order[types[varB]]:
        varA, varB = varB, varA
    typeA, typeB = types[varA], types[varB]

    if typeA == 'survival' and typeB == 'survival':
        row['note'] = 'both variables are survival outcomes'
        return row

    data = table[columnsOf(varA) + columnsOf(varB)].dropna()
    notes = []
    for var in (varA, varB):
        if types[var] == 'categorical':
            data, small = dropSmallGroups(data, var)
            if small:
                notes.append(f'dropped {var} groups with <{minGroupSize} patients: {", ".join(map(str, small))}')
    row['n'] = len(data)

    ## guardrails: enough patients, enough groups, something varies
    for var in (varA, varB):
        if types[var] == 'categorical' and data[var].nunique() < 2:
            notes.append(f'{var} has fewer than 2 groups with >={minGroupSize} patients')
        elif types[var] in ('continuous', 'ordinal') and data[var].nunique() < 2:
            notes.append(f'{var} has only one value')
    if typeA == 'survival' and data[survivalColumns[varA][1]].sum() == 0:
        notes.append('no events')
    if len(data) < minGroupSize:
        notes.append(f'fewer than {minGroupSize} patients')
    row['note'] = '; '.join(notes)
    if any(not note.startswith('dropped') for note in notes):
        return row

    if typeA == 'survival' and typeB == 'categorical':
        timeCol, eventCol = survivalColumns[varA]
        result = multivariate_logrank_test(data[timeCol], data[varB], data[eventCol])
        row.update(test='log-rank', statistic=result.test_statistic, p_value=result.p_value,
                   effect_measure='median months',
                   effect_detail=medianSurvival(data, timeCol, eventCol, varB))

    elif typeA == 'survival':
        timeCol, eventCol = survivalColumns[varA]
        model = data[[timeCol, eventCol, varB]].copy()
        unit = '1 step'
        if typeB == 'continuous':
            model[varB] = (model[varB] - model[varB].mean()) / model[varB].std()
            unit = '1 SD'
        cox = CoxPHFitter().fit(model, timeCol, eventCol)
        summary = cox.summary.loc[varB]
        row.update(test='Cox PH', statistic=summary['z'], p_value=summary['p'],
                   effect_measure=f'hazard ratio per {unit}', effect_size=summary['exp(coef)'],
                   effect_detail=f"95% CI {summary['exp(coef) lower 95%']:.2f}-{summary['exp(coef) upper 95%']:.2f}")

    elif typeA == 'categorical' and typeB == 'categorical':
        crosstab = pd.crosstab(data[varA], data[varB])
        chi2 = stats.chi2_contingency(crosstab, correction=False)[0]
        cramersV = np.sqrt(chi2 / (len(data) * (min(crosstab.shape) - 1)))
        if crosstab.shape == (2, 2):
            oddsRatio, p = stats.fisher_exact(crosstab)
            row.update(test='Fisher exact', statistic=oddsRatio, p_value=p, effect_detail=f'odds ratio {oddsRatio:.2f}')
        else:
            chi2, p, dof, expected = stats.chi2_contingency(crosstab)
            row.update(test='chi-square', statistic=chi2, p_value=p)
            if (expected < 5).any():
                row['note'] = '; '.join(filter(None, [row['note'], 'some expected counts <5']))
        row.update(effect_measure="Cramer's V", effect_size=cramersV)

    elif typeA == 'categorical':
        groups = {name: values[varB] for name, values in data.groupby(varA)}
        medians = '; '.join(f'{name}: {values.median():.2f}' for name, values in groups.items())
        if len(groups) == 2:
            (nameA, valuesA), (nameB, valuesB) = groups.items()
            u, p = stats.mannwhitneyu(valuesA, valuesB)
            rankBiserial = 1 - 2 * u / (len(valuesA) * len(valuesB))
            row.update(test='Mann-Whitney', statistic=u, p_value=p,
                       effect_measure=f'rank-biserial ({nameB} vs {nameA})', effect_size=rankBiserial)
        else:
            h, p = stats.kruskal(*groups.values())
            row.update(test='Kruskal-Wallis', statistic=h, p_value=p,
                       effect_measure='epsilon squared', effect_size=h / (len(data) - 1))
        row['effect_detail'] = f'median {varB}: {medians}'

    else:
        rho, p = stats.spearmanr(data[varA], data[varB])
        row.update(test='Spearman', statistic=rho, p_value=p, effect_measure='rho', effect_size=rho)

    return row

## Correlate every pair of variables (or only pairs involving `against`), then add
## Benjamini-Hochberg q-values because many tests are run
def correlateAll(table, types=variables, against=None):
    rows = []
    for varA, varB in combinations(types, 2):
        if {varA, varB} in sameSource:
            continue
        if against and varA not in against and varB not in against:
            continue
        rows.append(correlate(table, varA, varB, types))

    results = pd.DataFrame(rows)
    tested = results['p_value'].notna()
    results['q_value'] = np.nan
    if tested.any():
        results.loc[tested, 'q_value'] = stats.false_discovery_control(results.loc[tested, 'p_value'])

    return results.sort_values(['q_value', 'p_value'], na_position='last').reset_index(drop=True)

## Cox model with several variables at once, e.g. does expression still matter after age, sex and stage?
## Continuous variables are per 1 SD, ordinal per 1 step, categorical vs. their first group
def coxAnalysis(table, covariates, survival='overall_survival', types=variables):
    timeCol, eventCol = survivalColumns[survival]
    data = table[[timeCol, eventCol] + covariates].dropna()
    for var in covariates:
        if types[var] == 'continuous':
            data[var] = (data[var] - data[var].mean()) / data[var].std()
    categorical = [var for var in covariates if types[var] == 'categorical']
    data = pd.get_dummies(data, columns=categorical, drop_first=True, dtype=float)

    cox = CoxPHFitter().fit(data, timeCol, eventCol)
    summary = cox.summary[['exp(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%', 'p']]
    summary.columns = ['hazard_ratio', 'ci_lower', 'ci_upper', 'p_value']
    summary.index.name = 'variable'

    return summary, len(data), int(data[eventCol].sum())


## Shared look for every plot
def styleAxes(ax, title, xlabel, ylabel):
    ax.set_title(title, color=inkColor, fontsize=11, loc='left')
    ax.set_xlabel(xlabel, color=mutedInk)
    ax.set_ylabel(ylabel, color=mutedInk)
    ax.tick_params(colors=mutedInk)
    ax.grid(True, color=gridColor, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(gridColor)

## One plot for a pair of variables; the kind of plot follows the variable types
def plotPair(table, varA, varB, path, types=variables):
    order = {'survival': 0, 'categorical': 1, 'continuous': 2, 'ordinal': 2}
    if order[types[varA]] > order[types[varB]]:
        varA, varB = varB, varA
    typeA, typeB = types[varA], types[varB]
    data = table[columnsOf(varA) + columnsOf(varB)].dropna()
    title = f'{geneSymbol}: {varA} vs {varB}'
    fig, ax = plt.subplots(figsize=(7, 4.5))

    if typeA == 'survival':
        timeCol, eventCol = survivalColumns[varA]
        groupCol = varB
        if typeB != 'categorical':
            ## split continuous values at the median for the picture (the test stays continuous)
            groupCol = f'{varB} group'
            data[groupCol] = np.where(data[varB] >= data[varB].median(), 'High', 'Low')
        data, small = dropSmallGroups(data, groupCol)
        for color, (group, groupData) in zip(palette, data.groupby(groupCol)):
            kmf = KaplanMeierFitter().fit(groupData[timeCol], groupData[eventCol], label=f'{group} (n={len(groupData)})')
            kmf.plot_survival_function(ax=ax, color=color, linewidth=2, ci_alpha=0.12,
                                       show_censors=True, censor_styles={'ms': 5, 'marker': '|'})
        styleAxes(ax, title, 'Months', 'Survival probability')
        ax.set_ylim(0, 1.02)
        ax.legend(frameon=False, labelcolor=inkColor)

    elif typeA == 'categorical' and typeB == 'categorical':
        data, small = dropSmallGroups(data, varA)
        data, small = dropSmallGroups(data, varB)
        crosstab = pd.crosstab(data[varA], data[varB], normalize='index')
        bottom = np.zeros(len(crosstab))
        for color, column in zip(palette, crosstab.columns):
            ax.barh(crosstab.index.astype(str), crosstab[column], left=bottom, color=color,
                    edgecolor='white', linewidth=2, label=str(column), height=0.6)
            bottom += crosstab[column].to_numpy()
        styleAxes(ax, title, f'Fraction of patients ({varB})', varA)
        ax.legend(frameon=False, labelcolor=inkColor, bbox_to_anchor=(1, 1), loc='upper left')

    elif typeA == 'categorical':
        data, small = dropSmallGroups(data, varA)
        groups = sorted(data[varA].unique())
        values = [data.loc[data[varA] == group, varB] for group in groups]
        ax.boxplot(values, widths=0.5, showfliers=False,
                   medianprops={'color': inkColor, 'linewidth': 2}, boxprops={'color': mutedInk},
                   whiskerprops={'color': mutedInk}, capprops={'color': mutedInk})
        jitter = np.random.default_rng(0)
        for i, (color, groupValues) in enumerate(zip(palette, values), start=1):
            ax.scatter(i + jitter.uniform(-0.15, 0.15, len(groupValues)), groupValues, s=10, color=color, alpha=0.6)
        ax.set_xticks(range(1, len(groups) + 1), [f'{group}\n(n={len(v)})' for group, v in zip(groups, values)])
        styleAxes(ax, title, varA, varB)

    else:
        ax.scatter(data[varA], data[varB], s=12, color=palette[0], alpha=0.6)
        styleAxes(ax, title, varA, varB)

    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor='white')
    plt.close(fig)

## Grid of every variable pair, colored by -log10(q). Numbers are effect sizes, * = q < 0.05
def plotHeatmap(results, path, types=variables):
    names = list(types)
    strength = pd.DataFrame(np.nan, index=names, columns=names)
    labels = pd.DataFrame('', index=names, columns=names)
    for _, row in results.iterrows():
        if pd.isna(row['q_value']):
            text = 'skip'
            value = np.nan
        else:
            value = -np.log10(max(row['q_value'], 1e-10))
            text = f"{row['effect_size']:.2f}" if pd.notna(row['effect_size']) else ''
            text += '*' if row['q_value'] < 0.05 else ''
        for a, b in ((row['varA'], row['varB']), (row['varB'], row['varA'])):
            strength.loc[a, b] = value
            labels.loc[a, b] = text

    blues = LinearSegmentedColormap.from_list('blues', ['#cde2fb', '#86b6ef', '#2a78d6', '#184f95', '#0d366b'])
    blues.set_bad('#f0efec')
    fig, ax = plt.subplots(figsize=(10, 8.5))
    image = ax.imshow(strength.to_numpy(dtype=float), cmap=blues, vmin=0, vmax=max(2, np.nanmax(strength.to_numpy(dtype=float), initial=0)))
    ax.set_xticks(range(len(names)), names, rotation=45, ha='right', color=mutedInk)
    ax.set_yticks(range(len(names)), names, color=mutedInk)
    for i in range(len(names)):
        for j in range(len(names)):
            value = strength.iloc[i, j]
            ax.text(j, i, labels.iloc[i, j], ha='center', va='center', fontsize=7,
                    color='white' if pd.notna(value) and value > 1.2 else inkColor)
    for side in ax.spines.values():
        side.set_visible(False)
    colorbar = fig.colorbar(image, ax=ax, shrink=0.7)
    colorbar.set_label('-log10(q)   (1.3 = q of 0.05)', color=mutedInk)
    ax.set_title(f'{geneSymbol} ({studyId}): every variable pair\n'
                 'cell number = effect size (rho, HR, V, ...); * = q < 0.05; gray = not tested',
                 color=inkColor, fontsize=11, loc='left')
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor='white')
    plt.close(fig)

## Readable summary of the run
def writeSummary(path, table, results, cox, coxN, coxEvents):
    significant = results[results['q_value'] < 0.05]
    lines = [
        f'{geneSymbol} correlation analysis, {studyId}',
        f'Run: {datetime.now().isoformat(timespec="seconds")}',
        '',
        f'Patients (one primary tumor each): {len(table)}',
        f'  with copy number: {table["copy_number"].notna().sum()}',
        f'  with expression: {table["expression_z"].notna().sum()}',
        f'  sequenced: {table["mutated"].notna().sum()} ({(table["mutated"] == "Mutated").sum()} mutated)',
        f'  with overall survival: {table[["os_months", "os_event"]].notna().all(axis=1).sum()} '
        f'({int(table["os_event"].sum())} deaths)',
        '',
        f'Pairs tested: {results["p_value"].notna().sum()} (skipped: {results["p_value"].isna().sum()}, see note column)',
        f'Significant after Benjamini-Hochberg correction (q < 0.05): {len(significant)}',
    ]
    for _, row in significant.iterrows():
        lines.append(f"  {row['varA']} vs {row['varB']}: {row['test']}, q = {row['q_value']:.3g}, "
                     f"{row['effect_measure']} {row['effect_size']:.2f} {row['effect_detail']}".rstrip())
    lines += ['', f'Adjusted Cox model for overall survival (n = {coxN}, {coxEvents} deaths):', cox.round(3).to_string(), '',
              'Notes:',
              '- Exploratory, mostly unadjusted tests; correlation is not causation.',
              '- Survival is always tested with censoring (log-rank / Cox), never by correlating raw months.',
              '- Pairs between clinical variables (e.g. age vs stage) are included as context, not gene findings.',
              f'- Groups with fewer than {minGroupSize} patients are dropped from a test; see the note column.']
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return '\n'.join(lines)


def main():
    outputFolder = outputRoot / datetime.now().strftime('%Y%m%d_%H%M%S')
    (outputFolder / 'plots').mkdir(parents=True)

    table = buildPatientTable(*getAllData())
    table.to_csv(outputFolder / 'patient_table.csv', index=False)

    print('Correlating...', flush=True)
    results = correlateAll(table)
    results.to_csv(outputFolder / 'results.csv', index=False)

    cox, coxN, coxEvents = coxAnalysis(table, ['expression_z', 'age', 'sex', 'stage'])
    cox.to_csv(outputFolder / 'cox_adjusted.csv')

    print('Plotting...', flush=True)
    plotHeatmap(results, outputFolder / 'correlation_heatmap.png')
    significant = results.loc[results['q_value'] < 0.05, ['varA', 'varB']].itertuples(index=False, name=None)
    for varA, varB in dict.fromkeys(keyPlots + list(significant)):
        plotPair(table, varA, varB, outputFolder / 'plots' / f'{varA}__{varB}.png')

    print(writeSummary(outputFolder / 'summary.txt', table, results, cox, coxN, coxEvents))
    print(f'\nSaved to {outputFolder}')


if __name__ == "__main__":
    # getStudies()
    # getClinicalAttributes()
    main()
