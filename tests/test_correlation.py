"""Offline tests for the patient table and the correlation engine (no network)."""
import unittest

import numpy as np
import pandas as pd

from correlation import buildPatientTable, correlate, correlateAll, coxAnalysis, parseEvent, parseStage


def rawData():
    samples = [
        {'patientId': 'p1', 'sampleId': 's1-02', 'sampleType': 'Primary Solid Tumor'},
        {'patientId': 'p1', 'sampleId': 's1-01', 'sampleType': 'Primary Solid Tumor'},
        {'patientId': 'p2', 'sampleId': 's2-01', 'sampleType': 'Primary Solid Tumor'},
        {'patientId': 'p3', 'sampleId': 's3-06', 'sampleType': 'Metastatic'},
        {'patientId': 'p4', 'sampleId': 's4-01', 'sampleType': 'Primary Solid Tumor'},
    ]
    copyNumber = [{'sampleId': 's2-01', 'alteration': 1}, {'sampleId': 's1-02', 'alteration': -2},
                  {'sampleId': 's1-01', 'alteration': 0}]
    expressionZ = {'s1-01': -1.0, 's2-01': 2.0}
    expressionRaw = {'s1-01': 3.0, 's2-01': 15.0}
    sequencedIds = ['s1-01', 's2-01']
    mutations = [{'sampleId': 's2-01', 'mutationType': 'Missense_Mutation'}]
    clinical = {
        'p1': {'OS_MONTHS': '12', 'OS_STATUS': '1:DECEASED', 'AGE': '70', 'SEX': 'Male',
               'AJCC_PATHOLOGIC_TUMOR_STAGE': 'STAGE IIIA'},
        'p2': {'OS_MONTHS': '30', 'OS_STATUS': '0:LIVING', 'SEX': 'Female'},
        'p4': {'OS_MONTHS': '-1', 'OS_STATUS': 'DECEASED'},
    }
    return samples, copyNumber, expressionZ, expressionRaw, sequencedIds, mutations, clinical


## A synthetic cohort where expression drives survival and differs by sex
def syntheticCohort(n=200):
    rng = np.random.default_rng(1)
    expression = rng.normal(size=n)
    sex = np.where(np.arange(n) % 2 == 0, 'Female', 'Male')
    expression = expression + np.where(sex == 'Male', 1.0, 0.0)
    months = rng.exponential(60 / np.exp(0.8 * expression))
    return pd.DataFrame({
        'expression_z': expression,
        'age': rng.normal(65, 8, size=n),
        'sex': sex,
        'cna_group': rng.choice(['Deletion', 'Diploid', 'Gain/Amp'], size=n),
        'stage': rng.integers(1, 5, size=n).astype(float),
        'mutated': ['Mutated'] * 5 + ['Wild-type'] * (n - 5),
        'os_months': np.minimum(months, 100),
        'os_event': (months < 100).astype(float),
    })


testTypes = {'expression_z': 'continuous', 'age': 'continuous', 'sex': 'categorical',
             'cna_group': 'categorical', 'stage': 'ordinal', 'mutated': 'categorical',
             'overall_survival': 'survival'}


class TableTests(unittest.TestCase):
    def setUp(self):
        self.table = buildPatientTable(*rawData()).set_index('patientId')

    def test_one_primary_sample_per_patient(self):
        self.assertEqual(list(self.table.index), ['p1', 'p2', 'p4'])
        self.assertEqual(self.table.loc['p1', 'sampleId'], 's1-01')

    def test_joined_by_sample_id(self):
        self.assertEqual(self.table.loc['p1', 'copy_number'], 0)
        self.assertEqual(self.table.loc['p2', 'cna_group'], 'Gain/Amp')
        self.assertEqual(self.table.loc['p2', 'expression_z'], 2.0)
        self.assertEqual(self.table.loc['p2', 'expression_log2'], 4.0)

    def test_missing_copy_number_is_not_diploid(self):
        self.assertTrue(pd.isna(self.table.loc['p4', 'copy_number']))
        self.assertTrue(pd.isna(self.table.loc['p4', 'cna_group']))

    def test_wild_type_vs_not_sequenced(self):
        self.assertEqual(self.table.loc['p1', 'mutated'], 'Wild-type')
        self.assertEqual(self.table.loc['p2', 'mutated'], 'Mutated')
        self.assertEqual(self.table.loc['p2', 'mutation_type'], 'Missense_Mutation')
        self.assertTrue(pd.isna(self.table.loc['p4', 'mutated']))

    def test_clinical_parsing(self):
        self.assertEqual(self.table.loc['p1', 'os_event'], 1)
        self.assertEqual(self.table.loc['p2', 'os_event'], 0)
        self.assertTrue(pd.isna(self.table.loc['p4', 'os_event']))
        self.assertTrue(pd.isna(self.table.loc['p4', 'os_months']))
        self.assertEqual(self.table.loc['p1', 'stage'], 3)
        self.assertEqual(self.table.loc['p1', 'age'], 70)
        self.assertTrue(pd.isna(self.table.loc['p2', 'age']))

    def test_expression_group_split_at_median(self):
        self.assertEqual(self.table.loc['p2', 'expression_group'], 'High')
        self.assertEqual(self.table.loc['p1', 'expression_group'], 'Low')
        self.assertTrue(pd.isna(self.table.loc['p4', 'expression_group']))

    def test_parsers(self):
        self.assertEqual(parseEvent('1:PROGRESSION'), 1)
        self.assertEqual(parseEvent('0:CENSORED'), 0)
        self.assertTrue(np.isnan(parseEvent('[Not Available]')))
        self.assertEqual(parseStage('STAGE IV'), 4)
        self.assertEqual(parseStage('STAGE IB'), 1)
        self.assertTrue(np.isnan(parseStage('[Discrepancy]')))


class EngineTests(unittest.TestCase):
    def setUp(self):
        self.table = syntheticCohort()

    def test_test_is_picked_by_variable_types(self):
        expected = {
            ('expression_z', 'age'): 'Spearman',
            ('expression_z', 'sex'): 'Mann-Whitney',
            ('cna_group', 'expression_z'): 'Kruskal-Wallis',
            ('sex', 'cna_group'): 'chi-square',
            ('expression_z', 'overall_survival'): 'Cox PH',
            ('overall_survival', 'cna_group'): 'log-rank',
            ('stage', 'expression_z'): 'Spearman',
        }
        for (varA, varB), test in expected.items():
            with self.subTest(pair=(varA, varB)):
                self.assertEqual(correlate(self.table, varA, varB, testTypes)['test'], test)

    def test_real_signals_are_found(self):
        cox = correlate(self.table, 'expression_z', 'overall_survival', testTypes)
        self.assertLess(cox['p_value'], 0.001)
        self.assertGreater(cox['effect_size'], 1)
        self.assertLess(correlate(self.table, 'sex', 'expression_z', testTypes)['p_value'], 0.001)

    def test_small_groups_are_skipped(self):
        row = correlate(self.table, 'mutated', 'overall_survival', testTypes)
        self.assertTrue(np.isnan(row['p_value']))
        self.assertIn('fewer than 2 groups', row['note'])

    def test_no_events_is_skipped(self):
        table = self.table.assign(os_event=0.0)
        row = correlate(table, 'expression_z', 'overall_survival', testTypes)
        self.assertTrue(np.isnan(row['p_value']))
        self.assertIn('no events', row['note'])

    def test_missing_values_dropped_per_pair(self):
        table = self.table.copy()
        table.loc[:49, 'age'] = np.nan
        self.assertEqual(correlate(table, 'age', 'expression_z', testTypes)['n'], 150)
        self.assertEqual(correlate(table, 'sex', 'expression_z', testTypes)['n'], 200)

    def test_correlate_all_adds_q_values(self):
        results = correlateAll(self.table, testTypes)
        self.assertEqual(len(results), 21)
        tested = results['p_value'].notna()
        self.assertTrue((results.loc[tested, 'q_value'] >= results.loc[tested, 'p_value']).all())
        self.assertTrue(results.loc[~tested, 'q_value'].isna().all())

    def test_against_limits_pairs(self):
        results = correlateAll(self.table, testTypes, against=['overall_survival'])
        self.assertEqual(len(results), 6)
        self.assertTrue(((results['varA'] == 'overall_survival') | (results['varB'] == 'overall_survival')).all())

    def test_same_source_pairs_are_not_tested(self):
        types = {'copy_number': 'ordinal', 'cna_group': 'categorical'}
        table = pd.DataFrame({'copy_number': [0, 1] * 20, 'cna_group': ['Diploid', 'Gain/Amp'] * 20})
        self.assertEqual(len(correlateAll(table, types)), 0)

    def test_adjusted_cox(self):
        summary, n, events = coxAnalysis(self.table, ['expression_z', 'age', 'sex', 'stage'], types=testTypes)
        self.assertEqual(list(summary.index), ['expression_z', 'age', 'stage', 'sex_Male'])
        self.assertEqual(n, 200)
        self.assertLess(summary.loc['expression_z', 'p_value'], 0.001)


if __name__ == '__main__':
    unittest.main()
