import unittest
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import Models.DeltaBetaFilter


class TestDeltaBetaFilter(unittest.TestCase):   
    def setUp(self):
        data = pd.DataFrame(
            {
                'CpG': [1, 2, 3],
                '1': [0, 0.5, 1],
                '2': [0, 1, 0.5],
                '3': [0, 1, 0.5],
                '4': [0, 1, 0.5],
                '5': [0, 1, 0.5],
                '6': [1, 1, 0.5],
            }
        )
        dmp = pd.DataFrame(
            {
                'CpG': [1, 2, 3],
                'gene': ['A', 'A', 'B'],
                'feature': ['TSS200', 'TSS1500', 'none'],
                'adj.P.Val': [0.1, 0.04, 0.04],
            }
        )
        self.dBetaFilter = Models.DeltaBetaFilter.DeltaBetaFilter()
        self.data = data
        self.dmp = dmp
    
    def test_CalculateDeltaBeta(self):
        result = self.dBetaFilter.CalculateDeltaBeta(self.data, 1)
        self.assertEqual(result['dbeta'].tolist(), [0, 0.5, -0.5])
        self.assertEqual(result['avg'].tolist(), [0, 0.5, 1])

    def test_FilterDeltaBeta(self):
        dbeta = self.dBetaFilter.CalculateDeltaBeta(self.data, 1)
        result_t = self.dBetaFilter.FilterDeltaBeta(dbeta, self.dmp, True)
        result_f = self.dBetaFilter.FilterDeltaBeta(dbeta, self.dmp, False)
        self.assertEqual(result_t['CpG'].tolist(), [2])
        self.assertEqual(result_f['CpG'].tolist(), [2, 3])
    
    def test_DetermineDNAm(self):
        dbeta = self.dBetaFilter.CalculateDeltaBeta(self.data, 1)
        dbeta = self.dBetaFilter.FilterDeltaBeta(dbeta, self.dmp, False)
        result_hyper, result_hypo = self.dBetaFilter.DetermineDNAm(dbeta, 0.4, -0.4, 0.05)
        self.assertEqual(result_hyper['CpG'].tolist(), [2])
        self.assertEqual(result_hypo['CpG'].tolist(), [3])

    def test_DrawVolcanoPlot(self):
        dbeta = self.dBetaFilter.CalculateDeltaBeta(self.data, 1)
        dbeta = self.dBetaFilter.FilterDeltaBeta(dbeta, self.dmp, False)
        hyper, hypo = self.dBetaFilter.DetermineDNAm(dbeta, 0.4, -0.4, 0.05)
        result = self.dBetaFilter.DrawVolcanoPlot(dbeta, hyper, hypo)
        self.assertIsInstance(result, plt.Figure)

if __name__ == '__main__':
    unittest.main()