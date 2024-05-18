import unittest
import pandas as pd
import Models.CutpointCalculator
import Models.DeltaBetaFilter

class CutpointCalculatorTest(unittest.TestCase):
    def setUp(self):
        data = pd.DataFrame(
            {
                'CpG': [1, 2, 3],
                '1': [0, 0.5, 1],
                '2': [1, 1, 0.5],
                '3': [1, 1, 0.5],
            }
        )
        select_hyper = pd.DataFrame(
            {
                'CpG': [2]
            }
        )
        select_hypo = pd.DataFrame(
            {
                'CpG': [3]
            }
        )
        self.cutpointCul = Models.CutpointCalculator.CutpointCalculator()
        self.data = data
        self.select_hyper = select_hyper
        self.select_hypo = select_hypo
    
    def test_calculate_cutpoint(self):
        result_min = self.cutpointCul.CalculateCutpoint(self.data, self.select_hyper, 1, "hyper", "CpG", "min")
        result_mid = self.cutpointCul.CalculateCutpoint(self.data, self.select_hyper, 1, "hyper", "CpG", "mid")
        result_max = self.cutpointCul.CalculateCutpoint(self.data, self.select_hyper, 1, "hyper", "CpG", "max")
        self.assertEqual(result_min["cutpoint"].tolist(), [0.51])
        self.assertEqual(result_mid["cutpoint"].tolist(), [0.75])
        self.assertEqual(result_max["cutpoint"].tolist(), [0.99])
    
    def test_calculate_cutpoint_multi(self):
        result = self.cutpointCul.CalculateCutpoint(self.data, [self.select_hyper, self.select_hypo], 1, ["hyper", "hypo"], "CpG", "min")
        self.assertEqual(result["cutpoint"].tolist(), [0.51, 0.51])