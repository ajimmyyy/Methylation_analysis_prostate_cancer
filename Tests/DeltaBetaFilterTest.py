import unittest
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
        self.dBetaFilter = Models.DeltaBetaFilter.DeltaBetaFilter()
        self.data = data
    
    def test_CalculateDeltaBeta(self):
        result = self.dBetaFilter.CalculateDeltaBeta(self.data, 1)
        self.assertEqual(result['dbeta'].tolist(), [0, 0.5, -0.5])
        self.assertEqual(result['avg'].tolist(), [0, 0.5, 1])

if __name__ == '__main__':
    unittest.main()