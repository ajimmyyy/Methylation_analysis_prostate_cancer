import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import re

class CancerStageCalculator:
    def CancerStageSplit(self, allBetaDf: pd.DataFrame, sampleSheetDf: pd.DataFrame, cancerStage, stageColName = "clinical_T", sampleColName = "sample"):
        sampleSheetDf[stageColName].fillna('Unknown', inplace=True)

        if isinstance(cancerStage, str):
            cancerStage = [cancerStage]

        filteredSamples = sampleSheetDf[sampleSheetDf[stageColName].str.startswith(tuple(cancerStage))][sampleColName]
        filteredSamples = filteredSamples.astype(str)
        stageDf = allBetaDf.loc[:, allBetaDf.columns.isin(filteredSamples)]

        return stageDf