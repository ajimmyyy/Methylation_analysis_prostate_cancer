import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import re

class CancerStageCalculator:
    def MappingPatientTNM(self, cancerStageDf, IdDf):
        filtered_data = cancerStageDf[cancerStageDf['TCGA_barcode'].isin(IdDf['TCGA_barcode'])]
        df = pd.merge(IdDf, filtered_data, on = ["TCGA_barcode"], how = "outer")
        return df
    
    def MarkCancerStage(self, sampleSheet):
        tqdm.pandas(desc="find cancer stage")
        sampleSheet["cancer_stage"] = sampleSheet.progress_apply(self.__JudgeCancerStage, axis = 1)
        earlyPatient = sampleSheet[sampleSheet["cancer_stage"] == "early"].drop_duplicates()
        laterPatient = sampleSheet[sampleSheet["cancer_stage"] == "later"].drop_duplicates()

        return sampleSheet, earlyPatient, laterPatient
    
    def DrawDbetaDifference(self, dbataDifDf, binWidth, boundary, dbetaDifColumn = "dbate_dif"):
        fig = plt.figure(figsize=(8,6))
        std = dbataDifDf[dbetaDifColumn].std()
        mean = dbataDifDf[dbetaDifColumn].mean()

        hist, bins = np.histogram(dbataDifDf[dbetaDifColumn], bins = np.arange(-boundary, boundary + binWidth, binWidth))
        
        colors = ['gray' if (x < (mean - 2 * std) or x > (mean + 2 * std)) else 'blue' for x in np.arange(-boundary, boundary + binWidth, binWidth)]

        plt.bar(bins[:-1], hist, width=binWidth, align='edge', color = colors)
        plt.xlabel('delta bate difference')
        plt.ylabel('times frequency')

        return fig
    
    def __JudgeCancerStage(self, row):
        if row["Sample_Group"] == "Normal":
            return pd.Series({"cancer_stage": "None"})
        
        tumor_stage = str(row["tumor_stage"])
        node_stage = str(row["node_stage"])

        if  node_stage == "N0":
            if re.search(r"T1|T2", tumor_stage):
                return pd.Series({"cancer_stage": "early"})
            else:
                return pd.Series({"cancer_stage": "later"})
        else:
            return pd.Series({"cancer_stage": "later"})