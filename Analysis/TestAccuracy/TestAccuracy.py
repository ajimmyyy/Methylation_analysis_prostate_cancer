from configparser import ConfigParser
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, f1_score, recall_score
from MakeFile.FileSaver import FileSaver

class TestAccuracy:
    def __init__(self, betaDf: pd.DataFrame, cutpoints: pd.DataFrame, trueLabels: np.array):
        self.betaDf = betaDf
        self.cutpoints = cutpoints
        self.trueLabels = trueLabels
        self.predictionEach = {}
        self.predictions = None

    def PredictWithVoting(self):
        predictions = []
        featurePerformance = {feature: [] for feature in cutpoints['CpG']}
        for _, row in self.betaDf.T.iterrows():
            votes = 0
            for _, cutpoint_row in self.cutpoints.iterrows():
                feature = cutpoint_row['CpG']
                cutpoint = cutpoint_row['cutpoint']
                if row[feature] > cutpoint:
                    featurePerformance[feature].append(1)
                    votes += 1
                else:
                    featurePerformance[feature].append(0)
                    votes -= 1
            if votes > 0:
                predictions.append(1)
            else:
                predictions.append(0)

        self.predictions = np.array(predictions)
        self.predictionEach = featurePerformance
        return np.array(predictions)

    def DetailPredict(self):
        performanceDf = pd.DataFrame(self.predictionEach)
        performanceDf['True Label'] = self.trueLabels
        error = performanceDf.copy()
        for feature in error.columns:
            error[feature] = error[feature].astype(bool) ^ error['True Label'].astype(bool)
    
        plt.figure(figsize=(12, 8))
        sns.heatmap(performanceDf.T, annot=False, cmap='coolwarm', cbar=True, 
                mask=error.T, fmt=".0f", annot_kws={"size": 12},
                xticklabels=performanceDf.index.tolist() + ['True Label'], 
                yticklabels=performanceDf.columns)

        plt.title('Feature Prediction Performance for Each Sample')
        plt.xlabel('Sample Index')
        plt.ylabel('Feature')
        plt.show()

if __name__ == "__main__":
    configPath = "Analysis/TestAccuracy/config.ini"
    config = ConfigParser()
    config.read(configPath)

    betaDataDf = pd.read_csv(config["Paths"]["850K_DATA_PATH"], index_col=0)
    testGene = {
        'CpG': ['cg18759209', 'cg24530250', 'cg15229124'],
        'cutpoint': [0.29, 0.11, 0.29]
    }
    cutpoints = pd.DataFrame(testGene)
    betaDataDf = betaDataDf[betaDataDf.index.isin(cutpoints['CpG'])]
    trueLabels = np.array([0]*57 + [1]*57)

    test = TestAccuracy(betaDataDf, cutpoints, trueLabels)
    predictions = test.PredictWithVoting()

    accuracy = accuracy_score(trueLabels, predictions)
    f1 = f1_score(trueLabels, predictions)
    recall = recall_score(trueLabels, predictions)
    specificity = recall_score(trueLabels, predictions, pos_label=0)

    print(f'Accuracy: {accuracy}')
    print(f'F1-score: {f1}')
    print(f'Recall: {recall}')
    print(f'Specificity: {specificity}')
    test.DetailPredict()