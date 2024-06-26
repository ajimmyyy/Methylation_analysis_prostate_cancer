from configparser import ConfigParser
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import joblib
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.metrics import accuracy_score, f1_score, recall_score
from xgboost import plot_tree
from MakeFile.FileSaver import FileSaver
from RandomForest import TransformTrainData

if __name__ == "__main__":
    _configPath = "Analysis/Models/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    MODEL_PATH = _config["Paths"]["RANDOM_FOREST_TREE_PATH"]
    SAVE_PATH = _config["Paths"]["RANDOM_FOREST_IMPORTANCES_PATH"]

    # filter out the CpG sites
    _aucDf = pd.read_csv(_config["Paths"]["MEAN_HYPER_WARD_CHOOSE_PATH"], usecols=["CpG", "DNAm"])
    _aucDf = _aucDf[_aucDf['DNAm'] == "hyper"]
    keepFeature = _aucDf["CpG"].tolist()
    keepFeature.append("cancer")

    # read the testing data
    _df = pd.read_csv(_config["Paths"]["TEST_BETA_DATA_PATH"], index_col=0)
    _df = TransformTrainData(_df, 25)
    _df = _df[_df.columns.intersection(keepFeature)]
    _df = _df.iloc[1:]

    # read the testing data
    _testDf = pd.read_csv(_config["Paths"]["850K_DATA_PATH"], index_col=0)
    _testDf = TransformTrainData(_testDf, 57)
    _testDf = _testDf[_testDf.columns.intersection(keepFeature)]
    _testDf = _testDf.iloc[1:]
    _testDf = _testDf[_df.columns]

    # split the training, testing data into X and Y
    _testX = _testDf.drop(columns=["cancer"])
    _testY = _testDf["cancer"]

    _model = joblib.load(MODEL_PATH)

    trainPredicted = _model.predict(_testX)

    print("training results:")
    print("Accuracy: ", accuracy_score(_testY, trainPredicted))
    print("F1: ", f1_score(_testY, trainPredicted, average = "binary"))
    print("Recall: ", recall_score(_testY, trainPredicted, average = "binary"))
    print("Specificity: ", recall_score(_testY, trainPredicted, average = "binary", pos_label=0))

    _result = permutation_importance(_model, _testX, _testY, n_repeats=10, random_state=0)

    _importancesMean = _result.importances_mean
    _importancesStd = _result.importances_std
    _importance = _model.feature_importances_

    _sortedIndices = np.argsort(_importancesMean)
    _sortedScores = _importancesMean[_sortedIndices]

    data = {'CpG': _testX.columns, 'FeatureImportance': _importance, 'ImportanceMean': _importancesMean, 'ImportanceStd': _sortedScores}
    df = pd.DataFrame(data)
    FileSaver.SaveData(df, SAVE_PATH)

    plt.figure(figsize=(10, 6))
    plt.boxplot(_result.importances[_sortedIndices].T, vert=False,
            labels=np.array(_testX.columns)[_sortedIndices])
    plt.xlabel('Importance Score')
    plt.title('Permutation Importance')
    plt.show()
