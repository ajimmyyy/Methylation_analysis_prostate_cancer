from configparser import ConfigParser
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import joblib
from MakeFile.FileSaver import FileSaver
from RandomForest import TransformTrainData
from sklearn.metrics import accuracy_score, f1_score

if __name__ == "__main__":
    _configPath = "Analysis/Models/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    MODEL_PATH = _config["Paths"]["RANDOM_FOREST_TREE_PATH"]

    # # filter out the CpG sites
    # _aucDf = pd.read_csv(_config["Paths"]["AUC_GROUP_DATA_PATH"], usecols=["CpG", "DNAm"])
    # _aucDf = _aucDf[_aucDf['DNAm'] == "hyper"]
    # keepFeature = _aucDf["CpG"].tolist()
    # keepFeature.append("cancer")
    keepFeature = pd.read_csv(_config["Paths"]["RANDOM_FOREST_FEATURES_SELECTION_PATH"])["CpG"].tolist()

    # read the testing data
    _df = pd.read_csv(_config["Paths"]["TEST_BETA_DATA_PATH"], index_col=0)
    _df = TransformTrainData(_df, 25)
    _df = _df[_df.columns.intersection(keepFeature)]
    _df = _df.iloc[1:]

    # # read the testing data
    # _testDf = pd.read_csv(_config["Paths"]["450K_DATA_PATH"], index_col=0)
    # _testDf = TransformTrainData(_testDf, 10)
    # _testDf = _testDf[_testDf.columns.intersection(keepFeature)]
    # _testDf = _testDf.iloc[1:]
    # _testDf = _testDf[_df.columns]

    # # split the training, testing data into X and Y
    # _testX = _testDf.drop(columns=["cancer"])
    # _testY = _testDf["cancer"]

    _model = joblib.load(MODEL_PATH)

    # testPredicted = _model.predict(_testX)
    # accuracy = accuracy_score(_testY, testPredicted)
    # f1 = f1_score(_testY, testPredicted, average = "binary")
    # print(accuracy)
    # print("F1: ", f1)

    importance = _model.feature_importances_
    importance = pd.DataFrame({'CpG': _df.columns, 'Importance': importance})
    importance["Std"] = np.std([tree.feature_importances_ for tree in _model.estimators_], axis=0)
    importance = importance.sort_values(by="Importance", ascending=False)
    fig, ax = plt.subplots()
    x = range(importance.shape[0])
    y = importance["Importance"]
    yerr = importance["Std"]
    plt.bar(x, y, yerr=yerr, align="center")
    plt.show()
    FileSaver.SaveData(importance, "C:/Users/acer/Desktop/test/FeatureImportance.csv")
