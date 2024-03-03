from configparser import ConfigParser
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import joblib
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, f1_score
from sklearn.model_selection import GridSearchCV
from sklearn.tree import plot_tree
from collections import OrderedDict
from imblearn.over_sampling import RandomOverSampler, SMOTE, KMeansSMOTE, SVMSMOTE, ADASYN
from imblearn.pipeline import Pipeline
from MakeFile.FileSaver import FileSaver

def TransformTrainData(df, normalNum):
    target = [0] * normalNum + [1] * (len(df.columns) - normalNum)
    df.loc['cancer'] = target
    df = df.T
    return df

if __name__ == "__main__":
    _configPath = "Analysis/RandomForest/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    # filter out the CpG sites
    _aucDf = pd.read_csv(_config["Paths"]["AUC_GROUP_DATA_PATH"], usecols=["CpG", "DNAm"])
    _aucDf = _aucDf[_aucDf['DNAm'] == "hyper"]
    keepFeature = _aucDf["CpG"].tolist()
    keepFeature.append("cancer")


    # read the training data
    _trainDf = pd.read_csv(_config["Paths"]["TRAIN_BETA_DATA_PATH"], index_col=0)
    _trainDf = TransformTrainData(_trainDf, 50)
    _trainDf = _trainDf[_trainDf.columns.intersection(keepFeature)]
    _trainDf = _trainDf.iloc[1::2]


    # read the testing data
    _testDf = pd.read_csv(_config["Paths"]["TEST_BETA_DATA_PATH"], index_col=0)
    _testDf = TransformTrainData(_testDf, 50)
    _testDf = _testDf[_testDf.columns.intersection(keepFeature)]
    _testDf = _testDf.iloc[1::2]


    # split the training, testing data into X and Y
    _trainX = _trainDf.drop(columns=["cancer"])
    _trainY = _trainDf["cancer"]
    _testX = _testDf.drop(columns=["cancer"])
    _testY = _testDf["cancer"]


    # oversample the training data
    _trainX, _trainY = KMeansSMOTE().fit_resample(_trainX, _trainY) 


    # train the model
    # forest = RandomForestClassifier(n_estimators = 200, min_samples_leaf = 10, n_jobs=-1)
    forest = RandomForestClassifier(n_estimators = 100, max_features = 30, n_jobs = -1)
    forest.fit(_trainX, _trainY)


    trainPredicted = forest.predict(_trainX)
    accuracy = accuracy_score(_trainY, trainPredicted)
    f1 = f1_score(_trainY, trainPredicted, average = "binary")
    print(accuracy)
    print("F1: ", f1)


    # test the model
    testPredicted = forest.predict(_testX)
    accuracy = accuracy_score(_testY, testPredicted)
    f1 = f1_score(_testY, testPredicted, average = "binary")
    print(accuracy)
    print("F1: ", f1)

    importances = forest.feature_importances_
    forest_importances = pd.Series(importances,index=_trainX.columns).sort_values(ascending=True)
    forest_importances_top = forest_importances[-50:]

    fig, ax = plt.subplots()
    forest_importances_top.plot.barh(ax=ax)
    ax.set_title("Feature importances using MDI")
    ax.set_xlabel("Mean decrease in impurity")
    fig.tight_layout()
    plt.show()
    FileSaver.SavePlot(fig, _config["Paths"]["RANDOM_FOREST_IMPORTANCES_PATH"])

    # save the model
    treePath = _config.get('Paths', 'RANDOM_FOREST_TREE_PATH')
    joblib.dump(forest, treePath)


    # RANDOM_STATE = 123
    # ensemble_clfs = [
    #     (
    #         "RandomForestClassifier, max_features=5",
    #         RandomForestClassifier(
    #             warm_start=True,
    #             oob_score=True,
    #             max_features=5,
    #             random_state=RANDOM_STATE,
    #         ),
    #     ),
    #     (
    #         "RandomForestClassifier, max_features=20",
    #         RandomForestClassifier(
    #             warm_start=True,
    #             max_features=20,
    #             oob_score=True,
    #             random_state=RANDOM_STATE,
    #         ),
    #     ),
    #     (
    #         "RandomForestClassifier, max_features=30",
    #         RandomForestClassifier(
    #             warm_start=True,
    #             max_features=30,
    #             oob_score=True,
    #             random_state=RANDOM_STATE,
    #         ),
    #     ),(
    #         "RandomForestClassifier, max_features=40",
    #         RandomForestClassifier(
    #             warm_start=True,
    #             max_features=40,
    #             oob_score=True,
    #             random_state=RANDOM_STATE,
    #         ),
    #     ),
    # ]
    # error_rate = OrderedDict((label, []) for label, _ in ensemble_clfs)
    # min_estimators = 100
    # max_estimators = 150
    # for label, clf in ensemble_clfs:
    #     for i in range(min_estimators, max_estimators + 1, 5):
    #         clf.set_params(n_estimators=i)
    #         clf.fit(_trainX, _trainY)

    #         # Record the OOB error for each `n_estimators=i` setting.
    #         oob_error = 1 - clf.oob_score_
    #         error_rate[label].append((i, oob_error))

    # # Generate the "OOB error rate" vs. "n_estimators" plot.
    # for label, clf_err in error_rate.items():
    #     xs, ys = zip(*clf_err)
    #     plt.plot(xs, ys, label=label)

    # plt.xlim(min_estimators, max_estimators)
    # plt.ylim(0.01, 0.04)
    # plt.xlabel("n_estimators")
    # plt.ylabel("OOB error rate")
    # plt.legend(loc="upper right")
    # plt.show()