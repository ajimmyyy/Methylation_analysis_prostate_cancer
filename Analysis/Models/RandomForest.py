from configparser import ConfigParser
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import joblib
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, f1_score
from sklearn.model_selection import GridSearchCV
from sklearn.tree import plot_tree
from sklearn import tree
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
    _configPath = "Analysis/Models/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    # filter out the CpG sites
    # _aucDf = pd.read_csv(_config["Paths"]["AUC_GROUP_DATA_PATH"], usecols=["CpG", "DNAm"])
    # _aucDf = _aucDf[_aucDf['DNAm'] == "hyper"]
    # keepFeature = _aucDf["CpG"].tolist()
    # keepFeature.append("cancer")
    _featureDf = pd.read_csv(_config["Paths"]["LASSO_IMPORTANCES_PATH"], index_col=0)
    keepFeature = _featureDf.index.tolist()
    keepFeature = keepFeature[:46]
    keepFeature.append("cancer")

    # read the training data
    _trainDf = pd.read_csv(_config["Paths"]["TRAIN_BETA_DATA_PATH"], index_col=0)
    _trainDf = TransformTrainData(_trainDf, 25)
    _trainDf = _trainDf[_trainDf.columns.intersection(keepFeature)]
    _trainDf = _trainDf.iloc[1:]


    # read the testing data
    _testDf = pd.read_csv(_config["Paths"]["TEST_BETA_DATA_PATH"], index_col=0)
    _testDf = TransformTrainData(_testDf, 25)
    _testDf = _testDf[_testDf.columns.intersection(keepFeature)]
    _testDf = _testDf.iloc[1:]


    # split the training, testing data into X and Y
    _trainX = _trainDf.drop(columns=["cancer"])
    _trainY = _trainDf["cancer"]
    _testX = _testDf.drop(columns=["cancer"])
    _testY = _testDf["cancer"]


    # oversample the training data
    _trainX, _trainY = KMeansSMOTE().fit_resample(_trainX, _trainY) 

    # train the model
    # forest = RandomForestClassifier(n_estimators = 200, min_samples_leaf = 10, n_jobs=-1)
    forest = RandomForestClassifier(n_estimators = 50, n_jobs = -1)
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

    # save the model
    # treePath = _config.get('Paths', 'RANDOM_FOREST_TREE_PATH')
    # joblib.dump(forest, treePath)
    importances = forest.feature_importances_
    std = np.std([tree.feature_importances_ for tree in forest.estimators_], axis=0)
    forest_importances = pd.Series(importances, index=_testX.columns)
    fig, ax = plt.subplots()
    forest_importances.plot.bar(yerr=std, ax=ax)
    ax.set_title("Feature importances using MDI")
    ax.set_ylabel("Mean decrease in impurity")
    fig.tight_layout()
    plt.show()


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
    # max_estimators = 200
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
    # # plt.ylim(0.01, 0.04)
    # plt.xlabel("n_estimators")
    # plt.ylabel("OOB error rate")
    # plt.legend(loc="upper right")
    # plt.show()


    # pipeline = Pipeline([
    # ('oversampling', RandomOverSampler()),
    # ('clf', RandomForestClassifier(random_state=42))
    # ])

    # param_grid = {
    #     'oversampling': [RandomOverSampler(), SMOTE(), KMeansSMOTE(), SVMSMOTE(), ADASYN()],
    #     'clf__n_estimators': [50, 100, 200],
    #     'clf__max_features': [None, 10, 20, 30, 40 ,50]
    # }

    # grid_search = GridSearchCV(pipeline, param_grid=param_grid, cv=5, scoring='f1', n_jobs=-1)
    # grid_search.fit(_trainX, _trainY)

    # best_params = grid_search.best_params_
    # print("Best Parameters:", best_params)

    # results = grid_search.cv_results_
    # oversampling_names = [type(oversampler).__name__ for oversampler in param_grid['oversampling']]
    # param_combinations = [(oversampling_name, n_estimators, max_features) for oversampling_name in oversampling_names for n_estimators in param_grid['clf__n_estimators'] for max_features in param_grid['clf__max_features']]
    # f1_scores = results['mean_test_score']

    # plt.figure(figsize=(12, 8))
    # for i, oversampling_name in enumerate(oversampling_names):
    #     oversampling_indices = [idx for idx, (oversampling, _, _) in enumerate(param_combinations) if oversampling == oversampling_name]
    #     oversampling_f1_scores = [f1_scores[idx] for idx in oversampling_indices]
    #     plt.scatter(oversampling_indices, oversampling_f1_scores, label=oversampling_name)
    # plt.title('F1 Scores for Different Oversampling Methods')
    # plt.xlabel('Parameter Combinations')
    # plt.ylabel('F1 Score')
    # plt.xticks(np.arange(len(param_combinations)), [f'{oversampling}\n{n_estimators}\n{max_features}' for oversampling, n_estimators, max_features in param_combinations], rotation=45, ha='right')
    # plt.legend()
    # plt.tight_layout()
    # plt.show()