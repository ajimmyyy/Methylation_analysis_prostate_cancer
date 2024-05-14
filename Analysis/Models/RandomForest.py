from configparser import ConfigParser
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.feature_selection import RFECV
from imblearn.over_sampling import RandomOverSampler, SMOTE, KMeansSMOTE, SVMSMOTE, ADASYN
from imblearn.pipeline import Pipeline
import utils
from MakeFile.FileSaver import FileSaver

if __name__ == "__main__":
    _configPath = "Analysis/Models/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    # filter out the CpG sites
    _aucDf = pd.read_csv(_config["Paths"]["AUC_GROUP_DATA_PATH"])
    _aucDf = _aucDf[_aucDf['DNAm'] == "hyper"]
    keepFeature = _aucDf["CpG"].tolist()
    
    # read the training and testing data
    _trainDf = pd.read_csv(_config["Paths"]["TRAIN_BETA_DATA_PATH"], index_col=0)
    _testDf = pd.read_csv(_config["Paths"]["TEST_BETA_DATA_PATH"], index_col=0)
    
    # split the training, testing data into X and Y
    _trainX, _trainY, _testX, _testY = utils.SpliteTrainTest(_trainDf, _testDf, 25, keepFeature)

    # oversample the training data
    _trainX, _trainY = SVMSMOTE().fit_resample(_trainX, _trainY) 
    print(_trainY.value_counts())

    # train the model
    _rfModel = RandomForestClassifier(n_estimators = 300, n_jobs=-1)
    _rfecv = RFECV(estimator=_rfModel, min_features_to_select=40, step=1, cv=5, scoring='f1', n_jobs=-1)
    _rfecv.fit(_trainX, _trainY)

    # test the model
    utils.TestModelPerformance(_rfecv, _trainX, _trainY)
    print()
    utils.TestModelPerformance(_rfecv, _testX, _testY)

    # feature selection
    print("Optimal number of features : %d" % _rfecv.n_features_)
    print("Ranking of features : %s" % _rfecv.ranking_)
    scores = _rfecv.cv_results_['mean_test_score']
    stds = _rfecv.cv_results_['std_test_score']
    plt.figure()
    plt.title('RFECV')
    plt.xlabel('Number of features selected')
    plt.ylabel('Cross validation score (F1)')
    plt.plot(range(1, len(scores) + 1), scores, marker='o', linestyle='-')
    plt.fill_between(range(1, len(scores) + 1),
                    scores - stds,
                    scores + stds,
                    alpha=0.2)
    plt.show()

    feature_names = _trainX.columns
    selected_feature_names = [feature_names[i] for i in range(len(feature_names)) if _rfecv.support_[i]]
    results_df = _aucDf[_aucDf['CpG'].isin(selected_feature_names)]
    FileSaver.SaveData(results_df, _config["Paths"]["RANDOM_FOREST_FEATURES_SELECTION_PATH"])

    # importance = _rfModel.feature_importances_
    # importance = pd.DataFrame({'CpG': _trainX.columns, 'Importance': importance})
    # importance["Std"] = np.std([tree.feature_importances_ for tree in _rfModel.estimators_], axis=0)
    # importance = importance.sort_values(by="Importance", ascending=False)
    # fig, ax = plt.subplots()
    # x = range(importance.shape[0])
    # y = importance["Importance"]
    # yerr = importance["Std"]
    # plt.bar(x, y, yerr=yerr, align="center")
    # plt.show()
    # FileSaver.SaveDataframe(importance, _config["Paths"]["RANDOM_FOREST_IMPORTANCES_PATH"])

    # save the model
    # treePath = _config.get('Paths', 'RANDOM_FOREST_TREE_PATH')
    # joblib.dump(_rfModel, treePath)


    # # Parameter Tuning
    # pipeline = Pipeline([
    # ('oversampling', RandomOverSampler()),
    # ('clf', RandomForestClassifier(random_state=42))
    # ])

    # param_grid = {
    #     'oversampling': [KMeansSMOTE()],
    #     'oversampling__sampling_strategy': ['0.2', '0.3', '0.5', '0.7', 'auto'],
    #     'clf__n_estimators': [50, 100, 200],
    #     'clf__max_features': [2, 3, 4, 5, 10]
    # }

    # grid_search = GridSearchCV(pipeline, param_grid=param_grid, cv=10, scoring='f1', n_jobs=-1)
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