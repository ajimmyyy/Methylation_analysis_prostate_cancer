from configparser import ConfigParser
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import f1_score
from sklearn.feature_selection import RFECV
from xgboost import XGBClassifier
from imblearn.over_sampling import SVMSMOTE
import shap
import optuna
import utils
from MakeFile.FileSaver import FileSaver

def objective(trial):
    params = {
        "verbosity": 0,
        "objective": "binary:logistic",
        "tree_method": "exact",
        "booster": "gbtree",
        "lambda": trial.suggest_float("lambda", 1e-8, 1.0),
        "alpha": trial.suggest_float("alpha", 1e-8, 1.0),
        "subsample": trial.suggest_float("subsample", 0.2, 1.0),
        "colsample_bytree": trial.suggest_float("colsample_bytree", 0.2, 1.0),
        "max_depth" : trial.suggest_int("max_depth", 1, 9),
        "min_child_weight" : trial.suggest_int("min_child_weight", 2, 10),
        "eta" : trial.suggest_float("eta", 1e-8, 1.0)
    }
    xgboostModel = XGBClassifier(**params)
    predicted = xgboostModel.fit(_trainX, _trainY)
    f1 = f1_score(_testY, predicted, average = "binary")
    return f1

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
    _xgboostModel = XGBClassifier(
        n_estimators = 500, 
        max_depth = 4, 
        min_child_weight = 2, 
        subsample = 0.9,
        colsample_bytree = 0.8,
        reg_alpha= 0.1,
        reg_lambda= 1,
        learning_rate= 0.05,
        objective= 'binary:logistic'
    )
    eval_set = [(_testX, _testY)]
    _rfecv = RFECV(estimator=_xgboostModel, min_features_to_select=40, step=1, cv=5, scoring='f1', n_jobs=-1)
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
    FileSaver.SaveData(results_df, _config["Paths"]["XGBOOST_FEATURES_SELECTION_PATH"])


    # # Parameter optimization
    # study = optuna.create_study(direction="maximize")
    # study.optimize(objective, n_trials=100)

    # # Showing optimization results
    # print('Number of finished trials:', len(study.trials))
    # print('Best trial parameters:', study.best_trial.params)
    # print('Best score:', study.best_value)
