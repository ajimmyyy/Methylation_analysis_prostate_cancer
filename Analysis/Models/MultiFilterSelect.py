from configparser import ConfigParser
import numpy as np
import pandas as pd
from xgboost import XGBClassifier
from sklearn.metrics import f1_score, accuracy_score, recall_score
from sklearn.ensemble import RandomForestClassifier, VotingClassifier
from sklearn.linear_model import LogisticRegression 
from imblearn.over_sampling import SVMSMOTE
from sklearn.feature_selection import RFECV
from RandomForest import TransformTrainData
from MakeFile.FileSaver import FileSaver

# def FilterFeature(X: pd.DataFrame, y: pd.DataFrame, models, threshold = 0.1):
#     if not isinstance(models, list):
#         models = [models]

#     importances_list = []
#     for model in models:
#         keepNum = int(X.shape[1] * (1 - threshold))
#         rfecv = RFECV(estimator=model, min_features_to_select=keepNum, step=1, cv=5, scoring='f1', n_jobs=-1)
#         rfecv.fit(_trainX, _trainY)
#         importances_list.append(rfecv.support_)
        
#     combined_importances = np.ones(X.shape[1], dtype=bool)
#     for importance_mask in importances_list:
#         combined_importances = np.logical_or(combined_importances, importance_mask)

#     return X.iloc[:, combined_importances]

if __name__ == "__main__":
    _configPath = "Analysis/Models/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    # filter out the CpG sites
    _aucDf = pd.read_csv(_config["Paths"]["AUC_GROUP_DATA_PATH"])
    _aucDf = _aucDf[_aucDf['DNAm'] == "hyper"]
    keepFeature = _aucDf["CpG"].tolist()
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
    _trainX, _trainY = SVMSMOTE().fit_resample(_trainX, _trainY) 
    print(_trainY.value_counts())

    # filter the features
    voteClf = VotingClassifier(estimators=[
        ('rf', RandomForestClassifier(n_estimators=100)), 
        ('xgb', XGBClassifier(n_estimators=100)),
        ('lr', LogisticRegression(penalty='l1', C=1.0, solver='liblinear')),
        ],
        voting='soft',
        n_jobs=-1
    )
    voteClf.fit(_trainX, _trainY)
    
    trainPredicted = voteClf.predict(_trainX)
    accuracy = accuracy_score(_trainY, trainPredicted)
    f1 = f1_score(_trainY, trainPredicted, average = "binary")
    print(accuracy)
    print("F1: ", f1)
    
    trainPredicted = voteClf.predict(_testX)
    accuracy = accuracy_score(_testY, trainPredicted)
    f1 = f1_score(_testY, trainPredicted, average = "binary")
    print(accuracy)
    print("F1: ", f1)
    print("Recall: ", recall_score(_testY, trainPredicted, average = "binary"))
    print("Specificity: ", recall_score(_testY, trainPredicted, average = "binary", pos_label=0))