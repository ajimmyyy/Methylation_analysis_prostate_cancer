from configparser import ConfigParser
import time
import numpy as np
import pandas as pd
from xgboost import XGBClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFECV
from boruta import BorutaPy
from sklearn.model_selection import cross_val_score, StratifiedKFold
from imblearn.over_sampling import SVMSMOTE
import utils

RANDOM_SEED = 42

class ModelTrainer:
    def __init__(self, trainDf, testDf, normalCount, keepFeature):
        trainX, trainY, testX, testY = utils.SpliteTrainTest(trainDf, testDf, normalCount, keepFeature)
        self.trainX = trainX
        self.trainY = trainY
        self.testX = testX
        self.testY = testY.values
        self.selector = {}
        self.featureSelected = {}

    def TrainModel(self, models, method = 'RFECV', cv = 5, randomSeed = RANDOM_SEED):
        X, y = SVMSMOTE(random_state=randomSeed).fit_resample(self.trainX, self.trainY)
        results = None

        start_time = time.time()
        if method == 'RFECV':
            results = self.TrainModelRFECV(models, X, y, cv)
        elif method == 'Boruta':
            results = self.TrainModelBoruta(models, X, y, randomSeed)
        print(f"Training time: {time.time() - start_time}")

        return results

    def TrainModelRFECV(self, models, X, y, cv):
        results = {}
        for modelName, model in models.items():
            selector = RFECV(estimator=model, min_features_to_select=2, step=1, cv=cv, scoring='f1', n_jobs=-1)
            selector.fit(X, y)

            print(modelName)
            accuracy, f1, recall, specificity = utils.TestModelPerformance(selector, self.testX, self.testY)
            self.selector[modelName] = selector

            features = X.columns
            features = [features[i] for i in range(len(features)) if selector.support_[i]]
            self.featureSelected[modelName] = features

            results[modelName] = {
                "accuracy": accuracy,
                "f1": f1,
                "recall": recall,
                "specificity": specificity,
                "feature_num": len(features)
            }

        resultsDf = pd.DataFrame(results).T
        return resultsDf

    def TrainModelBoruta(self, models, X, y, randomSeed):
        np.int = np.int32
        np.float = np.float64
        np.bool = np.bool_
        trainX = X.values
        trainY = y.values
        results = {}

        for modelName, model in models.items():
            selector = BorutaPy(model, n_estimators='auto', verbose=2, random_state=randomSeed)
            selector.fit(trainX, trainY)
            self.selector[modelName] = selector

            features = X.columns
            features = [features[i] for i in range(len(features)) if selector.support_[i]]
            self.featureSelected[modelName] = features

            model.fit(X.loc[:, features], y)
            print(modelName)
            accuracy, f1, recall, specificity = utils.TestModelPerformance(model, self.testX.loc[:, features], self.testY)
            results[modelName] = {
                "accuracy": accuracy,
                "f1": f1,
                "recall": recall,
                "specificity": specificity,
                "feature_num": len(features)
            }

            
        resultsDf = pd.DataFrame(results).T
        return resultsDf

    def DrawSelectorPlot(self):
        for modelName, selector in self.selector.items():
            utils.DrawRFECVPlot(selector)

    def CrossFeatures(self):
        if self.featureSelected == {}:
            raise Exception("Please train the model first")
        featureSets = [set(features) for features in self.featureSelected.values()]
        commonFeatures = set.intersection(*featureSets) if featureSets else set()

        return commonFeatures

    def GetSelectedFeatures(self):
        return self.featureSelected.copy()
    
    def TestModelPerformance(self, models, cv):
        results = {}
        _cv = StratifiedKFold(n_splits=cv, shuffle=True, random_state=RANDOM_SEED)
        for modelName, model in models.items():
            X = self.trainX.loc[:, self.featureSelected[modelName]]
            scores = self.EvaluateModelCv(model, X, self.trainY, _cv)
            results[modelName] = scores
        return results

    def EvaluateModelCv(self, model, X, y, cv):
        scores = cross_val_score(model, X, y, cv=cv, scoring='f1', n_jobs=-1)
        return scores

def OverlapRatio(list1, list2):
    set1 = set(list1)
    set2 = set(list2)
    
    intersection = set1 & set2
    union = set1 | set2
    
    ratio_relative_to_first_list = len(intersection) / len(set1)
    ratio_relative_to_union = len(intersection) / len(union)
    
    return ratio_relative_to_first_list, ratio_relative_to_union

if __name__ == "__main__":
    configPath = "Analysis/Models/config.ini"
    config = ConfigParser()
    config.read(configPath)

    aucDf = pd.read_csv(config["Paths"]["AUC_GROUP_DATA_PATH"])
    aucDf = aucDf[aucDf['DNAm'] == "hyper"]
    keepFeature = aucDf["CpG"].tolist()

    # read the training and testing data
    trainDf = pd.read_csv(config["Paths"]["TRAIN_BETA_DATA_PATH"], index_col=0)
    testDf = pd.read_csv(config["Paths"]["TEST_BETA_DATA_PATH"], index_col=0)

    # train the model
    modelTrainer = ModelTrainer(trainDf, testDf, 25, keepFeature)
    models = {
        'RF': RandomForestClassifier(
            n_estimators=150, 
            max_depth=10, 
            max_leaf_nodes=5, 
            max_features='sqrt', 
            n_jobs=-1,
            random_state=RANDOM_SEED
        ),
        'XGB': XGBClassifier(
            n_estimators=150,
            subsample=0.7,
            max_depth=5,
            learning_rate=0.1,
            colsample_bytree=0.7,
            reg_alpha=1,
            n_jobs=-1,
            random_state=RANDOM_SEED
        )
    }

    testModels = {
        'RF': RandomForestClassifier( 
            n_jobs=-1,
            random_state=RANDOM_SEED
        ),
        'XGB': XGBClassifier(
            n_jobs=-1,
            random_state=RANDOM_SEED
        )
    }

    modelTrainer.TrainModel(models, method='Boruta')
    boruta_features = modelTrainer.GetSelectedFeatures()
    boruta_score = modelTrainer.TestModelPerformance(testModels, cv=10)
    
    modelTrainer.TrainModel(models, method='RFECV')
    rfecv_features = modelTrainer.GetSelectedFeatures()
    rfecv_score = modelTrainer.TestModelPerformance(testModels, cv=10)


    for modelName in boruta_features.keys():
        print(f"Model: {modelName}")
        print("RFECV")
        print(f"Features: {rfecv_features[modelName]}")
        print(f"Score: {rfecv_score[modelName]}")
        print(f'Average: {np.mean(rfecv_score[modelName])}')
        print("\n")
        print("Boruta")
        print(f"Features: {boruta_features[modelName]}")
        print(f"Score: {boruta_score[modelName]}")
        print(f'Average: {np.mean(boruta_score[modelName])}')
        print("\n")

    # crossFeatures = modelTrainer.CrossFeatures()
    # print(crossFeatures)