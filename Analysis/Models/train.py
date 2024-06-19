from configparser import ConfigParser
import pandas as pd
from xgboost import XGBClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFECV
from imblearn.over_sampling import SVMSMOTE
import utils
from MakeFile.FileSaver import FileSaver

class ModelTrainer:
    def __init__(self, trainDf, testDf, normalCount, keepFeature):
        trainX, trainY, testX, testY = utils.SpliteTrainTest(trainDf, testDf, normalCount, keepFeature)
        self.trainX = trainX
        self.trainY = trainY
        self.testX = testX
        self.testY = testY
        self.featureSelected = {}

    def TrainModel(self, models, cv = 5, randomSeed = 42):
        X, y = SVMSMOTE(random_state=randomSeed).fit_resample(self.trainX, self.trainY) 
        results = {}

        for modelName, model in models.items():
            rfecv = RFECV(estimator=model, min_features_to_select=2, step=1, cv=cv, scoring='f1', n_jobs=-1)
            rfecv.fit(X, y)

            # test the model
            print(f"Model: {modelName}")
            utils.TestModelPerformance(rfecv, X, y)
            accuracy, f1, recall, specificity = utils.TestModelPerformance(rfecv, self.testX, self.testY)
            utils.DrawRFECVPlot(rfecv)

            results[modelName] = {
                "accuracy": accuracy,
                "f1": f1,
                "recall": recall,
                "specificity": specificity,
            }

        features = X.columns
        features = [features[i] for i in range(len(features)) if rfecv.support_[i]]
        self.featureSelected = {modelName: features}

        resultsDf = pd.DataFrame(results).T
        return resultsDf

    def CrossFeatures(self):
        if self.featureSelected == {}:
            raise Exception("Please train the model first")
        featureSets = [set(features) for features in self.featureSelected.values()]
        commonFeatures = set.intersection(*featureSets) if featureSets else set()

        return commonFeatures

    def GetSelectedFeatures(self):
        return self.featureSelected

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
        "XGBoost": XGBClassifier(random_state=42),
        "RandomForest": RandomForestClassifier(random_state=42)
    }
    results = modelTrainer.TrainModel(models)
    print(results)

    crossFeatures = modelTrainer.CrossFeatures()
    print(crossFeatures)

