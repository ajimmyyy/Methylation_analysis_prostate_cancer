from configparser import ConfigParser
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
from sklearn.linear_model import LogisticRegression
from xgboost import XGBClassifier
from sklearn.ensemble import RandomForestClassifier
from imblearn.over_sampling import RandomOverSampler, SMOTE, ADASYN, SVMSMOTE, BorderlineSMOTE
from sklearn.model_selection import GridSearchCV
from imblearn.pipeline import Pipeline as ImbPipeline
from sklearn.metrics import recall_score, make_scorer
import utils

ROUND = 5

class ModelTrainer:

    def __init__(self, trainDf, testDf, normalCount, keepFeature):
        trainX, trainY, testX, testY = utils.SpliteTrainTest(trainDf, testDf, normalCount, keepFeature)
        self.trainX = trainX
        self.trainY = trainY
        self.testX = testX
        self.testY = testY
        self.results = None

    def TrainModel(self, oversamplers, models, scoring):
        results = {metric: [] for metric in scoring.keys()}
        results['model'] = []
        results['oversampler'] = []

        for oversampler_name, oversampler in oversamplers.items():
            for model_name, model in models.items():
                pipeline = ImbPipeline([
                    ('sampling', oversampler),
                    ('classification', model)
                ])
                
                grid_search = GridSearchCV(pipeline, param_grid = {}, scoring=scoring, refit='f1', cv=5, n_jobs=-1)
                grid_search.fit(self.trainX, self.trainY)
                
                best_model = grid_search.best_estimator_
                accuracy, f1, recall, specificity = utils.TestModelPerformance(best_model, self.testX, self.testY)
                
                results['model'].append(model_name)
                results['oversampler'].append(oversampler_name)
                results['accuracy'].append(accuracy)
                results['f1'].append(f1)
                results['recall'].append(recall)
                results['specificity'].append(specificity)

        results = pd.DataFrame(results)
        self.results = results

        return results

def DrawPlot(df):
    df['model_oversampler'] = df['model'] + ' ' + df['oversampler']
    df_melted = df.melt(id_vars=['model_oversampler'], value_vars=['accuracy', 'f1', 'recall', 'specificity'], var_name='metric', value_name='value')
    
    plt.figure(figsize=(14, 8))
    sns.scatterplot(data=df_melted, x='model_oversampler', y='value', hue='metric', style='metric', s=100)
    plt.xticks(rotation=90)
    plt.xlabel('Model and Oversampler Combination')
    plt.ylabel('Value')
    plt.title('Model and Oversampler Performance Metrics')
    plt.legend(title='Metric')

    plt.tight_layout()
    plt.show()

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

    # 定義評估指標
    scoring = {
        'accuracy': 'accuracy',
        'f1': 'f1',
        'recall': 'recall',
        'specificity': make_scorer(lambda y_true, y_pred: recall_score(y_true, y_pred, pos_label=0))
    }

    results = []

    for random_seed in np.random.randint(0, 1000, ROUND):
        oversamplings = {
            'Random': RandomOverSampler(random_state=random_seed),
            'SMOTE': SMOTE(random_state=random_seed),
            'SVMSMOTE': SVMSMOTE(random_state=random_seed),
            'ADASYN': ADASYN(random_state=random_seed),
            'BLSMOTE': BorderlineSMOTE(random_state=random_seed),
        }

        models = {
            'LR': LogisticRegression(
                C=10, 
                penalty='l1', 
                solver='saga',
                n_jobs=-1, 
                random_state=random_seed
            ),
            'RF': RandomForestClassifier(
                n_estimators=150, 
                max_depth=10, 
                max_leaf_nodes=5, 
                max_features='sqrt', 
                n_jobs=-1,
                random_state=random_seed
            ),
            'XGB': XGBClassifier(
                n_estimators=150,
                subsample=0.7,
                max_depth=5,
                learning_rate=0.1,
                colsample_bytree=0.7,
                reg_alpha=1,
                n_jobs=-1,
                random_state=random_seed
            )
        }

        modelTrainer = ModelTrainer(trainDf, testDf, 25, keepFeature)
        result = modelTrainer.TrainModel(oversamplings, models, scoring)
        results.append(result)

    results = pd.concat(results)
    average_results = results.groupby(['oversampler', 'model']).mean(numeric_only=True).reset_index()
    DrawPlot(average_results)