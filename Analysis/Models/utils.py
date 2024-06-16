import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import accuracy_score, f1_score, recall_score

def TransformTrainData(df: pd.DataFrame, normalNum: int):
    target = [0] * normalNum + [1] * (len(df.columns) - normalNum)
    df.loc['cancer'] = target
    df = df.T
    return df

def SpliteTrainTest(trainDf: pd.DataFrame, testDf: pd.DataFrame, normalCount: int, keepFeature: list):
    if "cancer" not in keepFeature:
        keepFeature.append("cancer")

    trainDf = TransformTrainData(trainDf, normalCount)
    trainDf = trainDf[trainDf.columns.intersection(keepFeature)]
    trainDf = trainDf.iloc[1:]

    testDf = TransformTrainData(testDf, 25)
    testDf = testDf[testDf.columns.intersection(keepFeature)]
    testDf = testDf.iloc[1:]
    testDf = testDf[testDf.columns]

    trainX = trainDf.drop(columns=["cancer"])
    trainY = trainDf["cancer"]
    testX = testDf.drop(columns=["cancer"])
    testY = testDf["cancer"]

    return trainX, trainY, testX, testY

def TestModelPerformance(model, testX, testY):
    testPredicted = model.predict(testX)
    accuracy = accuracy_score(testY, testPredicted)
    f1 = f1_score(testY, testPredicted, average = "binary")
    recall = recall_score(testY, testPredicted, average = "binary")
    specificity = recall_score(testY, testPredicted, average = "binary", pos_label=0)

    print("Accuracy: ", accuracy)
    print("F1: ", f1)
    print("Recall: ", recall)
    print("Specificity: ", specificity)
    print()

    return accuracy, f1, recall, specificity

def CorrelationHeatmap(train):
    correlations = train.corr()

    fig, ax = plt.subplots(figsize=(10,10))
    sns.heatmap(correlations, vmax=1.0, center=0, fmt='.2f', cmap="rainbow",
                square=True, linewidths=.5, annot=True, cbar_kws={"shrink": .70}
                )
    plt.show()
