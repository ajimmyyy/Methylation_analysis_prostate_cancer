import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial import distance
from scipy.cluster import hierarchy
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

class GosemsimCalculator:
    def CalculateDistanceMatrixGeometricMean():
        return
    # ClusterGoterm(self, gosemsimDf, distanceMethod = 'average', clusterNumRange = range(5, 51))
    # 對距離資料分群
    # Parameters:
    # gosemsimDf: DataFrame，擁有兩兩距離的資料
    # distanceMethod: String，距離算法，詳見https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
    # clusterNumRange: Rrange，分群數量範圍
    # Return:
    # DataFrame
    def ClusterGoterm(self, gosemsimDf, distanceMethod = 'average', clusterNumRange = range(5, 51)):
        scaler = StandardScaler()
        scaledData = scaler.fit_transform(gosemsimDf)
        rowMatrix = hierarchy.linkage(scaledData, method = distanceMethod)

        scores = self.__selectClusterNum(scaledData, rowMatrix, clusterNumRange)

        bestClusterNum = clusterNumRange[scores.index(max(scores))]
        bestSilhouetteScore = max(scores)

        clusters = hierarchy.fcluster(rowMatrix, t = bestClusterNum, criterion = 'maxclust')
        data = {
            "gene": gosemsimDf.columns.tolist(),
            "cluster": clusters
        }
        df = pd.DataFrame(data)

        silhouetteFig = plt.figure(figsize=(20, 40))
        plt.plot(clusterNumRange, scores, marker='o')
        plt.title('Silhouette Score vs Number of Clusters')
        plt.xlabel('Number of Clusters')
        plt.ylabel('Silhouette Score')

        return df, silhouetteFig

    # DrawHeatmap(self, gosemsimDf, distanceMethod = 'average')
    # 畫出熱圖
    # Parameters:
    # gosemsimDf: DataFrame，擁有兩兩距離的資料
    # distanceMethod: String，距離算法，詳見https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
    # Return:
    # Figure，熱圖
    def DrawHeatmap(self, gosemsimDf, distanceMethod = 'average'):
        fig = plt.figure(figsize=(20, 20))
        row_linkage = hierarchy.linkage(gosemsimDf, method = distanceMethod)
        dn = hierarchy.dendrogram(row_linkage, labels = gosemsimDf.columns.tolist())
        return fig
    
    def __selectClusterNum(self, distanceNatrix, linkMatrix, clusterNumRange):
        silhouetteScores = []

        for numClusters in clusterNumRange:
            clusters = hierarchy.fcluster(linkMatrix, t = numClusters, criterion = 'maxclust')
            silhouetteAvg = silhouette_score(distanceNatrix, clusters)    
            silhouetteScores.append(silhouetteAvg)
        
        return silhouetteScores