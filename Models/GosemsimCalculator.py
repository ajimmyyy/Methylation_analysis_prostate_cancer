import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import fcluster
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn_extra.cluster import KMedoids
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from collections import Counter

class GosemsimCalculator:
    # ClusterHierarchy(self, gosemsimDf, distanceMethod = 'average', clusterNumRange = range(5, 51))
    # 對距離資料分群
    # Parameters:
    # gosemsimDf: DataFrame，擁有兩兩相似度的資料
    # distanceMethod: String，距離算法，詳見https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
    # clusterNumRange: Rrange，分群數量範圍
    # Return:
    # DataFrame
    def ClusterHierarchy(self, gosemsimDf, distanceMethod = 'average', clusterNumRange = range(5, 51)):
        scaler = StandardScaler()
        scaledData = scaler.fit_transform(gosemsimDf)
        rowMatrix = hierarchy.linkage(scaledData, method = distanceMethod)

        silhouetteScores = []
        for numClusters in clusterNumRange:
            clusters = hierarchy.fcluster(rowMatrix, t = numClusters, criterion = 'maxclust')
            silhouetteAvg = silhouette_score(scaledData, clusters)    
            silhouetteScores.append(silhouetteAvg)

        bestClusterNum = clusterNumRange[silhouetteScores.index(max(silhouetteScores))]

        clusters = hierarchy.fcluster(rowMatrix, t = bestClusterNum, criterion = 'maxclust')
        data = {
            "gene": gosemsimDf.columns.tolist(),
            "cluster": clusters
        }
        df = pd.DataFrame(data)

        silhouetteFig = plt.figure(figsize=(20, 40))
        plt.plot(clusterNumRange, silhouetteScores, marker='o')
        plt.title('Silhouette Score vs Number of Clusters')
        plt.xlabel('Number of Clusters')
        plt.ylabel('Silhouette Score')
        plt.annotate(f'Best Cluster Number: {bestClusterNum}', 
            xy=(bestClusterNum, max(silhouetteScores)),
            xytext=(10, -30),
            textcoords='offset points',
            arrowprops=dict(facecolor='black', arrowstyle='wedge,tail_width=0.7', lw=1.5),
            fontsize=12,
            color='black',
            weight='bold')

        return df, silhouetteFig
    
    # ClusterNearestNeighbors(self, gosemsimDf)
    # 對距離資料分群
    # Parameters:
    # gosemsimDf: DataFrame，擁有兩兩相似度的資料
    # Return:
    def ClusterNeighborJioning(self, gosemsimDf, clusterNum = 5):
        distanceMatrix = 1 - gosemsimDf.abs()   
        matrixArray = distanceMatrix.to_numpy().tolist()
        labels = distanceMatrix.columns.tolist()

        n = len(matrixArray)
        distanceArray = [[0] * (i + 1) for i in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                distanceArray[j][i] = matrixArray[i][j]

        distMatrix = DistanceMatrix(labels, distanceArray)

        constructor = DistanceTreeConstructor()
        NJTree = constructor.nj(distMatrix)

        clusters = self.__TreeFCluster(NJTree, clusterNum)

        fig = plt.figure(figsize=(20, 20), dpi=100)
        ax = fig.add_subplot(1, 1, 1)
        Phylo.draw(NJTree, axes=ax, do_show=False)

        return clusters, fig

    # ClusterKMedoids(self, gosemsimDf, distanceMethod = 'average', clusterNumRange = range(5, 51))
    # 對距離資料分群
    # Parameters:
    # gosemsimDf: DataFrame，擁有兩兩相似度的資料
    # distanceMethod: String，距離算法
    # clusterNumRange: Rrange，分群數量範圍
    # Return:
    # DataFrame
    def ClusterKMedoids(self, gosemsimDf, distanceMethod = 'euclidean', clusterNumRange = range(5, 51)):
        distanceMatrix = (1 - gosemsimDf.abs()).to_numpy()

        silhouetteScores = []
        for numClusters in clusterNumRange:
            kmedoids = KMedoids(n_clusters=numClusters, metric=distanceMethod, method="pam").fit(distanceMatrix)
            labels = kmedoids.labels_
            silhouetteAvg = silhouette_score(distanceMatrix, labels)
            silhouetteScores.append(silhouetteAvg)

        bestClusterNum = clusterNumRange[silhouetteScores.index(max(silhouetteScores))]
        kmedoids = KMedoids(n_clusters=bestClusterNum, metric=distanceMethod, method="pam").fit(distanceMatrix)
        data = {
            "gene": gosemsimDf.columns.tolist(),
            "cluster": kmedoids.labels_
        }
        df = pd.DataFrame(data)

        silhouetteFig = plt.figure(figsize=(20, 40))
        plt.plot(clusterNumRange, silhouetteScores, marker='o')
        plt.title('Silhouette Score vs Number of Clusters')
        plt.xlabel('Number of Clusters')
        plt.ylabel('Silhouette Score')
        plt.annotate(f'Best Cluster Number: {bestClusterNum}', 
             xy=(bestClusterNum, max(silhouetteScores)),
             xytext=(10, -30),
             textcoords='offset points',
             arrowprops=dict(facecolor='black', arrowstyle='wedge,tail_width=0.7', lw=1.5),
             fontsize=12,
             color='black',
             weight='bold')

        return df, silhouetteFig

    # DrawHierarchy(self, gosemsimDf, distanceMethod = 'average')
    # 畫出熱圖
    # Parameters:
    # gosemsimDf: DataFrame，擁有兩兩距離的資料
    # distanceMethod: String，距離算法，詳見https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
    # Return:
    # Figure，熱圖
    def DrawHierarchy(self, gosemsimDf, distanceMethod = 'average', num_clusters=None):
        fig = plt.figure(figsize=(20, 20))
        row_linkage = hierarchy.linkage(gosemsimDf, method = distanceMethod)

        if num_clusters is not None:
            clusters = fcluster(row_linkage, num_clusters, criterion='maxclust')
            max_d = max(row_linkage[:, 2])
            color_threshold = row_linkage[-(num_clusters-1), 2]
        else:
            color_threshold = None
            
        dn = hierarchy.dendrogram(row_linkage, labels = gosemsimDf.columns.tolist(), color_threshold=color_threshold)
        return fig
    
    def __TreeFCluster(self, tree, minCount):
        depths = tree.depths(True)
        min_depth = min(depth for depth, count in Counter(depths.values()).items() if count >= minCount)
        groups = [key for key, value in depths.items() if value == min_depth]

        clusters = {
            "gene": [],
            "cluster": []
        }
        
        for index, clade in enumerate(groups):
            leaves_below_clade = self.__GetLeavesBelowNode(clade)
            clusters["gene"].extend(leaves_below_clade)
            clusters["cluster"].extend([index] * len(leaves_below_clade))

        return pd.DataFrame(clusters)
    
    def __GetLeavesBelowNode(self, node):
        leaves = []
        if node.is_terminal():
            leaves.append(node.name)
        else:
            for child in node.clades:
                leaves.extend(self.__GetLeavesBelowNode(child))
        return leaves