import pandas as pd
from scipy.cluster import hierarchy
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

def select_cluster_num(distance_matrix, link_matrix, cluster_range):
    silhouette_scores = []

    for num_clusters in cluster_range:
        clusters = hierarchy.fcluster(link_matrix, t = num_clusters, criterion = 'maxclust')
        silhouette_avg = silhouette_score(distance_matrix, clusters)    
        silhouette_scores.append(silhouette_avg)
    
    return silhouette_scores

if __name__ == "__main__":
    fn_gosemsim = "Data/110820004/Data-gosemsim/geometric_mean/mgeneSim_geometric_mean_hyper.csv"
    # fn_gosemsim = "Data/110820004/Data-gosemsim/MF/mgeneSim_hyper.csv"
    # fn_gosemsim = "Data/110820004/Data-gosemsim/BP/mgeneSim_hyper.csv"
    # fn_gosemsim = "Data/110820004/Data-gosemsim/CC/mgeneSim_hyper.csv"
    fn_group = "Data/110820004/Data-gosemsim/geometric_mean/goterm_group_hyper.csv"
    fn_heatmap = "Data/110820004/Data-gosemsim/geometric_mean/goterm_heatmap_hyper.png"

    df_go = pd.read_csv(fn_gosemsim, index_col=0)

    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(df_go)

    row_matrix = hierarchy.linkage(scaled_data, method = 'complete')
    dendrogram = hierarchy.dendrogram(row_matrix, no_labels=True)

    range = range(5, 51)
    scores = select_cluster_num(scaled_data, row_matrix, range)

    best_num_clusters = range[scores.index(max(scores))]
    best_silhouette_score = max(scores)

    clusters = hierarchy.fcluster(row_matrix, t = best_num_clusters, criterion = 'maxclust')
    data = {
        "gene": df_go.columns.tolist(),
        "cluster": clusters
    }
    df = pd.DataFrame(data)
    df.to_csv(fn_group, sep=',', encoding='utf-8', index=False)
    plt.savefig(fn_heatmap)

    print(f"Best Number of Clusters: {best_num_clusters}")
    print(f"Best Silhouette Score: {best_silhouette_score}")

    plt.figure(figsize=(8, 6))
    plt.plot(range, scores, marker='o')
    plt.title('Silhouette Score vs Number of Clusters')
    plt.xlabel('Number of Clusters')
    plt.ylabel('Silhouette Score')
    plt.show()