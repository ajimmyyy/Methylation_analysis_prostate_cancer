import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

# 
#測試單一位點
#  
fn_beta = "C:/Users/acer/Desktop/Data-origin/train/all_beta_normalized.csv"
normal_num = 50
test_cpg = "cg12161228"
cutpoint = 0.17

# 取得位點
data_beta = pd.read_csv(fn_beta)
row = data_beta[data_beta[data_beta.columns[0]] == test_cpg]

# 區分normal, tumor
transform_row = row.to_numpy()[0]
normal_beta = transform_row[1:normal_num + 1:2]
tumor_beta = transform_row[normal_num + 1::2]

# 混淆矩陣
TN = np.sum(normal_beta < cutpoint)
FP = np.sum(normal_beta > cutpoint)
FN = np.sum(tumor_beta < cutpoint)
TP = np.sum(tumor_beta > cutpoint)

Recall = TP/(TP+FN)
Precision = TP/(TP+FP)
Accuracy = (TP+TN)/(TP+FP+FN+TN)
F1 = 2 * Precision * Recall / (Precision + Recall)

# 輸出
hist_n, bin_edges = np.histogram(normal_beta, bins=np.arange(0, 1.01, 0.01))
hist_t, bin_edges = np.histogram(tumor_beta, bins=np.arange(0, 1.01, 0.01))
fig, axs = plt.subplots(2, 1, figsize=(8, 10))

axs[0].bar(bin_edges[:-1], hist_n, width=0.01, color = "blue", alpha=0.5)
axs[0].bar(bin_edges[:-1], hist_t, width=0.01, color = "red", alpha=0.5)
axs[0].set_title('Bar Chart of Data')
axs[0].set_xlabel('Value')
axs[0].set_ylabel('Count')

axs[1].boxplot([normal_beta, tumor_beta])
axs[1].set_xticks([1,2],['normal' ,'tumor'])
axs[1].set_title('Boxplot of Data')

plt.tight_layout()
plt.show()

print("Recall: ", Recall)
print("Precision: ", Precision)
print("Accuracy: ", Accuracy)
print("F1: ", F1)
print("Median dif: ", (np.median(tumor_beta) - np.median(normal_beta)) / np.median(tumor_beta))