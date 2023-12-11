import pandas as pd

if __name__ == "__main__":
    fn_group = "Data/110820004/Data-gosemsim/geometric_mean/goterm_group_hyper.csv"
    fn_i = "Data/110820004/Data-ROC_AUC/comorbidity_group_auc.csv"
    fn_max_dbeta = "Data/110820004/Data-gosemsim/result/goterm_group_max_dbate.csv"
    fn_max_auc = "Data/110820004/Data-gosemsim/result/goterm_group_max_auc.csv"
    fn_max_F1 = "Data/110820004/Data-gosemsim/result/goterm_group_max_F1.csv"

    df_group = pd.read_csv(fn_group)
    df = pd.read_csv(fn_i)

    merged_df = pd.merge(df, df_group, on='gene')
    idx_max_dbeta = merged_df.groupby('cluster')['dbeta'].idxmax()
    idx_max_auc = merged_df.groupby('cluster')['auc'].idxmax()
    idx_max_F1 = merged_df.groupby('cluster')['F1'].idxmax()

    merged_df.loc[idx_max_dbeta].to_csv(fn_max_dbeta, sep=',', encoding='utf-8', index=False)
    merged_df.loc[idx_max_auc].to_csv(fn_max_auc, sep=',', encoding='utf-8', index=False)
    merged_df.loc[idx_max_F1].to_csv(fn_max_F1, sep=',', encoding='utf-8', index=False)
