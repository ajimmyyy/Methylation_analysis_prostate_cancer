import pandas as pd

def tidy_GO_term(DataFrame):
    df = DataFrame.assign(Genes=DataFrame['Genes'].str.split(', ')).explode('Genes')
    df = df.groupby('Genes')['Term'].agg(list).reset_index()
    df.columns.values[0] = 'gene'
    return df

if __name__ == "__main__":
    fn_BP = 'Data/110820004/Data-GO_term/BP.txt'
    fn_CC = 'Data/110820004/Data-GO_term/CC.txt'
    fn_MF = 'Data/110820004/Data-GO_term/MF.txt'
    fn_hyperhypo = "Data/110820004/Data-ROC_AUC/comorbidity_group_auc.csv"
    fn_o = 'Data/110820004/Data-GO_term/comorbidity_group_go_term.csv'

    df_BP = pd.read_csv(fn_BP, sep='\t', lineterminator='\n', usecols=['Term', 'Genes'])
    df_BP['Term'] = df_BP['Term'].str[:10]
    df_CC = pd.read_csv(fn_CC, sep='\t', lineterminator='\n', usecols=['Term', 'Genes'])
    df_CC['Term'] = df_CC['Term'].str[:10]
    df_MF = pd.read_csv(fn_MF, sep='\t', lineterminator='\n', usecols=['Term', 'Genes'])
    df_MF['Term'] = df_MF['Term'].str[:10]

    df_hyerhypo = pd.read_csv(fn_hyperhypo, usecols=["gene", "DNAm"])


    df_BP = tidy_GO_term(df_BP)
    df_BP.columns.values[1] = 'BP'
    df_CC = tidy_GO_term(df_CC)
    df_CC.columns.values[1] = 'CC'
    df_MF = tidy_GO_term(df_MF)
    df_MF.columns.values[1] = 'MF'
    df_merge = df_BP.merge(df_CC,on='gene', how = "outer").merge(df_MF,on='gene', how = "outer")

    df_merge = df_merge.merge(df_hyerhypo, on = "gene", how = "outer")

    df_merge.to_csv(fn_o, sep=',', encoding='utf-8', index=False)


