from configparser import ConfigParser
import pandas as pd
from MakeFile.FileSaver import FileSaver

def RemoveUnimportantGenes(df, pathwaysCol):
    sorted_df = df.iloc[df[pathwaysCol].apply(lambda x: len(x) if isinstance(x, list) else 0).argsort()]

    for idx, row in sorted_df.iterrows():
        covered_pathways = set()
        for other_idx, other_row in df.iterrows():
            if idx != other_idx and isinstance(other_row[pathwaysCol], list):
                covered_pathways.update(other_row[pathwaysCol])

        if isinstance(row[pathwaysCol], list) and set(row[pathwaysCol]).issubset(covered_pathways):
            df.drop(idx, inplace=True)

    return df

if __name__ == "__main__":
    _configPath = "Analysis/Pathway/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _pathwayDf = pd.read_csv(_config["Paths"]["PATHWAY_DATA_PATH"])
    _featureDf = pd.read_csv(_config["Paths"]["FEATURE_DATA_PATH"])

    _terms = {}

    for index, row in _pathwayDf.iterrows():
        genes = row['Genes'].split(', ')
        term = row['Term'].split(':')[0]

        for gene in genes:
            if gene in _terms:
                _terms[gene].append(term)
            else:
                _terms[gene] = [term]

    data = []
    for gene, terms in _terms.items():
        data.append({'gene': gene, 'pathways': terms})

    _df = pd.DataFrame(data)
    _df = pd.merge(_featureDf, _df, on='gene', how='left', suffixes=('', ''))
    _df = RemoveUnimportantGenes(_df, 'pathways')
    
    FileSaver.SaveData(RemoveUnimportantGenes(_df, 'pathways'), _config["Paths"]["PATHWAY_FILTER_DATA_PATH"])
