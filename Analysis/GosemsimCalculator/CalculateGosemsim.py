from configparser import ConfigParser
import pandas as pd
from GOntoSim import GOntoSim

if __name__ == '__main__':
    _configPath = "Analysis/GosemsimCalculator/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _goTermDf = pd.read_csv(_config["Paths"]["HYPER_GO_TERM_PATH"])
    _goTermDf = _goTermDf.dropna()
    _goTermDf['GO term accession'] = _goTermDf['GO term accession'].apply(GOntoSim.goterm_replaced_by)

    _mfDf = _goTermDf[_goTermDf['GO domain'] == 'molecular_function']
    _goTermList = _mfDf.groupby('Gene name')['GO term accession'].apply(list).reset_index().values.tolist()

    _uniqueGoterms = set(_goTermDf.loc[:, 'GO term accession'].tolist())
    _uniqueGoterms = list(_uniqueGoterms)

    print(_uniqueGoterms)

    # Wang method
    S_values = [(x, GOntoSim.Semantic_Value(x, 'wang')) for x in _uniqueGoterms]
    S_values = dict(S_values)

    _similarityDf = GOntoSim.Similarity_Matrix(_goTermList, "wang", S_values)
    _similarityDf = pd.DataFrame(_similarityDf, index=_mfDf['Gene name'].unique(), columns=_mfDf['Gene name'].unique())

    # # GontoSim method
    # S_values = [(x, GOntoSim.Semantic_Value(x, 'Baseline_LCA_avg')) for x in _uniqueGoterms]
    # S_values = dict(S_values)

    # _similarityDf = GOntoSim.Similarity_Matrix(_goTermList, "GOntoSim", S_values)
    # _similarityDf = pd.DataFrame(_similarityDf, index=_mfDf['Gene name'].unique(), columns=_mfDf['Gene name'].unique())
    # _similarityDf.to_csv(_config["Paths"]["MF_GONTO_SIMILARITY_HYPER_PATH"])

