from configparser import ConfigParser
import pandas as pd
import networkx as nx
import functools
from pygosemsim import download, graph, similarity, term_set
from MakeFile.FileSaver import FileSaver

if __name__ == '__main__':
    _configPath = "Analysis/GosemsimCalculator/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _goTermDf = pd.read_csv(_config["Paths"]["HYPER_GO_TERM_PATH"])
    _goTermDf = _goTermDf.dropna()

    _grouped = _goTermDf.groupby(['Gene name', 'GO domain'])
    _goTermDict = {}

    for (gene, go_domain), group_df in _grouped:
        if gene not in _goTermDict:
            _goTermDict[gene] = {'BP': [], 'CC': [], 'MF': []}
        if go_domain == 'biological_process':
            _goTermDict[gene]['BP'].extend(group_df['GO term accession'].tolist())
        elif go_domain == 'cellular_component':
            _goTermDict[gene]['CC'].extend(group_df['GO term accession'].tolist())
        elif go_domain == 'molecular_function':
            _goTermDict[gene]['MF'].extend(group_df['GO term accession'].tolist())

    try:
        download.obo("go-basic")
        download.gaf("goa_human")
    except:
        pass

    G = graph.from_resource("go-basic")
    similarity.precalc_lower_bounds(G)
    sf = functools.partial(term_set.sim_func, G, similarity.wang)

    _similarityDf = pd.DataFrame(columns=_goTermDict.keys(), index=_goTermDict.keys())

    for gene1 in _goTermDict:
        for gene2 in _goTermDict:
            similarity = term_set.sim_bma(_goTermDict[gene1]["MF"], _goTermDict[gene2]["MF"], sf)
            _similarityDf.loc[gene1, gene2] = similarity

    _similarityDf.fillna(0, inplace=True)
    _similarityDf.to_csv(_config["Paths"]["MF_WANG_SIMILARITY_HYPER_PATH"])
