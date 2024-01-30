from configparser import ConfigParser
import pandas as pd
import csv
from biomart import BiomartServer

if __name__ == '__main__':
    _configPath = "Analysis/GosemsimCalculator/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _aucDf = pd.read_csv(_config["Paths"]["AUC_GROUP_DATA_PATH"])
    _geneList = _aucDf[_aucDf['DNAm'] == 'hyper']['gene'].tolist()

    server = BiomartServer("http://www.ensembl.org/biomart")
    interpro = server.datasets['hsapiens_gene_ensembl']
    response = interpro.search({
        'filters': {
            'wikigene_name': _geneList
        },
        'attributes': [
            "external_gene_name",
            "name_1006"
        ]
    })