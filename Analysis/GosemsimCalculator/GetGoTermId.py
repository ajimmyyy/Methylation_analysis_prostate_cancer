from configparser import ConfigParser
import pandas as pd
import csv
from biomart import BiomartServer
from MakeFile.FileSaver import FileSaver

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
            "go_id",
            "namespace_1003"
        ]
    }, header = 1)

    data = []
    for entry in response.iter_lines():
        data.append(entry.decode('utf-8').split('\t'))
    header = data[0]
    data = data[1:]
    data = [sublist for sublist in data if all(item != '' for item in sublist)]
    df = pd.DataFrame(data, columns=header)
    FileSaver.SaveData(df, _config["Paths"]["HYPER_GO_TERM_PATH"])