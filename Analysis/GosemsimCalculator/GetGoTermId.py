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
    _geneList = _aucDf[_aucDf['DNAm'] == 'hypo']['gene'].tolist()

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
    df = pd.DataFrame(data, columns=header)
    FileSaver.SaveDataframe(df, _config["Paths"]["HYPO_GO_TERM_PATH"])

    # with open(_config["Paths"]["HYPO_GO_TERM_PATH"], 'w', newline='') as csvfile:
    #     csvwriter = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    #     for line in response.iter_lines():
    #         line = line.decode('utf-8')
    #         line_parts = line.split("\t")
    #         if all(line_part.strip() for line_part in line_parts):
    #             csvwriter.writerow(line_parts)