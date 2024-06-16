import hydra
from omegaconf import DictConfig
import pandas as pd
from src.DeltaBetaFilter import DeltaBetaFilter
from src.utils import SaveData

@hydra.main(config_path="../configs", config_name="config")
def main(cfg: DictConfig):
    _betaDataDf = pd.read_csv(cfg.row_data.train.beta_data_path) 
    _betaDataDf.columns.values[0] = "CpG"

    _dbateFilter = DeltaBetaFilter()
    _dfOut = _dbateFilter.CalculateDeltaBeta(_betaDataDf, 25)

    SaveData(_dfOut, cfg.Paths.DBETA_DATA_PATH)

if __name__ == "__main__":
    _configPath = "Analysis/DeltaBetaFilter/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    _betaDataDf = pd.read_csv(_config["Paths"]["BETA_DATA_PATH"]) 
    _betaDataDf.columns.values[0] = "CpG"

    _dbateFilter = DeltaBetaFilter()
    _dfOut = _dbateFilter.CalculateDeltaBeta(_betaDataDf, 25)

    SaveData(_dfOut, _config["Paths"]["DBETA_DATA_PATH"])