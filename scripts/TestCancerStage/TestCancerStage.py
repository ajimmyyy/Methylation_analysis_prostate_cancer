from configparser import ConfigParser
import pandas as pd
import matplotlib.pyplot as plt
from Models.DeltaBetaFilter import DeltaBetaFilter
from MakeFile.FileSaver import FileSaver

def DrawDiffBeta(earlyDf,laterDf):
    earlyDf.set_index(earlyDf.columns[0], inplace=True)
    laterDf.set_index(laterDf.columns[0], inplace=True)
    earlyDf = earlyDf.T
    laterDf = laterDf.T
    earlyDf = earlyDf[75:]
    laterDf = laterDf[75:]

    fig, ax = plt.subplots()
    box_width = 0.35
    positions1 = range(len(earlyDf.columns))
    positions2 = [x + box_width for x in positions1]

    earlyDf.boxplot(ax=ax, positions=positions1, widths=box_width, showfliers=False, patch_artist=True, color='skyblue')
    laterDf.boxplot(ax=ax, positions=positions2, widths=box_width, showfliers=False, patch_artist=True, color='salmon')

    ax.set_xticks([p + box_width / 2 for p in positions1])
    ax.set_xticklabels(earlyDf.columns.tolist())
    ax.set_title('Early vs Late Stage Cancer')
    ax.legend()

    plt.show()

def CalculateDeltaBeta(earlyDf, laterDf):
    dbetaFilter = DeltaBetaFilter()
    earlyDbetaDf = dbetaFilter.CalculateDeltaBeta(earlyDf, 75)
    laterDbetaDf = dbetaFilter.CalculateDeltaBeta(laterDf, 75)

    mergedDf = pd.merge(earlyDbetaDf, laterDbetaDf, left_index=True, right_index=True, suffixes=('_early', '_later'))
    mergedDf["dbeta_difference"] = (mergedDf["dbeta_later"] - mergedDf["dbeta_early"]).abs()

    max_dbeta_diff = mergedDf["dbeta_difference"].max()
    min_dbeta_diff = mergedDf["dbeta_difference"].min()
    average_dbeta_diff = mergedDf["dbeta_difference"].mean()

    print("max dbeta diff:", max_dbeta_diff)
    print("mean dbeta diff:", average_dbeta_diff)
    print("min dbeta diff:", min_dbeta_diff)

    return mergedDf


if __name__ == "__main__":
    _configPath = "Analysis/TestCancerStage/config.ini"
    _config = ConfigParser()
    _config.read(_configPath)

    filterGene = ["GAPDH", "ACTB", "MYH14", "PRHOXNB"]
    _dmpDf = pd.read_csv(_config["Paths"]["DMP_DATA"], usecols=["CpG", "gene"])
    _earlyDf = pd.read_csv(_config["Paths"]["EARLY_STAGE_DATA"])
    _laterDf = pd.read_csv(_config["Paths"]["LATE_STAGE_DATA"])

    # df = CalculateDeltaBeta(_earlyDf, _laterDf)
    # FileSaver.SaveData(df, "C:/Users/acer/Desktop/test/test.csv")

    for gene in filterGene:
        filerCpG = _dmpDf[_dmpDf["gene"] == gene]["CpG"].tolist()
        earlyDf = _earlyDf[_earlyDf["CpG"].isin(filerCpG)]
        laterDf = _laterDf[_laterDf["CpG"].isin(filerCpG)]

        print("Gene:", gene)
        CalculateDeltaBeta(earlyDf, laterDf)
        DrawDiffBeta(earlyDf,laterDf)








