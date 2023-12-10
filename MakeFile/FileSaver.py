import pandas as pd
import matplotlib.pyplot as plt
import pickle

class FileSaver:
    @staticmethod
    def SaveDataframe(dataframe, filename):
        dataframe.to_csv(filename, index=False)

    @staticmethod
    def SavePlot(plot, filename):
        plot.savefig(filename)

    @staticmethod
    def SaveList(dataList, filename):
        with open(filename, 'wb') as file:
            pickle.dump(dataList, file)