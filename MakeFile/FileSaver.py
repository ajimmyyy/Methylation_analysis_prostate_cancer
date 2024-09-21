import pandas as pd
import matplotlib.pyplot as plt
import pickle

class FileSaver:
    @staticmethod
    def SaveData(data, filename, index=False):
        if isinstance(data, pd.DataFrame):
            data.to_csv(filename, index=index)
        if isinstance(data, plt.Figure or plt.Axes):
            data.figure.savefig(filename)
        if isinstance(data, list):
            with open(filename, 'wb') as file:
                pickle.dump(data, file)