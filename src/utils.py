import pandas as pd
import matplotlib.pyplot as plt
import pickle

def SaveData(data, filename):
    if isinstance(data, pd.DataFrame):
        data.to_csv(filename, index=False)
    if isinstance(data, plt.Figure or plt.Axes):
        data.figure.savefig(filename)
    if isinstance(data, list):
        with open(filename, 'wb') as file:
            pickle.dump(data, file)