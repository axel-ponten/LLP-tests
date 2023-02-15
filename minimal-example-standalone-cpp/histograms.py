import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#df = pd.read_csv("LLP_dataframe_small_cs.csv")
df = pd.read_csv("LLP_dataframe_2.csv")


for key in df.keys():
    plt.figure()
    plt.hist(df[key], bins=20, label=str(len(df[key])))
    plt.title(key)
    plt.legend()

plt.show()
