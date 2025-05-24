import pandas as pd

df = pd.read_pickle("data_with_msa.pkl")
print(df.iloc[0,2])
print("==========================")
print(df.iloc[1,2])