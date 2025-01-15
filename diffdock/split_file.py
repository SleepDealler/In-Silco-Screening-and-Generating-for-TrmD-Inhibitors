import pandas as pd

df = pd.read_csv('all_samples.csv')

for guy in range(8):
    if guy < 7:
        df2 = df.iloc[guy * 1000: (guy + 1) * 1000]
    else:
        df2 = df.iloc[guy * 1000:]
    df2.to_csv(f'./splittt/siema{guy}.csv')
