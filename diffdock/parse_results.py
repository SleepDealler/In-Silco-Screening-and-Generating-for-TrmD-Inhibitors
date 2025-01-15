from IPython.display import display
from glob import glob
import pandas as pd
import re

results = []
for folder in glob('CHEMBL*'):
    for file in glob(f'./{folder}/rank1_*'):
        file = str(file)
        results.append((file[2:-9].split('/')[0], file[-9:-4]))

df = pd.DataFrame(results, columns=['id', 'confidence'])
df.confidence = df.confidence.apply(lambda text: re.findall(r"[-+]?\d*\.?\d+|\d+", text)[0])
df.confidence = df.confidence.apply(float)

df = df.sort_values(by='confidence', ascending=False)

df.to_csv('results.csv', index=False)