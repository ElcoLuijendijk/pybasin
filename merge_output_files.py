import os
import pandas as pd

__author__ = 'elco'

folder = '/home/elco/python_scripts/pybasin/model_output/MB/final_results_23mar2016_2stage_cooling'

files = os.listdir(folder)

csv_files = [os.path.join(folder, file) for file in files if file[-4:] == '.csv']

df = pd.read_csv(csv_files[0])
dfo = df
for csv_file in csv_files[1:]:
    dfn = pd.read_csv(csv_file)

    df = pd.concat([df, dfn])

df = df[dfo.columns]
df = df.dropna(subset=['well'])

fn = os.path.join(folder, 'model_results_merged.csv')
df.to_csv(fn, index=False)

print 'done'
