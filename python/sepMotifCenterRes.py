# create subtables by motif name

import pandas as pd
import argparse
import os



parser = argparse.ArgumentParser(description = 'arguments')
parser.add_argument('--input', '-i', type = str, help = 'path of merged motif scanning result file', required = True)
parser.add_argument('--output', '-o', type = str, help = 'directory of motif subtables', required = True)
args = parser.parse_args()

try:
    centerPath = args.input
    outputDir = args.output
except Exception as e:
    print(e)

if not os.path.exists(outputDir):
    os.makedirs(outputDir, exist_ok=True)

mm9_motifCenter = pd.read_csv(centerPath, header = None)
for key, group in mm9_motifCenter.groupby(mm9_motifCenter.iloc[:, 1]):
    pre = "_".join(key.split('.')[:-1])
    group.to_csv(f'{outputDir}/{pre}.txt', columns=[0,2,3], index = False, header = None, sep = '\t')