import sys, os
import pandas as pd
import numpy as np
from peakdetect import peakdetect

input = '/home/pejo/Scrivania/nept_pattern/profiles/profile1.csv'
output = "/home/pejo/Scrivania/prova.csv"
thresh = sys.argv[3]

assert isinstance(input, str)
assert os.path.isfile(input)
assert isinstance(output, str)
assert isinstance(thresh, float)
assert 0 < thresh < 1

# reads profile file
prof = pd.read_csv(input, header=None)
assert len(prof.columns) == 1
assert prof[0].dtype == "float"


occ = np.asarray(prof[0])
occ[occ < thresh] = 0


# Retrieve the nucleosome coordinates
pos = np.where(occ > 0)[0]
idx = np.append(True, np.diff(pos) > 1)
st = pos[idx] + 1

pos = np.where(occ > 0)[0]
idx = np.append(np.diff(pos) > 1, True)
en = pos[idx] + 139
ln = en - st

peaks = peakdetect(occ, lookahead=40)
indexesMax = np.array([peaks[0][el][0] for el in range(len(peaks[0]))])
indexesMax = indexesMax[np.array([np.where(indexesMax > st[i])[0][0] for i in range(st.size-1)])]
peak = np.zeros(st.size-1)
peak[np.where( (indexesMax > st[:-1]) & (indexesMax < en[:-1]) & (ln[:-1] < 150))] = \
    indexesMax[(indexesMax > st[:-1]) & (indexesMax < en[:-1]) & (ln[:-1] < 150)]


pd.DataFrame({'start': st[:-1], 'end': en[:-1], "length": ln[:-1], "peak": peak}).to_csv(output)