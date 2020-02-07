import sys, os
import numpy as np
import pandas as pd
from peakdetect import peakdetect

# reads profile file
print("Checking input file")
path = os.path.dirname(sys.argv[1])

input = "/home/pejo/Scrivania/nept_pattern/profiles/profile1.csv"
output = "/home/pejo/phased.csv"

prof = pd.read_csv(input, header=None, )
occ = np.asarray(prof[0])
peaks = peakdetect(occ, lookahead=50)

indexesMax = np.array([peaks[0][el][0] for el in range(len(peaks[0]))])
indexesMin = np.array([peaks[1][el][0] for el in range(len(peaks[1]))])
coord = np.sort(np.concatenate((indexesMax, indexesMin)))

# retrieves distance between peaks
rr = coord[1:] - coord[:-1]

# sets the width to 0 if it is smaller than 70bp or bigger than 90bp
rr[(rr < 70) | rr > 90] = 0
#prova = [el if 70 <= el <= 90 else 0 for el in rr]

# keeps only those elements with a width different from 0
index = np.arange(len(rr))
rrind = np.where(rr != 0)[0]
    #[el for el in index if prova[el] != 0]
#prova = None

st = np.array(coord[rrind[0]])
en = np.array(coord[rrind[0] + 1])

for el in zip(rrind[:-1], rrind[1:]):
    if el[1] - el[0] == 1:
        en[-1] = coord[el[1] + 1]
    else:
        st = np.append(st, coord[el[1]])
        en = np.append(en, coord[el[1] + 1])

pd.DataFrame({'start': st, 'end': en, "length": en - st}).to_csv(output)

