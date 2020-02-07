import sys, os  # , #peakutils
import numpy as np
import pandas as pd
from peakdetect import peakdetect
import time

# reads profile file
input = sys.argv[1]
output = sys.argv[2]

prof = pd.read_csv(input, header=None)

occ = np.asarray(prof[0])

# finds indexes of peaks above the threshold (dummy, very low value)
# and with at least 50bp distance ideal to remove little bumps
peaks = peakdetect(occ, lookahead=50)

indexesMax = np.array([peaks[0][el][0] for el in range(len(peaks[0]))])
indexesMin = np.array([peaks[1][el][0] for el in range(len(peaks[1]))])
valsMax = np.array([peaks[0][el][1] for el in range(len(peaks[0]))])
valsMin = np.array([peaks[1][el][1] for el in range(len(peaks[1]))])
st = indexesMax - 75
en = indexesMax + 75

#remove bumps

# evaluations
olap_r = st[1:] > en[:-1] # parte dal primo
olap_l = np.append(olap_r[1:], True)
#st_crsp = st[1:-1] - indexesMin[:-1]
#en_crsp = indexesMin - en[:-1]
l_h_check = valsMax[1:] / 2 > valsMin #parte dal primo minimo e secondo massimo
r_h_check = valsMax[:-1] / 2 > valsMin  #parte dal secondo


fcond = np.where(olap_r[:-1] & olap_l[:-1] & l_h_check[:-1] & r_h_check[1:])[0] + 1

wdst = st[fcond]
wden = en[fcond]
wdln = wden - wdst
#np.where(np.diff(fcond) > 1)[0]

fzst = np.append(1, st[fcond[np.where(np.diff(fcond) -1 > 1)[0]]+1])
fzen = np.append(en[fcond[np.where(np.diff(fcond) -1 > 1)[0]]-1][0],
                 en[fcond[np.where(np.diff(fcond) -1 > 1)[0]+1]-1])
fzln = fzen - fzst

nucst = np.sort(np.append(wdst, fzst))
nucen = np.sort(np.append(wden, fzen))
nucln = nucen - nucst

print("Writing output file")
pd.DataFrame({'start': nucst, 'end': nucen, "length": nucln}).to_csv(output)

