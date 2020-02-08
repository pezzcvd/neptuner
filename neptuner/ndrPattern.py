import sys, os
import pandas as pd
import numpy as np
from peakdetect import peakdetect
from plusminus1 import plusminus1

annotation = '/home/pejo/Scrivania/nept_pattern/annotation.csv'
profile = '/home/pejo/Scrivania/nept_pattern/profiles/profile1.csv'
nucleosome = "/home/pejo/neptuner_prj/neptuner/neptuner/nucl1.csv"
chr_name = "ChrI_A_nidulans_FGSC_A4"
#output


annot = pd.read_csv(annotation)
annot = annot.loc[(annot.seqnames) == chr_name,]
annot = annot.sort_values(by=["start", "end"])

prof = pd.read_csv(profile, header=None)
nucl = pd.read_csv(nucleosome)

peaks = peakdetect(prof, lookahead=50)
maxPos = np.array([peaks[0][c][0] for c in range(len(peaks[0]))])
minPos = np.array([peaks[1][c][0] for c in range(len(peaks[1]))])

prom = annot

prom.loc[(prom.strand == "+"), "end"] = prom.loc[(prom.strand == "+"), "start"] + 100
prom.loc[(prom.strand == "+"), "start"] = prom.loc[(prom.strand == "+"), "start"] - 500

prom.loc[(prom.strand == "-"), "start"] = prom.loc[(prom.strand == "-"), "end"] - 100
prom.loc[(prom.strand == "-"), "end"] = prom.loc[(prom.strand == "-"), "end"] + 500
prom["width"] = prom["end"] - prom["start"]

nucInt = np.array([pd.Interval(a,b) for a, b in zip(nucl["start"], nucl["end"])])
promInt = np.array([pd.Interval(a,b) for a, b in zip(prom["start"], prom["end"])])

#getpm1 = lambda x: plusminus1(x, prom, nucl, promInt, nucInt, maxPos)
res = [plusminus1(i, prom, nucl, promInt, nucInt, maxPos) for i in range(1,10)] #range(1,prom.shape[0])]
pd.DataFrame({'name': np.array([prom.iloc[i]["Name"] for i in range(1,10)]), #range(1,prom.shape[0])]),
              'plus1st': np.array([res[c][0] for c in range(len(res))]),
              'plus1en': np.array([res[c][1] for c in range(len(res))]),
              'minus1st': np.array([res[c][2] for c in range(len(res))]),
              'minus1en': np.array([res[c][3] for c in range(len(res))]),
              "ndr": np.maximum(np.array([res[c][0] for c in range(len(res))]) -
                            np.array([res[c][3] for c in range(len(res))]),
                            np.array([res[c][2] for c in range(len(res))]) -
                            np.array([res[c][1] for c in range(len(res))]))})

#for p in promInt: per ogni promotore
#for n in nucInt:
 #   p.overlaps(n)

# strand +
olaps = np.where(np.array([promInt[1].overlaps(n) for n in nucInt]))[0]
nnuc = nucl.iloc[olaps].shape[0] - 1
candidates = nucl.iloc[olaps[nnuc] - 2:olaps[nnuc] + 2]
#confronto tra distanze di picchi (corrente e precedente)
# per stabilire nucleosoma+1
idx = np.where(np.max(np.diff(np.intersect1d(maxPos[maxPos < np.max(candidates["end"])],
                                             maxPos[maxPos > np.min(candidates["start"])])))) #maxend and minstart
minus1 = candidates.iloc[idx[0][0]]
plus1 = candidates.iloc[idx[0][0] + 1]

#strand -
olaps = np.where(np.array([promInt[4].overlaps(n) for n in nucInt]))[0]
candidates = nucl.iloc[olaps[0] - 1:olaps[0] + 3]
idx = np.where(np.max(np.diff(np.intersect1d(maxPos[maxPos < np.max(candidates["end"])],
                                             maxPos[maxPos > np.min(candidates["start"])])))) #maxend and minstart
plus1 = candidates.iloc[idx[0][0]]
minus1 = candidates.iloc[idx[0][0] + 1]


np.max(np.array([res[c][0] for c in range(len(res))]) - np.array([res[c][3] for c in range(len(res))]),
                            np.array([res[c][2] for c in range(len(res))]) -
                            np.array([res[c][1] for c in range(len(res))]))