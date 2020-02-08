import numpy as np


def plusminus1(i, promoters, nucleosomes, pIntervals, nIntervals, maxCoord):
    olaps = np.where(np.array([pIntervals[i].overlaps(n) for n in nIntervals]))[0]
    if promoters.iloc[i]["strand"] == "+":
        nnuc = nucleosomes.iloc[olaps].shape[0] - 1
        candidates = nucleosomes.iloc[olaps[nnuc] - 2:olaps[nnuc] + 2]
        idx = np.where(np.max(np.diff(np.intersect1d(maxCoord[maxCoord < np.max(candidates["end"])],
                                                     maxCoord[maxCoord > np.min(
                                                         candidates["start"])]))))  # maxend and minstart
        minus1 = candidates.iloc[idx[0][0]]
        plus1 = candidates.iloc[idx[0][0] + 1]
    else:
        candidates = nucleosomes.iloc[olaps[0] - 1:olaps[0] + 3]
        idx = np.where(np.max(np.diff(np.intersect1d(maxCoord[maxCoord < np.max(candidates["end"])],
                                                     maxCoord[maxCoord > np.min(
                                                         candidates["start"])]))))  # maxend and minstart
        plus1 = candidates.iloc[idx[0][0]]
        minus1 = candidates.iloc[idx[0][0] + 1]
    return plus1["start"], plus1["end"], minus1["start"], minus1["end"]
