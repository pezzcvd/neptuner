import multiprocessing as mp
import pandas as pd
import numpy as np
from scipy.stats import poisson

import sys, os, math


def readstack(o, ss, ee):
    o[ss:ee+1] += 1


def getpvals(st, en):
    lng = max(en)
    old = np.array([0] * lng)
    [readstack(old, sta, end) for sta, end in np.nditer([st, en])][-1]

    expect = poisson.rvs(np.mean(old), size=10000)
    pvals = poisson.pmf(old, np.mean(old))

    pd.DataFrame({'pval': pvals}).to_csv(opt+"_pvals_example.csv",
                                         index=False, header=False)
    return

getpvals(bf["start"], bf["end"])


bedfile = pd.read_csv("/home/pejo/neptuner_prj/validation/neurospora/mn1.bed", sep='\t', header=None)
st = bedfile[1]
en = bedfile[2]

#st = st[:1000]
#en = en[:1000]
lng = max(en)
old = np.array([0]*lng)
[readstack(old, sta, end) for sta, end in np.nditer([st, en])][-1]

expect = poisson.rvs(np.mean(old), size=10000)
pvals = poisson.pmf(old, np.mean(old))

pd.DataFrame({'pval': pvals}).to_csv("/home/pejo/neptuner_prj/validation/neurospora/mnase1_pvals_example.csv", index=False, header=False)




#####
for i in range(10000):
    s = bedfile.iloc[i]["start"]
    e = bedfile.iloc[i]["end"]
    old[s:e] += 1

fin = old[:10000]




import numpy as np
expect = np.random.poisson(np.mean(fin), 10000)
