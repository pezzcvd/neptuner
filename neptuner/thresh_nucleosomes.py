import sys, os
import pandas as pd
import numpy as np

###-MAIN FUNCTION-###
def thresh_nucleosomes(ipt, opt, thr):
    ## Identify nucleosome coordinates from profile, exploting a threshold value

    # Assertions on parameters
    # Input is a string and an existing file
    assert isinstance(ipt, str)
    assert os.path.isfile(ipt)
    # Output is a string
    assert isinstance(opt, str)
    # thresh is a float in the range [0, 1]
    assert isinstance(thr, float)
    assert 0 < thr < 1

    # reads profile file
    prof = pd.read_csv(ipt, header=None)
    assert prof.shape[1] == 4
    assert prof[0].unique().size == 1
    assert isinstance(prof[0][0], str)
    # Second column (start) is integer
    occ = np.array(prof[1])
    assert all([isinstance(occ[i], np.int64) for i in range(prof.shape[0])])
    # Third column (end) is integer
    occ = np.array(prof[2])
    assert all([isinstance(occ[i], np.int64) for i in range(prof.shape[0])])
    # Third column (end) is integer
    occ = np.array(prof[3])
    assert all([isinstance(occ[i], np.float64) for i in range(prof.shape[0])])

    # Thresholding
    occ[occ < thr] = 0
    # Retrieve coordinates
    nucleosomes = coordinates(occ)
    # Write output file
    pd.DataFrame({'start': nucleosomes[0], 'end': nucleosomes[1], "length": nucleosomes[2]}).to_csv(opt)
    return

###-INTERNAL FUNCTIONS-###

def coordinates(pr):
    # Retrieve the signal above the threshold coordinates
    pos = np.where(pr > 0)[0]
    idx = np.append(True, np.diff(pr) > 1)
    st = pos[idx] + 1
    en = pos[idx]
    ln = en - st
    return np.array([st, en, ln])

###-END FUNCTIONS-###

#input = '/home/pejo/Scrivania/chr1.bedgraph'
#output = "/home/pejo/Scrivania/prova.csv"
# User parameters
input = sys.argv[1]
output = sys.argv[2]
thresh = sys.argv[3]
# Call main function
thresh_nucleosomes(input, output, thresh)