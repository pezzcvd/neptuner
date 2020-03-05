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
    # Extend intervals to 150bp
    nucleosomes = padding(nucleosomes[0], nucleosomes[1], nucleosomes[2])
    # Set condition
    olap_r = nucleosomes[0][1:] > nucleosomes[1][:-1]
    olap_l = np.append(olap_r[1:], True)
    condition = np.where(olap_r & olap_l)[0] + 1
    wd_nucs = well_defined(nucleosomes[0], nucleosomes[1], condition)
    fz_nucs = fuzzy(nucleosomes[0], nucleosomes[1], condition)
    # Combining
    nucleosomes = nucleosome_coordinates(wd_nucs, fz_nucs)

    # Write output file
    pd.DataFrame({'start': nucleosomes[0], 'end': nucleosomes[1], "length": nucleosomes[2]}).to_csv(opt)
    return


###-INTERNAL FUNCTIONS-###


def coordinates(pr):
    # Retrieve the signal above the threshold coordinates
    pos = np.where(pr > 0)[0]
    idx = np.append(0, np.where(np.diff(pos) > 1)[0])
    st = pos[idx[:-1] + 1]
    en = pos[idx[1:]]
    ln = en - st
    return np.array([st, en, ln])


def padding(s, e, le):
    idx = np.where(le < 150)[0]
    rest = np.zeros(le.size)
    rest[idx] = np.ceil((150 - le[idx])/2)
    s = s - rest
    e = e + rest
    return np.array(s, e, le)



def well_defined(s, e, cnd):
    # Identify well-defined peaks coordinates

    wdst = s[cnd]
    wden = e[cnd]
    wdln = wden - wdst
    return np.array([wdst, wden, wdln])


def fuzzy(s, e, cnd):
    # Identify fuzzy peaks coordinates

    fzst = s[cnd[np.where(np.diff(cnd) - 1 > 1)[0]] + 1]
    fzen = e[cnd[np.where(np.diff(cnd) - 1 > 1)[0] + 1] - 1]  # to fix one number in append
    fzln = fzen - fzst
    return np.array([fzst, fzen, fzln])


def nucleosome_coordinates(wd, fz):
    # List all nucleosomes coordinates

    nucst = np.append(wd[0], fz[0])
    nucst = np.sort(nucst)

    nucen = np.append(wd[1], fz[1])
    nucen = np.sort(nucen)
    nucln = nucen - nucst
    return np.array([nucst, nucen, nucln])


###-END FUNCTIONS-###


#input = '/home/pejo/Scrivania/chr1.bedgraph'
#output = "/home/pejo/Scrivania/prova.csv"
# User parameters
input = sys.argv[1]
output = sys.argv[2]
thresh = sys.argv[3]
# Call main function
thresh_nucleosomes(input, output, thresh)
