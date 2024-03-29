import sys, os
import pandas as pd
import numpy as np


###-MAIN FUNCTION-###
# USAGE: pyhton3 nucleosomes.py input_noextension <output.txt>


def nucleosomes(ipt, opt):
    ## Identify nucleosome coordinates from profile, including well-defined and fuzzy nucleosomes.

    # Assertions on parameters
    # Input is a string and a bedgraph file
    assert isinstance(ipt, str)
    #assert os.path.isfile(ipt)
    #assert ipt.endswith(".bedgraph")
    # Output is a bedgraph file
    assert isinstance(opt, str)
    assert opt.endswith(".txt")

    # Read the input file
    prof = pd.read_csv(ipt+".bedgraph", header=None)
    assert prof.shape[1] == 4
    assert prof[0].unique().size == 1
    #assert isinstance(prof[0][0], str)
    # Second column (start) is integer
    occ = np.array(prof[1])
    assert all([isinstance(occ[i], np.int64) for i in range(prof.shape[0])])
    # Third column (end) is integer
    occ = np.array(prof[2])
    assert all([isinstance(occ[i], np.int64) for i in range(prof.shape[0])])
    # Third column (end) is integer
    occ = np.array(prof[3])
    assert all([isinstance(occ[i], np.float64) for i in range(prof.shape[0])])
    occ[occ < 0.05] = 0 #max(occ) - (max(occ) * 0.95)

    # finds indexes of peaks above the threshold (dummy, very low value)
    # and with at least 50bp distance ideal to remove little bumps
    peaks = peakdetect(occ, lookahead=50)

    # Get maximum positions
    indexesmax = np.array([peaks[0][el][0] for el in range(len(peaks[0]))])
    indexesmin = np.array([peaks[1][el][0] for el in range(len(peaks[1]))])

    # fixlengths
    if indexesmax[0] < indexesmin[0]:
        indexesmax = indexesmax[1:]
    if indexesmax[-1] > indexesmin[-1]:
        indexesmax = indexesmax[:-1]
#    # Pvals filtering

    pvals = pd.read_csv(ipt + "_pvals.csv", header=None)
    pvals = np.array(pvals[0])
    sig = np.where(pvals[indexesmax] < 1e-3)[0]
    lst = sig[-1]
    #sig1 = sig[:-1]
    indexesmax = indexesmax[sig]
    indminbef = indexesmin[sig]
    indminaft = np.append(indexesmin[sig][1:], indexesmin[lst + 1])


    # Get maximum values
    valsmax = np.array([peaks[0][el][1] for el in range(len(peaks[0]))])
    valsmax = valsmax[sig]
    # Get minimum values
    valsmin = np.array([peaks[1][el][1] for el in range(len(peaks[1]))])
    valsminbef = valsmin[sig]
    valsminaft = np.append(valsmin[sig][1:], valsmin[lst + 1])


    #if indexesmax[0] > indexesmin[0]:
        #valsmin = valsmin[1:]
    #if indexesmax[-1] < indexesmin[-1]:
        #valsmin = valsmin[:-1]

    # Establish nucleosome intervals [max - 75, max + 75]
    st = indexesmax - 75
    en = indexesmax + 75

    fcond = nucleosome_conditions(st, en, valsmax, valsminbef, valsminaft)
    well_defined_peaks = well_defined(st, en, fcond)
    fuzzy_peaks = fuzzy(st, en, fcond)
    rest_peaks = rest(st, en, fcond)
    nucleos = nucleosome_coordinates(well_defined_peaks, fuzzy_peaks, rest_peaks)

    # Writing output file
    pd.DataFrame({'start': nucleos[0], 'end': nucleos[1], "length": nucleos[2]}).to_csv(opt, index=False)
    return


###-INTERNAL FUNCTIONS-###


def peakdetect(y_axis, x_axis=None, lookahead=500, delta=0):
    """
    Converted from/based on a MATLAB script at http://billauer.co.il/peakdet.html

    Algorithm for detecting local maximas and minmias in a signal.
    Discovers peaks by searching for values which are surrounded by lower
    or larger values for maximas and minimas respectively

    keyword arguments:
    y_axis -- A list containg the signal over which to find peaks
    x_axis -- A x-axis whose values correspond to the 'y_axis' list and is used
        in the return to specify the postion of the peaks. If omitted the index
        of the y_axis is used. (default: None)
    lookahead -- (optional) distance to look ahead from a peak candidate to
        determine if it is the actual peak (default: 500)
        '(sample / period) / f' where '4 >= f >= 1.25' might be a good value
    delta -- (optional) this specifies a minimum difference between a peak and
        the following points, before a peak may be considered a peak. Useful
        to hinder the algorithm from picking up false peaks towards to end of
        the signal. To work well delta should be set to 'delta >= RMSnoise * 5'.
        (default: 0)
            Delta function causes a 20% decrease in speed, when omitted
            Correctly used it can double the speed of the algorithm

    return -- two lists [maxtab, mintab] containing the positive and negative
        peaks respectively. Each cell of the lists contains a tupple of:
        (position, peak_value)
        to get the average peak value do 'np.mean(maxtab, 0)[1]' on the results
    """
    maxtab = []
    mintab = []
    dump = []  # Used to pop the first hit which always if false

    length = len(y_axis)
    if x_axis is None:
        x_axis = range(length)

    # perform some checks
    if length != len(x_axis):
        raise ValueError("Input vectors y_axis and x_axis must have same length")
    if lookahead < 1:
        raise ValueError("Lookahead must be above '1' in value")
    if not (np.isscalar(delta) and delta >= 0):
        raise ValueError("delta must be a positive number")

    # needs to be a numpy array
    y_axis = np.asarray(y_axis)

    # maxima and minima candidates are temporarily stored in
    # mx and mn respectively
    mn, mx = np.Inf, -np.Inf

    # Only detect peak if there is 'lookahead' amount of points after it
    for index, (x, y) in enumerate(zip(x_axis[:-lookahead], y_axis[:-lookahead])):
        if y > mx:
            mx = y
            mxpos = x
        if y < mn:
            mn = y
            mnpos = x

        ####look for max####
        if y < mx - delta and mx != np.Inf:
            # Maxima peak candidate found
            # look ahead in signal to ensure that this is a peak and not jitter
            if y_axis[index:index + lookahead].max() < mx:
                maxtab.append((mxpos, mx))
                dump.append(True)
                # set algorithm to only find minima now
                mx = np.Inf
                mn = np.Inf

        ####look for min####
        if y > mn + delta and mn != -np.Inf:
            # Minima peak candidate found
            # look ahead in signal to ensure that this is a peak and not jitter
            if y_axis[index:index + lookahead].min() > mn:
                mintab.append((mnpos, mn))
                dump.append(False)
                # set algorithm to only find maxima now
                mn = -np.Inf
                mx = -np.Inf

    # Remove the false hit on the first value of the y_axis
    try:
        if dump[0]:
            maxtab.pop(0)
            # print "pop max"
        else:
            mintab.pop(0)
            # print "pop min"
        del dump
    except IndexError:
        # no peaks were found, should the function return empty lists?
        pass

    return maxtab, mintab


def nucleosome_conditions(s, e, ma, mib, mia):
    # Selects indexes of well-defined nucleosomes, setting some conditions:
    # nuclesomes must not overlap and peak height must be at least twice as high as the sorrounding minima

    # Setting conditions
    # Overlapping nucleosomes on the right
    olap = s[1:] > e[:-1]  # parte dal -primo
    # Overlapping nucleosomes on the left
    #olap_l = np.append(olap_r[1:], True)
    # Current maximum is higher than two times the minimumm on the left
    l_h_check = ma / 2 > mib  # parte dal primo minimo e secondo massimo
    # Current maximum is higher than two times the minimumm on the right
    r_h_check = ma / 2 > mia  # parte dal secondo
    # Applying conditions to select well-defined peaks
    return np.where(olap & l_h_check[:-1] & r_h_check[1:])[0] + 1


def well_defined(s, e, cnd):
    # Identify well-defined peaks coordinates

    wdst = s[cnd]
    wden = e[cnd]
    wdln = wden - wdst
    return np.array([wdst, wden, wdln])


def fuzzy(s, e, cnd):
    # Identify fuzzy peaks coordinates

    fzst = s[cnd[np.where(np.diff(cnd) - 1 > 1)[0]] + 1]
    fzen = e[cnd[np.where(np.diff(cnd) - 1 > 1)[0] + 1] - 1]
    fzln = fzen - fzst
    return np.array([fzst, fzen, fzln])


def rest(s, e, cnd):
    # Identify rest of peaks (special case)

    rest = cnd[np.where(np.diff(cnd) - 1 == 1)[0]] + 1
    restst = s[rest]
    resten = e[rest]
    restln = resten - restst
    return np.array([restst, resten, restln])


def nucleosome_coordinates(wd, fz, rs):
    # List all nucleosomes coordinates

    nucst = np.append(wd[0], fz[0])
    nucst = np.append(nucst, rs[0])
    nucst = np.sort(nucst)

    nucen = np.append(wd[1], fz[1])
    nucen = np.append(nucen, rs[1])
    nucen = np.sort(nucen)
    nucln = nucen - nucst
    return np.array([nucst, nucen, nucln])


###-END FUNCTIONS-###


# User parameters
input = sys.argv[1]
output = sys.argv[2]
# Call main function
nucleosomes(input, output)
