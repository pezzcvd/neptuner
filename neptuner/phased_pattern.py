import sys
import numpy as np
import pandas as pd


###-MAIN FUNCTION-###


def phased_pattern(ipt, opt):
    ## Identify phased patterns, arrays of nucleosomes that lie at regular distance from each other

    # Read the input file
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

    coord = extrema(occ)
    rrind = phased_indexes(coord)
    phased_r = populate(coord, rrind)

    # Write output file
    pd.DataFrame({'start': phased_r[0], 'end': phased_r[1], "length": phased_r[1] - phased_r[0]}).to_csv(opt, index=False)
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


def extrema(pr):
    # Finds maxima and minima and orders them
    peaks = peakdetect(pr, lookahead=50)
    indexesmax = np.array([peaks[0][el][0] for el in range(len(peaks[0]))])
    indexesmin = np.array([peaks[1][el][0] for el in range(len(peaks[1]))])
    extrema = np.sort(np.concatenate((indexesmax, indexesmin)))
    return extrema


def phased_indexes(cd):
    # retrieves distance between peaks
    rr = cd[1:] - cd[:-1]
    # sets the width to 0 if it is smaller than 70bp or bigger than 90bp
    rr[(rr < 70) | rr > 90] = 0
    # keeps only those elements with a width different from 0
    rridx = np.where(rr != 0)[0]
    return rridx


def populate(cd, rri):
    # Sets first element coordinates of phased regions
    st = np.array(cd[rri[0]])
    en = np.array(cd[rri[0] + 1])
    # Populates the two vectors
    # Scans the indexes of regular regions,
    # identifies clusters of consecutive indexes
    # and takes the correspondent genomic coordinates as coordinates of phased regions
    for el in zip(rri[:-1], rri[1:]):
        if el[1] - el[0] == 1:
            en[-1] = cd[el[1] + 1]
        else:
            st = np.append(st, cd[el[1]])
            en = np.append(en, cd[el[1] + 1])
    return np.array([st, en])


###-END FUNCTIONS-###


# User parameters
input = sys.argv[1]
output = sys.argv[2]
# Call main function
phased_pattern(input, output)
