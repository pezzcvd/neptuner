import os, sys
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression


###-MAIN FUNCTION-###
## USAGE: python3 ndr_pattern.py <annotation.csv> <profile.bedgraph> <nucleosomes.txt> chromosome_name <output.txt>

def ndr_pattern(annotation, profile, nucleosome, chr_name, outp):
    # Input controls
    input_controls(annotation, profile, nucleosome, chr_name, outp)

    # Read tables for the different parameters and setup
    # Read input files
    annot = pd.read_csv(annotation)
    # Select one chromosome and subselect annotation file
    annot = annot.loc[annot.seqnames == chr_name, ]
    # Sort genes
    annot = annot.sort_values(by=["start", "end"])
    # Add transcription starting site information
    annot["tss"] = np.sort(np.append(annot["start"][(annot['strand'] == "+")],
                                     annot["end"][(annot['strand'] == "-")]))
    # Filter out uORF genes
    names = annot["Name"].to_numpy()
    names = names.astype(str)
    namesid = np.where(np.char.endswith(names, "-uORF"))[0]
    annot = annot.reset_index(drop=True)
    annot = annot.drop(namesid, axis=0)
    # Reset row numbers
    annot = annot.reset_index(drop=True)
    prof = pd.read_csv(profile, header=None)
    prof = np.array(prof[3])
    nucl = pd.read_csv(nucleosome)

    # Operations
    # Peak calculation, it has to be outside the function!!!
    peaks = peakdetect(prof, lookahead=50)

    # Maxima positions
    maxima = extrema(peaks)

    # Fitting a linear model to genes profile extrema and retrieving slope
#    slope = np.array([lm_intercept(annot.iloc[r], peaks) for r in range(2, annot.shape[0])])
    # identifying promoter regions
    promoters = promoter(annot)
    # Set intervals for nucleosomes and promoters
    nucint = np.array([pd.Interval(a, b) for a, b in zip(nucl["start"], nucl["end"])])
    promint = np.array([pd.Interval(a, b) for a, b in zip(promoters["start"], promoters["end"])])
    # Calculate info for plus1 and minus peaks
    pm_peaks = [plusminus1(r, promoters, nucl, promint, nucint, maxima) for r in range(3, promoters.shape[0])]

    new_annot = extended_annotation(annot, pm_peaks, promoters)
    naidx = np.where(pd.isnull(new_annot["ndr_pattern"]))[0]
    new_annot.loc[naidx, "ndr_pattern"] = False

    # DIRE PATTERN SI/NO
    new_annot.to_csv(outp, index=False)
    return


###-INTERNAL FUNCTIONS-###


def input_controls(a, p, n, c, o):
    # Assertions on parameters
    # Annotation is a string and an existing .csv file
    assert isinstance(a, str)
    assert os.path.isfile(a)
    assert a.endswith(".csv")
    # Profile is a string and an existing file
    assert isinstance(p, str)
    assert os.path.isfile(p)
    # Nucleosome is a string and an existing file
    assert isinstance(n, str)
    assert os.path.isfile(n)
    # Chr_name is a string
    assert isinstance(c, str)
    # Output is a string
    assert isinstance(o, str)

    # Read annotation file
    anno = pd.read_csv(a)
    # Annotation file has at least six fields:
    # seqnames, start, end, width, strand, Name
    assert anno.shape[1] >= 6
    # First column (chromosome) is integer
    sq = np.array(anno["seqnames"])
    assert all([isinstance(sq[i], str) for i in range(anno.shape[0])])
    # Second column (start) is integer
    st = np.array(anno["start"])
    assert all([isinstance(st[i], np.int64) for i in range(anno.shape[0])])
    # Third column (end) is integer
    en = np.array(anno["end"])
    assert all([isinstance(en[i], np.int64) for i in range(anno.shape[0])])
    # Fourth column (length) is integer
    ln = np.array(anno["width"])
    assert all([isinstance(ln[i], np.int64) for i in range(anno.shape[0])])
    # Fifth column (strand) is integer
    std = np.array(anno["strand"])
    assert all([std[i] in ("+", "-") for i in range(anno.shape[0])])
    # Sixth column (gene_id) is integer
    gid = np.array(anno["Name"])
    assert all([isinstance(gid[i], str) for i in range(anno.shape[0])])

    # Read the profile file
    pro = pd.read_csv(p, header=None)
    # Profile file has at four fields
    assert pro.shape[1] == 4
    # First column (chromosome), only one chromosome present and is a string
    assert pro[0].unique().size == 1
    assert isinstance(pro[0][0], str)
    # Second column (start) is integer
    occ = np.array(pro[1])
    assert all([isinstance(occ[i], np.int64) for i in range(pro.shape[0])])
    # Third column (end) is integer
    occ = np.array(pro[2])
    assert all([isinstance(occ[i], np.int64) for i in range(pro.shape[0])])
    # Fourth column (score) is integer
    occ = np.array(pro[3])
    assert all([isinstance(occ[i], np.float64) for i in range(pro.shape[0])])

    # Further controls on nucleosomes file
    nuc = pd.read_csv(n)
    # Nucleosomes file has at three fields
    assert nuc.shape[1] == 3
    # First column (start) is integer
    st = np.array(nuc["start"])
    assert all([isinstance(st[i], np.int64) for i in range(nuc.shape[0])])
    # Second column (end) is integer
    en = np.array(nuc["end"])
    assert all([isinstance(en[i], np.int64) for i in range(nuc.shape[0])])
    # Third column (length) is integer
    ln = np.array(nuc["length"])
    assert all([isinstance(ln[i], np.int64) for i in range(nuc.shape[0])])

    # Further controls on chr_name
    assert c in anno["seqnames"].unique()
    return


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


def extrema(pe):
    # Idenfity extrema
    toolow = np.where(np.array([pe[0][c][1] < 0.1 for c in range(len(pe[0]))]))
    maxpos = np.array([pe[0][c][0] for c in range(len(pe[0]))])
    maxpos = np.delete(maxpos, toolow)
    return maxpos


def lm_intercept(an, pe):
    # Maxima calculation, positions and values
    maxp = np.array([pe[0][c][0] for c in range(len(pe[0]))])
    toolow = np.where(np.array([pe[0][c][1] < 0.1 for c in range(len(pe[0]))]))[0]
    idx1 = np.intersect1d(np.where(maxp > an[1])[0], np.where(maxp < an[2])[0])
    idx1 = np.delete(idx1, toolow)
    maxp = maxp[idx1]
    maxv = np.array([pe[0][c][1] for c in range(len(pe[0]))])
    maxv = maxv[idx1]

    # Minima calculation, positions and values
    minp = np.array([pe[1][c][0] for c in range(len(pe[1]))])
    toolow = np.where(np.array([pe[0][c][1] < 0.1 for c in range(len(pe[0]))]))
    idx2 = np.intersect1d(np.where(minp > an[1])[0], np.where(minp < an[2])[0])
    idx2 = np.delete(idx2, toolow)
    minp = minp[idx2]
    minv = np.array([pe[1][c][1] for c in range(len(pe[1]))])
    minv = minv[idx2]

    # Joining positions and values vectors to set up x and y for the linear regression
    x = np.append(maxp, minp).reshape((-1, 1))
    y = np.append(maxv, minv)

    try:
        # Model fitting
        model = LinearRegression().fit(x, y)
        # If the strand is negative we multiply the slope by -1
        if an[4] == "-":
            model.intercept_ = -model.intercept_

        return model.intercept_
    except ValueError:
        return np.nan


def promoter(a):
    prom = a.copy()

    prom.loc[(prom.strand == "+"), "end"] = prom.loc[(prom.strand == "+"), "start"] + 100
    prom.loc[(prom.strand == "+"), "start"] = prom.loc[(prom.strand == "+"), "start"] - 500

    prom.loc[(prom.strand == "-"), "start"] = prom.loc[(prom.strand == "-"), "end"] - 100
    prom.loc[(prom.strand == "-"), "end"] = prom.loc[(prom.strand == "-"), "end"] + 500
    prom["width"] = prom["end"] - prom["start"]
    return prom


def find_candidates(strand, n, olp, mc):
    if strand == "+":
        # Number of nucleosomes in the overlap
        nnuc = n.iloc[olp].shape[0] - 1
        # Candidates are last three nucleosomes in the overlap region plus the following one
        cand = n.iloc[olp[nnuc] - 2:olp[nnuc] + 2]
        # Takes position of max distance between peaks in the promoter region
        # Index of minus1 peaks
        idx = int(np.argmax(np.diff(np.intersect1d(mc[mc < np.max(cand["end"])],
                                                   mc[mc > np.min(cand["start"])]))))  # maxend and minstart
        # cddm = candidate max
        # Identify candidates in a window
        cddm = np.intersect1d(mc[mc < np.max(cand["end"])],
                              mc[mc > np.min(cand["start"])])[idx]
        cddm = pd.Interval(cddm, cddm)
        cddp = np.array([pd.Interval(a, b) for a, b in zip(cand["start"], cand["end"])])
    else:
        # Define candidates
        cand = n.iloc[olp[0] - 1:olp[0] + 3]
        if np.intersect1d(mc[mc < np.max(cand["end"])],
                          mc[mc > np.min(cand["start"])]).size == 1:
            cddm = np.intersect1d(mc[mc < np.max(cand["end"])],
                                  mc[mc > np.min(cand["start"])])[0]
        else:
            idx = int(np.argmax(np.diff(np.intersect1d(mc[mc < np.max(cand["end"])],
                                                       mc[mc > np.min(cand["start"])]))))  # maxend and minstart
            cddm = np.intersect1d(mc[mc < np.max(cand["end"])],
                                  mc[mc > np.min(cand["start"])])[idx]
        cddm = pd.Interval(cddm, cddm)
        cddp = np.array([pd.Interval(a, b) for a, b in zip(cand["start"], cand["end"])])
    return np.array([cand, cddm, cddp])


def peaks_coordinates(strand, cand, candmax, candpos):
    if strand == "+":
        # Saving plus1 and minus1 coordinates
        # If the real maximum is among the candidates
        if np.where(np.array([i.overlaps(candmax) for i in candpos]))[0][0] >= cand.shape[0] - 1:
            minus1 = cand.iloc[cand.shape[0] - 2]
            plus1 = cand.iloc[cand.shape[0] - 1]
        # Otherwise
        else:
            minus1 = cand.iloc[np.where(np.array([i.overlaps(candmax) for i in candpos]))[0][0]]
            # following peak is plus1
            plus1 = cand.iloc[np.where(np.array([i.overlaps(candmax) for i in candpos]))[0][0] + 1]
    else:
        if np.where(np.array([i.overlaps(candmax) for i in candpos]))[0][0] >= cand.shape[0] - 1:
            plus1 = cand.iloc[cand.shape[0] - 2]
            minus1 = cand.iloc[cand.shape[0] - 1]
        # Otherwise
        else:
            plus1 = cand.iloc[np.where(np.array([i.overlaps(candmax) for i in candpos]))[0][0]]
            minus1 = cand.iloc[np.where(np.array([i.overlaps(candmax) for i in candpos]))[0][0] + 1]
    return plus1["start"], plus1["end"], minus1["start"], minus1["end"]


def plusminus1(i, p, nu, pi, ni, mx):
    try:
        # positions of nucleosomes that overlap the promoter region
        olaps = np.where(np.array([pi[i].overlaps(n) for n in ni]))[0]
        # If there are overlaps
        if olaps.shape[0] > 0:
            candidates = find_candidates(p.iloc[i]["strand"], nu, olaps, mx)
            plmi1 = peaks_coordinates(p.iloc[i]["strand"], candidates[0], candidates[1], candidates[2])
            return plmi1
        else:
            return np.nan, np.nan, np.nan, np.nan
    except:
        return np.nan, np.nan, np.nan, np.nan


def extended_annotation(a, r, p):
    newinfo = pd.DataFrame({'plus1st': np.array([r[c][0] for c in range(len(r))]),
                            'plus1en': np.array([r[c][1] for c in range(len(r))]),
                            'minus1st': np.array([r[c][2] for c in range(len(r))]),
                            'minus1en': np.array([r[c][3] for c in range(len(r))]),
                            "ndr": np.maximum(np.array([r[c][0] for c in range(len(r))]) -
                                              np.array([r[c][3] for c in range(len(r))]),
                                              np.array([r[c][2] for c in range(len(r))]) -
                                              np.array([r[c][1] for c in range(len(r))]))})
    newinfo["ndr_pattern"] = (newinfo["plus1en"] - newinfo["plus1st"] == 150) & \
                             (newinfo["minus1en"] - newinfo["minus1st"] == 150) & \
                             (newinfo["ndr"] > 40)
    a["pStart"] = p["start"]
    a["pEnd"] = p["end"]

    finalannot = a.iloc[1:, ]
    finalannot = finalannot.reset_index(drop=True)
    finalannot = finalannot.join(newinfo)
    return finalannot


###-END FUNCTIONS-###


annotatio = sys.argv[1]
profil = sys.argv[2]
nucleosom = sys.argv[3]
chr_nam = sys.argv[4]
output = sys.argv[5]

ndr_pattern(annotatio, profil, nucleosom, chr_nam, output)
