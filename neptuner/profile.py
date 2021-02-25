import multiprocessing as mp
import pandas as pd
import numpy as np
from scipy.stats import poisson
import sys, os, math

###-MAIN FUNCTION-###
# USAGE: python3 profile.py <input.bed> outfile_noextension


def profile(ipt, opt):
    ## Calculates the nucleosome profile from bed file containing read coordinates

    # Assertions on parameters
    # Input is a string and a bed file
    assert isinstance(ipt, str)
    assert os.path.isfile(ipt)
    assert ipt.endswith(".bed")
    # Output is a bedgraph file
    assert isinstance(opt, str)
    #assert opt.endswith(".bedgraph")

    # Number of cpus
    ncpu = mp.cpu_count()
    # Bedfile loading and preprocessing
    global bedfile
    bedfile = bedfile_preproc(ipt, opt)
    chro = bedfile["chromosome"][0]
    # Identifying highest coordinate
    fin = np.max(bedfile["end"])
    # final array declaration
    global old
    old = np.zeros(fin + 1)
    # Rows per subfile (used in the subfile function)
    global rps
    rps = math.ceil(bedfile.shape[0] // ncpu)

    # Parallel profile computation
    pool = mp.Pool(ncpu)  # mp.cpu_count())
    prof = pool.map(subfiles, [p for p in range(ncpu)])
    pool.close()
    prof = np.sum(prof, axis=0)

    # Normalization
#    norm_prof = normalize(prof)

    # Write output
#    pd.DataFrame({'chr': np.array(chro), 'start': np.arange(norm_prof.size) + 1,
#                  'end': np.arange(norm_prof.size) + 2, "score": norm_prof}).to_csv(opt, index=False, header=False)
    pd.DataFrame({'chr': np.array(chro), 'start': np.arange(prof.size) + 1,
                  'end': np.arange(prof.size) + 2, "score": prof}).to_csv(opt+".bedgraph", index=False, header=False)
    print(max(prof))
#    print(max(norm_prof))
    return


###-INTERNAL FUNCTIONS-###


def readstack(o, ss, ee):
    o[ss:ee+1] += 1


def getpvals(st, en, outp):
    lng = max(en)
    old = np.array([0] * lng)
    [readstack(old, sta, end) for sta, end in np.nditer([st, en])][-1]

    expect = poisson.rvs(np.mean(old), size=10000)
    pvals = poisson.pmf(old, np.mean(old))

    pd.DataFrame({'pval': pvals}).to_csv(outp+"_pvals.csv",
                                         index=False, header=False)
    return


def bedfile_preproc(bfi, outpu):
    # It loads the csv file given by the user, name the fields and filters reads shorter than 100bp or
    # longer than 200bp. After the preprocessing the file is ready to be used for the profile calculation.

    # Reading and setting up input file
    bf = pd.read_csv(bfi, sep='\t', header=None)
    # Additional controls
    # First column (chromosome) is a string and all are the same
    assert bf.shape[1] >= 3
    assert bf[0].unique().size == 1
    #assert isinstance(bf[0][0], str)
    # Second column (start) is integer
    a = np.array(bf[1])
    assert all([isinstance(a[i], np.int64) for i in range(bf.shape[0])])
    # Third column (end) is integer
    a = np.array(bf[2])
    assert all([isinstance(a[i], np.int64) for i in range(bf.shape[0])])

    # Rename columns
    bf = bf.iloc[:, 0:3]
    bf.rename(columns={0: "chromosome", 1: "start", 2: "end"}, inplace=True)
    bf["length"] = bf["end"] - bf["start"]
    # Filtering reads shorter than 100bp or longer than 200bp
    #bf = bf[(bf["length"] < 200) & (bf["length"] > 100)]
    bf = bf.reset_index(drop=True)

    getpvals(bf["start"], bf["end"], outpu)

    return bf


def subfiles(i):
    # The bedfile is divided in a number of subfiles which are going to be processed indipendently.
    # i is the part of the original file, it goes from 0 to the number of used cpus minus one.
    # For each of the subfiles a partial result array is created and the is going to be stored the
    # result of the triangle kernel calculation.
    # N.B. The function takes only one parameter and refers to the global environment "rps" value
    # to correctly use the pool.map function from the multiparallel library.

    global bedfile
    global rps
    global old
    # Identifies the subfile
    subf = bedfile.iloc[i * rps: rps * (i + 1)]
    # Defines the subfile result file initialized with zeros
    new = np.zeros(np.max(subf["end"]) - np.min(subf["start"]) + 1)
    st = np.array(subf["start"])
    start = max(np.min(subf["start"]) - 1, 0)
    en = np.array(subf["end"])
    # Index of the smallest coordinate for this subfile
    # Used to correctly integrate the result chunks in the correct position of the final result
    count = np.min(subf["start"])
    # Apply the triangle kernel to the subfile
    new = [triangle(new, s, e, count) for s, e in np.nditer([st, en])][-1]
    # Integrate the result chunk in the complete result vector at the corresponding position
    old[np.min(subf["start"]):np.max(subf["end"])+1] = old[np.min(subf["start"]):np.max(subf["end"])+1] + new
    return old


def triangle(tot, st, en, co):
    # Calculate the triangle kernel for one read/line of the bedfile.
    # The triangle is isoscele and has area of 1

    base = en - st + 1
    height = int(math.floor(base / 2))
    curtrg = np.zeros(base)
    # Lambda functions for even and odd length cases
    even1 = lambda x: (4 * (x - st)) / (base ** 2 - 1)
    odd1 = lambda x: (4 * (x - st)) / base ** 2
    # Even case
    if 2 * height < base:
        # First half
        curtrg[:height] = even1(np.arange(st, st + height))
        # Middle
        curtrg[height] = 2 / base
        # Second half
        curtrg[height + 1:base] = np.flip(curtrg[:height])
    # Odd length
    if 2 * height == base:
        # First half
        curtrg[:height] = odd1(np.arange(st, st + height))
        # Second half
        curtrg[height:base] = np.flip(curtrg[:height])
    # Update the resulting subvector
    tot[(st - co):(en + 1 - co)] = tot[(st - co):(en + 1 - co)] + curtrg
    return tot


def normalize(pr):
    # Normalize the final profile so that is in the range [0, 1]
    return np.divide(pr, max(pr))


###-END FUNCTIONS-###


# User parameteres
input = sys.argv[1]
output = sys.argv[2]
# Call main function
profile(input, output)
