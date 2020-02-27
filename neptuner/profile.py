import multiprocessing as mp
import pandas as pd
import numpy as np
import sys, os, math


def bedfile_preproc(bfi):
    # It loads the csv file given by the user, name the fields and filters reads shorter than 100bp or
    # longer than 200bp. After the preprocessing the file is ready to be used for the profile calculation.

    # Reading and setting up input file
    bf = pd.read_csv(bfi, sep='\t', header=None)
    # Additional controls
    # First column (chromosome) is a string and all are the same
    assert bf.shape[1] == 3
    assert bf[0].unique().size == 1
    assert isinstance(bf[0][0], str)
    # Second column (start) is integer
    a = np.array(bf[1])
    assert all([isinstance(a[i], np.int64) for i in range(bf.shape[0])])
    # Third column (end) is integer
    a = np.array(bf[2])
    assert all([isinstance(a[i], np.int64) for i in range(bf.shape[0])])

    # Rename columns
    bf.rename(columns={0: "chromosome", 1: "start", 2: "end"}, inplace=True)
    bf["length"] = bf["end"] - bf["start"]
    # Filtering reads shorter than 100bp or longer than 200bp
    bf = bf[(bf["length"] < 200) & (bf["length"] > 100)]
    bf = bf.reset_index(drop=True)
    return bf


def subfiles(i):
    # The bedfile is divided in a number of subfiles which are going to be processed indipendently.
    # i is the part of the original file, it goes from 0 to the number of used cpus minus one.
    # For each of the subfiles a partial result array is created and the is going to be stored the
    # result of the triangle kernel calculation.
    # N.B. The function takes only one parameter and refers to the global environment "rps" value
    # to correctly use the pool.map function from the multiparallel library.

    # Identifies the subfile
    subf = bedfile.iloc[i * rps: rps * (i + 1)]
    # Defines the subfile result file initialized with zeros
    new = np.zeros(np.max(subf["end"]) - np.min(subf["start"]) + 1)
    st = np.array(subf["start"])
    en = np.array(subf["end"])
    # Index of the smallest coordinate for this subfile
    # Used to correctly integrate the result chunks in the correct position of the final result
    count = np.min(subf["start"])
    # Apply the triangle kernel to the subfile
    new = [triangle(new, s, e, count) for s, e in np.nditer([st, en])][-1]
    # Integrate the result chunk in the complete result vector at the corresponding position
    old[np.min(subf["start"]) - 1:np.max(subf["end"])] = old[np.min(subf["start"] - 1):np.max(subf["end"])] + new
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


#tsti = time.time()
#
#bfi='/home/pejo/Scrivania/nept_pattern/correct_files/chr1.bed'


# User parameteres
input = sys.argv[1]
output = sys.argv[2]

# Assertions on parameters
# Input is a string and a bed file
assert isinstance(input, str)
assert os.path.isfile(input)
assert input.endswith(".bed")
# Output is a bedgraph file
assert isinstance(output, str)
assert output.endswith(".bedgraph")

# Number of cpus
ncpu = mp.cpu_count()
# Bedfile loading and preprocessing
bedfile = bedfile_preproc(input)
chr = bedfile["chromosome"][0]
# Identifying highest coordinate
fin = np.max(bedfile["end"])
# final array declaration
old = np.zeros(fin + 1)
# Rows per subfile (used in the subfile function)
rps = math.ceil(bedfile.shape[0] / ncpu)

# Parallel profile computation
pool = mp.Pool(ncpu)  # mp.cpu_count())
prof = pool.map(subfiles, [p for p in range(ncpu)])
pool.close()
prof = np.sum(prof, axis=0)

# Normalization
norm_prof = normalize(prof)

# Write output
pd.DataFrame({'chr': np.array(chr), 'start': np.arange(norm_prof.size) + 1,
              'end': np.arange(norm_prof.size) + 2, "score": norm_prof}).to_csv(output, index=False, header=False)

