import pandas as pd
import sys
import pdb

def find_intersection_range(d, i, N):
    end_val = d.iloc[i]["chrEnd"]
    j = i+1
    while (j < N and d.iloc[j]["chrStart"] <= end_val):
        j = j + 1

    return (j)

def filter_nonintersecting_bed_file(d, output_bed_file):

    merged_indices = []
    start_indices = []
    end_indices = []

    all_chr = d["chr"].unique().tolist()
    for chr in all_chr:
        cd = d.loc[d["chr"] == chr]
        si = 0
        N = len(cd)
        while (si < N):
            print([chr, si, N])
            ei = find_intersection_range(cd, si, N)
            #q_vals = cd.iloc[si:ei]["q"].idxmax()
            # take the first one
            merged_indices.append(cd.index[si])
            start_indices.append(si)
            end_indices.append(ei)
            si = ei

    d_filtered = d.iloc[merged_indices]
    d_filtered.to_csv(output_bed_file,
                      header=None, sep="\t", index=False)


d = pd.read_csv(sys.argv[1], header=None, sep="\t",
                names=["chr", "chrStart", "chrEnd", "peakName",
                       "score", "strand", "sequence", "reverseC"])
filter_nonintersecting_bed_file(d, sys.argv[2])
