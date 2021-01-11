import pandas as pd
import numpy as np
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
            print([chr, si, N, cd.iloc[si]["chrStart"]])
            ei = find_intersection_range(cd, si, N)
            sub_cd_q = cd.iloc[si:ei]["q"]
            q_vals = sub_cd_q.idxmax()
            q_vals_i = si + np.argmax(sub_cd_q)
            #if cd.iloc[q_vals]["chrStart"] == 4414281:
            #    pdb.set_trace()
            #if cd.iloc[q_vals]["chrStart"] == 4414682:
            #    pdb.set_trace()
            merged_indices.append(q_vals)
            start_indices.append(si)
            end_indices.append(ei)
            si = find_intersection_range(cd, q_vals_i, N)

    d_filtered = d.iloc[merged_indices]
    d_filtered.to_csv(output_bed_file,
                      header=None, sep="\t", index=False)


d = pd.read_csv(sys.argv[1], header=None, sep="\t",
                names=["chr", "chrStart", "chrEnd", "q"])
filter_nonintersecting_bed_file(d, sys.argv[2])
    

