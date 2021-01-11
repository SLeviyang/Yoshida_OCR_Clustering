import pandas as pd
import os
import pdb
import numpy as np
import matplotlib.pyplot as plt

# collects summits into +/- 250 bp windows across all bed files in directory
def merge_bed_file(dir, output_bed_file):
    files = os.listdir(dir)
    paths = [dir + f for f in files]

    if (not os.path.isfile(output_bed_file)):
        os.environ["files"] = " ".join(paths)
        # chr, chrStart, chrEnd, summit, q
        os.system("cat $files | cut -f 1,2,3,9,10 > merged.txt")
        os.system("awk '{print $1 \"\t\" $2+$5-250 \"\t\" $2+$5+250 \"\t\" $4}' merged.txt  > summit_merged.txt")

        os.environ["sm"] = "summit_merged.txt"
        os.environ["obd"] = output_bed_file

        os.system("sort -k1,1 -k2,2n -k4nr $sm > $obd")
        os.system("rm summit_merged.txt merged.txt")
    else:
        print("merged and sorted file already exists")


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
            q_vals = cd.iloc[si:ei]["q"].idxmax()
            merged_indices.append(q_vals)
            start_indices.append(si)
            end_indices.append(ei)
            si = ei

    d_filtered = d.iloc[merged_indices]
    d_filtered.to_csv(output_bed_file,
                      header=None, sep="\t", index=False)

def construct_OCR_matrix(master_peaks_bed_file,
                         dir):

    files = os.listdir(dir)
    paths = [dir + f for f in files]
    os.environ["master_bed"] = master_peaks_bed_file
    os.environ["PATH"] = os.environ["PATH"] + ":/opt/local/bin"
    os.system("bedtools sort -i $master_bed > master_sorted.bed")

    OCR_vals = []
    for p in paths:
        print(p)
        os.environ["bed"] = p
        os.system("cut -f 1,2,3 $bed > sample_temp.bed")
        os.system("bedtools sort -i sample_temp.bed > sample_temp_sorted.bed")
        os.system("bedtools closest -a master_sorted.bed -b sample_temp_sorted.bed -d -t first | cut -f 1,2,3,4,8 > sample_closest.bed")
        # need to konw how many fiels in each line
        d = pd.read_csv("sample_closest.bed", header=None, sep="\t",
                        names=["chr", "chrStart", "chrEnd", "q", "distance"])
        vals = (d["distance"] == 0)*1
        OCR_vals.append(vals.tolist())

    m = np.array(OCR_vals)

    return m
