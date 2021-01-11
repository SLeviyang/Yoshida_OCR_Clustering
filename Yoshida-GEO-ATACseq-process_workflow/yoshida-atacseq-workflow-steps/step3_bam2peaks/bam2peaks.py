import numpy as np
import pandas as pd
import os
import sys

full_df=pd.read_csv('../GEO_info.csv')
# for some reason "orcGFP" cell type bam files were not created, remove
full_df = full_df[full_df["cell_name"] != "orcGFP"]

all_cell_name=full_df["cell_name"].tolist()
all_cell_name=set(all_cell_name)

# check that all bam files have been created, otherwise repeat step 2!
os.system("aws s3 ls s3://yoshida-atacseq/bam/ | tr -s ' ' | cut -d ' ' -f 4 > existing_bam.txt")
existing_bam = set()
if os.path.getsize("existing_bam.txt") > 0:
    existing_bam1 = pd.read_csv('existing_bam.txt', header=None)[0].astype("string")
    for i in range(len(existing_bam1)):
        z = existing_bam1[i]
        existing_bam.add(z[0:(len(z)-4)])

all_SRR = set(full_df["SRR"].tolist())
missing_SRR = all_SRR.difference(existing_bam)
if (len(missing_SRR) != 0):
    sys.exit("missing bam files!  cannot proceed to peak construction!")


# generate list of existing cell folders
os.system("aws s3 ls s3://yoshida-atacseq/peaks/ | tr -s ' ' | cut -d ' ' -f 3 > existing_peaks.txt")
existing_peaks = set()
if os.path.getsize("existing_peaks.txt") > 0:
    existing_peaks1 = pd.read_csv('existing_peaks.txt', header=None)[0].astype("string")
    for i in range(len(existing_peaks1)):
        z = existing_peaks1[i]
        existing_peaks.add(z[0:(len(z)-1)])


active_cell_names = list(all_cell_name.difference(existing_peaks))
already_cell_names = all_cell_name.intersection(existing_peaks)

# this file is not used below
pd.Series(active_cell_names).to_csv("missing_cell_types.txt", header=False, index=False)

# iterate over active_cell_names
nm = len(active_cell_names)

if nm != 0:
  for cell_name in active_cell_names:
    accession = full_df[full_df["cell_name"]==cell_name]["SRR"]
    groupName = pd.Series(cell_name)
    accession.to_csv("accessions.txt", header=False, index=False)
    groupName.to_csv("groupName.txt", header=False, index=False)
#
    print("groupName: " + cell_name)
    print("sending the following accessions to cromwell...")
    print(accession.tolist())
#
    os.system("sudo java -jar ../../cromwell.jar run -i bam2peaks_input.txt bam2peaks.wdl")
    print("deleting cromwell execution directory to save space")
    os.system("sudo rm -r *-exec* *-workflow*")
