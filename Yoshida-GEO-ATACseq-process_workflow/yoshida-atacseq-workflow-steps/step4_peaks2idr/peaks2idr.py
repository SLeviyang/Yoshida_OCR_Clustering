import numpy as np
import pandas as pd
import os

GEO_d=pd.read_csv('../GEO_info.csv')
cell_names = set(GEO_d["cell_name"].tolist())

# generate list of cell types with peaks in s3
os.system("aws s3 ls s3://yoshida-atacseq/peaks/ | tr -s ' ' | cut -d ' ' -f 3 > existing_peaks.txt")
if os.path.getsize("existing_peaks.txt") > 0:
    existing_peaks = pd.read_csv('existing_peaks.txt', header=None)[0].astype("string")
else:
    existing_peaks = pd.Series([], dtype="string")

# generate list of cell types with peaks in s3
os.system("aws s3 ls s3://yoshida-atacseq/idr/ | tr -s ' ' | cut -d ' ' -f 4 > existing_idr.txt")
if os.path.getsize("existing_idr.txt") > 0:
    existing_idr = pd.read_csv('existing_idr.txt', header=None)[0].astype("string")
else:
    existing_idr = pd.Series([], dtype="string")

active_cell_names = []
for cc in cell_names:
    if any(existing_idr.str.contains(cc, regex=False)):
      print("idr already exist for " + cc)
    else:
      if not any(existing_peaks.str.contains(cc, regex=False)):
        print("peaks do not exist for " + cc)
      else:
        active_cell_names.append(cc)

active_cell_names.sort()
pd.Series(active_cell_names, dtype="string").to_csv("missing_cell_types.txt", header=False, index=False)


# # iterate over active SRR
nm = len(active_cell_names)

for i in range(nm):
    d = pd.Series(active_cell_names[i], dtype="string")
    d.to_csv("cell_type.txt", header=False, index=False)
#
    print("sending the following cell type to cromwell...")
    print(active_cell_names[i])
#
    os.system("sudo java -jar ../../cromwell.jar run -i peaks2idr_inputs.txt peaks2idr.wdl")
    print("deleting cromwell execution directory to save space")
    os.system("sudo rm -r *-exec* *-workflow*")
