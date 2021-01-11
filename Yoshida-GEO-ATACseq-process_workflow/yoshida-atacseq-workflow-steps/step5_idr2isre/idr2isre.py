import numpy as np
import pandas as pd
import os
import pdb

GEO_d=pd.read_csv('../GEO_info.csv')
cell_names = set(GEO_d["cell_name"].tolist())

# generate list of cell types with idr in s3
os.system("aws s3 ls s3://yoshida-atacseq/idr/ | tr -s ' ' | cut -d ' ' -f 4 > existing_idr.txt")
if os.path.getsize("existing_idr.txt") > 0:
    existing_idr = pd.read_csv('existing_idr.txt', header=None)[0].astype("string")
else:
    existing_idr = pd.Series([], dtype="string")

# generate list of cell types with ISRE in s3
os.system("aws s3 ls s3://yoshida-atacseq/ISRE/ | tr -s ' ' | cut -d ' ' -f 4 > existing_ISRE.txt")
if os.path.getsize("existing_ISRE.txt") > 0:
    existing_ISRE = pd.read_csv('existing_ISRE.txt', header=None)[0].astype("string")
else:
    existing_ISRE = pd.Series([], dtype="string")

active_cell_names = []
for cc in cell_names:
    if any(existing_ISRE.str.contains(cc, regex=False)):
      print("ISRE OCR already exist for " + cc)
    else:
      if not any(existing_idr.str.contains(cc, regex=False)):
        print("idr do not exist for " + cc)
      else:
        active_cell_names.append(cc)

active_cell_names.sort()
pd.Series(active_cell_names, dtype="string").to_csv("missing_cell_types.txt", header=False, index=False)


# # iterate over active SRR
nm = len(active_cell_names)

batch_size = 6
start_i = np.arange(0, nm, batch_size)
end_i = [min(i, nm) for i in (start_i + batch_size)]
nb = len(start_i)

for i in range(nb):
    cc = active_cell_names[start_i[i]:end_i[i]]
    d = pd.Series(cc, dtype="string")
    d.to_csv("cell_type.txt", header=False, index=False)
#
    print("sending the following cell type to cromwell...")
    print(cc)
#
    os.system("sudo java -jar ../../cromwell.jar run idr2isre.wdl")
    print("deleting cromwell execution directory to save space")
    os.system("sudo rm -r *-exec* *-workflow*")
