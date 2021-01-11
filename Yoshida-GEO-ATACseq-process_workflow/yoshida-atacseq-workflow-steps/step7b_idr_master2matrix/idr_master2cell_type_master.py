import numpy as np
import pandas as pd
import os
import pdb

GEO_d=pd.read_csv('../GEO_info.csv')
cell_names = set(GEO_d["cell_name"].tolist())

# generate list of cell types with master in s3
os.system("aws s3 ls s3://yoshida-atacseq/master2/ | tr -s ' ' | cut -d ' ' -f 4 > existing_master.txt")
bad_e = ["nothing"]
if os.path.getsize("existing_master.txt") > 0:
    existing_masters = pd.read_csv('existing_master.txt', header=None)[0].astype("string")
    # some existing masters may be bad, likely because of an AWS download issue.  Check for those
    for e in existing_masters:
        os.environ["e"] = e
        os.system("aws s3 cp s3://yoshida-atacseq/master2/$e $e")
        d = pd.read_csv(e, header=None, sep="\t",
            names=["chr", "chrStart", "chrEnd", "peakName",
                   "score", "strand", "sequence", "reverseC", "distance"])
        if sum(d["distance"]) == 0:
            bad_e.append(e)
        os.system("rm $e")
else:
    existing_masters = pd.Series([], dtype="string")

bad_masters = pd.Series(bad_e, dtype="string")

active_cell_names = []
for cc in cell_names:
    if any(bad_masters.str.contains(cc, regex=False)):
       active_cell_names.append(cc)
    elif any(existing_masters.str.contains(cc, regex=False)):
        print("master already exist for " + cc)
    else:
        active_cell_names.append(cc)

# remove orcGFP (no bam files)
active_cell_names.remove('orcGFP')

active_cell_names.sort()
pd.Series(active_cell_names, dtype="string").to_csv("missing_cell_types.txt", header=False, index=False)


# # iterate over active SRR
nm = len(active_cell_names)

for i in range(nm):
    d = pd.Series(active_cell_names[i], dtype="string")
    d.to_csv("cell_type.txt", header=False, index=False)

    GEO_d.loc[GEO_d["cell_name"] == active_cell_names[i]]["SRR"].to_csv("accessions.txt", header=False, index=False)
#
    print("sending the following cell type to cromwell...")
    print(active_cell_names[i])
#
    os.system("sudo java -jar ../../cromwell.jar run -i isre_maser2cell_type_master_inputs.txt isre_master2cell_type_master.wdl")
    print("deleting cromwell execution directory to save space")
    os.system("sudo rm -r *-exec* *-workflow*")
