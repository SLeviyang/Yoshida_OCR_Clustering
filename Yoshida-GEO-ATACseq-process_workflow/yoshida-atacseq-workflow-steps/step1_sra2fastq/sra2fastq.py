import numpy as np
import pandas as pd
import os

# how many SRA to download concurrently
ACCESSION_BATCH_SIZE = 5

all_SRR=pd.read_csv('../GEO_info.csv')["SRR"].tolist()

# generate list of existing fastq files on s3
os.system("aws s3 ls s3://yoshida-atacseq/fastq/ | tr -s ' ' | cut -d ' ' -f 4 > existing_fastq.txt")
if os.path.getsize("existing_fastq.txt") > 0:
    existing_fastq = pd.read_csv('existing_fastq.txt', header=None)[0].astype("string")
else:
    existing_fastq = pd.Series([], dtype="string")

os.system("rm existing_fastq.txt")

# find all the SRR that have not been processed
missing_SRR = []
for cSRR in all_SRR:
    if not any(existing_fastq.str.contains(cSRR)):
        missing_SRR.append(cSRR)
    else:
        print("SRR already in s3:")
        print(cSRR)

# iterate over missing SRR
nm = len(missing_SRR)
if nm != 0:
  starts = np.arange(0, nm, ACCESSION_BATCH_SIZE)
  ends = np.append(np.delete(starts, 0), nm)
#
  nbatch=len(starts)
  for i in range(nbatch):
    cs = starts[i]
    ce = ends[i]
    d = pd.Series(missing_SRR[cs:ce], dtype="string")
    d.to_csv("accessions.txt", header=False, index=False)
#
    print("sending the following accessions to cromwell...")
    print(missing_SRR[cs:ce])
#
    os.system("java -jar ../cromwell.jar run -i sra2fastq_inputs.txt sra2fastq.wdl")
    print("deleting cromwell execution directory to save space")
    os.system("sudo rm -r *-exec* *-workflow*")
