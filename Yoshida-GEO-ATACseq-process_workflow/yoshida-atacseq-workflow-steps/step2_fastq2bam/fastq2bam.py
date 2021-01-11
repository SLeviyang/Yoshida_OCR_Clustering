import numpy as np
import pandas as pd
import os

# how many fastq to process concurrently
ACCESSION_BATCH_SIZE = 5

all_SRR=pd.read_csv('../GEO_info.csv')["SRR"].tolist()

# generate list of existing fastq files on s3
os.system("aws s3 ls s3://yoshida-atacseq/fastq/ | tr -s ' ' | cut -d ' ' -f 4 > existing_fastq.txt")
if os.path.getsize("existing_fastq.txt") > 0:
    existing_fastq = pd.read_csv('existing_fastq.txt', header=None)[0].astype("string")
else:
    existing_fastq = pd.Series([], dtype="string")

# generate list of existing bam files on s3
os.system("aws s3 ls s3://yoshida-atacseq/bam/ | tr -s ' ' | cut -d ' ' -f 4 > existing_bam.txt")
if os.path.getsize("existing_bam.txt") > 0:
    existing_bam = pd.read_csv('existing_bam.txt', header=None)[0].astype("string")
else:
    existing_bam = pd.Series([], dtype="string")

#os.system("rm existing_fastq.txt existing_bam.txt")

# find accessions with no bam file but with both fastq
active_SRR = []
for cSRR in all_SRR:
    fq1 = cSRR + "_1.fastq"
    fq2 = cSRR + "_2.fastq"
    if any(existing_fastq.str.match(fq1)) and any(existing_fastq.str.match(fq2)):
      if not any(existing_bam.str.contains(cSRR)):
          active_SRR.append(cSRR)
      else:
          print("bam already exists for " + cSRR)
    else:
        print("fastq not present for " + cSRR)

# iterate over active SRR
nm = len(active_SRR)
if nm != 0:
  starts = np.arange(0, nm, ACCESSION_BATCH_SIZE)
  ends = np.append(np.delete(starts, 0), nm)
#
  nbatch=len(starts)
  for i in range(1):
    cs = starts[i]
    ce = ends[i]
    d = pd.Series(active_SRR[cs:ce], dtype="string")
    d.to_csv("accessions.txt", header=False, index=False)
#
    print("sending the following accessions to cromwell...")
    print(active_SRR[cs:ce])
#
    os.system("sudo java -jar ../../cromwell.jar run -i fastq2bam_inputs.txt fastq2bam.wdl")
    print("deleting cromwell execution directory to save space")
    os.system("sudo rm -r *-exec* *-workflow*")
