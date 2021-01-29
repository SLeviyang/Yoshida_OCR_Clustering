import pandas as pd
import numpy as np
import re
import pdb

# find positions in l2 that match the elements of l1, 
# -1 if not present
def match(l1, l2):
    d2 = pd.DataFrame({'x':l2, 
                       'index':np.arange(len(l2))})
    d1 = pd.DataFrame({'x':l1})
    
    merged = d1.merge(d2, how="left")
    for i in range(len(merged)):
        if np.isnan(merged.at[i,"index"]):
            merged.at[i,"index"] = -1
    
    return merged["index"].to_list()

def rint():
    return str(np.round(1E13*np.random.uniform()).astype(np.int64))

# Computes nucleotide frequency for each sequence in a list of sequences
def nucleotide_frequency_table(seqs, nucs=["A", "C", "G", "T", 
                                           "GC", "CG"]):
    
    chars = np.array([list(s) for s in seqs])
    nuc_table = np.zeros([len(chars),len(nucs)])
    for i in range(len(nucs)):
      nuc_table[:,i] = [len(re.findall(nucs[i], s))/len(s) for s in seqs]
       
    tb = pd.DataFrame(nuc_table, columns=nucs)
    return tb