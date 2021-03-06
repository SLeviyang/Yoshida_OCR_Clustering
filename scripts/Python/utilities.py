import pandas as pd
import numpy as np
import igraph as ig
import re
import pdb
import os

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

def ontology_analysis(genes):
      
    
    gene_file = "temp_genes_" + str(np.random.uniform()) + ".csv"
    out_file = "temp_genes_out_" + str(np.random.uniform()) + ".csv"
    
    pd.Series(genes).to_csv(gene_file, sep=",", 
                            index=False, header=False)
    os.environ["gf"] = gene_file
    os.environ["of"] = out_file
    os.system("/usr/local/bin/Rscript ../R/ontology.R $gf $of")
      
    z =  pd.read_csv(out_file, sep=",")
    z = z[z["term_size"] < 500]
    os.remove(gene_file)
    os.remove(out_file)
    return z

def find_root(g):
     # check that g is a tree
    parents = [len(g.neighbors(x, mode="IN")) for x in g.vs]
    if not set(parents) == set([0,1]):
        False, -1
    root_ind = np.where(np.array(parents)==0)[0]
    if not len(root_ind) == 1: 
        False, -1
        
    return True, root_ind[0]