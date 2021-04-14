import workflow
import TF as tf
import manuscript_master as mm

import sys
import pdb
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pylab as plt
import configuration as conf

def create_M_construction_table():
    idrs = [0.01, 0.05, 0.10, 0.15]
    tb = []
    for idr in idrs:
      w = workflow.workflow("idr", edge_FDR=0.001,
                              idr_FDR=idr)
      m = w.load_matrix()
      m = m[np.sum(m, 1) > 0,]
      total = len(m)
      cs = np.sum(np.sum(m, 1) <= 2)/total
      
      tb.append({'idr':idr,
                 'loci':total,
                 'cell specific':cs})
      
    tb = pd.DataFrame(tb)
        
    outfile = mm.TABLE_DIR + "M_construction.csv" 
    tb.to_csv(outfile, index=False)
    
    return tb
        
        