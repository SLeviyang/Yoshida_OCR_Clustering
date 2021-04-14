import workflow
import manuscript_master as mm
import tree_cluster 

import sys
import pdb
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pylab as plt
import configuration as conf

def make_column_cluster_tree_figure():
     w = workflow.workflow("idr", 0.001, 0.01)
     
     fig, (ax1, ax2) = plt.subplots(1, 2)
     
     tc = w.form_tree_cluster_from_fits(3)
     tc.treeplot(target=ax1)
     
     tc = w.form_tree_cluster_from_fits(8)
     tc.treeplot(target=ax2)
     
     outfile = mm.FIGURE_DIR + "column_cluster_trees" + ".jpeg"
     plt.savefig(outfile, bbox_inches='tight')
     
     return None
 
def make_column_cluster_var_figure():
     w = workflow.workflow("idr", 0.001, 0.01)
     m = w.load_m_list()
     m = m[1:len(m)]
     
     m_tot = np.concatenate(m)
     m_mean = np.array([np.mean(mm) for mm in m])
     m_var = np.array([np.var(mm) for mm in m])
     
     R = m_tot.shape[0]
     C = m_tot.shape[1]
     
     v = [] 
     for k in range(1,13):
         
         if k > 1:
           tc = w.form_tree_cluster_from_fits(k=k)
           a = tc.get_assignments()
         else:
           a = np.zeros(C)
         var_cc = 0
         var_wr = 0
         var_wc = 0
         f_tot = 0
         for r in range(len(m)):
           for cl in range(k):
          
               Mij = m[r][:,a==cl]
               Mrc = np.mean(Mij)
               Mc = np.mean(Mij, axis=0)
               
               f = np.size(Mij)/(R*C)
               f_tot = f_tot + f
               
               var_cc = var_cc + f*(Mrc - m_mean)**2
               var_wr = var_wr + f*np.mean((Mij - Mc)**2)
               var_wc = var_wc + f*np.mean((Mc - Mrc)**2)
               
         v.append({'k':k,
                   'cc':var_cc/m_var,
                      'wr':var_wr/m_var,
                      'wc':var_wc/m_var,
                      'check':(var_cc + var_wr + var_wc)/m_var})
     
               
     v = pd.DataFrame(v)
     return (v)
                         
         
   