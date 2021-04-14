import workflow
import TF as tf
import pdb
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pylab as plt

FIGURE_DIR = "../../writeup/manuscript_v1/figures/"

# helper functions

def make_TF_barplots(w,
                     TF,
                     min_q=5,
                     outfile=None):

  m = w.load_matrix()
  c = tf.TF_ChIP(w)
  tb = c.load_ChIP_table()
  tb = tb[(tb["TF"] == TF) & tb.active]
  b_list = []
  for i in tb.index:
    print(i)
    b = c.intersect_bed(i, min_q=min_q)
    b.insert(b.shape[1], "cell_type", tb["cell_type"][i])
    b.insert(b.shape[1], "dataset_index", i)
    b_list = b_list + [b]
    
  b = pd.concat(b_list)
  #b = b[b["cluster"] <= 20]
  
  cell_types = np.array(w.get_cell_types())
  d_list = []
  for i in b["dataset_index"].unique():
    print(i)
    cb = b[b["dataset_index"] == i]
    ct = cb["cell_type"].iloc[0]
    counts = cb["cluster"].value_counts()
    d = pd.DataFrame({"cluster":[str(x) for x in counts.index],
                    "count":counts})
    d.insert(d.shape[1], "cell_type", ct)
    d.insert(d.shape[1], "percent",
             d["count"]/len(cb))
    
    # now counts of local
    cbloc = cb[(cb["nopen"] <= 5) & (cb["cluster"] == -1)]
    counts1 = np.sum(m[cbloc["row"],:], 0)
    d1 = pd.DataFrame({"cluster":cell_types,
                    "count":counts1})
    d1.insert(d1.shape[1], "cell_type", ct)
    d1.insert(d1.shape[1], "percent",
             d1["count"]/len(cb))
    
    # add on total loc
    d2 = pd.DataFrame({"cluster":"cell specific",
                       "count":len(cbloc),
                       "cell_type":ct,
                       "percent":len(cbloc)/len(cb)},
                      index=[0])
                      
    
    d_list.append(d)
    d_list.append(d1)
    d_list.append(d2)
    
  d = pd.concat(d_list)
  d = d[d["percent"] > 0.05]
  d = d[d["cluster"] != "-1.0"]
  
  
  sns.barplot(x = "cluster", y="percent", 
              hue="cell_type",
              data=d)
  plt.xticks(rotation=90)
  
  if outfile is not None:
    plt.savefig(outfile, bbox_inches='tight')
    
  return d

##############################################
# figures
  
def make_significant_TF_figure(w):
  #w = workflow.workflow("idr", 0.001, 0.01)
  t = tf.TF_motif(w)
  e = t.enrichment_matrix(sig=True, cutoff=0.25)
  e = e.loc[:,np.sum(e > 0, 0) > 0]
  sns.clustermap(e > 0, row_cluster=False)
  
  outfile = FIGURE_DIR + "TF_" + "sig" + ".jpeg"
  plt.savefig(outfile, bbox_inches='tight')
  

def make_EOMES_figure():
    
  w = workflow.workflow("idr", 0.001, 0.01)
  outfile = FIGURE_DIR + "TF_" + "EOMES" + ".jpeg"
  d = make_TF_barplots(w, TF="EOMES", min_q=5,
                   outfile=outfile)
  
  return d
  
def make_PAX5_figure():
    
  w = workflow.workflow("idr", 0.001, 0.01)
  outfile = FIGURE_DIR + "TF_" + "PAX5" + ".jpeg"
  make_TF_barplots(w, TF="PAX5", min_q=5,
                   outfile=outfile)



