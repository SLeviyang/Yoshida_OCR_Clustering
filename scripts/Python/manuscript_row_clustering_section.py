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


def create_row_clustering_information():
    edge_FDR = [0.01, 0.001, 0.0001]
    idr_FDR = [0.01, 0.05, 0.10, 0.15]
    in_tb = pd.DataFrame([[x,y] for x in idr_FDR for y in edge_FDR],
                    columns=["idr_FDR", "edge_FDR"]).sort_values("idr_FDR")
    
    tb = []
    for i in range(len(in_tb)):
      
      print([i, "of", len(in_tb)])
  
      w = workflow.workflow("idr", 
                            edge_FDR=in_tb["edge_FDR"].iloc[i],
                            idr_FDR=in_tb["idr_FDR"].iloc[i])
      m = w.load_matrix()
      m = m[np.sum(m, 1) > 2,:]
      total = len(m)
    
      c = w.load_clusters()
      base_c = np.sum(c["cluster"] <= 19)
      total_c = len(c)
    
      glob = np.sum(np.sum(m, 1) >= (m.shape[1]-2))
    
      tb.append({'total':total,
            'Louvain':total_c,
            'Louvain_p':total_c/total,
            'clusters':base_c,
            'clusters_p':base_c/total,
            'global':glob/total})
      
     
   
    tb_out = pd.concat([in_tb, pd.DataFrame(tb)], axis=1)
    
    outfile = mm.TABLE_DIR + "row_cluster_information.csv" 
    tb_out.to_csv(outfile, index=False)
    return tb_out

def make_row_cluster_table():
    w = workflow.workflow("idr", 0.001, 0.01)
    rcv = w.get_row_cluster_view()
    tb = rcv.get_cluster_info(20).loc[:,["cluster",
                                         "size",
                                        "promoters",
                                        "enhancers"]]
    tb["promoters"] = np.round(100*tb["promoters"])/100
    tb["enhancers"] = np.round(100*tb["enhancers"])/100
    
    outfile = mm.TABLE_DIR + "row_cluster_table.tex" 
    tb.to_latex(outfile, index=False)
    return tb

def make_row_cluster_matrix_heatmap():
    w = workflow.workflow("idr", 0.001, 0.01)
    rcv = w.get_row_cluster_view()
    m_list = rcv.form_m_list(min_cluster_size=30)
    ncluster = len(m_list)
    
    means = []
    indices = []
    for index in range(ncluster):
      cm = rcv.form_cluster_matrix(index)
      cm_mean = np.mean(cm, 0)
      means.append(cm_mean.tolist())
      indices.append(index)
        
    # arrange by column cluster
    ct = mm.load_cell_type_master_table()
    w_ct = np.array(w.get_cell_types())
    if not np.all(ct["cell_type"] == w_ct):
        sys.exit("mismatch between master table and workflow!")
    ct = ct.merge(pd.DataFrame({'cell_type':w.get_cell_types(),
                                'w_index':range(w_ct.shape[0])}), 
                  on="cell_type")
    col_order = ct.sort_values("cluster")["w_index"]
    means_m = np.array(means)[:,col_order]
    tb = pd.DataFrame(means_m)
    tb.index = indices
    sns.heatmap(tb, cmap="viridis")
    
    outfile = mm.FIGURE_DIR + "cluster_matrix" + ".jpeg"
    plt.savefig(outfile, bbox_inches='tight')
    
    return means_m

def make_bicluster_matrix_heatmap():
    w = workflow.workflow("idr", 0.001, 0.01)
    a = w.biclusters_to_response_control(k=12)
    
     # arrange by column cluster
    ct = mm.load_cell_type_master_table()
    w_ct = np.array(w.get_cell_types())
    if not np.all(ct["cell_type"] == w_ct):
        sys.exit("mismatch between master table and workflow!")
    ct = ct.merge(pd.DataFrame({'cell_type':w.get_cell_types(),
                                'w_index':range(w_ct.shape[0])}), 
                  on="cell_type")
    col_order = ct.sort_values("cluster")["w_index"]
    means_m = np.array(a)[:,col_order]
    tb = pd.DataFrame(means_m)
    sns.heatmap(tb, cmap="viridis")
    
    outfile = mm.FIGURE_DIR + "cluster_matrix_01" + ".jpeg"
    plt.savefig(outfile, bbox_inches='tight')
    

    
    
