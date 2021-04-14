import workflow
import TF as tf
import pdb
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pylab as plt
import configuration as conf

FIGURE_DIR = "../../writeup/manuscript_v1/figures/"
TABLE_DIR = "../../writeup/manuscript_v1/tables/"

def load_cell_type_master_table():
    return pd.read_csv(TABLE_DIR + "cell_type_master_table.csv")
    

def make_cell_type_master_table():
    w = workflow.workflow("idr", 0.001, 0.01)
    tc = w.form_tree_cluster_from_fits(k=8)
    a = tc.get_assignments()
    
    colors = ["red", "green","blue", "yellow",
              "pink", "aqua", "grey", "brown"]
    name = ["Stem+", "MMP", "MF", "ILC", "NK", "T", "B", "B.pl"]
    a_colors = [colors[i] for i in a]
    a_names = [name[i] for i in a]
  
    tb = pd.DataFrame({'cell_type':w.get_cell_types(),
                       'cluster':a+1,
                       'color':a_colors,
                       'cluster_name':a_names})
    
    # load the short cell type names
    tb_cn = pd.read_csv(conf.INPUT_DATA + "short_cell_type_names.csv",
                        sep=",")
    tb = tb.merge(tb_cn, on="cell_type")
    tb.to_csv(TABLE_DIR + "cell_type_master_table.csv", index=False)
    
    return tb
    
    
    