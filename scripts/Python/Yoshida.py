import configuration as conf
import utilities
import tree_cluster as tc
import peak_clusters as pk

import pandas as pd
import os
import igraph as ig
import numpy as np
import pdb

# Class to create and access data provided by Yoshida et al.
#
# Input:
# @param raw_data_file the csv file provided by Yoshida et al as their 
# Table S2A.
#
# Output:
# @param matrix_file path to file containing signal (rows are loci, columns
#  are cell types)
# @param bed_file path to a bed file describing locations of loci on mm10
class Yoshida_data:
    
    raw_data_file = conf.INPUT_DATA + "Yoshida_raw_download.csv"
    
    yoshida_dir = conf.DATA_DIR + "Yoshida/"
    matrix_file = yoshida_dir + "Yoshida_matrix.csv"
    bed_file = yoshida_dir + "Yoshida.bed"
    
    def __init__(self):
        if not os.path.isdir(self.yoshida_dir):
            os.mkdir(self.yoshida_dir)
        if not os.path.isfile(self.matrix_file):
            print("creating Yoshida matrix...")
            self.create_matrix()
        if not os.path.isfile(self.bed_file):
            print("creating Yoshida bed...")
            self.create_bed()
    
    def create_bed(self):
        d = pd.read_csv(self.raw_data_file, sep=",")
        valid_chr = ["chr" + str(i) for i in range(20)] + ["chrX"]
        d = d[d["chrom"].isin(valid_chr)]
      
        s = d["Summit"].to_numpy()
        chrom = d["chrom"]
        p = d["_-log10_bestPvalue"]
        
        t = pd.DataFrame({'chr':chrom,
                           'chrStart':s,
                           'chrEnd':s,
                           'p':p})
        t.to_csv(self.bed_file, sep="\t", header=None,
                 index=False)
        
    def load_bed(self):
        d = pd.read_csv(self.bed_file, sep="\t", header=None)
        return d
        
    def create_matrix(self):
        d = pd.read_csv(self.raw_data_file, sep=",")
        d = d.iloc[:,8:99]
        d.to_csv(self.matrix_file, sep=",", index=False)
        
    def load_matrix(self, as_dataframe=False):
        d = pd.read_csv(self.matrix_file, sep=",")
        if not as_dataframe:
          return d.to_numpy()
        else:
          return d
    
    def load_raw_data(self):
        d = pd.read_csv(self.raw_data_file, sep=",")
        return d
        
    
# Creates the igraph Graph object for the Yoshida tree.
#
# Constructor parameters
#   @param just_root should the graph just contain descendants of "LTHSC.34-.BM
#   @param edgelist_file csv file containing edgelist, columns node1, node2
#
# Object fields
#   @field g the igraph Graph object
class Yoshida_tree:
    
    edgelist_file = conf.INPUT_DATA + "Yoshida_edge_list.csv"
   
    def load_igraph(self, just_root=True):
        e = self.load_edgelist()
        e.to_csv("temp_Yoshida_graph.csv", sep=" ", index=False, header=False)
        g = ig.Graph.Read_Ncol("temp_Yoshida_graph.csv", directed=True)
        
        if just_root:
            root = "LTHSC.34-.BM"
            descendants = g.neighborhood(root, 
                                         order=1000, mode="out")
            g = g.induced_subgraph(descendants)
                
        os.remove("temp_Yoshida_graph.csv")
        
        return g
    
    def load_edgelist(self):
        e = pd.read_csv(self.edgelist_file, sep=",")
        e = e.loc[e["node1"].notna() & e["node2"].notna()]
        
        return e
    
 
    
    
    
        
        
        
        
        
    
  
            
        
        
        
        
        