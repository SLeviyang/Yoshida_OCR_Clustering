import configuration as conf
import peak_clusters as pk
import utilities
import Yoshida

import os, sys
import pdb
import pandas as pd
import numpy as np
import glmnet_python 
import sklearn.preprocessing as skpre
import sklearn.linear_model as sklin

import matplotlib.pyplot as plt
import itertools as it
import seaborn as sns

from Bio.SeqUtils import GC
import qnorm


# TF motif score matrices are created by a call to
# create_motif_score_matrices() in the R script
# create_motif_score_matrices.R.  I use R since it has quick
# connection to the chromVAR motifs used in Yoshida et al.
class TF:
    
    motif_dir = conf.DATA_DIR + "TF_motifs/"
    raw_motif_dir = motif_dir + "raw_matrices_distal/"
    scaled_motif_dir = motif_dir + "scaled_matrices/"
    
    nucs = ["A", "C", "G", "T", "GC", "CG"]
     
    def __init__(self, max_cluster_number=20):
      self.ncluster = max_cluster_number + 1
      if not os.path.isdir(self.motif_dir) or not os.path.isdir(self.raw_motif_dir):
        sys.exit("TF score matrices have not been created.  Run R function!")
        
      for i in range(max_cluster_number+1):
          file = self.get_raw_matrix_filename(i, null=False) 
          nfile = self.get_raw_matrix_filename(i, null=True)
          if not os.path.isfile(file):
            sys.exit(["motif score matrix missing for cluster", i])
          if not os.path.isfile(nfile):
            sys.exit(["null motif score matrix missing for cluster", i])
            
      if not os.path.isdir(self.scaled_motif_dir):
          os.mkdir(self.scaled_motif_dir)
          self.create_genomics_tables(max_cluster_number)
          self.create_scaled_score_matrices(max_cluster_number)
          
    
    
    # I can't use the genomics table from peaks_cluster because the window
    # size of the TF search (201) is different than the window size of the 
    # loci (501).
    def create_genomics_tables(self, max_cluster_number):
  
      for i in range(max_cluster_number+1):
        print(["loading sequences for cluster", i])
        in_file = self.get_sequences_filename(i) 
        out_file = self.get_genomics_table_filename(i)
       
        seqs = pd.read_csv(in_file).iloc[:,0].tolist()
        nuc_tb = utilities.nucleotide_frequency_table(seqs, nucs=self.nucs)
        
        nuc_tb.to_csv(out_file, sep=",", index=False)
      
    # quantile normalize motif scores, normalizing across
    # motifs and across all clustersß
    def create_scaled_score_matrices(self, max_cluster_number):
        
      m_list = []
      null_m_list = []
      clusters = []

      # read in unscaled data created by R script
      for i in range(max_cluster_number+1):
        print(["loading scores for cluster", i])
        in_file = self.get_raw_matrix_filename(i, null=False) 
        in_nfile = self.get_raw_matrix_filename(i, null=True) 
         
        m = pd.read_csv(in_file, sep=",")
        nm = pd.read_csv(in_nfile, sep=",")
        
        if not len(m) == len(nm):
            sys.exit("sampled and null score matrices should have same dimensions")
        
        if i == 1:
            motif_names = m.columns
        m_list.append(m.to_numpy())
        null_m_list.append(nm.to_numpy())
        clusters = clusters + np.repeat(i, len(m)).tolist()
      
      m_pre = np.concatenate(m_list, axis=0)
      null_m_pre = np.concatenate(null_m_list, axis=0)
      
      m = qnorm.quantile_normalize(m_pre)
      null_m = qnorm.quantile_normalize(null_m_pre)
      #m = skpre.scale(m_pre)
      #null_m = skpre.scale(null_m_pre)
      
      tb = pd.DataFrame(m, columns=motif_names)
      tb.insert(0, "cluster", clusters)
      null_tb = pd.DataFrame(null_m, columns=motif_names)
      null_tb.insert(0, "cluster", clusters)
  
      # save scaled data, split by cluster
      for i, ctb in tb.groupby("cluster"):
        print(["writing scaled scores for cluster", i])  
        out_file = self.get_scale_matrix_filename(i, null=False)
        ctb.drop("cluster", axis=1).to_csv(out_file, sep=",", index=False)
    
      # save null scaled data, split by cluster
      for i, ctb in null_tb.groupby("cluster"):
        print(["writing null scaled scores for cluster", i])  
        out_file = self.get_scale_matrix_filename(i, null=True)
        ctb.drop("cluster", axis=1).to_csv(out_file, sep=",", index=False)

    ############### input filenames
    def get_raw_matrix_filename(self, i, null=False):
        if not null:
          out_file = self.raw_motif_dir + "motif_scores_cluster_" \
                        + str(i) + ".csv"
        else: 
          out_file = self.raw_motif_dir + "motif_scores_cluster_" \
                      + str(i) + "_permuted.csv"
        return out_file
    
    def get_sequences_filename(self, i):
        out_file = self.raw_motif_dir + "motif_scores_cluster_" \
                      + str(i) + "_sequences.csv"
        return out_file
    
    ############### output filenames
    def get_genomics_table_filename(self, i):
        out_file = self.scaled_motif_dir + "scaled_motif_scores_cluster_" \
                        + str(i) + "_genomics.csv"
                                    
        return (out_file)               
    
    def get_scale_matrix_filename(self, i, null=False):
        if not null:
          out_file = self.scaled_motif_dir + "scaled_motif_scores_cluster_" \
                        + str(i) + ".csv"
        else: 
          out_file = self.scaled_motif_dir + "scaled_null_motif_scores_cluster_" \
                      + str(i) + ".csv"
        return out_file
    
    ###############################  
    def get_motif_names(self):
        tb = pd.read_csv(self.get_scale_matrix_filename(0), sep=",",
                         nrows=1, header=None)
        motifs = tb.to_numpy().tolist()[0]
        return motifs
    
    def get_Yoshida_motif_names(self):
        return Yoshida.Yoshida_data().get_significant_TF()
        
    def load_raw_score_matrix(self, cluster_index, null=False):
      tb = pd.read_csv(self.get_raw_matrix_filename(cluster_index, 
                                                      null),
                       sep=",")
      m = tb.to_numpy()
      return m
        
    def load_score_matrix(self, cluster_index, null=False):
      tb = pd.read_csv(self.get_scale_matrix_filename(cluster_index, 
                                                      null),
                       sep=",")
      m = tb.to_numpy()
      return m
  
    def load_genomic_table(self, index):
        tb = pd.read_csv(self.get_genomics_table_filename(index),
                       sep=",")
            
        return tb
    
    def raw_FDR(self, FDR=.01):
        call_list = []
        q = 1 - FDR
        for index in range(self.ncluster):
          print(["calculating calls for cluster", index])
          m = self.load_raw_score_matrix(index, null=False)
          mn = self.load_raw_score_matrix(index, null=True)
          nc = mn.shape[1]
      
          cutoffs = [np.quantile(mn[:,i], q=q) for i in range(nc)]
          calls = [np.mean(m[:,i] > cutoffs[i]) for i in range(nc)]
          call_list.append(calls)
          
        tb = pd.DataFrame(call_list, 
                          columns = self.get_motif_names())
        tb = tb.loc[:,np.max(tb, 0) > .2]
        return (tb)
  
    def load_scores_for_motif(self, motif, 
                              ncluster=None,
                              null=False):
        if ncluster is None:
          ncluster = self.ncluster
        all_motifs = np.array(self.get_motif_names() )
        mindex = np.where(all_motifs == motif)[0]
        m_list =[]
        for i in range(ncluster):
          print(["loading cluster", i])
          cm = self.load_score_matrix(i, null=null)
          cmf = cm[:,mindex].flatten()
          
          m_list.append(cmf[~np.isnan(cmf)].tolist())
          
        return m_list
  
    def heatmap(self, use_Yoshida=True, null=False):
        mean_list = []
        all_TF = self.get_motif_names()
        if use_Yoshida:
            int_TF = self.get_Yoshida_motif_names()
        else:
            int_TF = all_TF  
        int_TF = [s.upper() for s in int_TF]
        all_TF = [s.upper() for s in all_TF]
        
        active_cols = [x in int_TF for x in all_TF]
        active_TF = [x for x in all_TF if x in int_TF]
      
        for i in range(self.ncluster):
            m = self.load_score_matrix(i, null=null)[:,active_cols]
            mean_list.append(np.mean(m, 0) > 2)
            
        m_TF = np.array(mean_list)
        sns.heatmap(m_TF)
        return (active_TF)
                             
            
       
        
        
        
        
        
        