import configuration as conf
import os
import pdb
import pandas as pd
import numpy as np
import glmnet_python 

import sklearn.linear_model as sklin




class TF:
    
    motif_dir = conf.DATA_DIR + "TF_motifs/"
    
    def __init__(self, max_cluster_number=14):
      m_list = []
      null_m_list = []
      for i in range(max_cluster_number+1):
        file = self.motif_dir + "motif_scores_cluster_" + str(i) + ".csv"
        nfile = self.motif_dir + "motif_scores_cluster_" + str(i) + "_permuted.csv"
         
        m = pd.read_csv(file, sep=",")
        nm = pd.read_csv(nfile, sep=",")
        
        if i == 0:
            motif_names = m.columns
        m_list.append(m.to_numpy())
        null_m_list.append(nm.to_numpy())
        
      self.motif_names = motif_names
      self.m_list = m_list
      self.null_m_list = null_m_list
      
    def get_motif_names(self):
        return self.motif_names
        
    def load_score_matrix(self, index, null=False):
      if null:
          m = self.null_m_list[index]
      else:
          m = self.m_list[index]
          
      return m
  
    def find_within_cluster_enriched(self, index1, index2):
        m = self.load_score_matrix(index1, null=False)
        nm = self.load_score_matrix(index2, null=False)
        motif_names = self.get_motif_names()
        
        y = np.concatenate((np.repeat(1, m.shape[0]), 
                         np.repeat(0, nm.shape[0])))
        X = np.concatenate((m, nm), axis=0)
        
        #X = X[0:100,1:50]
        #y = np.concatenate((np.repeat(1, 50), 
         #                np.repeat(0, 50)))*1.0
      
        f =  sklin.LogisticRegressionCV(cv=5, 
                                        random_state=0,
                                        max_iter=1000,
                                        penalty="l1")
        f.fit(X, y)
        
        pdb.set_trace()
       
        
        
        
        
        
        