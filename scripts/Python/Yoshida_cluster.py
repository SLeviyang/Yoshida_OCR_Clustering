import configuration as conf
import utilities
import tree_cluster as tc
import peak_clusters as pk
import Yoshida as yosh

import pandas as pd
import os
import igraph as ig
import numpy as np
import pdb

   
class Yoshida_cluster:
    
    def  __init__(self, min_peak_cluster_size=100):
        
        g = yosh.Yoshida_tree(just_root=True).g
        
        p = pk.peak_clusters()
        cluster_sizes = p.get_cluster_sizes(min_peak_cluster_size)
        uclusters = cluster_sizes["cluster"]
        m_list = [p.form_cluster_matrix(ii) for ii in uclusters]
    
        column_names = p.get_matrix_column_names()
        
        shared = set(column_names).intersection(g.vs["name"])
        gsub = g.induced_subgraph(shared)
        
        ind = utilities.match(gsub.vs["name"], column_names)
        
        new_m = [mm[:,ind] for mm in m_list]
        new_g = gsub
        
        self.g = new_g
        self.m = new_m
        
    def create_tree_cluster(self, K):
        return tc.tree_cluster(self.g, self.m, K=K)
        
        self.gc = tc.tree_cluster(self.g)
        
    def cluster_variance(self, assignments):
        m = self.m
        ncell_types = m[0].shape[1]
        K = len(np.unique(assignments))
        a = assignments
        
        var = np.zeros(len(m))
        for j in range(len(m)):
          row_clust_size = m[j].shape[0]*m[j].shape[1]
          for k in range(K):
            cell_type_means = np.mean(m[j][:,a==k], axis=0)
            var[j] = var[j] + np.var(cell_type_means)
          var[j] = var[j]/K
            
        return {'row_cluster_var':var,
                'total_var':np.mean(var)}
        
        
    def fitK(self, K, num_runs=3):
        
        gc = self.create_tree_cluster(K=K)
        gc.fit(num_runs)
        
        a = gc.get_assignments()
        pa = self.cluster_variance(a)
        
        return pa
    
    # determine fit if each cell type is a separate cluster
    def fitKinf(self):
        m = self.m
        ncell_types = m[0].shape[1]
        a = np.arange(ncell_types)
        pa = self.cluster_variance(a)
        
        return pa
    
      # determine fit if each cell type is a separate cluster
    def fitK1(self):
        m = self.m
        ncell_types = m[0].shape[1]
        a = np.zeros(ncell_types)
        pa = self.cluster_variance(a)
        
        return pa
        
    def fit(self, Kmin=2, Kmax=6, num_runs=3):
        
        Kvals = np.arange(Kmin, Kmax+1)
        pa = []
        
        for i in range(len(Kvals)):
            f = self.fitK(K=Kvals[i], num_runs=num_runs)
            pa.append(f["total_var"])           
        
        maxvar = self.fitK1()["total_var"]
        tb = pd.DataFrame({'K':Kvals,
                            'var':np.array(pa)/maxvar})
            
        return tb
    
    
    

