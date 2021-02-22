import pdb
import numpy as np
import igraph as ig
import sys
import pandas as pd
import time
from sklearn.cluster import AgglomerativeClustering

'''
cell_cluster class:
    @params:
        K = # clusters
        m_list = list of row cluster dataframes
    initialize_assignments
        pass in assignments or default to hierarchical clustering
    optimize
        K-means algorithm
            initialize assignments & mediods
            assign each cell type to closest mediod
            recompute mediods
            repeat till next position has equal residual or 1000 iterations
'''

class cell_cluster:
    
    def __init__(self, K, m_list):
        self.K = K
        self.m = m_list
        
        ncols = list(set([m.shape[1] for m in m_list]))
        if not len(ncols) == 1:
            sys.exit("matrices in m_list must have same number of columns")
       
        # checks on K
        if K < 2 or K > ncols[0]:
            sys.exit("K out of range")
        
        self.assignments = None
        self.mediods = None 
        
    def initialize_assignments(self, assignments = []):
        if len(assignments)>0:
            if not self.K == len(set(assignments)):
                pdb.set_trace()
                sys.exit("number of components does not equal number of clusters") 
            self.assignments = assignments
        else:
            model = AgglomerativeClustering(n_clusters=self.K,linkage='complete')
            df_list = [pd.DataFrame(m) for m in self.m]
            M = pd.concat(df_list)
            clustering = model.fit(M.T)
            self.assignments = clustering.labels_
        
        self.update_mediods()
                
    def update_mediods(self):
        K = self.K
        m = self.m
        r = range(len(m))
        a = self.assignments
        
        cluster_means = [[np.mean(m[i][:,a==k]) for k in range(K)] for i in r]
        self.mediods = np.array(cluster_means)
        
    def compute_residual2(self):
        K = self.K
        mediods = self.mediods 
        a = self.assignments
        m = self.m
        
        lm = len(m)
        col_ind = [np.where(a == k)[0] for k in range(self.K)]
        ss2 = 0
        
        for j in range(lm):
          for k in range(K):
            res2 = (m[j][:,col_ind[k]] - mediods[j,k])**2
            ss2 = ss2 + np.sum(res2)
            
        return ss2
                    
    def assignCelltoMediod(self):
        mediods = self.mediods 
        a = self.assignments
        m = self.m
        k = self.K
        
        assignments = np.empty(len(a),dtype=int)
        for cell in range(len(a)):
            res = [0]*k
            for j in range(len(m)):
                for cluster in range(k):
                    res[cluster] += np.sum((m[j][:,cell] - mediods[j,cluster])**2)
            assignments[cell] = np.argmin(res)
        self.assignments = assignments    
    
    # single starting point optimization
    def optimize(self, assignments = None, max_iter = 1000):
        start = time.time()
        
        print("beginning cell cluster optimization")
        if not assignments is None:
            self.initialize_components(assignments = assignments)
            
        if self.assignments is None:
            sys.exit("assignments must be passed or clustering initialized.")
            
        previous_loss = self.compute_residual2()
        iteration = 1
        while iteration<=max_iter:
           print(["epoch", iteration, previous_loss])
           self.assignCelltoMediod()
           self.update_mediods()
           current_loss = self.compute_residual2()
           
           if current_loss > previous_loss:
                pdb.set_trace()
                sys.exit("unexpected state!")
           elif current_loss < previous_loss:
                previous_loss = current_loss
                iteration = iteration + 1
           else:
                break
            
        end = time.time()
        print(["optimization time", end - start])
    
    def get_assignments(self):
        return self.assignments
    
    def get_mediods(self):
        return self.mediods
    
    def treeplot(self,g,savepath='',
                 vertex_label=False, 
                 m_index=None):
        a = self.assignments
        vs = {}
        
        if a is not None:
            pal = ig.drawing.colors.ClusterColoringPalette(self.K)
            vs["vertex_color"] = pal.get_many(a)
                  
        vs["bbox"] = (1200, 1000)
        if m_index is None:
            vs["vertex_size"] = 20
        else:
            vs["vertex_size"] = 30*np.mean(self.m[m_index], 0) + 1
        vs["vertex_label_size"] = 20
        if vertex_label:
            vs["vertex_label"] = g.vs['name']
        else:
            vs["vertex_label"] = [str(i)  for i in range(g.vcount())]
        vs["vertex_label_dist"] = 1.5
           
        layout = g.layout_reingold_tilford(mode="all")
     
        if savepath == '':
            pl = ig.plot(g, layout=layout, **vs)
            pl.show()
        else:
            ig.plot(g,savepath,layout=layout, **vs)