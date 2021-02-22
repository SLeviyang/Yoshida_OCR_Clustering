import configuration as conf
import Yoshida
import utilities

import numpy as np
import pandas as pd
import pdb
import os
from sknetwork.clustering import Louvain
import igraph as ig
from scipy.sparse import csr_matrix
from scipy.stats import hypergeom
import matplotlib.pyplot as plt
import itertools as it
import seaborn as sns

# Create row clusters from a binary master
#
# @param m a binary matrix
# @param FDR the FDR with which edges are called between rows
# @param homo_cutoff only consider rows with more than this many 1's
# and 0's.  Rows with < homo_cutoff 1's are ignored.  Rows with >
# m.shape[1] - homo_cutoff 1's are grouped together in a single 
# cluster
class row_cluster: 
   
    def __init__(self, m, FDR, homo_cutoff=2):
        self.np = null_row_cluster()
        # binary matrix
        self.m = m.astype("float")
        self.FDR = FDR
        
        self.edges = self.compute_edges(FDR=FDR)
        self.clusters = self.compute_clusters(self.edges, 
                                              homo_cutoff=homo_cutoff)
        
    # computation methods
    def compute_edges_between_blocks(self, n1, n2, FDR):
        
        print(["starting computation for edges", n1, n2])
  
        m = self.m
        nc = m.shape[1]
        
        ind1 = (np.where(np.sum(m, 1) == n1))[0]
        mm1 = m[ind1,:]
        N1 = len(mm1)
        
        ind2 = (np.where(np.sum(m, 1) == n2))[0]
        mm2 = m[ind2,:]
        N2 = len(mm2)
        
        dist_m = np.dot(mm1, np.transpose(1 - mm2)) \
                 + np.dot(1-mm1, np.transpose(mm2))
        if n1 == n2:
          np.fill_diagonal(dist_m, nc)
    
        # get the sampled distance counts 
        dc = np.zeros([N1, nc+1])
        for k in range(N1):
          cc = np.unique(dist_m[k,:], return_counts=True)
          ind = np.int32(cc[0])
          dc[k,ind] = dc[k,ind] + cc[1]
           
        null_dist_pmf = null_row_cluster().load_distribution(n2)[n1,:]
        null_dc = N2*null_dist_pmf
        null_cumsum = np.cumsum(null_dc)

        cutoffs = np.zeros(N1)
        print("identifying significant edges...")
        for i in range(N1):
            sample_cumsum = np.cumsum(dc[i,:])
            pass_FDR = np.where((null_cumsum <= sample_cumsum*FDR) &
                                (sample_cumsum > 0))[0]
            if len(pass_FDR) == 0:
                cutoffs[i] = -1
            else:
                cutoffs[i] = np.max(pass_FDR)
                
        e1 = []
        e2 = []
        print("gathering significant edges...")
        for i in range(N1):
            ce2 = np.where(dist_m[i,:] <= cutoffs[i])[0]
            if not len(ce2) == 0:
                e1 = e1 + np.repeat(i, len(ce2)).tolist()
                e2 = e2 + ce2.tolist()
                
        print(["found", len(e1), "edges between blocks",
                            n1, n2])  
            
       
        
        if len(e1) == 0:
            return None
            
        edge_d = pd.DataFrame({'e1':ind1[e1], 'e2':ind2[e2]})
        if  n1 == n2:
           edge_d = edge_d.loc[edge_d["e1"] < edge_d["e2"]]
       
        return edge_d
    
    
    def compute_edges(self, FDR):
        nc = self.m.shape[1]
        
        #d_cutoff = null_peak_clusters().read_cutoff_file()
        #d_cutoff = d_cutoff[d_cutoff["rs1"] <= d_cutoff["rs2"]]
        
        edges = []
        for n1 in range(3,nc-3):
            for delta in [0,1]:
                if n1+delta <= nc-2:
                  n2 = n1 + delta
                  ce = self.compute_edges_between_blocks(n1, n2, FDR)
                  if not ce is None:
                     print(["found", len(ce), "edges between blocks",
                            n1, n2])
                     edges.append(ce)
                  else:
                     print(["found", 0, "edges between blocks",
                            n1, n2]) 
        
        all_edges = pd.concat(edges)
        return all_edges
       
    
    def louvain(self, edges):
        
        unique_edges = set(edges.iloc[:,0].to_list() + 
                           edges.iloc[:,1].to_list())
        unique_edges = sorted(unique_edges)
        nv = len(unique_edges)
        ud = pd.DataFrame({'e1':unique_edges})
        ud = ud.reset_index(drop=False)
        
        e1_list = edges.merge(ud, how="left")["index"].to_list()
        ud = ud.rename({'e1':'e2'}, axis=1)
        e2_list = edges.merge(ud, how="left")["index"].to_list()
        
        values = np.repeat(1, len(edges))
       
        m = csr_matrix((values, (e1_list, e2_list)), 
                       shape=(nv, nv))
        
        louvain = Louvain(verbose=True)
        labels = louvain.fit_transform(m)
        
        cluster_d = pd.DataFrame({'row':unique_edges,
                          'cluster':labels}).sort_values(["cluster","row"])
        
        return cluster_d
        
    
    def compute_clusters(self, edges, homo_cutoff=2):
        
        cluster_d = self.louvain(edges) 
    
        # get rows that are mostly open
        m = self.m
        homo_open = np.where(np.sum(m, 1) >= (m.shape[1] - homo_cutoff))[0]
        d_home = pd.DataFrame({'row':homo_open,'cluster':-1})
        
        cluster_d.loc[~cluster_d["row"].isin(homo_open)]
        
        joint = d_home.append(cluster_d)
        joint["cluster"] = joint["cluster"] + 1
        joint = joint.sort_values("cluster")
          
        return joint
        



            
# Constructs null distribution for row clustering
# dist(k, n, N) = probability of k overlaps between two rows with n and
# N 1's, respectively.
#
# This class computes the value dist(k,n,N) over all k,n,N from 0,1,2,...,78
# The two-d array for dist(-,n,-) is saved to a file for each value of n.
class null_row_cluster:
    null_peak_cluster_directory = conf.DATA_DIR + "null_row_cluster/"
   
    def __init__(self, max_val=78):
      
        self.max_val = max_val
        
        if not os.path.isdir(self.null_peak_cluster_directory):
            os.mkdir(self.null_peak_cluster_directory)   
            self.create_distance_distribution_matrices()
        
    # distribution of hamming distance between two rows with n and N 1's,
    # respectively
    #    
    # hypergeometric(M, n, N) gives number of 1's that agree
    # notation follows python hypergeom parameters
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.hypergeom.html
    # M = number of cell types
    # n = number of 1's in a row
    # N = number of 1's in other row
    # k = number of shared 1's between the two rows
    # hamming distance = (n - k) + (N-k)
    def compute_distance_distribution(self, M, n, N):
        allk = range(M+1)
        pmf = hypergeom.pmf(allk, M, n, N)
        
        dist = np.zeros(M+1)
        for k in range(np.min([M, n+N])):
            hamming_dist = np.int32((n-k) + (N-k))
            if hamming_dist > M:
                hamming_dist = M
            if hamming_dist < 0:
                hamming_dist = 0
            dist[hamming_dist] = dist[hamming_dist] + pmf[k]
            
        return dist
    
    # compute dist(-,n,-) over all n
    def compute_distance_distribution_matrix(self, N):
        M = self.max_val
        distm = np.zeros([M,M+1])
        for n in range(M):
            distm[n,:] = self.compute_distance_distribution(M, n, N)
            
        return distm
            
    def create_distance_distribution_matrices(self):
        for N in range(self.max_val+1):
            print(["creating distance matrix", N])
            d = self.compute_distance_distribution_matrix(N)
          
            cls = ["dist" + str(i) for i in range(d.shape[1])]
            tb = pd.DataFrame(d, columns=cls)
            
            outf = self.null_peak_cluster_directory \
                        + "distance_distribution_" + str(N) + ".csv"
            tb.to_csv(outf, sep=",", index=False)
            
    def load_distribution(self, N):
         outf = self.null_peak_cluster_directory \
                        + "distance_distribution_" + str(N) + ".csv"
         return pd.read_csv(outf, sep=",").to_numpy()
    
            
    
    
    
  