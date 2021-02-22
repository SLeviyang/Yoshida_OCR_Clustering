import os,sys
import pdb
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import igraph as ig
import scipy.stats as stats
import seaborn as sns
from sklearn.cluster import AgglomerativeClustering as AggCl
from scipy.cluster.hierarchy import dendrogram

import workflow
import Yoshida
import utilities

class TF_motif_row_cluster:
    
  def __init__(self, 
               enhancer=True,
               edge_FDR=0.001, idr_FDR=0.01):

    w = workflow.workflow("idr", edge_FDR, idr_FDR)
    rcv = w.get_row_cluster_view()
    self.rcv = rcv
    
    cluster_TF_quantile = []
    cluster_TF_means = []
    for cluster in range(20):
   
    
      print(["computing scores for cluster", cluster])
      # TF table (loci by TF)
      scores = w.load_TF_score_matrix(cluster, null=False)
      nscores = w.load_TF_score_matrix(cluster, null=True)
    
      genomics = rcv.form_cluster_genomics(cluster)
      genomics.index = scores.index
      if enhancer:
        scores = scores[np.abs(genomics["dTSS"]) > 3000]
        nscores = nscores[np.abs(genomics["dTSS"]) > 3000]
    
      nTF = scores.shape[1]
      nloci = scores.shape[0]
      TF = scores.columns
    
      q = []
      for i in range(nTF):
          cquant = np.quantile(nscores.iloc[:,i], q=0.95)
          q.append(np.mean(scores.iloc[:,i] > cquant))
          
      cluster_TF_quantile.append(np.array(q)*(np.array(q) > .15))
      cluster_TF_means.append(np.mean(scores, axis=0))
      
    q_tb = pd.DataFrame(cluster_TF_quantile, columns = TF)
    q_tb = q_tb.loc[:,np.sum(q_tb, axis=0) > 0]
    m_tb = pd.DataFrame(cluster_TF_means, columns = TF)
    
    agg = AggCl(n_clusters=q_tb.shape[1],linkage="complete")
    model = agg.fit(np.transpose(q_tb))
    
    self.clust = model
    
    sns.clustermap(np.transpose(q_tb), col_cluster=False)
    self.q_tb = q_tb
    self.m_tb = m_tb
                
    
class TF_expression_row_cluster:
    
  def __init__(self, cluster,
               edge_FDR=0.001, idr_FDR=0.01):

    w = workflow.workflow("idr", edge_FDR, idr_FDR)
    rcv = w.get_row_cluster_view()
    
    
    # atacseq stuff   
    atacseq = pd.DataFrame(rcv.form_cluster_matrix(cluster), 
                           columns = rcv.get_cell_types())
    atacseq_cell_types = rcv.get_cell_types()    
    
    # RNAseq info (table is genes by cell types)
    r = Yoshida.Yoshida_RNAseq()
    rnaseq_tb = r.load_quantiled()
    genes = rnaseq_tb["gene"]
    rnaseq = rnaseq_tb.drop("gene", axis=1)
    rnaseq.index = genes
    rnaseq_cell_types = rnaseq.columns
    
    # restrict RNAseq to TF
    scores = w.load_TF_score_matrix(cluster, null=False)
    TF = scores.columns.to_list()
    rnaseq = rnaseq.iloc[[g in TF for g in genes],:]
    
    
    
    
    
    
    # restrict rnaseq and atacseq to joint cell types
    joint_cell_types = list(set(atacseq_cell_types).intersection(rnaseq_cell_types))
        
    # put everything in terms of Yoshida tree
    g = Yoshida.Yoshida_tree().load_igraph()
    active_v = [i for i,ct in enumerate(g.vs["name"]) if ct in joint_cell_types]
    self.g =  g.induced_subgraph(active_v)
        
    self.rnaseq = rnaseq.loc[:,self.g.vs["name"]]
    self.atacseq = atacseq.loc[:,self.g.vs["name"]]
    self.OCR = np.mean(self.atacseq.to_numpy(), axis=0)
    self.TF = rnaseq.index
    self.rcv = rcv
    
    cc = []
    r = self.rnaseq.to_numpy()
    TF = self.TF
    for i in range(len(TF)):
        cc.append(np.corrcoef(self.OCR, r[i,:])[0,1])
        
    cc = np.array(cc)
    q = pd.DataFrame({'TF':TF,
                      'corr':np.abs(cc),
                      'dir':cc/np.abs(cc)})
    self.q = q.sort_values("corr", ascending=False)
    
    #     # plot the mean OCR of cluster (TF is None) or
#     # the expression of a TF (TF is not None)
  def treeplot(self, TF = None):
        if TF is not None:
          vals = 2*self.rnaseq.loc[TF,:].to_numpy()
        else:
          OCR = self.OCR
          scaled_OCR = (OCR - np.min(OCR))/(np.max(OCR) - np.min(OCR))
          vals = 25*scaled_OCR
        
        g = self.g
        
        vs = {}
        vs["bbox"] = (1200, 1000)
        vs["vertex_size"] = vals
        vs["vertex_label_size"] = 20
        vs["vertex_label"] = [str(i)  for i in range(g.vcount())]
        vs["vertex_label_dist"] = 1.5
           
        layout = g.layout_reingold_tilford(mode="all")
     
        pl = ig.plot(g, layout=layout, **vs)
        pl.show()
  
    
    
  
    
          
          
       
    
    


# Analyze the distribution of motif scores for a collection
# of sequences
# @param w a workflow object

# @param
# class TF_cluster:
    
#     nucs = ["A", "C", "G", "T", "GC", "CG"]
     
#     def __init__(self, cluster,
#                  edge_FDR=0.001, idr_FDR=0.01):
        
#         w = workflow.workflow("idr", edge_FDR, idr_FDR)
        
#         # TF table (loci by TF)
#         self.scores = w.load_TF_score_matrix(cluster, null=False)
#         self.nscores = w.load_TF_score_matrix(cluster, null=True)
        
#         # RNAseq info (table is genes by cell types)
#         r = Yoshida.Yoshida_RNAseq()
#         rnaseq_tb = r.load_quantiled()
#         genes = rnaseq_tb["gene"]
#         self.rnaseq = rnaseq_tb.drop("gene", axis=1)
#         self.rnaseq.index = genes
#         rnaseq_cell_types = self.rnaseq.columns
        
#         # ATACseq info (table is loci by cell types)
#         self.atacseq = pd.DataFrame(w.load_m_list()[cluster],
#                                     columns = w.get_cell_types()) 

#         # restrict atacseq results to cell types in rnaseq
#         # restrict rnaseq and scores to shared TF/genes
#         self.atacseq = self.atacseq.loc[:,rnaseq_cell_types]
      
#         joint = list(set(self.scores.columns).intersection(self.rnaseq.index))
#         self.scores = self.scores.loc[:,joint]
#         self.nscores = self.nscores.loc[:,joint]
#         self.rnaseq = self.rnaseq.loc[joint,:]
        
#         # put everything in terms of Yoshida tree
#         g = Yoshida.Yoshida_tree().load_igraph()
#         active_v = [i for i,ct in enumerate(g.vs["name"]) if ct in rnaseq_cell_types]
#         self.g =  g.induced_subgraph(active_v)
        
#         self.rnaseq = self.rnaseq.loc[:,self.g.vs["name"]]
#         self.atacseq = self.atacseq.loc[:,self.g.vs["name"]]
                                             
#         self.calculate_TF_quantiles()
#         self.calculate_OCR_accessibility()
        
#         plt.scatter(self.acc,self.TF_p)
        
#     def calculate_TF_quantiles(self):
#         s = self.scores.to_numpy()
#         ns = self.nscores.to_numpy()
        
#         ns_q = np.quantile(ns, q=0.95, axis=0)
#         TF_p = [np.mean(s[:,i] > ns_q[i]) for i in range(s.shape[1])]
        
#         self.TF_p = np.array(TF_p)
        
#     def calculate_OCR_accessibility(self):
#         r = self.rnaseq.to_numpy()
#         a = np.mean(self.atacseq.to_numpy(), axis=0)
#         acc = [np.corrcoef(a, r[i,:])[0,1] for i in range(r.shape[0])]
        
#         self.acc = acc
        
#     def get_cell_types(self):
#         return self.g.vs["name"]
        
#     # plot the mean OCR of cluster (TF is None) or
#     # the expression of a TF (TF is not None)
#     def treeplot(self, TF=None):
#         if TF is not None:
#             vals = 2*self.rnaseq.loc[TF,:].to_numpy()
#         else:
#             vals = 30*np.mean(self.atacseq.to_numpy(), axis=0)
        
#         print(np.max(vals))
#         g = self.g
        
#         vs = {}
#         vs["bbox"] = (1200, 1000)
#         vs["vertex_size"] = vals
#         vs["vertex_label_size"] = 20
#         vs["vertex_label"] = [str(i)  for i in range(g.vcount())]
#         vs["vertex_label_dist"] = 1.5
           
#         layout = g.layout_reingold_tilford(mode="all")
     
#         pl = ig.plot(g, layout=layout, **vs)
#         pl.show()
        
        
   