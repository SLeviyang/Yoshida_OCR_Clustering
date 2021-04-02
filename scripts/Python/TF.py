import os,sys
import pdb
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import igraph as ig
import scipy.stats as stats
import seaborn as sns
import copy
from sklearn.cluster import AgglomerativeClustering as AggCl
from scipy.cluster.hierarchy import dendrogram
import scipy.stats as stats
import sklearn.preprocessing as skpre
from sklearn.cluster import KMeans

from scipy.sparse import csr_matrix
from sknetwork.clustering import Louvain


import workflow
import Yoshida
import utilities
import configuration as conf


class TF_ChIP:
    
  GEO_table = conf.INPUT_DATA + "ChIP_GEO.csv"
  TF_ChIP_dir = conf.DATA_DIR + "ChIP/"
  aws_peak_dir = "yoshida-chipseq/peaks/"
  
  def __init__(self, w):
    self.w = w
    if not os.path.isdir(self.TF_ChIP_dir):
      os.mkdir(self.TF_ChIP_dir)
      
    tb = self.load_ChIP_table()
    for i in range(len(tb)):
      GEO = tb["GEO"][i]
      bed_file = self.get_bedfile(i)
      if not os.path.isfile(bed_file):
        os.environ["bf"] = bed_file
        os.environ["if"] = "s3://" + self.aws_peak_dir + GEO + ".bed" 
        os.system("/usr/local/bin/aws s3 cp $if $bf")
      
      
  def load_ChIP_table(self):
      tb = pd.read_csv(self.GEO_table)
      return tb
  
  def get_bedfile(self, i):
    tb = self.load_ChIP_table()
    
    GEO = tb["GEO"][i]
    cell_type = tb["cell_type"][i]
    authors = tb["authors"][i]
    
    bed_file = self.TF_ChIP_dir + authors + "_" + cell_type + "_" + GEO + ".bed"
    return bed_file

  def load_bed(self, i):
      return pd.read_csv(self.get_bedfile(i),
                         sep="\t",
                         header=None,
                         names=["chr", "chrStart", "chrEnd",
                                  "name", "score", "strand", 
                                  "signalValue", "p", "q", "peak"])
  
  def intersect_bed(self, i,
                    enhancer=True,
                    min_q=10):
      
      m = self.w.load_matrix()
      msum = np.sum(m, 1)
      
      # restrict TF bed file to summits
      TFbedfile = "temp_TF_" + str(np.random.uniform()) + ".bed"
      TF_b = self.load_bed(i)
      summit = TF_b["chrStart"] + TF_b["peak"]
      in_b = pd.DataFrame({"chr":TF_b["chr"],
                           "chrStart":summit,
                           "chrEnd":summit,
                           "q":TF_b["q"]})
      in_b.to_csv(TFbedfile, sep="\t", header=None, index=False)
      
      
      sortfile = "temp_sort_" + str(np.random.uniform()) + ".bed"
      outfile = "temp_" + str(np.random.uniform()) + ".bed"
      os.environ["sf"] = sortfile
      os.environ["TFbf"] = TFbedfile
      os.environ["mbf"] = self.w.bed_file
      os.environ["btools"] = conf.BEDTOOLS_PATH
      os.environ["of"] = outfile
      
      os.system('$btools sort -i $TFbf > $sf')
      os.system('$btools closest -d -t first -a $mbf -b $sf > $of') 
      b = pd.read_csv(outfile, header=None, sep="\t")
   
      b = b.iloc[:,[0,1,2,6,7]]
      b.columns = ["chr", "chrStart", "chrEnd", "-log10q", "dist"]
     
      # add cluster and dTSS info
      clusters =self. w.load_clusters()
      dTSS = self.w.load_genomics()["dTSS"]
      peak_clusters = np.zeros(len(b)) - 1
      peak_clusters[clusters["row"]] = clusters["cluster"]
      b.insert(b.shape[1], "dTSS", dTSS)
      b.insert(b.shape[1], "cluster", peak_clusters)
      b.insert(b.shape[1], "row", range(len(b)))
      b.insert(b.shape[1], "nopen", msum)
      
      os.remove(sortfile)
      os.remove(outfile)
      os.remove(TFbedfile)
      
      b = b[(b["dist"] == 0) & (b["-log10q"] > min_q)]
      if enhancer:
          b = b[(b["dTSS"] > 3000) & (b["dTSS"] < 1E5)]
      else:
          b = b[b["dTSS"] < 500]
      return b
  
  
      
      
      
   
# used in the first portion of the TF section of the manuscript
class TF_motif:
    
  z_prefix = conf.DATA_DIR + "Yoshida/chromVar_z_"
  cell_type_prefix = conf.DATA_DIR + "Yoshida/chromVar_cell_type"
  
  # M, X, E, B are in the notation of chromVar online methods
  def __init__(self,  w,  
               enhancer=True, 
               just_clusters=True,
               k=5, nperm=50):
        
    self.w = w
    self.rcv = w.get_row_cluster_view()
    self.X = w.load_matrix()
    self.M = w.load_TF_matrix()
    
    self.sigTF = Yoshida.Yoshida_TF().get_TF(significant=True)
    # restrict to significant TF
    self.bed = w.load_bed()
    self.cell_types = w.get_cell_types()
    
    genomics = w.load_genomics()
    if enhancer:
        genomics = genomics[genomics["dTSS"] >= 3000]
    else:
        genomics = genomics[genomics["dTSS"] <= 500]
    g_index = genomics.index

    if just_clusters:
      nclust = len(w.load_m_list())
      clusters = w.load_clusters()
      c_index = clusters[clusters["cluster"] < nclust]["row"]
    else:
      c_index = np.arange(len(self.X))
      
    self.peak_clusters = np.zeros(len(self.X)) - 1
    self.peak_clusters[clusters["row"]] = clusters["cluster"]
      
    active_index = list(set(c_index).intersection(g_index))
    # keep only enhancer or promoter peaks
    self.genomics = w.load_genomics().iloc[active_index,:]
    self.X = self.X[active_index,:]
    self.M = self.M.iloc[active_index,:]
    self.bed = self.bed.iloc[active_index,:]
    self.peak_clusters = self.peak_clusters[active_index]
    
    # remove peaks that are all 0's as that crashes chromVar
    active_peaks = np.sum(self.X, 1) != 0
    self.genomics = self.genomics.loc[active_index,:]
    self.X = self.X[active_peaks,]
    self.M = self.M.loc[active_peaks,]
    self.bed = self.bed.loc[active_peaks,]
    self.peak_clusters = self.peak_clusters[active_peaks]
    
    if not os.path.isfile(self.z_prefix + ".csv"):
      self.generate_chromVar_z_cluster_file(k=12, nperm=nperm)
    if not os.path.isfile(self.cell_type_prefix + ".csv"):
      self.generate_chromVar_cell_type_file(nperm=nperm)
   
    self.z_across = pd.read_csv(self.z_prefix  + ".csv",
                                sep=",")
    self.znull_across = pd.read_csv(self.z_prefix  + "_null.csv",
                                sep=",")
    self.r_frac = pd.read_csv(self.z_prefix  + "_fraction.csv",
                                sep=",")
    
    self.ct = pd.read_csv(self.cell_type_prefix + ".csv",
                                sep=",")
    self.ct.index = self.cell_types
    self.ct_null = pd.read_csv(self.cell_type_prefix + "_null.csv",
                                sep=",")
    self.ct_null.index = self.cell_types
    self.ct_frac = pd.read_csv(self.cell_type_prefix  + "_fraction.csv",
                                sep=",")
    self.ct_frac.index = self.cell_types
  
    return None

  # constructor methods
 
  def generate_chromVar_cell_type_file(self,
                                       nperm=50):

    X = self.X
    M = self.M.to_numpy()
      
    ct = np.array(self.cell_types)
    z = []
    r_f = []
    znull_max = []

    for i in range(len(ct)):
    
      print(["computing z values for cell type", ct[i]])
      response_cell_types = (ct == ct[i])
      control_cell_types = (ct != ct[i])
    
      # response loci are the loci that are open in the cell type
      response_peaks = self.X[:,i] == 1
      control_peaks = self.X[:,i] == 0
    
      zinfo = self.compute_chromVar_z(M, X,
                                      response_cell_types,
                                      control_cell_types,
                                      response_peaks,
                                      control_peaks,
                                      nperm=nperm)
      z.append(zinfo["z"])
      r_f.append(zinfo["response_frac"])
      znull_max.append(zinfo["null_max"])
      
    z = pd.DataFrame(z, columns=self.M.columns)
    znull = pd.DataFrame(znull_max, columns=self.M.columns)
    r_f = pd.DataFrame(r_f, columns=self.M.columns)
    

    z.to_csv(self.cell_type_prefix + ".csv", sep=",", index=False)
    znull.to_csv(self.cell_type_prefix + "_null.csv", sep=",", 
                 index=False)
    r_f.to_csv(self.cell_type_prefix + "_fraction.csv", sep=",", 
                 index=False)

  def generate_chromVar_z_cluster_file(self, 
                                       k,
                                       nperm=50):

    X = self.X
    M = self.M.to_numpy()
      
    a = self.w.biclusters_to_response_control(k=k)
    z = []
    r_f = []
    znull_max = []

    for i in range(len(a)):
    
      print(["computing z values for cluster", i])
      response_cell_types = (a[i,:] == 1)
      control_cell_types = (a[i,:] == 0)
    
    
      response_peaks = self.peak_clusters == i
      control_peaks = self.peak_clusters != i
    
      zinfo = self.compute_chromVar_z(M, X,
                                      response_cell_types,
                                      control_cell_types,
                                      response_peaks,
                                      control_peaks,
                                      nperm=nperm)
      z.append(zinfo["z"])
      r_f.append(zinfo["response_frac"])
      znull_max.append(zinfo["null_max"])
      
    z = pd.DataFrame(z, columns=self.M.columns)
    znull = pd.DataFrame(znull_max, columns=self.M.columns)
    r_f = pd.DataFrame(r_f, columns=self.M.columns)
    
    z.to_csv(self.z_prefix +  ".csv", sep=",", index=False)
    znull.to_csv(self.z_prefix + "_null.csv", sep=",", 
                 index=False)
    r_f.to_csv(self.z_prefix + "_fraction.csv", sep=",", 
                 index=False)
   
    
 
    
        
  # see online methods of chromVar, Schep 2017 Nature Methods,
  def compute_chromVar_z(self, 
                         M, X,
                         response_cell_types,
                         control_cell_types,
                         response_peaks,
                         control_peaks,
                         nperm=50):  
  
    Xr = np.sum(X[response_peaks,:][:,response_cell_types], axis=1)
    Mr = M[response_peaks,:]
    
    Xc = np.sum(X[control_peaks,:][:,control_cell_types], axis=1)
    Mc = M[control_peaks,:]   
    nTF = self.M.shape[1]
    
    # fraction of open loci that are open and have motif
    r_N_open_and_motif = np.transpose(Mr) @ Xr
    r_N_open = np.sum(Xr)
    r_f = r_N_open_and_motif/r_N_open
    
    c_N_open_and_motif = np.transpose(Mc) @ Xc
    c_N_open = np.sum(Xc)
    c_f = c_N_open_and_motif/c_N_open
    
    t_f = (r_N_open_and_motif + c_N_open_and_motif)/(r_N_open + c_N_open)
    
    z = (r_f - c_f)/t_f
    
    zn_list = []
    r_list = []
    c_list = []
    
    M_joint = np.concatenate([Mr, Mc], axis=0)
    X_joint = np.concatenate([Xr, Xc])
   
    for i in range(nperm):
        print(["permutation", i])
        rowp = np.random.permutation(len(M))
        M_null = M_joint[rowp,:]
        X_null = X_joint[rowp]
        
        Xrn = X_null[response_peaks]
        Mrn = M_null[response_peaks,:]
    
        Xcn = X_null[control_peaks]
        Mcn = M_null[control_peaks,:]   
        
        r = (np.transpose(Mrn) @ Xrn)/np.sum(Xrn)
        c = (np.transpose(Mcn) @ Xcn)/np.sum(Xcn)
       
        r_list.append(r)
        c_list.append(z)
        
        zn_list.append((r-c)/t_f) 
        
    zmax = np.max(zn_list, axis=0)
    
    return {'z':z,
            'response_frac':r_f,
            'control_frac':c_f,
            'total_frac':t_f,
            'null_max':zmax}

  # methods to analyze enrichment over all TF using chromVar
  # methods
  def enrichment_matrix(self, 
                        sig=True,
                        cutoff=0.25):
    d = self.score_df(sig=sig, cell_types=False)
    d = d[(d["across"] > cutoff) & d["across_sig"]]
    
    TF = sorted(d["TF"].unique().tolist())
    nclust = np.max(d["cluster"] + 1)
    m = pd.DataFrame(np.zeros([nclust, len(TF)]),
                     columns=TF)
    for i in range(len(d)):
        cluster = d["cluster"].iloc[i]
        cTF = d["TF"].iloc[i]
        m.loc[m.index[cluster], cTF] = d["across"].iloc[i]
        
    return m

  def score_df(self, sig=True,
                  cell_types=False):
      
      if sig:
        TF = self.sigTF
      else:
        TF = self.M.columns
        
      if cell_types:
        z_across = self.ct
        z_null = self.ct_null
        frac = self.ct_frac
        factors = self.cell_types
      else:
        z_across = self.z_across
        z_null = self.znull_across
        frac = self.r_frac
        factors = np.arange(len(z_across))
      
      za = np.array(z_across.loc[:,TF])
      za_n = np.array(z_null.loc[:,TF])
      
      za_sig = za > za_n
      #zw_sig = zw > zw_n
     
      f = np.array(frac.loc[:,TF])
    
      d_list = []
      for i in range(len(za)):
          d_list.append(pd.DataFrame({'TF':TF,
                         'cluster':np.repeat(factors[i], len(TF)),
                        'across':za[i,:],
                        'across_sig':za_sig[i,:],
                       # 'within':zw[i,:],
                       # 'within_sig':zw_sig[i,:],
                        'frac':f[i,:]}))
           
      d = pd.concat(d_list, axis=0)
  
      return d
  
  # get methods
  
  def get_TF(self, sig=False):
    if sig:
      return self.sigTF
    else:
      return self.M.columns
      
  
      
# Handles RNAseq over TF
# @param w workflow object
# @param significant should only the significant TF used in Yoshida
# be considered, or all TF used in Yoshida?
class TF_RNAseq:
    
  DE_file = conf.DATA_DIR + "Yoshida/Yoshida_TF_RNAseq_DE.csv"
  
  def __init__(self, w):

    # RNAseq info (table is genes by cell types)
    rnaseq = Yoshida.Yoshida_RNAseq().load_quantiled_RNAseq()
    rnaseq.index = rnaseq["gene"]
    rnaseq = rnaseq.drop("gene", axis=1)
    
    # ATACseq
    m_list = w.load_m_list()
    atacseq = pd.DataFrame(np.array([np.mean(m, axis=0) for m in m_list]),
                           columns=w.get_cell_types())
    
    # limit to TF
    TF = Yoshida.Yoshida_TF().get_TF(significant=False)
    
    # limit and order by Yoshida tree cell tyeps
    g = Yoshida.Yoshida_tree().load_igraph()
    all_ct = g.vs["name"]
    r_ct = rnaseq.columns
    a_ct = atacseq.columns
    
    # joint genes, cell types
    joint_ct = list(set(all_ct).intersection(r_ct))
    joint_ct = list(set(joint_ct).intersection(a_ct))
    joint_genes = list(set(rnaseq.index).intersection(TF))
    
    self.g = g.induced_subgraph(joint_ct)
    self.rnaseq = rnaseq.loc[joint_genes,self.g.vs["name"]]
    self.atacseq = atacseq.loc[:,self.g.vs["name"]]
    self.num_clusters = len(self.atacseq)
    
    self.atacseq_rc = self.cluster_atacseq_into_response_control()
    if not os.path.isfile(self.DE_file):
      self.DE = self.compute_DE()
      self.DE.to_csv(self.DE_file, sep=",", index=False)
    else:
      self.DE = pd.read_csv(self.DE_file, sep=",")
      
    
    
  def cluster_atacseq_into_response_control(self):
    atacseq = self.atacseq
  
    atac_01 = []
    for i in range(len(atacseq)):
      X = atacseq.iloc[i,:]
      X = X.to_numpy().reshape([len(X), 1])
      kmeans = KMeans(n_clusters=2, random_state=0).fit(X)
      labels = kmeans.labels_
      if np.mean(X[labels==1,:]) < np.mean(X[labels==0,:]):
        labels = 1 - labels
   
      atac_01.append(labels)
      
    atac_01 = np.array(atac_01)
    return atac_01

  def compute_DE(self):
      y = Yoshida.Yoshida_RNAseq()
      
      arc = self.atacseq_rc
      expr = copy.deepcopy(self.rnaseq)
      expr = 2**expr
      expr.insert(0, "gene", expr.index)
      expr_file = "expression_file_temp" + str(np.random.uniform()*1E10) + ".csv"
      expr.to_csv(expr_file, sep=",", index=False)
      
      DE_list = []
      rnaseq = self.rnaseq.to_numpy()
      for i in range(len(arc)):
        contrasts = arc[i,:]
        tb = y.DE(contrasts, expr_file)
        tb.insert(0, "cluster", i)
        
        # compute separation
        r_mins = np.quantile(rnaseq[:,contrasts==1], 
                             q=.25,
                             axis=1)
        c_maxs = np.quantile(rnaseq[:,contrasts==0],
                             q=.75,
                             axis=1)
    
        tb.insert(tb.shape[1], "sep", r_mins - c_maxs)
        DE_list.append(tb)
        
      
      os.remove(expr_file)
        
      return pd.concat(DE_list)
  
  def significant_DE(self, logFC=2, pval=0.001):
      DE = self.DE
      sig_DE = DE[(np.abs(DE["logFC"]) > logFC) & (DE["adj.P.Val"] < pval)]
      
      sig_TF = sig_DE["gene"].unique()
      nc = self.num_clusters
      sig_m = np.zeros([len(sig_TF), nc])
      for i in range(len(sig_TF)):
        TF = sig_TF[i]
        active_sig_DE = sig_DE[sig_DE["gene"]==TF]
        sig_m[i,active_sig_DE["cluster"]] = active_sig_DE["logFC"]
        
      sig_tb = pd.DataFrame(sig_m, index=sig_TF)
      
      return sig_tb
         
 
  def scatter(self, TF, cluster):
      TF = self.rnaseq.loc[TF,:]
      OCR = self.atacseq.iloc[cluster,:]
      colors = np.repeat("black", len(TF))
      rc = self.atacseq_rc[cluster,:]
      colors[rc==1] = "red"
      jitter = 0.075 * np.random.uniform(size=len(TF))
      plt.scatter(TF, OCR + jitter, c=colors)
      plt.xlabel("rnaseq")
      plt.ylabel("atacseq")
      
      return {'response':TF[rc==1],
              'control':TF[rc==0]}
      
  # show TF expression or atacseq
  def treeplot(self, TF = None, cluster=None):
        if TF is not None:
          vals = 2*self.rnaseq.loc[TF,:].to_numpy()
        else:
          OCR = self.atacseq.iloc[cluster,:]
          vals = 25*OCR
        
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
  
    
    
  
    
          