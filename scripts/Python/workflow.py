import configuration as conf
import Yoshida
import utilities
import row_clusters
import mm10
import optimize_tree_cluster as otc
import tree_cluster

import sys, os
import numpy as np
import pandas as pd
import pdb
import matplotlib.pyplot as plt
import igraph as ig
import seaborn as sns
from pandarallel import pandarallel as pll
from sklearn.cluster import KMeans


def make_all_workflows(ncore=4):
  edge_FDR = [0.01, 0.001, 0.0001]
  idr_FDR = [0.01, 0.05, 0.10, 0.15]
  tb = pd.DataFrame([[x,y] for x in idr_FDR for y in edge_FDR],
                    columns=["idr_FDR", "edge_FDR"]).sort_values("idr_FDR")
  
  # make idr workflows
  for i in range(len(tb)):
      workflow("idr", 
               edge_FDR=tb["edge_FDR"].iloc[i],
               idr_FDR=tb["idr_FDR"].iloc[i])
  
  # count_cutoffs = [2,5,10]
  # tb = pd.DataFrame([[x,y] for x in count_cutoffs for y in edge_FDR],
  #                   columns=["cutoff", "edge_FDR"]).sort_values("cutoff")
  
  # tb.apply(lambda s : workflow("count", 
  #                              edge_FDR=s["edge_FDR"],
  #                              count_cutoff=s["cutoff"]),
  #          axis=1)


# For a particular binary matrix and bed file, run through the
# workflow of row clustering
#
# @param binary_matrix_source One of "count" or "idr"
# @param edge_FDR FDR for calling edges in clustering
# @param idr_FDR float giving FDR for idr.  Only relevant for "idr" source
# @param count_cutoff scaler giving count cutoff.  Only relevant for "count"
# @param n_starting_points number of starting points to test for tree clustering,
# This number of random and smart (see optimize_tree_cluster) starting points are
# used
# @param max_K columns cluster for k=2,3,..,max_K


class workflow:
    
    # contains TF motifs
    #meme_file = conf.INPUT_DATA + "JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt"
    meme_file = conf.INPUT_DATA + "meme.txt"
    def __init__(self, 
                 binary_matrix_source,
                 edge_FDR,
                 idr_FDR = None,
                 count_cutoff = None,
                 n_starting_points = 10,
                 max_K=12,
                 ncore=4):
        
        # check if genomics should be implemented
        if conf.BEDTOOLS_PATH == "":
            self.include_genomics = False
        else:
            self.include_genomics = True
        self.edge_FDR = edge_FDR
        
        if not binary_matrix_source in ["count", "idr"]:
            sys.exit("binary_matrix_source muast be count or idr")
            
        if binary_matrix_source == "count" and (count_cutoff is None):
            sys.exit("count_cutoff must be specified for count source")
            
        if binary_matrix_source == "idr" and (idr_FDR is None):
            sys.exit("idr_FDR must be specified for idr source")
            
        if binary_matrix_source == "count":
            tag = "Yoshida_count_matrix_cutoff_" + str(count_cutoff) 
        else:
            tag = "Yoshida_idr_matrix_idr_" + str(idr_FDR)
        tag = tag + "_edge_FDR_" + str(edge_FDR)
            
        print(["building workflow", tag])
            
        self.output_dir = conf.DATA_DIR + tag + "/"
        if not os.path.isdir(conf.DATA_DIR):
            os.mkdir(conf.DATA_DIR)
        if not os.path.isdir(self.output_dir):
            os.mkdir(self.output_dir)
            
        self.g = Yoshida.Yoshida_tree().load_igraph()
        self.m_list = None
        self.ncore = ncore
            
        # OUTPUT FILES
        # starting data for the workflow
        self.matrix_file = self.output_dir + "matrix.csv"
        self.bed_file = self.output_dir  + "loci.bed"
        
        # genomics info
        self.genomics_file = self.output_dir + "genomics.csv"
        self.sequences_file = self.output_dir  + "sequences.csv"
        
        # row clustering files
        self.edge_file = self.output_dir  + "edges.csv"
        self.cluster_file = self.output_dir + "clusters.csv"
        # tree cluster files
        self.starting_points_file = self.output_dir + "starting_points.csv"
        self.fits_file = self.output_dir + "tree_cluster_fits.csv"
        
        # run the workflow
        
        # 1a. transfer the binary matrix and bed file
        if ((not os.path.isfile(self.matrix_file)) or 
            (not os.path.isfile(self.bed_file))):
          print("creating binary matrix and bed files...")
          self.create_workflow_inputs(binary_matrix_source, 
                                      idr_FDR, count_cutoff)
        
        # 1b. if genomics is needed, use bed file to create
        # a genomics table, the loci sequences
        if self.include_genomics:
            if not os.path.isfile(self.sequences_file):
              print("creating sequence file...")
              self.create_sequences()
            if not os.path.isfile(self.genomics_file):
              print("creating genomics file...")
              self.create_genomics()
              
   
        # 2. apply row clustering to the binary matrix
        if ((not os.path.isfile(self.edge_file)) or 
            (not os.path.isfile(self.cluster_file))):
          print("clustering rows of binary matrix...")
          self.cluster_rows(edge_FDR)
          
        self.m_list = self.get_row_cluster_view().form_m_list(min_cluster_size=30)
          
        # 3. create starting points for column clustering
        if not os.path.isfile(self.starting_points_file):
            self.create_starting_points(N=n_starting_points,
                                        start_k=2,
                                        end_k=max_K,
                                        g=self.g, 
                                        m_list=self.m_list)
            
        # 4. apply column clustering optimizations using starting
        # points
        if not os.path.isfile(self.fits_file):
            self.create_fits(g=self.g,
                             m_list=self.m_list)
            
    # constructor methods
        
    def create_workflow_inputs(self, binary_matrix_source, 
                               idr_FDR, count_cutoff):
        if binary_matrix_source == "count":
            y = Yoshida.Yoshida_ATACseq_counts()
            tb = y.load_count_table()
            b = y.load_bed()
            
            m = (tb.to_numpy() > count_cutoff)*np.int32(1)
        else:
            y = Yoshida.Yoshida_ATACseq_idr(idr_FDR)
            tb = y.load_idr_table()
            b = y.load_bed()
            
            m = tb.to_numpy()
            
            
        # select relevant cell types
        m_cell_types = tb.columns   
        cell_types = self.get_cell_types()
        missing_ct = set(cell_types).difference(m_cell_types)
        if len(missing_ct) > 0:
            print(missing_ct)
            sys.exit("some cell types missing in data")
            
        y_ind = utilities.match(cell_types, m_cell_types)
        tb_out = pd.DataFrame(m[:,y_ind], 
                              columns = cell_types)
        
        b.to_csv(self.bed_file, sep="\t", index=False, header=False)
        tb_out.to_csv(self.matrix_file, sep=",", index=False)
        
    def create_sequences(self):
      tb = self.load_bed()[["chr", "chrStart", "chrEnd"]]
      tempfile = "temp_master_peaks" + str(np.round(1E10 * np.random.uniform())) + ".bed"
     
      tb.to_csv(tempfile, sep="\t", header=None, index=False)
      seqs = mm10.mm10_genome().sequences(tempfile)
      
      os.remove(tempfile)
      pd.DataFrame({'seq':seqs}).to_csv(self.sequences_file, header=None, index=False)
     
  
    def create_genomics(self):
      seqs = self.load_sequences()
      TSSfile = "temp_TSS" + str(np.round(1E10 * np.random.uniform())) + ".bed"
      outfile = "temp_out" + str(np.round(1E10 * np.random.uniform())) + ".bed"
      
      temp_bed = "temp_loci" + str(np.round(1E10 * np.random.uniform()))  \
                     + ".bed"
      self.load_bed().to_csv(temp_bed, sep="\t", header=None,
                             index=False)
      mm10.mm10_biomart().to_bed(TSSfile, padding=0)
      
      os.environ["btools"] = conf.BEDTOOLS_PATH
      os.environ["tf"] = TSSfile
      os.environ["of"] = outfile
      os.environ["peakbed"] = temp_bed
      
      os.system('$btools closest -a $peakbed -b $tf -t first -d > $of')
      
      # chr1	3119692	3120192	chr1	3671498	3671498	XKR4	551306
      bedtools_tb = pd.read_csv(outfile, sep="\t", header=None,
                       names=["chr", "chrStart", "chrEnd",
                              "x", "y", "z", "gene", "dTSS"])
      
      # debug check
      peak_tb = pd.read_csv(temp_bed, sep="\t", header=None)
      if not len(peak_tb) == len(bedtools_tb):
          sys.exit("problem with bedtools closest script!")
      bedtools_tb = bedtools_tb[["dTSS", "gene"]]
          
      # nuc freqs
      nucs = ["A", "C", "G", "T", "GC", "CG"]
      nuc_tb = utilities.nucleotide_frequency_table(seqs, nucs=nucs)
    
      tb = pd.concat([bedtools_tb, nuc_tb], axis=1)
      tb.to_csv(self.genomics_file, sep=",", index=False)
      
      os.remove(TSSfile)
      os.remove(outfile)
      os.remove(temp_bed)
      
  
        
    def cluster_rows(self, edge_FDR):
        m = self.load_matrix()
        rc = row_clusters.row_cluster(m, FDR=edge_FDR)
        rc.edges.to_csv(self.edge_file, sep=",", index=False)
        rc.clusters.to_csv(self.cluster_file, sep=",", index=False)
        
    def create_starting_points(self, N, start_k, end_k,
                               g, m_list):      
        tb_list = []
        for k in range(start_k, end_k+1):
            print(["generating starting points for K=", k])
            oc = otc.optimize_tree_cluster(k, g, m_list, 
                                           nCPU=self.ncore)
            oc.generate_starting_points(n_random=N, n_smart=N)
            ctb = oc.get_starting_points(as_dataFrame=True)
            ctb.insert(ctb.shape[1], "K", k)
            tb_list.append(ctb)
            
        tb = pd.concat(tb_list, axis=0)
        tb.to_csv(self.starting_points_file, sep=",", index=False)
        
    def create_fits(self, g, m_list):
        tb = self.load_starting_points().groupby("K")
        fit_tb_list = []
        for k,ctb in tb:
            temp_file = "temp_sp_" + str(np.random.uniform()) + ".csv"
            ctb = ctb.drop("K", axis=1)
            ctb.to_csv(temp_file, sep=",", index=False)
            oc = otc.optimize_tree_cluster(k, g, m_list,
                                           starting_points_file=temp_file,
                                           nCPU=self.ncore)
            # check
            if not ctb.shape[0] == len(oc.get_starting_points()):
                sys.exit("this should never happen!")
            oc.optimize()
            ctb = oc.get_fits(as_dataFrame=True)
            ctb.insert(ctb.shape[1], "K", k)
            fit_tb_list.append(ctb)
        
        fit_tb = pd.concat(fit_tb_list, axis=0)
        fit_tb.to_csv(self.fits_file, sep=",", index=False)
        
    # accessor methods
    def load_matrix(self):
        return pd.read_csv(self.matrix_file, sep=",").to_numpy()
        
    def load_bed(self):
        return pd.read_csv(self.bed_file, sep="\t", header=None,
                           names=["chr", "chrStart", "chrEnd"])
    
    def load_genomics(self):
        if self.include_genomics:
          return pd.read_csv(self.genomics_file, sep=",")
        else:
          return None
    
    def load_sequences(self):
        if self.include_genomics:
          m =  pd.read_csv(self.sequences_file, sep=",")
          return m.iloc[:,0].tolist()
        else:
          return None
    
    def get_cell_types(self):
        return Yoshida.Yoshida_tree().get_cell_types()
    
    def load_clusters(self):
        return pd.read_csv(self.cluster_file, sep=",")
    
    def load_starting_points(self):
        return pd.read_csv(self.starting_points_file, sep=",")
    
    def load_fits(self):
        return pd.read_csv(self.fits_file, sep=",")
    
    def load_tree(self):
        return self.g
    
    # input matrix split into list of matrices corresponding to row clusters
    def load_m_list(self):
        return self.m_list
    
    # returns a tree cluster object
    def form_tree_cluster_from_fits(self, k, index=None):
        tb = self.load_fits()
        tb = tb[tb["K"] == k]
        
        if index is None:
            index = np.argmin(tb["residual"])
        tb = tb.drop(["residual","generated","K"], axis=1)
        a = tb.iloc[index].to_numpy()
        
        g = self.load_tree()
        m_list = self.load_m_list()
        tc = tree_cluster.tree_cluster(k, g, m_list)
        tc.initialize_components(assignments=a)
        
        return tc
        
    def load_TF_matrix(self):
        
        motif_file = self.output_dir + "motif.csv"
       
        if not os.path.isfile(motif_file):
            print(["motif match file has not been created",
                   "see script/R/create_TF_matrices.R"])
      
            
        m = pd.read_csv(motif_file, sep=",")
        
        return m
    
    # for each row cluster, split columns by on/off. 
    # If k=0 then each column is indpendently assigned,
    # otherwise the tree-based column clusters are assigned
    def biclusters_to_response_control(self, 
                                      m_list=None,
                                      k=0):
      rcv = self.get_row_cluster_view()
      atacseq = rcv.column_means_of_clusters(m_list=m_list)
     
      atac_01 = []
      for i in range(len(atacseq)):
        X = atacseq[i,:]
        X = X.reshape([len(X), 1])
        kmeans = KMeans(n_clusters=2, random_state=0).fit(X)
        labels = kmeans.labels_
        if np.mean(X[labels==1,:]) < np.mean(X[labels==0,:]):
          labels = 1 - labels
   
        atac_01.append(labels)
      
      atac_01 = np.array(atac_01)
      if k > 0:
        tc = self.form_tree_cluster_from_fits(k=k)
        assignments = tc.get_assignments()
        for ck in range(k):
          print(ck)
          cols = np.where(assignments == ck)[0]
          clust_votes = np.mean(atac_01[:,cols],1)
          for cc in cols:
            atac_01[:,cc] = 1 * (clust_votes > .5)
        
      return atac_01
        

    
    def get_row_cluster_view(self):
        return row_cluster_view(self.load_matrix(),
                                self.load_clusters(),
                                self.load_tree(),
                                self.load_genomics(),
                                self.load_sequences(),
                                include_genomics=self.include_genomics)
        
    
        
# Class allows for viewing row clusters with different
# anotations
class row_cluster_view:
    
    def __init__(self, m, clusters, tree, genomics, seqs,
                 include_genomics):
        # DataFrame with columns row and cluster
        self.clusters = clusters
        self.m = m
        self.g = tree
        
        self.genomics = genomics
        self.sequences = seqs
        self.include_genomics = include_genomics
    
    def load_clusters(self):
        return self.clusters
    
    def load_matrix(self, restrict_to_clusters=False):
        if restrict_to_clusters:
            return self.m[self.load_clusters()["row"],:]
        else:
          return self.m
    
    def get_cell_types(self):
        return self.g.vs["name"]
             
    def get_cluster_info(self, num_clusters=None):
        if num_clusters is None:
            num_clusters = np.max(self.clusters["cluster"]) + 1
        info = []
        for i in range(num_clusters):
            m = self.form_cluster_matrix(i)
            if self.include_genomics:
              g = self.form_cluster_genomics(i)
              ci = {'cluster':i,
                  'size':m.shape[0],
                  'promoters':np.mean(g["dTSS"] < 500),
                  'enhancers':np.mean(g["dTSS"] > 3000),
                  'fraction_on':np.mean(m)}
            else:
                ci = {'cluster':i,
                  'size':m.shape[0],
                  'fraction_on':np.mean(m)}
            info.append(ci)
        
        return pd.DataFrame(info)
    
    def form_m_list(self, min_cluster_size=30):
        ci = self.get_cluster_info()
        ci = ci[ci["size"] >= min_cluster_size]
        
        active_clusters = ci["cluster"].tolist()
        m_list = []
        for i in active_clusters:
            m_list.append(self.form_cluster_matrix(i))
            
        return m_list
    
    def form_cluster_matrix(self, index):
        
        clusters = self.load_clusters()
        rows = clusters.loc[clusters["cluster"]==index]["row"]
      
        m = self.m[rows,:]
        return m
    
    def form_cluster_sequences(self, index):
        clusters = self.load_clusters()
        rows = clusters.loc[clusters["cluster"]==index]["row"]
        
        all_seqs = self.sequences
        seqs = [all_seqs[i] for i in rows]
        
        return seqs
    
    def form_cluster_genomics(self, index):
        clusters = self.load_clusters()
        rows = clusters.loc[clusters["cluster"]==index]["row"]
        
        g = self.genomics
        g = g.iloc[rows]
        
        return g
    
    def column_means_of_clusters(self, m_list=None):
      if m_list is None:
        m_list = self.form_m_list()
      
      atacseq = np.array([np.mean(m, axis=0) for m in m_list])
  
      return atacseq
    
    
    # plot column means
    def scatterplot(self, index):      
          m = self.form_cluster_matrix(index)
          m_cell_types = self.get_matrix_column_names()
          y_cell_types = Yoshida.Yoshida_tree().load_igraph().vs["name"]
        
          m = m[:,utilities.match(y_cell_types, m_cell_types)]
          
      
          vc = np.mean(m, 0)
          plt.xlim(0, m.shape[1]+1)
          plt.ylim(0, 1)
          plt.title("Cluster" + str(index))

          plt.scatter(range(len(vc)), vc)
    
    def treeplot(self, index, m=None):
        if m is None:
          m = self.form_cluster_matrix(index)
        m_cell_types = self.get_cell_types()
        nct = len(m_cell_types)
        
        g = Yoshida.Yoshida_tree().load_igraph()
        
        mm = np.mean(m, axis=0)
        val_d = {m_cell_types[i]: mm[i] for i in range(nct)}
        
        vals = np.array([val_d[ct] for ct in g.vs["name"]])
        #plt.ylim(0, 1)
        #plt.scatter(range(len(vals)), vals)
        
        vs = {}
        vs["bbox"] = (1200, 1000)
        vs["vertex_size"] = 30*vals
        vs["vertex_label_size"] = 20
        vs["vertex_label"] = [str(i)  for i in range(g.vcount())]
        vs["vertex_label_dist"] = 1.5
           
        layout = g.layout_reingold_tilford(mode="all")
     
        pl = ig.plot(g, layout=layout, **vs)
        pl.show()
        
    def heatmap(self, ncluster=None, min_cluster_size=40):
        if ncluster is None:
          ncluster = np.max(self.clusters["cluster"]) + 1
        means = []
        indices = []
        for index in range(ncluster):
            cm = self.form_cluster_matrix(index)
            if len(cm) > min_cluster_size:
              cm_mean = np.mean(cm, 0)
              means.append(cm_mean.tolist())
              indices.append(index)
        
        means_m = np.array(means)
        tb = pd.DataFrame(means_m)
        tb.index = indices
        sns.heatmap(tb)
      

    
        
            
            
        
            
        
        
        