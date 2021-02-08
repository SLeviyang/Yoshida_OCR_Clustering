import configuration as conf
import Yoshida
import utilities
import row_clusters
import mm10

import sys, os
import numpy as np
import pandas as pd
import pdb
import matplotlib.pyplot as plt
import igraph as ig
import seaborn as sns

def make_all_workflows():
    FDR = [0.01, 0.05, 0.10, 0.15]
    for f in FDR:
        workflow("idr", FDR=f)
        
    cutoffs = [2,5,10]
    for cc in cutoffs:
        workflow("count", count_cutoff=cc)


# For a particular binary matrix and bed file, run through the
# workflow of row clustering
#
# @param binary_matrix_source One of "count" or "idr"
# @param FDR float giving FDR for idr.  Only relevant for "idr" source
# @param count_cutoff scaler giving count cutoff.  Only relevant for "count"

class workflow:
    
    def __init__(self, binary_matrix_source,
                 FDR = None,
                 count_cutoff = None):
        
        # check if genomics should be implemented
        if conf.BEDTOOLS_PATH == "":
            self.include_genomics = False
        else:
            self.include_genomics = True
        
        if not binary_matrix_source in ["count", "idr"]:
            sys.exit("binary_matrix_source muast be count or idr")
            
        if binary_matrix_source == "count" and (count_cutoff is None):
            sys.exit("count_cutoff must be specified for count source")
            
        if binary_matrix_source == "idr" and (FDR is None):
            sys.exit("FDR must be specified for idr source")
            
        if binary_matrix_source == "count":
            tag = "Yoshida_count_matrix_cutoff_" + str(count_cutoff) 
        else:
            tag = "Yoshida_idr_matrix_FDR_" + str(round(100*FDR))
            
        print(["building workflow", tag])
            
        self.output_dir = conf.DATA_DIR + tag + "/"
        if not os.path.isdir(conf.DATA_DIR):
            os.mkdir(conf.DATA_DIR)
        if not os.path.isdir(self.output_dir):
            os.mkdir(self.output_dir)
            
        # output files
        self.matrix_file = self.output_dir + tag + "_matrix.csv"
        self.bed_file = self.output_dir + tag + "_loci.bed"
        self.edge_file = self.output_dir + tag + "_edges.csv"
        self.cluster_file = self.output_dir + tag + "_clusters.csv"
        self.genomics_file = self.output_dir + tag + "_genomics.csv"
        self.sequences_file = self.output_dir + tag + "_sequences.csv"
        
        # run the workflow
        
        # 1a. transfer the binary matrix and bed file
        if ((not os.path.isfile(self.matrix_file)) or 
            (not os.path.isfile(self.bed_file))):
          print("creating binary matrix and bed files...")
          self.create_workflow_inputs(binary_matrix_source, 
                                      FDR, count_cutoff)
        
        # 1b. if genomics is needed, use bed file to create
        # a genomics table and the loci sequences
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
          self.cluster_rows()
          
        
    def create_workflow_inputs(self, binary_matrix_source, 
                               FDR, count_cutoff):
        if binary_matrix_source == "count":
            y = Yoshida.Yoshida_ATACseq_counts()
            tb = y.load_count_table()
            b = y.load_bed()
            
            m = (tb.to_numpy() > count_cutoff)*np.int32(1)
        else:
            y = Yoshida.Yoshida_ATACseq_idr(FDR)
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
        
        b.to_csv(self.bed_file, sep="\t", index=False)
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
    
        
        
    def cluster_rows(self):
        m = self.load_matrix()
        rc = row_clusters.row_cluster(m)
        rc.edges.to_csv(self.edge_file, sep=",", index=False)
        rc.clusters.to_csv(self.cluster_file, sep=",", index=False)
            
    def load_matrix(self):
        return pd.read_csv(self.matrix_file, sep=",").to_numpy()
        
    def load_bed(self):
        return pd.read_csv(self.bed_file, sep="\t",
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
    
    def load_tree(self):
        return Yoshida.Yoshida_tree().load_igraph()
    
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
      

    
        
            
            
        
            
        
        
        