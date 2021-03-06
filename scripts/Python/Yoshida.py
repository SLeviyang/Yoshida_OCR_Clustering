import configuration as conf
import mm10 

import pandas as pd
import os, sys
import igraph as ig
import numpy as np
import matplotlib.pyplot as plt
import pdb
import qnorm


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
    root = "LTHSC.34-.BM"
   
    def load_igraph(self, just_root=True):
        e = self.load_edgelist()
        e.to_csv("temp_Yoshida_graph.csv", sep=" ", index=False, header=False)
        g = ig.Graph.Read_Ncol("temp_Yoshida_graph.csv", directed=True)
        
        if just_root:
            descendants = g.neighborhood(self.root, 
                                         order=1000, mode="out")
            g = g.induced_subgraph(descendants)
                
        os.remove("temp_Yoshida_graph.csv")
        
        return g
    
    def get_root_vertex(self, return_name=False):
        if return_name:
          return self.root
        else:
          all_names = np.array(self.load_igraph().vs["name"])
          return np.where(all_names == self.root)[0][0]
    
    def get_cell_types(self):
        return self.load_igraph().vs["name"]
    
    def load_edgelist(self):
        e = pd.read_csv(self.edgelist_file, sep=",")
        e = e.loc[e["node1"].notna() & e["node2"].notna()]
        
        return e

# Classes to access the  idr peaks.
class Yoshida_ATACseq_idr:
    
 
  def __init__(self, FDR):
     
        if not FDR in [0.01, 0.05, 0.10, 0.15]:
            sys.exit("This FDR is not available")
            
        FDR_s = str(round(100*FDR))
        if len(FDR_s) == 1:
            FDR_s = "0" + FDR_s
        FDR_s = "FDR_" + FDR_s
            
        # input files (downloaded from S3)
        id = conf.INPUT_DATA + "Yoshida_idr/"
        self.input_bed_file = id + "master_idr_nonoverlapping_" + FDR_s + ".bed"
        self.input_matrix_file = id + "master_idr_matrix_" + FDR_s + ".txt"
        
 # labels for columns of input matrix.  
  def get_input_cell_types(self):
      d = pd.read_csv(self.input_matrix_file, sep=" ", nrows=1, header=None)
      cn = d.to_numpy().tolist()[0]
      # remove .bed
      cn = [s.replace(".bed", "") for s in cn]
      
      conv = self.GEO_to_Yoshida_conversions()
      for i in range(len(conv)):
          ind = cn.index(conv.at[i,"GEO"])
          cn[ind] = conv.at[i,"Yoshida"]
          
      return cn
  
  @staticmethod
  def GEO_to_Yoshida_conversions():
      # note: orc is not included in the Yoshida matrix
      Yoshida = ['FRC.SLN', 
                 'MF.LP.SI', 
                 'MPP3.48+.BM', 
                 'MPP4.135+.BM', 
                 'T8.IEL.LCMV.d7.SI',
                 "orc"]
      GEO = ['FRC.CD140a+.Madcam-.CD35-.SLN',
             'MF.SI.LP',
             'MMP3.48+.BM',
             'MMP4.135+.BM',
             'T8.IEL.LCMV.d7.Gut',
             'orc']
      conversion = pd.DataFrame({'Yoshida':Yoshida,
                                 'GEO':GEO})
      
      return  conversion
        
  def load_idr_table(self):
      tb = pd.read_csv(self.input_matrix_file, sep=" ")
      tb.columns = self.get_input_cell_types()
      
      return tb
  
  def load_bed(self):
      return pd.read_csv(self.input_bed_file, sep="\t", header=None,
                         names=["chr", "chrStart", "chrEnd", "score"]).iloc[:,0:3]
        
 
    
        

# Class to create and access ATACseq data provided by Yoshida et al.
#
# Input:
# @param raw_data_file the csv file provided by Yoshida et al as their 
# Table S2A.
#
# Output:
# @param matrix_file path to file containing signal (rows are loci, columns
#  are cell types).  
# @param bed_file path to a bed file describing locations of loci on mm10
class Yoshida_ATACseq_counts:
    
    raw_data_file = conf.INPUT_DATA + "Yoshida_Table_S2.csv"
    
    # def __init__(self):
    #     if not os.path.isdir(self.yoshida_dir):
    #         os.mkdir(self.yoshida_dir)
    #     if not os.path.isfile(self.matrix_file):
    #         print("creating Yoshida matrix...")
    #         self.create_matrix()
    #     if not os.path.isfile(self.bed_file):
    #         print("creating Yoshida bed...")
    #         self.create_bed()
            
    
    def load_bed(self):
        d = self.load_raw_data()
      
        s = d["Summit"].to_numpy()
        chrom = d["chrom"]
        #p = d["_-log10_bestPvalue"]
        
        t = pd.DataFrame({'chr':chrom,
                           'chrStart':s-125,
                           'chrEnd':s+125})
        return t
     
    def load_count_table(self):
        d = self.load_raw_data()
        d = d.iloc[:,8:99]
        
        return d
    
    def load_raw_data(self):
        d = pd.read_csv(self.raw_data_file, sep=",")
        
        valid_chr = ["chr" + str(i) for i in range(20)] + ["chrX"]
        d = d[d["chrom"].isin(valid_chr)]
        return d
    
# Transcription factors used in Yoshida et al
class Yoshida_TF:
    
    TF_file = conf.INPUT_DATA + "Yoshida_Table_S5.csv"
    
    def get_TF(self, significant=False):
        tb = pd.read_csv(self.TF_file)
        if significant:
          tf = tb[tb["Significant"]]["TF"].tolist()
        else:
          tf = tb["TF"].tolist()
          
        return [z.upper() for z in tf]
    
# Processes RNAseq data 
# Limits to cell types in Yoshida tree and averages counts across replicates
#
# The following RNAseq data.frames are created (genes by cell types).
# In all cases counts are given in log2
# "raw" (count): raw counts restricted to cell types in Yoshida_tree
# "averaged" (count):  raw counts averaged across cell type replicates
# "quantiled" (log2): averaged counts with gene counts quantile normalized 
# across cell types
    
        
# as samples.  There are multiple samples for each cell type.
# Cell types restricted to cell types in ATACseq data
class Yoshida_RNAseq:
    GEO_RNAseq = conf.INPUT_DATA + "GSE109125_Gene_count_table_GENCODE_vM25.csv"
    
    RNAseq_dir = conf.DATA_DIR + "Yoshida/"
    RNAseq_file_prefix = RNAseq_dir + "Yoshida_RNAseq_"
    
    def __init__(self):
      self.y_ATACseq = Yoshida_ATACseq_counts()
      if not os.path.isdir(self.RNAseq_dir):
         os.mkdir(self.RNAseq_dir)
            
      if not os.path.isfile(self.RNAseq_file_prefix + "raw.csv"):
        self.create_raw_RNAseq()  
      if not os.path.isfile(self.RNAseq_file_prefix + "averaged.csv"):
        self.create_averaged_RNAseq()  
      if not os.path.isfile(self.RNAseq_file_prefix + "quantiled.csv"):
        self.create_quantiled_RNAseq()  
      #if not os.path.isfile(self.RNAseq_file_prefix + "zoom.csv"):
      #  self.create_zoom_RNAseq() 
      self.g = Yoshida_tree().load_igraph()
        
    ## constructor methods
        
    def create_raw_RNAseq(self):
        tb = pd.read_csv(self.GEO_RNAseq, sep=",")
        
        tb_genes = [g.upper() for g in tb["GeneSymbol"]]
        # limit to genes in mm10 biomart since I have TSS for those
        mm10_tb = mm10.mm10_biomart().load()
        mm10_genes = [g.upper() for g in mm10_tb["gene_name"]]
       
        
        # limit columns to cell types in Yoshida tree
        col_names = (tb.columns.tolist()).copy()
        present_in_Yoshida = np.repeat(False, len(col_names))
        # skip GeneSymbol
        present_in_Yoshida[0] = True
        for i in range(1,len(col_names)):
            # remove sample number, e.g #2
            cc = self.convert_RNAseq_to_Yoshida(col_names[i])
            if not cc is None:
                present_in_Yoshida[i] = True
                col_names[i] = cc
                
        tb.columns = col_names
        tb = tb.loc[:,present_in_Yoshida]
        tb.insert(0, "gene", tb_genes)
   
        tb = tb.drop(["GeneSymbol"], axis=1)       
        tb = tb[tb["gene"].isin(mm10_genes)]
                     
        
        tb.to_csv(self.RNAseq_file_prefix + "raw.csv", 
                  sep=",", index=False)
        
    # converts a cell type name in RNAseq GEO datasets to the
    # corresponding cell type name in Yoshida_tree
    def convert_RNAseq_to_Yoshida(self, raw_name):
      # last two characters are sample number, e.g #1
      name = raw_name[:-2]
      if name in Yoshida_tree().get_cell_types():
          return name
      
      # conversions from RNAseq name to Yoshida name
      RNAseq = ['MF.SI.LP',
             'MMP3.48+.BM',
             'MMP4.135+.BM',
             'T8.IEL.LCMV.d7.Gut']
      Yoshida = ['MF.LP.SI',
                'MPP3.48+.BM',
                'MPP4.135+.BM',
                'T8.IEL.LCMV.d7.SI']
    
      conversion = pd.DataFrame({'Yoshida':Yoshida,
                                 'RNAseq':RNAseq})
      if name in RNAseq:
          print([name, "found in converted RNAseq"])
          index = RNAseq.index(name)
          return conversion.at[index,"Yoshida"]
      else:
          print([name, "not found in Yoshida"])
          return None
      
    ## accessor methods
      
    # computes mean of cell types replicates and
    # filter out silent genes
    def create_averaged_RNAseq(self, active_gene_cutoff=0):
        m = self.load_raw_RNAseq()
        if active_gene_cutoff > 0:
            max_m = np.max(m.iloc[:,1:m.shape[1]], 1)
            m = m.loc[max_m >= active_gene_cutoff]
  
        exp_list = []
        genes = m["gene"]
        y_cell_types = Yoshida_tree().load_igraph().vs["name"]
        cell_types = list(set(y_cell_types).intersection(m.columns))
        for ct in cell_types:
            ind = np.where(m.columns == ct)[0]
            m_ct = m.iloc[:,ind]
            m_ct_mean = np.mean(m_ct, 1)
            exp_list.append(m_ct_mean)
             
        m = pd.DataFrame(np.transpose(np.array(exp_list)),
                         columns=cell_types)
        m.insert(0, "gene", genes.tolist())
            
        m.to_csv(self.RNAseq_file_prefix + "averaged.csv", 
                 sep=",", index=False)
                   
        return m
    
    def create_quantiled_RNAseq(self):
        tb = self.load_averaged_RNAseq()
        m = tb.drop("gene", axis=1).to_numpy()
        
        mq = qnorm.quantile_normalize(np.log2(1+m), axis=1)
        tb_out = pd.DataFrame(mq)
        tb_out.insert(0, "gene", tb["gene"].tolist())
        tb_out.set_axis(tb.columns.tolist(), axis=1, inplace=True)
        
        tb_out.to_csv(self.RNAseq_file_prefix + "quantiled.csv", 
                      sep=",", index=False)
        
    
    # applies a variance stabilizing transform across the genes.
    # See Yoshida methods "OCR Variance Component Analysis".  They
    # apply the Anscombe transformation.  To estimate the dispersion
    # parameter (phi in Yoshida and k in Ascombe) I used the R
    # package varistran on the average RNAseq file, see dispersion.R
    def create_zoom_RNAseq(self):
        
      os.environ["rf"] = self.transform_R_script
      os.system('/usr/local/bin/Rscript $rf')
        
    
    def get_cell_types(self):
        d = pd.read_csv(self.RNAseq_file_prefix + "averaged.csv", 
                        sep=",", nrows=1,
                        header=None)
        
        return d.iloc[0,1:d.shape[1]].tolist()
    
    ###  accessor method
    def load_raw_RNAseq(self):
        return pd.read_csv(self.RNAseq_file_prefix + "raw.csv",
                           sep=",")
    
    def load_averaged_RNAseq(self):
        return pd.read_csv(self.RNAseq_file_prefix + "averaged.csv",
                           sep=",")
    
    def load_quantiled_RNAseq(self):
        return pd.read_csv(self.RNAseq_file_prefix + "quantiled.csv", 
                           sep=",")
    
    # use limma-zoom to perform a differential expression
    # analysis
    # constrasts should be a binary array where 1 is
    # response and 0 is control
    @staticmethod
    def DE(constrasts,
           expressions_file,
           DE_R_script="../R/DE.R"):
        
      c_file = "contrasts_" + str(np.random.uniform()) + ".txt"
      o_file = "DE_output_" + str(np.random.uniform()) + ".csv"
      
      pd.DataFrame({'category':constrasts}).to_csv(c_file,
                                                   index=None)
      
      os.environ["c_file"] = c_file
      os.environ["e_file"] = expressions_file
      os.environ["o_file"] = o_file
      os.environ["rp"] = DE_R_script
      
      os.system('$rp $e_file $c_file $o_file')
      
      d =  pd.read_csv(o_file, sep=",")
     
      os.remove(c_file)
      os.remove(o_file)
     
      
      return d
  
      
        
  
 
