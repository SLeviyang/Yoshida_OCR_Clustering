import configuration as conf
import Yoshida

import pandas as pd
import os
import sys
import pdb
import numpy as np
import mm10 

# Classes to access the master idr matrix
class master_peaks:
    
  master_peaks_dir = conf.DATA_DIR + "master_peaks/"
  table_file = conf.INPUT_DATA + "master_peak_table.bed"
  input_matrix_file = conf.INPUT_DATA + "master_peak_matrix.txt"
  matrix_file = master_peaks_dir + "master_peak_matrix.csv"
  sequence_file = master_peaks_dir + "master_peak_sequences.txt"
  genomic_file = master_peaks_dir + "master_peak_genomic_info.csv"
  
  def __init__(self):
        if not os.path.isdir(self.master_peaks_dir):
            os.mkdir(self.master_peaks_dir)
        if not os.path.isfile(self.input_matrix_file):
            sys.exit("input peak matrix is missing!")
        if not os.path.isfile(self.matrix_file):
                self.create_matrix()
        if not os.path.isfile(self.table_file):
            sys.exit("peak bed file is missing!")
        if not os.path.isfile(self.sequence_file):
            if not os.path.isfile(conf.BEDTOOLS_PATH):
              print("bedtools path not found. skipping OCR sequence creation")
            else:
              print("creating sequences...")
              self.create_sequences()
        if not os.path.isfile(self.genomic_file):
            if not os.path.isfile(conf.BEDTOOLS_PATH):
              print("bedtools path not found. skipping genomic information creation")
            else:
              print("creating genomic information file...")
              self.create_genomic_information()
        
          
  def create_matrix(self):
      pdb.set_trace()
      input_matrix_cell_types = self.get_input_cell_types()
      y_cell_types = Yoshida.Yoshida_tree().load_igraph().vs["name"]
      input_matrix = pd.read_csv(self.input_matrix_file, sep=" ").to_numpy()
      
      active_cols = [x in y_cell_types for x in input_matrix_cell_types]
      active_colnames = [x for x in input_matrix_cell_types if x in y_cell_types]
      m = input_matrix[:,active_cols]
      
      d = pd.DataFrame(m)
      d.columns = active_colnames
      
      d.to_csv(self.matrix_file, sep=",", index=False)
      
      
  def create_sequences(self):
      tb = self.load_table()[["chr", "chrStart", "chrEnd"]]
      tempfile = "temp_master_peaks" + str(np.round(1E10 * np.random.uniform())) + ".bed"
     
      tb.to_csv(tempfile, sep="\t", header=None, index=False)
      seqs = mm10.mm10_genome().sequences(tempfile)
      
      os.remove(tempfile)
      pd.DataFrame({'seq':seqs}).to_csv(self.sequence_file, header=None, index=False)
      
      return seqs 
  
  def create_genomic_information(self):
      TSSfile = "temp_TSS" + str(np.round(1E10 * np.random.uniform())) + ".bed"
      outfile = "temp_out" + str(np.round(1E10 * np.random.uniform())) + ".bed"
      
      mm10.mm10_biomart().to_bed(TSSfile, padding=0)
      os.environ["btools"] = conf.BEDTOOLS_PATH
      os.environ["peakbed"] = self.table_file
      
      os.environ["tf"] = TSSfile
      os.environ["of"] = outfile
      
      os.system('$btools closest -a $peakbed -b $tf -t first -d > $of')
      tb = pd.read_csv(outfile, sep="\t", header=None,
                       names=["chr", "chrStart", "chrEnd",
                              "score", "x", "y", "z", "dTSS"])
      
      # debug check
      peak_tb = pd.read_csv(self.table_file, sep="\t", header=None)
      if not len(peak_tb) == len(tb):
          sys.exit("problem with bedtools closest script!")
          
      tb = tb["dTSS"]
      tb.to_csv(self.genomic_file, sep=",", index=False)
      
      os.remove(TSSfile)
      os.remove(outfile)
      
   
  def get_cell_types(self):
      d = pd.read_csv(self.matrix_file, sep=",", nrows=1, header=None)
      cn = d.to_numpy().tolist()[0]
      
      return cn
    
  # there is a descrepancy between a few cell type names in the Yoshida provided csv
  # file/article and in the cell type names I produce in my s3 bucket because of typos
  # in the GEO.csv file (see yoshida-atacseq workflow).  Below I fix those by hand.
 
  # labels for columns of input matrix
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
  
  def load_sequences(self):
      if os.path.isfile(self.sequence_file):
        d = pd.read_csv(self.sequence_file, header=None)
        return d[0].to_list()
      else:
        print("sequence file not created!")
        return None
    
  def load_genomic_information(self):
      if os.path.isfile(self.genomic_file):
        d = pd.read_csv(self.genomic_file, sep=",")
        return d
      else:
        print("genomic file not created!")
        return None
      
  def load_matrix(self):
      d = pd.read_csv(self.matrix_file, sep=",")
      d = d.to_numpy()
      return d
      
  def load_input_matrix(self):
      d = pd.read_csv(self.input_matrix_file, sep=" ")
      d = d.to_numpy()
      return d
  
  def load_table(self):
      d = pd.read_csv(self.table_file, sep="\t", header=None,
                      names=["chr", "chrStart", "chrEnd", "score"])
      return d
    
    

