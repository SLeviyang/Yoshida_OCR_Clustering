import configuration as conf

from pybiomart import Dataset
import numpy as np
import pandas as pd
import os
import sys
import pdb

class mm10_biomart:
    #Downloads, saves, and returns biomart annotations for mm10
    mm10_dir = conf.DATA_DIR + "mm10/"
    csvfile_path = mm10_dir + "mm10_biomart.csv"
    
    def __init__(self):
        if not os.path.isdir(self.mm10_dir):
            os.mkdir(self.mm10_dir)
        if not os.path.isfile(self.csvfile_path):
            print("creating local biomart tables...")
            self.create_biomart_tables()
    
    def create_biomart_tables(self):
    
        dataset = Dataset(name='mmusculus_gene_ensembl',
                  host='http://www.ensembl.org')
        
        att = ['ensembl_gene_id',
               "ensembl_transcript_id",
               "refseq_mrna",
               "chromosome_name", 
               "transcription_start_site",
               "strand", 
               "transcript_count", 
               "external_gene_name",
               "transcript_biotype"]

        gene_data = dataset.query(attributes=att)
        gene_data.columns = ["gene_id", "transcript_id",
                       "refseq_id", "chr", 
                       "TSS", "strand", "transcript_count",
                       "gene_name", "transcript_type"]
  
        # filter for protein coding genes and core chr
        gene_data = gene_data.loc[gene_data["transcript_type"]=="protein_coding"]
        core_chr = [str(x) for x in np.arange(1, 20)] + ["X"]
        gene_data = gene_data.loc[[c in core_chr for c in gene_data["chr"]]]
  
        # tidy up
        gene_data["gene_name"] = gene_data["gene_name"].str.upper()
        gene_data["chr"] = ["chr" + s for s in gene_data["chr"]]
  
        gene_data.to_csv(self.csvfile_path, sep=",", index=False)
  
        return gene_data

    def load(self):
        d = pd.read_csv(self.csvfile_path, sep=",")
        return d
    
    def to_bed(self, outfile, padding=0,
               additional_columns=["gene_name"]):
        d = self.load()
        
        b = pd.DataFrame({'chr':d["chr"],
                          'chrStart':d["TSS"] - padding,
                          'chrEnd':d["TSS"] + padding})
    
        if additional_columns is not None:
            for ac in additional_columns:
              b.insert(b.shape[1], ac, d[ac])
        
        b = b.sort_values(["chr", "chrStart"])
        b.to_csv(outfile, sep="\t", index=False, header=None)

class mm10_genome:
    
    g_path = conf.MM10_FASTA
    gi_path = conf.MM10_FASTA + ".fai"
   
    def __init__(self):
        if not os.path.isfile(self.g_path):
            sys.exit("mm10 fasta file does not exist")
        if not os.path.isfile(self.gi_path):
            sys.exit("mm10 fasta file index does not exist")
        
    def sequences(self, bed_file):
     
        os.environ["btools"] = conf.BEDTOOLS_PATH
        os.environ["bf"] = bed_file
        os.environ["ff"] = self.g_path
        
        tempfile = "temp" + str(np.round(1E10 * np.random.uniform())) + ".bed"
        trashfile = "trash" + str(np.round(1E10 * np.random.uniform())) + ".bed"
        os.environ["tf"] = tempfile
        os.environ["trf"] = trashfile
        os.system('$btools getfasta -fi $ff -bed $bf -bedOut -fo $tf > $trf')
        seqs = pd.read_csv(tempfile, header=None, names=["seqs"])["seqs"].to_list()
        os.remove(tempfile)
        os.remove(trashfile)
        
        return seqs
        
        
        
