import configuration as conf
import peak_clusters as pc
import gene_clusters as gc
import master_peaks as mp
import utilities
import Yoshida

import igraph as ig
import os
import sys
import pdb
import pandas as pd
import numpy as np
import seaborn as sns
import sklearn.preprocessing as skpre
import scipy.stats as stats
import matplotlib.pyplot as plt


# for each gene determine the number of OCR within each 
# peak_cluster split by enhancer and promoter
class genes_and_peaks:
    
    gene_directory = conf.DATA_DIR + "gene_and_peaks/"
    gene_file = gene_directory + "gene_and_peaks_table.csv"
    
    def __init__(self, peak_FDR=1E-3, gene_FDR=1E-6):
        self.peak_FDR = peak_FDR
        self.gene_FDR = gene_FDR
    
        if not os.path.isdir(self.gene_directory):
            os.mkdir(self.gene_directory)
        if not os.path.isfile(self.gene_file):
            self.create_gene_table()
            
        self.table = pd.read_csv(self.gene_file, sep=",")
        
    def create_gene_table(self):
        
        # form locus, gene, TSS table
        genomic = mp.master_peaks().load_genomic_information()
        genomic.insert(0, "locus", genomic.index)
        genomic = genomic[["locus", "gene", "dTSS"]]
        
        gclust = gc.gene_clusters(self.gene_FDR).load_clusters(True)
        gclust = gclust.rename(columns={'cluster':'gene_cluster'})
        gclust = gclust.drop("row", axis=1)
        
        pclust = pc.peak_clusters(self.peak_FDR).load_clusters()
        pclust = pclust.rename(columns={'cluster':'peak_cluster',
                                        'row':'locus'})
     
        genomic = genomic.merge(gclust, on="gene", how="left")
        no_gene_cluster = np.isnan(genomic["gene_cluster"])
        genomic.loc[no_gene_cluster,'gene_cluster'] = -1
        genomic["gene_cluster"] = genomic["gene_cluster"].astype('int32')
        
        genomic = genomic.merge(pclust, on="locus", how="left")
        no_peak_cluster = np.isnan(genomic["peak_cluster"])
        genomic.loc[no_peak_cluster,'peak_cluster'] = -1
        genomic["peak_cluster"] = genomic["peak_cluster"].astype('int32')
        
        genomic.to_csv(self.gene_file, sep=",", index=False)
        
    def get_table(self):
        return self.table
    
    def form_cluster_table(self, 
                           start_dTSS=0,
                           end_dTSS=1E10,
                           sum_across_peaks=True,
                           max_gene_cluster=14,
                           max_peak_cluster=20,
                           sig=.05,
                           merge_peak_clusters_within_genes=False):
        tb = self.get_table()
        tb = tb[(tb["dTSS"] >= start_dTSS) &
                (tb["dTSS"] <= end_dTSS)]
        tb = tb[(tb["gene_cluster"] <= max_gene_cluster) &
                (tb["peak_cluster"] <= max_peak_cluster)]
        
        tb = tb[tb["peak_cluster"] >= 0]
        tb = tb[tb["gene_cluster"] >= 0]
        
        print(["size before merge", tb.shape])
        if merge_peak_clusters_within_genes:
            tb.drop_duplicates(["gene", "peak_cluster"], inplace=True)
        print(["size after merge", tb.shape])
        
        ng_clusters = max_gene_cluster + 1
        np_clusters = max_peak_cluster + 1
        
        ngenes = np.zeros(ng_clusters)
        npeaks = np.zeros(np_clusters)
        
        for i in range(ng_clusters):
            genes_in_cluster = tb[tb["gene_cluster"]==i]["gene"].unique()
            ngenes[i] = genes_in_cluster.shape[0]
        for i in range(np_clusters):
            npeaks[i] = len(tb[tb["peak_cluster"]==i])
            
        total_ngenes = np.sum(ngenes)
        total_npeaks = np.sum(npeaks)
               
        # clusters range from -1 to maxs
        m = np.ones([max_gene_cluster+1, max_peak_cluster+1])
        counts = np.zeros([max_gene_cluster+1, max_peak_cluster+1])
        
        for names, gm in tb.groupby(["gene_cluster", "peak_cluster"]):
            gene_ind = names[0] 
            peak_ind = names[1] 
            
            # number of peaks in cell
            # number of total peaks
            if sum_across_peaks:
              successes = len(gm)
              trials = npeaks[peak_ind]
              p = ngenes[gene_ind]/total_ngenes
            else:
              successes = len(gm)
              trials = tb[tb["gene_cluster"]==gene_ind].shape[0]
              p = npeaks[peak_ind]/total_npeaks
            
            pval = stats.binom_test(successes, 
                                    n=trials, 
                                    p=p, 
                                    alternative="greater")
            
            m[gene_ind, peak_ind] = pval
            counts[gene_ind, peak_ind] = successes
           
        
        if sum_across_peaks:
            cutoff = sig/m.shape[0]
        else:
            cutoff = sig/m.shape[1]
        m_single = np.zeros([max_gene_cluster+1, max_peak_cluster+1])
        m_single[m < cutoff] = 1
        
        sns.heatmap(m_single, 
                    annot=np.round(-10*np.log10(1E-9 + m)))
        #sns.heatmap(counts, annot=True)
    

 