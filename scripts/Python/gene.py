import configuration as conf
import peak_clusters as pc
import master_peaks as mp
import utilities
import Yoshida

import igraph as ig
import os
import sys
import pdb
import pandas as pd
import numpy as np

class gene_OCR:
    
    gene_OCR_directory = conf.DATA_DIR + "genes/"
    
    gene_OCR_TSS_file = conf.DATA_DIR + "genes/gene_TSS_OCRs.csv"
    # distal enhancers (terminology from Yoshida)
    gene_OCR_DE_file = conf.DATA_DIR + "genes/gene_DE_OCRs.csv"
    
    def __init__(self):
        self.master_m = mp.master_peaks().load_matrix()
        self.master_g = mp.master_peaks().load_genomic_information()
        self.clusters = pc.peak_clusters().load_clusters()
        
        if not os.path.isdir(self.gene_OCR_directory):
            os.mkdir(self.gene_OCR_directory)
        if not os.path.isfile(self.gene_OCR_TSS_file):
            tb = self.create_gene_table(0, 1000)
            tb.to_csv(self.gene_OCR_TSS_file, sep=",", index=False)
        if not os.path.isfile(self.gene_OCR_DE_file):
            tb = self.create_gene_table(3000, 1E5)
            tb.to_csv(self.gene_OCR_DE_file, sep=",", index=False)
        
        
    def create_gene_table(self, 
                          dTSS_start, dTSS_end):
        
        max_cluster_number = np.max(self.clusters["cluster"])
        mast = mp.master_peaks()
        m = mast.load_matrix()
        g = mast.load_genomic_information()
        nloci = m.shape[0]
        
        p = pc.peak_clusters().load_clusters()
        p = p.loc[p["cluster"] <= max_cluster_number]
        
        # cluster of each peak, -1 means no clusters
        assignments = np.repeat(-1, nloci)
        pg = p.groupby("cluster")
        for cluster,pgc in pg:
            assignments[pgc["row"]] = cluster
        g.insert(len(g.columns), "cluster", assignments)
        
        all_genes = sorted(g["gene"].unique().tolist())
        
        g = g[(g["dTSS"] >= dTSS_start) & (g["dTSS"] <= dTSS_end)]
        cluster_counts = {'gene':all_genes}
        gg = g.groupby("cluster")
        for cluster,ggc in gg:
            if cluster == -1:
                cluster_name = "cNone"
            else:
                cluster_name = "c" + str(cluster)
            counts = np.repeat(0, len(all_genes))
            sc = ggc["gene"].value_counts()
            ind = utilities.match(sc.index.tolist(), all_genes)
            # debug
            if -1 in ind:
                sys.exit("unidentified gene!")
            counts[ind] = sc.tolist()
            cluster_counts[cluster_name] = counts
         
        return pd.DataFrame(cluster_counts)
    
    def get_gene_loci(self, gene):
        g = self.master_g
        g = g[g["gene"] == gene]
        if len(g) == 0:
            return None
        
        clusters_ind = utilities.match(g.index, self.clusters["row"])
        gene_clusters = []
        for ind in clusters_ind:
            if ind == -1:
                gene_clusters.append(-1)
            else:
                gene_clusters.append(self.clusters.at[ind,"cluster"])
        
     
        cm = np.sum(self.master_m[g.index,:], 1)
        
        tb = pd.DataFrame({'index':range(len(g.index)),
                          'locus':g.index,
                          'dTSS':g["dTSS"],
                          'nOCR':cm,
                          'cluster':gene_clusters,
                          'gene':g["gene"]})
        
        return tb
        
    def load_gene_table(self):
        tb_TSS = pd.read_csv(self.gene_OCR_TSS_file, sep=",")
        tb_TSS.insert(1, "type", "TSS")
        
        tb_DE = pd.read_csv(self.gene_OCR_DE_file, sep=",")
        tb_DE.insert(1, "type", "DE")
        
        return tb_TSS.append(tb_DE)
    
    # show locus state (0/1) across all cell types
    # index gives the number of the locus within the gene, not
    # the locus value (row) itself
    def plottree(self, gene, index):
      g = Yoshida.Yoshida_tree().load_igraph()
      y_ct = g.vs["name"]
      ct = mp.master_peaks().get_cell_types()
      
      if not set(y_ct) == set(ct):
            sys.exit("tree and matrix cell types do not match!")
      map_ct = utilities.match(y_ct, ct)
     
      locus = self.get_gene_loci(gene).index[index]
      OCR = self.master_m[locus,map_ct]
      
      vs = {}
      vs["bbox"] = (1200, 1000)
      vs["vertex_size"] = 25*OCR + 5
      vs["vertex_label_size"] = 20
      vs["vertex_label"] = [str(i)  for i in range(len(y_ct))]
      vs["vertex_label_dist"] = 1.5
           
      layout = g.layout_reingold_tilford(mode="all")
     
      pl = ig.plot(g, layout=layout, **vs)
      pl.show()
      
      
      
    
    def test(self, n1, randomize=False):
        master = mp.master_peaks()
        m = master.load_matrix()
        m = m[np.sum(m,1)==n1,:]
        nc = m.shape[1]
        
        if randomize:
            for i in range(len(m)):
                m[i,:] = m[i,np.random.permutation(nc)]
        
        ct = master.get_cell_types()
        joint = []
        for i in range(len(m)):
            ind = np.where(m[i,:]==1)[0]
            ct_ind = [ct[j] for j in ind]
         
            joint.append("___".join(ct_ind))
            
        g = pd.Series(joint)
        return g.value_counts().sort_values(ascending=False)
      
   
        
  