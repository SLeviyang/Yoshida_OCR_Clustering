#This file contains the code for tables 1 and 2
#Note Table 2 is separated into functions for each row

import workflow
import numpy as np
import pandas as pd

#Table 1: Percent of Cell Specific Loci by FDR
def compareIDRFDRs(idr_fdr_list = [0.1,0.05,0.01], edge_fdr=0.001):
    pct_cell_specific = pd.DataFrame(np.zeros((1,len(idr_fdr_list))),columns = idr_fdr_list)
    
    #Loop through each idr fdr
    for idr_fdr in idr_fdr_list:
        w = workflow.workflow("idr", edge_FDR=edge_fdr, idr_FDR=idr_fdr, ncore=4)
        m = w.load_matrix()
        
        #Remove all rows with all 0s
        m = m[np.sum(m,axis=1)>0,:]
        
        #100 * num rows with 1 or 2 / num rows
        pct_cell_specific.at[0,idr_fdr] = 100 * m[np.sum(m,axis=1)<=2,:].shape[0] / m.shape[0]
        
    return pct_cell_specific

#Table 2 - row 1: Percent of Loci with an edge
def compareEdgeFDR_PctLociClustered(edge_fdr_list = [0.1,0.05,0.01,0.001,0.0001], idr_fdr=0.01):
    pct_loci_clustered = pd.DataFrame(np.zeros((1,len(edge_fdr_list))),columns = edge_fdr_list)
    
    #Loop through each edge fdr
    for edge_fdr in edge_fdr_list:
        w = workflow.workflow("idr", edge_FDR=edge_fdr, idr_FDR=idr_fdr, ncore=4)
        m = w.load_matrix()
        
        #Loci that are not cell specific
        m = m[np.sum(m,axis=1)>2,:]
        c = w.load_clusters() 
        
        #Num loci with assignment / num non-cell specific loci
        pct_loci_clustered.at[0,edge_fdr] = 100 * c.shape[0] / m.shape[0]
        
    return pct_loci_clustered

#Table 2 - row 2: Number of Clusters larger than 40 loci
def compareEdgeFDR_LargeClusters(edge_fdr_list = [0.1,0.05,0.01,0.001,0.0001], idr_fdr=0.01):
    large_clusters = pd.DataFrame(np.zeros((1,len(edge_fdr_list))),columns = edge_fdr_list)
    
    #Loop through each edge fdr
    for edge_fdr in edge_fdr_list:
        w = workflow.workflow("idr", edge_FDR=edge_fdr, idr_FDR=idr_fdr, ncore=4)
        c = w.load_clusters() 
        
        #Num clusters larger than 40 loci
        large_clusters.at[0,edge_fdr] = (c.groupby('cluster').count()>40).sum()
        
    return large_clusters  
    
#Table 2 - row 3: Percent of Loci in the top 20 clusters
def compareEdgeFDR_PctLociTop20(edge_fdr_list = [0.1,0.05,0.01,0.001,0.0001], idr_fdr=0.01):
    pct_loci_top20 = pd.DataFrame(np.zeros((1,len(edge_fdr_list))),columns = edge_fdr_list)
    
    #Loop through each edge fdr
    for edge_fdr in edge_fdr_list:
        
        w = workflow.workflow("idr", edge_FDR=edge_fdr, idr_FDR=idr_fdr, ncore=4)
        c = w.load_clusters()
        
        #Keep top 20 largest clusters
        top20 = c.groupby('cluster').count().sort_values(by=['row'],ascending=False).head(20)
        
        #Num loci in top 20 clusters / num loci with cluster assignments
        pct_loci_top20.at[0,edge_fdr] = 100 * top20.sum() / c.shape[0]
        
    return pct_loci_top20