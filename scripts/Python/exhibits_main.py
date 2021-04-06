#This code creates the exhibits and tables
#by importing and calling the functions in
#tables.py and exhibits.py

import os
os.chdir(r'D:\leviyang\github\code')
import workflow
import cell_clusters
import exhibits
import tables
import copy
import pandas as pd
import numpy as np
from eval_clusters import eval_clusters

if __name__ == '__main__':
    
    '''
    Code for tables 1 and 2
    '''
    
    #Table 1 dataframe - Table manually created in LaTex
    IDR_pct_cell_specific_rows = tables.compareIDRFDRs()
    
    #Table 2 dataframes - Table manually created in LaTex
    Edge_pct_loci_clustered = tables.compareEdgeFDR_PctLociClustered()
    Edge_large_clusters = tables.compareEdgeFDR_LargeClusters()
    Edge_top20 = tables.compareEdgeFDR_PctLociTop20()
    
    '''
    Code for Figures 3, 4, 8, 9, and 10
    '''
    
    #Set parameters to edge fdr = 0.001 and idr fdr = 0.01
    w = workflow.workflow("idr", edge_FDR=0.001, idr_FDR=0.01, ncore=4)
    
    #Names of each cell type
    cell_types = w.get_cell_types()
    
    #M list for this set of parameters
    m_list = w.load_m_list()
    
    #Get tree - used in plotting the assignments
    g = w.load_tree()
    
    #k is the number of column clusters
    for k in [4,7,8]:
        
        #Get assignments for k clusters based on our methodology
        fit = w.form_tree_cluster_from_fits(k).assignments
        
        #Converts each bicluster to two mean values
        _, optimized_twomeans_df, _, _ = eval_clusters.characterizeClusters(m_list,fit)  
        
        #Heatmap where each bicluster is assigned the mean of the bicluster
        exhibits.clusterMapMeans(copy.deepcopy(m_list), fit, cell_types, k)
        
        #Heatmap where each bicluster is assigned to two mean values
        exhibits.clusterMap2Means(copy.deepcopy(m_list), optimized_twomeans_df, fit, cell_types, k)
        
        #Heatmap where each bicluster is assigned 0 or 1
        exhibits.clusterMapBinary(copy.deepcopy(m_list), optimized_twomeans_df, fit, cell_types, k)
                
        cc = cell_clusters.cell_cluster(k,m_list)
        cc.initialize_assignments(fit)
        
        #Treeplot using our algorithm - respecting the tree
        cc.treeplot(g, 'D:/leviyang/github/data_output/treeplot_optimized_wrt_tree_' + str(k) + '.png')
        
        #Initialize assignments using hierarchical clustering
        cc.initialize_assignments()        
        
        #Treeplot using hierarchical clustering
        cc.treeplot(g, 'D:/leviyang/github/data_output/treeplot_hierarchical_' + str(k) + '.png')
        
        #Optimize without regard for the tree
        cc.optimize()
        
        #Treeplot using K-means clustering
        cc.treeplot(g, 'D:/leviyang/github/data_output/treeplot_optimized_no_tree_' + str(k) + '.png')
        
    
    '''
    Code for Figures 5, 6, and 7
    '''
    
    #Set parameters to edge fdr = 0.001 and idr fdr = 0.01
    w = workflow.workflow("idr", edge_FDR=0.001, idr_FDR=0.01, ncore=4)
    
    #Names of each cell type
    cell_types = w.get_cell_types()
    
    #M list for this set of parameters
    m_list = w.load_m_list()    
    
    #DataFrames used containing the series to plot
    pct_var_across = pd.DataFrame(np.zeros((11,3)), index = range(2,13), columns = ['Hier', 'Optimized No Tree', 'Optimized Tree'])
    pct_info = pd.DataFrame(np.zeros((11,3)), index = range(2,13), columns = ['Hier', 'Optimized No Tree', 'Optimized Tree'])
    pct_info_approx = pd.DataFrame(np.zeros((11,3)), index = range(2,13), columns = ['K Means', 'Binary', 'Means'])
    
    #k is the number of column clusters
    for k in range(2,13):
        
        #Get assignments for k clusters based on our methodology
        fit = w.form_tree_cluster_from_fits(k).assignments
                
        #Percent of variance across clusters using our methodology
        pct_var_across_optimized_tree, _, _       = eval_clusters.evaluateClusters(m_list,fit)
        
        #Percent of variance retained using our methodology when approximated using the three different methods
        _, _, pct_info_optimized_tree_twomeans = eval_clusters.characterizeClustersKMeans(m_list,fit)
        _, _, pct_info_optimized_tree_means    = eval_clusters.characterizeClustersMeans(m_list,fit)
        _, _, pct_info_optimized_tree_binary   = eval_clusters.characterizeClustersBinary(m_list,fit)
        
        #Add percent retained variance for this k - 3 approximation methods - 1 clustering method
        pct_info_approx.at[k,'K Means'] = pct_info_optimized_tree_twomeans
        pct_info_approx.at[k,'Binary'] = pct_info_optimized_tree_binary
        pct_info_approx.at[k,'Means'] = pct_info_optimized_tree_means
        
        #Hierarchical clustering
        cc = cell_clusters.cell_cluster(k,m_list)
        cc.initialize_assignments()
        
        #Percent of variance across clusters using hierarchical clustering
        pct_var_across_hier, _, _  = eval_clusters.evaluateClusters(m_list,cc.assignments)
        
        #Percent of variance retained using high-low approximation
        _, _,pct_info_hier = eval_clusters.characterizeClustersBinary(m_list,cc.assignments)
        
        #K-means clustering
        cc.optimize()
        
        #Percent of variance across clusters using k-means clustering
        pct_var_across_optimized, _, _  = eval_clusters.evaluateClusters(m_list,cc.assignments)
        
        #Percent of varaiance retained using high-low approximation
        _, _, pct_info_optimized = eval_clusters.characterizeClustersBinary(m_list,cc.assignments)
        
        #Add percent variance across clusters for this k - 3 clustering methods
        pct_var_across.at[k,'Hier'] = pct_var_across_hier
        pct_var_across.at[k,'Optimized No Tree'] = pct_var_across_optimized
        pct_var_across.at[k,'Optimized Tree'] = pct_var_across_optimized_tree    
 
        #Add percent retained variance for this k - 1 approximation method - 3 clustering methods
        pct_info.at[k,'Hier'] = pct_info_hier
        pct_info.at[k,'Optimized No Tree'] = pct_info_optimized
        pct_info.at[k,'Optimized Tree'] = pct_info_optimized_tree_binary
        
    #Plot Figure 5
    exhibits.pctEnergyCharacterizedAproxx(pct_info_approx)
    
    #Plot Figure 6
    exhibits.pctVarAcrossClusters(pct_var_across)
    
    #Plot Figure 7
    exhibits.pctEnergyCharacterized(pct_info)