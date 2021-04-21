import workflow
import manuscript_master as mm
import tree_cluster 
import cell_clusters

import sys
import pdb
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pylab as plt
import matplotlib as mpl
import configuration as conf
from sklearn.cluster import KMeans


def make_column_cluster_tree_figure():
     w = workflow.workflow("idr", 0.001, 0.01)
     
     fig, (ax1, ax2) = plt.subplots(1, 2)
     
     tc = w.form_tree_cluster_from_fits(3)
     tc.treeplot(target=ax1)
     
     tc = w.form_tree_cluster_from_fits(8)
     tc.treeplot(target=ax2)
     
     outfile = mm.FIGURE_DIR + "column_cluster_trees" + ".jpeg"
     plt.savefig(outfile, bbox_inches='tight')
     
     return None
 
def make_column_cluster_var_figure():
     w = workflow.workflow("idr", 0.001, 0.01)
     m = w.load_m_list()
     m = m[1:len(m)]
     
     m_tot = np.concatenate(m)
     m_mean = np.array([np.mean(mm) for mm in m])
     m_var = np.array([np.var(mm) for mm in m])
     
     R = m_tot.shape[0]
     C = m_tot.shape[1]
     
     v = [] 
     for k in range(1,13):
         
         if k > 1:
           tc = w.form_tree_cluster_from_fits(k=k)
           a = tc.get_assignments()
         else:
           a = np.zeros(C)
         var_cc = 0
         var_wr = 0
         var_wc = 0
         f_tot = 0
         for r in range(len(m)):
           for cl in range(k):
          
               Mij = m[r][:,a==cl]
               Mrc = np.mean(Mij)
               Mc = np.mean(Mij, axis=0)
               
               f = np.size(Mij)/(R*C)
               f_tot = f_tot + f
               
               var_cc = var_cc + f*(Mrc - m_mean)**2
               var_wr = var_wr + f*np.mean((Mij - Mc)**2)
               var_wc = var_wc + f*np.mean((Mc - Mrc)**2)
               
         v.append({'k':k,
                   'cc':var_cc/m_var,
                      'wr':var_wr/m_var,
                      'wc':var_wc/m_var,
                      'check':(var_cc + var_wr + var_wc)/m_var})
     
               
     v = pd.DataFrame(v)
     return (v)
 
def make_var_across_clusters_lineplot():
    
    def pctVarianceAcrossClusters(m_list, assignments):
        
        #Number of cell type clusters
        K = len(set(assignments))
        
        #Initialize all variances to 0
        total_var_rowclusters = 0
        total_var_acrosscolclusters = 0
        
        #Group cell types by cluster
        cell_cluster_grouped = [np.where(assignments == k)[0] for k in range(K)]

        
        #Loops through locus clusters
        for m in m_list:
            
            #Mean of entire row cluster
            m_bar=m.mean()
            
            #Variance across the cell type means
            total_var_rowclusters+=np.var(np.mean(m, axis=0))
            
            #Loop through cell type clusters
            for k in range(K):
                
                #Scale variance by cluster size
                scale = len(cell_cluster_grouped[k])/len(assignments)
                
                #Variance of cell type means across clusters
                total_var_acrosscolclusters += scale*(np.mean(m[:,cell_cluster_grouped[k]]) - m_bar)**2
                
        return 100*total_var_acrosscolclusters/total_var_rowclusters
    
        
    def make_lineplot(df):
        
        #Font to match latex
        mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']},size=10)
        
        #Define figure size
        fig, ax = plt.subplots(figsize=(6,4))
        
        #Add series to the plot
        plt.plot(df.index,df['Hierarchical Clustering'], linestyle = '-', marker = 'o', color = 'black', label = 'Hierarchical Clustering')
        plt.plot(df.index,df['K-Means Clustering'], linestyle = '-', marker = 'o',color = 'blue', label = 'K-Means Clustering')
        plt.plot(df.index,df['Tree Clustering'], linestyle = '-', marker = 'o',color = 'red', label = 'Tree Clustering')
        
        #Legend location and axis labels
        plt.legend(loc='upper left',frameon=False)
        plt.ylabel('Percent of Variance Across Clusters')
        plt.xlabel('Number of Cell Type Clusters')
        
        #Remove top bar from graph
        ax.tick_params(labelright=True,right=True)
        ax.spines['top'].set_visible(False)
        
        #Set x-axis ticks
        xticks = np.arange(2,13)
        locator = mpl.ticker.FixedLocator(xticks)
        ax.xaxis.set_major_locator(locator)
        
        #Set axis limits
        ax.set_xlim(1,13)
        ax.set_ylim(0,100)
        
        #Save image
        fig.savefig(conf.DATA_DIR + '/pctvaracrossclusters.png', dpi=400, bbox_inches='tight', pad_inches=0.01)
      
    
    
    #Set parameters to edge fdr = 0.001 and idr fdr = 0.01
    w = workflow.workflow("idr", edge_FDR=0.001, idr_FDR=0.01, ncore=4)
    
    #M list for this set of parameters
    m_list = w.load_m_list()    
    
    #DataFrames used containing the series to plot
    pct_var_across = pd.DataFrame(np.zeros((11,3)), index = range(2,13), columns = ['Hierarchical Clustering', 'K-Means Clustering', 'Tree Clustering'])
    
    #k is the number of column clusters
    for k in range(2,13):
        
        #Get assignments for k clusters based on our methodology
        fit = w.form_tree_cluster_from_fits(k).assignments
                
        #Percent of variance across clusters using our methodology
        pct_var_across_optimized_tree = pctVarianceAcrossClusters(m_list,fit)
        
        #Hierarchical clustering
        cc = cell_clusters.cell_cluster(k,m_list)
        cc.initialize_assignments()
        
        #Percent of variance across clusters using hierarchical clustering
        pct_var_across_hier  = pctVarianceAcrossClusters(m_list,cc.assignments)
        
        #K-means clustering
        cc.optimize()
        
        #Percent of variance across clusters using k-means clustering
        pct_var_across_optimized  = pctVarianceAcrossClusters(m_list,cc.assignments)
        
        #Add percent variance across clusters for this k - 3 clustering methods
        pct_var_across.at[k,'Hierarchical Clustering'] = pct_var_across_hier
        pct_var_across.at[k,'K-Means Clustering'] = pct_var_across_optimized
        pct_var_across.at[k,'Tree Clustering'] = pct_var_across_optimized_tree    
        
    #Plot Figure
    make_lineplot(pct_var_across)
    
#Percent of variance retained by the different approximations  
def make_pct_var_retained_aproxx_lineplot():
    
    def pctVarRetainedMeansApprox(m_list, assignments):
        
        #Number of cell type clusters
        K = len(set(assignments))
        
        #Initialize all variances to 0
        total_var_rowclusters = 0
        total_var_characterized_rowclusters_within = 0
        
        #Group cell types by cluster
        cell_cluster_grouped = [np.where(assignments == k)[0] for k in range(K)]
        
        #Loop through locus clusters
        for m in m_list:
            
            #Mean of each cell type cluster in this locus cluster
            m_cluster_mean = np.array([m[:,cell_cluster_grouped[k]].mean() for k in range(K)])
            
            #Variance of cell type means in this locus cluster
            total_var_rowclusters+=np.var(np.mean(m, axis=0))
            
            #Loop through cell type clusters
            for k in range(K):
                
                #Variance within cell type cluster
                total_var_characterized_rowclusters_within += (1/len(assignments))*np.square(np.mean(m[:,cell_cluster_grouped[k]],axis=0) - m_cluster_mean[k]).sum()

        return 100*(1-total_var_characterized_rowclusters_within/total_var_rowclusters) 
    
    def pctVarRetainedKMeansApprox(m_list, assignments):
        
        #Set number of means in K-means algorithm to 2
        kmeans = KMeans(n_clusters = 2)
        
        #Number of cell type clusters
        K = len(set(assignments))
        
        #Initialize all variances to 0
        total_var_rowclusters = 0
        total_var_characterized_rowclusters_within = 0
        
        #Group cell types by cluster
        cell_cluster_grouped = [np.where(assignments == k)[0] for k in range(K)]
        
        #Loop through locus clusters
        for m in m_list:
            
            #Mean of each cell type in this locus cluster
            m_cluster_mean = np.array([m[:,cell_cluster_grouped[k]].mean() for k in range(K)])
            
            #Run K-means
            kmeans.fit(m_cluster_mean.reshape(-1,1))
            
            #Variance of cell type means in this locus cluster
            total_var_rowclusters+=np.var(np.mean(m, axis=0))
            
            #Loop through cell type clusters
            for k in range(K):
                
                #Mean assigned to this cell type cluster
                cluster_mean = kmeans.cluster_centers_[kmeans.labels_[k]]
                
                #Variance within this cell type cluster
                total_var_characterized_rowclusters_within += (1/len(assignments))*np.square(np.mean(m[:,cell_cluster_grouped[k]],axis=0) - cluster_mean).sum()
        
        return 100*(1-total_var_characterized_rowclusters_within/total_var_rowclusters)    
        
    def make_lineplot(df):
      
        #Font to match latex
        mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']},size=10)
        
        #Define figure size
        fig, ax = plt.subplots(figsize=(6,4))
        
        #Add series to the plot
        plt.plot(df.index,df['Means'], linestyle = '-', marker = 'o', color = 'black', label = 'Means approximation')
        plt.plot(df.index,df['High/Low'], linestyle = '-', marker = 'o',color = 'red', label = 'High/Low approximation')
        
        #Legend location and axis labels
        plt.legend(loc = 'upper left', frameon=False)
        plt.ylabel('Percent of Across Column Variance Captured')
        plt.xlabel('Number of Cell Type Clusters')
        
        #Remove top bar from graph
        ax.tick_params(labelright=True,right=True)
        ax.spines['top'].set_visible(False)
        
        #Set x-axis ticks
        xticks = np.arange(2,13)
        locator = mpl.ticker.FixedLocator(xticks)
        ax.xaxis.set_major_locator(locator)
        
        #Set axis limits
        ax.set_xlim(1,13)
        ax.set_ylim(0,100)
        
        #Save image
        fig.savefig(conf.DATA_DIR + '/pctenergycharacterized_approx.png', dpi=400, bbox_inches='tight', pad_inches=0.01)
    
    #Set parameters to edge fdr = 0.001 and idr fdr = 0.01
    w = workflow.workflow("idr", edge_FDR=0.001, idr_FDR=0.01, ncore=4)
    
    #M list for this set of parameters
    m_list = w.load_m_list()    
    
    #DataFrames used containing the series to plot
    pct_info_approx = pd.DataFrame(np.zeros((11,2)), index = range(2,13), columns = ['High/Low', 'Means'])
    
    #k is the number of column clusters
    for k in range(2,13):
        
        #Get assignments for k clusters based on our methodology
        fit = w.form_tree_cluster_from_fits(k).assignments
        
        #Add percent retained variance for this k - 3 approximation methods - 1 clustering method
        pct_info_approx.at[k,'High/Low'] = pctVarRetainedKMeansApprox(m_list,fit)
        pct_info_approx.at[k,'Means'] = pctVarRetainedMeansApprox(m_list,fit)   
        
    #Plot Figure
    make_lineplot(pct_info_approx)
    
    
#Percent of variance retained by the three different clustering methods - using the high/low approximation
def make_pct_var_retained_method_lineplot():
    
    def pctVarRetainedKMeansApprox(m_list, assignments):
        
        #Set number of means in K-means algorithm to 2
        kmeans = KMeans(n_clusters = 2)
        
        #Number of cell type clusters
        K = len(set(assignments))
        
        #Initialize all variances to 0
        total_var_rowclusters = 0
        total_var_characterized_rowclusters_within = 0
        
        #Group cell types by cluster
        cell_cluster_grouped = [np.where(assignments == k)[0] for k in range(K)]
        
        #Loop through locus clusters
        for m in m_list:
            
            #Mean of each cell type in this locus cluster
            m_cluster_mean = np.array([m[:,cell_cluster_grouped[k]].mean() for k in range(K)])
            
            #Run K-means
            kmeans.fit(m_cluster_mean.reshape(-1,1))
            
            #Variance of cell type means in this locus cluster
            total_var_rowclusters+=np.var(np.mean(m, axis=0))
            
            #Loop through cell type clusters
            for k in range(K):
                
                #Mean assigned to this cell type cluster
                cluster_mean = kmeans.cluster_centers_[kmeans.labels_[k]]
                
                #Variance within this cell type cluster
                total_var_characterized_rowclusters_within += (1/len(assignments))*np.square(np.mean(m[:,cell_cluster_grouped[k]],axis=0) - cluster_mean).sum()
        
        return 100*(1-total_var_characterized_rowclusters_within/total_var_rowclusters)    
        
    def make_lineplot(df):
        
        #Font to match latex
        mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']},size=10)
        
        #Define figure size
        fig, ax = plt.subplots(figsize=(6,4))
        
        #Add series to the plot
        plt.plot(df.index,df['Hierarchical Clustering'], linestyle = '-', marker = 'o', color = 'black', label = 'Hierarchical Clustering')
        plt.plot(df.index,df['K-Means Clustering'], linestyle = '-', marker = 'o',color = 'blue', label = 'K-Means Clustering')
        plt.plot(df.index,df['Tree Clustering'], linestyle = '-', marker = 'o',color = 'red', label = 'Tree Clustering')
        
        #Legend location and axis labels
        plt.legend(loc = 'upper left', frameon=False)
        plt.ylabel('Percent of Across Column Variance Captured')
        plt.xlabel('Number of Cell Type Clusters')
        
        #Remove top bar from graph
        ax.tick_params(labelright=True,right=True)
        ax.spines['top'].set_visible(False)
        
        #Set x-axis ticks
        xticks = np.arange(2,13)
        locator = mpl.ticker.FixedLocator(xticks)
        ax.xaxis.set_major_locator(locator)
        
        #Set axis labels
        ax.set_xlim(1,13)
        ax.set_ylim(0,100)
        
        #Save image
        fig.savefig(conf.DATA_DIR + '/pctenergycharacterized_methods.png', dpi=400, bbox_inches='tight', pad_inches=0.01)
    
    #Set parameters to edge fdr = 0.001 and idr fdr = 0.01
    w = workflow.workflow("idr", edge_FDR=0.001, idr_FDR=0.01, ncore=4)
    
    #M list for this set of parameters
    m_list = w.load_m_list()    
    
    #DataFrames used containing the series to plot
    pct_info_methods = pd.DataFrame(np.zeros((11,3)), index = range(2,13), columns = ['Hierarchical Clustering', 'K-Means Clustering', 'Tree Clustering'])
    
    #k is the number of column clusters
    for k in range(2,13):
        
        #Get assignments for k clusters based on our methodology
        fit = w.form_tree_cluster_from_fits(k).assignments
                
        #Percent of variance across clusters using our methodology
        pct_info_optimized_tree = pctVarRetainedKMeansApprox(m_list,fit)
        
        #Hierarchical clustering
        cc = cell_clusters.cell_cluster(k,m_list)
        cc.initialize_assignments()
        
        #Percent of variance across clusters using hierarchical clustering
        pct_info_hier  = pctVarRetainedKMeansApprox(m_list,cc.assignments)
        
        #K-means clustering
        cc.optimize()
        
        #Percent of variance across clusters using k-means clustering
        pct_info_optimized  = pctVarRetainedKMeansApprox(m_list,cc.assignments)
        
        #Add percent variance across clusters for this k - 3 clustering methods
        pct_info_methods.at[k,'Hierarchical Clustering'] = pct_info_hier
        pct_info_methods.at[k,'K-Means Clustering'] = pct_info_optimized
        pct_info_methods.at[k,'Tree Clustering'] = pct_info_optimized_tree   

    #Plot Figure
    make_lineplot(pct_info_methods)