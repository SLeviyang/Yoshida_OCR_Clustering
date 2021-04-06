#This code creates an evalute clusters class with all static member functions
#It is used as supporting functions for the exhibits in exhibits.py

from sklearn.cluster import KMeans
import numpy as np
import pandas as pd


class eval_clusters:
    
    #Amount missing from total varaince is col variance
    
    #m_list is a list of row clusters, assignments is an array denoting to which cluster each cell type is assigned
    
    #Determines the percent of variance attributed to cell type clusters,
    #percent of variance attributed to locus clusters, and percent of
    #variance in the columns of the original matrix
    @staticmethod
    def evaluateClusters(m_list, assignments):
        
        #Number of cell type clusters
        K = len(set(assignments))
        
        #Initialize all variances to 0
        total_var = 0
        total_var_rowclusters = 0
        total_var_withincolclusters = 0
        total_var_acrosscolclusters = 0
        
        #Group cell types by cluster
        cell_cluster_grouped = [np.where(assignments == k)[0] for k in range(K)]
        
        #Tracks percent of variance in columns by each row cluster
        pct_var_incols = []
        
        #Loops through locus clusters
        for m in m_list:
            
            #Mean of entire row cluster
            m_bar=m.mean()
            
            #Variance of entire row cluster
            total_var+=np.var(m)
            
            #Variance across the cell type means
            total_var_rowclusters+=np.var(np.mean(m, axis=0))
            
            #Tracks variance in each locus cluster
            pct_var_incols.append(100*np.var(np.mean(m, axis=0))/total_var)
            
            #Loop through cell type clusters
            for k in range(K):
                
                #Scale variance by cluster size
                scale = len(cell_cluster_grouped[k])/len(assignments)
                
                #Variance of cell type means within a cluster
                total_var_withincolclusters += scale*np.var(np.mean(m[:,cell_cluster_grouped[k]],axis=0))
                
                #Variance of cell type means across clusters
                total_var_acrosscolclusters += scale*(np.mean(m[:,cell_cluster_grouped[k]]) - m_bar)**2
                
        return 100*total_var_acrosscolclusters/total_var_rowclusters, 100*total_var_rowclusters/total_var, pct_var_incols   
    
    
    #Determines the percent of variance retained by assigning each bicluster the mean of its elements
    @staticmethod
    def characterizeClustersMeans(m_list,assignments):      
        
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

        #Does not return a new characterized version of the data
        return None, None, 100*(1-total_var_characterized_rowclusters_within/total_var_rowclusters) 
    
    #Determines the percent of variance retained by assigning each bicluster one of two means based on K-means    
    @staticmethod
    def characterizeClustersKMeans(m_list,assignments):      
        
        #Set number of means in K-means algorithm to 2
        kmeans = KMeans(n_clusters = 2)
        
        #Number of cell type clusters
        K = len(set(assignments))
        
        #New m_list which will replace elements with their new approximation
        characterized_m_list = []
        
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
            
            #Assign each cell type the mean it is closest to
            characterized_m_list.append([kmeans.cluster_centers_[label][0] for label in kmeans.labels_])
            
            #Variance of cell type means in this locus cluster
            total_var_rowclusters+=np.var(np.mean(m, axis=0))
            
            #Loop through cell type clusters
            for k in range(K):
                
                #Mean assigned to this cell type cluster
                cluster_mean = kmeans.cluster_centers_[kmeans.labels_[k]]
                
                #Variance within this cell type cluster
                total_var_characterized_rowclusters_within += (1/len(assignments))*np.square(np.mean(m[:,cell_cluster_grouped[k]],axis=0) - cluster_mean).sum()
            
        #Create dataframe version of characterized_m_list    
        df = pd.concat([pd.DataFrame(cm).T for cm in characterized_m_list])
        df.index = range(df.shape[0])
        
        return characterized_m_list, df, 100*(1-total_var_characterized_rowclusters_within/total_var_rowclusters)      
    
    #Determines the percent of variance retained by assigning each bicluster high or low 
    @staticmethod
    def characterizeClustersBinary(m_list,assignments):     
        
        #Convert row to elements of 0 or 1
        def toBinary(arr):
            
            #If both means are greater than .8 then assign open to all biclusters with this row cluster
            if arr.min()>.8:
                return np.array([1,1])
            
            #If both means are less than .2 then assign closed to all biclusters with this row cluster
            if arr.max() <.2:
                return np.array([0,0])  
            
            #Otherwise assign open to biclusters with the higher mean and closed to lower
            max_val = arr.max()
            to_return = []
            for val in arr:
                if val==max_val:
                    to_return.append(1)
                else:
                    to_return.append(0)
            return to_return
        
        #Set number of means in K-means algorithm to 2
        kmeans = KMeans(n_clusters = 2)
        
        #Number of cell type clusters
        K = len(set(assignments))
        
        #New m_list which will replace elements with their new approximation
        characterized_m_list = []
        
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
            
            #Assign each cell type the mean it is closest to
            characterized_m_list.append([kmeans.cluster_centers_[label][0] for label in kmeans.labels_])
            
            #Variance of cell type means in this locus cluster
            total_var_rowclusters+=np.var(np.mean(m, axis=0))
            
            #Loop through cell type clusters
            for k in range(K):
                
                #Mean assigned to this cell type cluster - converted to binary
                cluster_mean = toBinary(kmeans.cluster_centers_)[kmeans.labels_[k]]
                
                #Variance within this cell type cluster
                total_var_characterized_rowclusters_within += (1/len(assignments))*np.square(np.mean(m[:,cell_cluster_grouped[k]],axis=0) - cluster_mean).sum()
            
        #Create dataframe version of characterized_m_list  
        df = pd.concat([pd.DataFrame(cm).T for cm in characterized_m_list])
        df.index = range(df.shape[0])
        
        return characterized_m_list, df, 100*(1-total_var_characterized_rowclusters_within/total_var_rowclusters)      
    