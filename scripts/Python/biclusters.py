import pandas as pd
import numpy as np
import math
import random
from scipy.special import betaln
import os
import master_peaks as mp
import configuration as conf
import peak_clusters

class biclusters:
    clusters = peak_clusters.peak_clusters().clusters
    m = np.delete(mp.master_peaks().load_matrix().astype("float"),clusters['row'],0)
    cell_types = mp.master_peaks().get_cell_types()
    bicluster_directory = conf.DATA_DIR + "bicluster_unassigned/"
    bicluster_file = bicluster_directory + "biclusters.xlsx"
    bicluster_analysis_file = bicluster_directory + "bicluster_analysis.xlsx"
    
    def __init__(self):
        if not os.path.isdir(self.bicluster_directory):
            os.mkdir(self.bicluster_directory)
    
    def load_biclusters(self,num_biclusters=1):       
        filepath = self.bicluster_file
        if not os.path.isfile(filepath):
            self.LAS(num_biclusters,save=True)
        self.biclusters = []
        sheet_names = pd.ExcelFile(filepath).sheet_names
        for sheet in sheet_names:
            self.biclusters.append(pd.read_excel(filepath,index_col=0))
        
    '''
    Find a submatrix U of X that approximately maximizes the score function S
    Subtract the average of U from each of its elements in X
    Return to search
    Stop when S(U) falls below a threshold or a user-defined number of submatrices are produced
    '''
    def LAS(self, num_biclusters = 1, search_iterations = 1000, save=False):
        def scoreFunction(U, X):
            def binomln(n,k):
                #Assumes binom(n,k) >=0
                return -betaln(1+n-k,1+k)-math.log(n+1)
            m,n = X.shape
            k,l = U.shape
            log_m_k = binomln(m,k)
            log_n_l = binomln(n,l)
            z = -U.values.mean()*np.sqrt(k*l)
            return -log_m_k - log_n_l -(np.log(.5) - .717*z - .416*z**2)
        
        def search(X, prev_search, row_cap):
            m,n = X.shape
            while True:
                k = random.randint(1,min(m,row_cap))
                l = random.randint(1,int(n))
                if (k,l) not in prev_search:
                    break
            B_old = []
            A_old = []
            A = [1]*k
            B = random.sample(range(0,n), l)
            
            while set(A)!=set(A_old) and set(B)!=set(B_old):
                A_old = A
                B_old = B
                X_B = X.iloc[:,B]
                A = X_B.sum(axis=1).nlargest(k).index.tolist()
                X_A = X.loc[A,:]
                B = X_A.sum(axis=0).reset_index(drop=True).nlargest(l).index.tolist()
                
            return X_A.iloc[:,B], k ,l

        original_m = pd.DataFrame(self.m)
        X = original_m.copy()
        X = (X-X.mean())/X.std()
        X = X/abs(X)*np.log(1+ np.abs(X))
        row_cap = int(X.shape[0]/2)
        biclusters = []
        iteration = 1
        while len(biclusters) < num_biclusters:        
            print('Computing bicluster ', str(iteration))
            max_score = -np.inf
            prev_search=[]
            for i in range(search_iterations):
                U,k,l = search(X, prev_search,row_cap)
                prev_search.append((k,l))
                latest_score = scoreFunction(U,X)
                if latest_score > max_score:
                    U_star = U
                    max_score = latest_score
            X.drop(U_star.index,axis=0,inplace=True)
            X = (X-X.mean())/X.std()
            biclusters.append(original_m.loc[U_star.index,U_star.columns])    
            iteration+=1
    
        if save:
            column_mapping = dict(zip(np.arange(len(self.cell_types)),self.cell_types))
            writer = pd.ExcelWriter(self.bicluster_file)
            analysis_writer = pd.ExcelWriter(self.bicluster_analysis_file)
            i=0
            for bicluster in biclusters:
                outside_cluster = original_m.loc[bicluster.index,[col for col in original_m.columns if col not in bicluster.columns]]
                bicluster.columns = [column_mapping.get(col) for col in bicluster.columns]
                outside_cluster.columns = [column_mapping.get(col) for col in outside_cluster.columns]
                bicluster.to_excel(writer,'cluster_'+str(i))
                pd.DataFrame(bicluster.mean(),columns=['in_cluster_mean']).to_excel(analysis_writer,'cluster_'+str(i), startcol=0)
                pd.DataFrame(outside_cluster.mean(),columns=['out_cluster_mean']).to_excel(analysis_writer,'cluster_'+str(i), startcol=3)
                i+=1
            writer.save()
            analysis_writer.save()
        
        self.biclusters = biclusters