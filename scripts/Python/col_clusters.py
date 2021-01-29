'''
column_clusters class (parent class)
    Constructor parameters:
        k: (int). Number of clusters.
        min_size: (int). Minimum size row cluster to use.
        num_rand: (int). Number of random starting positions to make.
        
initial_column_clusters class (inherits column_clusters)
    Constuctor parameters:
        k: (int). Number of clusters.
        min_size: (int). Minimum size row cluster to use.
        num_rand: (int). Number of random starting positions to make.

    Summary:
        Defines possible initial cut edges.
        Solves for the residual of all combinations of cut edges with k clusters
        Solves for the residual of num_rand random initial cut edges
        Saves excel file with assignments (map), cut edges, residual, and how the observation was generated (random versus selected)
        
optimal_column_clusters class (inherits column_clusters)  
    Constuctor parameters:
        k: (int). Number of clusters.
        min_size: (int). Minimum size row cluster to use.
        num_rand: (int). Number of random starting positions to make.
        n: (int). Smallest n residuals to use as initial cut edges.
        num_iters: (int). Number of times to optimize each initial starting position.
        
    Key Function:
        load_column_clusters(self, overwrite=False, plot=False):
            overwrite: (boolean). Rerun and overwrite previously existing files.
            plot: (boolean). Plot the optimal tree_cluster (minimum residual).
            
        Summary:
            Loads both initial and optimal column clusters.
            Saves output if there are no prexisting files or overwrite == True.
        
    Summary:
         Uses output from initial column clusters (either loads existing or runs automatically)
         Selects n initial positions with the smallest n residuals
         Creates tree_cluster object from each initial positions and solves for optimal clustering
         Saves excel file with initial assignments (map), initial cut edges, initial residual, how the observation was generated,
         optimal assignments (map), optimal cut edges, optimal residual.
         
Note:
    Both initial_column_clusters and optimal_column_clusters have functions that use multiprocessing. Multiprocessing is demanding 
    on the computer while the functions are being called but greatly improve computation speed.
    
-------------------------------------------------------------------------------------------------------     
        
Example code that saves initial and optimal clusters:
    occ = optimal_column_clusters(k=5, min_size=200, num_rand=1000, n = 40, num_iters = 1)
    occ.load_column_clusters(overwrite=False, plot=False)
'''

import os
import pandas as pd
import networkx as nx
import numpy as np
import configuration as conf
import Yoshida
import peak_clusters
import master_peaks as mp
import multiprocessing
import itertools
import ast
from pandarallel import pandarallel
from col_cluster_suppfunctions import clusterScores
from tree_cluster import tree_cluster

class column_clusters:
    col_cluster_directory = conf.DATA_DIR + 'col_clusters/'
    
    def __init__(self,k,min_size,num_rand):
        self.num_clusters = k
        self.g = Yoshida.Yoshida_tree().load_igraph()
        self.rootnode = self.g.vs['name'][0]
        self.clusters = peak_clusters.peak_clusters().clusters
        self.num_rand = num_rand
        self.min_size = min_size
        master_peaks = mp.master_peaks()
        self.m = pd.DataFrame(master_peaks.load_matrix().astype('float'),columns=master_peaks.get_cell_types()).loc[self.clusters.groupby('cluster').filter(lambda x: len(x) > min_size)['row']]
        if not os.path.isdir(self.col_cluster_directory):
            os.mkdir(self.col_cluster_directory)

    def create_rowclusterM(self):
        m_wclusterassign = pd.merge(self.m,self.clusters,left_index=True,right_on='row',how='left')
        m_wclusterassign.set_index('row',inplace=True)
        M = []
        for cluster in m_wclusterassign.cluster.unique():
            M.append(m_wclusterassign[m_wclusterassign['cluster']==cluster].drop('cluster',axis=1).values)
        return M
            
class initial_column_clusters(column_clusters):
    
    def __init__(self,k,min_size,num_rand):
        column_clusters.__init__(self,k,min_size,num_rand)
        self.col_cluster_directory = self.col_cluster_directory+'initial_clusters/'
        self.col_cluster_file = self.col_cluster_directory+str(self.num_clusters) + '_initialcolclusters.xlsx'
        nxG = nx.from_edgelist(self.g.get_edgelist(),nx.DiGraph())
        self.nxG_map = dict(zip(nxG.nodes(),self.g.vs['name']))
        self.nxG = nx.relabel_nodes(nxG,self.nxG_map)
        self.results = pd.DataFrame(data = None, columns = ['initial_cluster_assign','initial_cut_edges','initial_cluster_residual','generated'])
        if not os.path.isdir(self.col_cluster_directory):
            os.mkdir(self.col_cluster_directory)

    def add_rand(self):
        treecluster = tree_cluster(self.g,self.create_rowclusterM(),self.num_clusters)
        for n in range(self.num_rand):
            treecluster.initialize_components()
            newrow = pd.DataFrame([treecluster.compute_residual2()],columns=['initial_cluster_residual'])
            newrow['initial_cluster_assign'] = [dict(zip(treecluster.g.vs['name'],treecluster.get_assignments()))]
            newrow['initial_cut_edges'] = [treecluster.cut_edges]
            newrow['generated'] = 'random'
            self.results = self.results.append(newrow[['initial_cluster_assign','initial_cut_edges','initial_cluster_residual','generated']])
    
    def load_column_clusters(self,overwrite=False):
        self.get_qualified_cut_vertices()
        if (os.path.isfile(self.col_cluster_file)) & (overwrite==False):
            df = pd.read_excel(self.col_cluster_file)
            df['initial_cluster_assign'] = [ast.literal_eval(val) for val in df['initial_cluster_assign'].values]
            df['initial_cut_edges'] = [ast.literal_eval(val) for val in df['initial_cut_edges'].values]
            self.results = df
        else:
            self.score_cut_combinations()
            self.add_rand()
            self.save_results()
            
    def get_qualified_cut_vertices(self):
        def criteria(G,node):
            num_direct_inedges = len(G.in_edges(node))
            num_indirect_outedges = len(nx.nodes(nx.dfs_tree(G, node)))-1
            if num_direct_inedges==0:
                return True
            parent = list(G.predecessors(node))[0]
            num_direct_outedges_parent = len(G.out_edges(parent))
            if (num_indirect_outedges>1) & (num_direct_outedges_parent > 1):
                return True
            return False
        G = self.nxG
        qualified_cut_vertices = []
        for node in G.nodes():
            if criteria(G,node):
                qualified_cut_vertices.append(node)
        self.qualified_cut_vertices = qualified_cut_vertices
        
    def score_cut_combinations(self,save=False):
        def convertToDict(clustering):
            returndict = {}
            for cluster in range(len(clustering)):
                for node in clustering[cluster]:
                    returndict.update({node:cluster})
            return returndict
        cut_vertices = self.qualified_cut_vertices.copy()
        cut_vertices.remove(self.rootnode)
        inv_nxG_map = {v: k for k, v in self.nxG_map.items()}
        combinations = [list(combination) for combination in itertools.combinations(cut_vertices,self.num_clusters-1)]
        cut_edges = [[[(u,inv_nxG_map.get(v))] for v in comb for u in self.g.predecessors(inv_nxG_map.get(v))] for comb in combinations]
        combinations = [[self.rootnode] + comb for comb in combinations]
        m_wclusterassign = pd.merge(self.m,self.clusters,left_index=True,right_on='row',how='left')
        m_wclusterassign.set_index('row',inplace=True)
        m_wclusterassign.index = m_wclusterassign.index.rename('index')
        p = multiprocessing.Pool(processes=multiprocessing.cpu_count()-1)
        results = p.starmap(clusterScores,zip(combinations,itertools.repeat(self.nxG),itertools.repeat(m_wclusterassign)))   
        results = pd.DataFrame(results, columns = ['initial_cluster_residual', 'initial_cluster_assign'])
        results['generated'] = 'selected'
        results['initial_cluster_assign'] = results['initial_cluster_assign'].map(lambda x: convertToDict(x))
        results['initial_cut_edges'] = cut_edges
        self.results = self.results.append(results[['initial_cluster_residual','initial_cut_edges','initial_cluster_assign','generated']],ignore_index=True)
        
    def save_results(self):
        self.results.to_excel(self.col_cluster_file,index=False)
        
class optimal_column_clusters(column_clusters):
    
    def __init__(self,k,min_size,num_rand, n = 100, num_iters = 1):
        column_clusters.__init__(self,k,min_size,num_rand)
        self.col_cluster_directory = self.col_cluster_directory + 'optimal_clusters/'
        self.col_cluster_file = self.col_cluster_directory+str(self.num_clusters) + '_optimalcolclusters.xlsx'
        self.num_iters = num_iters
        self.n = n
        if not os.path.isdir(self.col_cluster_directory):
            os.mkdir(self.col_cluster_directory)

    def load_column_clusters(self, overwrite=False, plot=False):
        if (os.path.isfile(self.col_cluster_file)) & (overwrite==False):
            df = pd.read_excel(self.col_cluster_file)
            df['initial_cluster_assign'] = [ast.literal_eval(val) for val in df['initial_cluster_assign'].values]
            df['initial_cut_edges'] = [ast.literal_eval(val) for val in df['initial_cut_edges'].values]
            self.optimal_results = df
        else:
            icc = initial_column_clusters(self.num_clusters,self.min_size,self.num_rand)
            icc.load_column_clusters(overwrite=overwrite)
            self.icc = icc
            self.optimal_col_cluster()
            self.save_results()
        if plot:
            self.plot_optimal()

    def optimal_col_cluster(self):    
        def apply_optimal(row):
            treecluster = tree_cluster(self.g,self.create_rowclusterM(),self.num_clusters)
            min_residual = np.inf
            for interation in range(self.num_iters):
                treecluster.initialize_components(initial_cutedges=row['initial_cut_edges'])
                treecluster.fit()
                residual = treecluster.compute_residual2()
                if residual<min_residual:
                    assignments = dict(zip(treecluster.g.vs['name'],treecluster.get_assignments()))
                    cut_edges = treecluster.cut_edges
                    min_residual = residual
            return assignments, cut_edges, min_residual
        
        initial = self.icc.results
        initial.sort_values(by=['initial_cluster_residual'],ascending=True,inplace=True)
        initial=initial[:self.n]
        
        pandarallel.initialize()
        
        initial[['optimal_cluster_assign','optimal_cut_edges','optimal_residual']] = initial.parallel_apply(apply_optimal,axis=1, result_type="expand")
        
        self.optimal_results = initial
       
    def save_results(self):
        self.optimal_results.to_excel(self.col_cluster_file,index=False)
            
    def plot_optimal(self):
        results = self.optimal_results
        results.sort_values(by=['optimal_residual'],ascending=True,inplace=True)
        optimal_tree = tree_cluster(self.g,self.create_rowclusterM(),self.num_clusters)
        optimal_tree.initialize_components(initial_cutedges=results.iloc[0]['optimal_cut_edges'])
        optimal_tree.treeplot(savepath=self.col_cluster_directory + str(self.num_clusters) + '_optimal_tree.pdf')