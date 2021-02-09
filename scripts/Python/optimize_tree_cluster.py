import os, sys
import pandas as pd
import networkx as nx
import numpy as np
import multiprocessing
import itertools
import ast
import pdb

import tree_cluster
import utilities
import Yoshida
import peak_clusters
import master_peaks as mp
import configuration as conf
import utilities

class initial_clusters:
    
    col_cluster_directory = conf.DATA_DIR + 'col_clusters/'
    
    def __init__(self, k, g, peak_cluster, 
                 num_row_clusters=20, 
                 nCPU=4):
        
        if not os.path.isdir(self.col_cluster_directory):
            os.mkdir(self.col_cluster_directory)
            
        self.col_cluster_file = self.col_cluster_directory + str(k) + '_initial.csv'
            
        self.nCPU = nCPU
        self.num_clusters = k
        
        self.g = g
        
        # form matrices for each row cluster, matching columns to vertex names
        g_names = self.g.vs["name"]
        pc_names = peak_cluster.get_cell_types()
        if not set(g_names) == set(pc_names):
            sys.exit("tree vertex names and peak cluster cell types must match")
        match_ind = utilities.match(g_names, pc_names)
        m_list = []
        for i in range(num_row_clusters+1):
            m = peak_cluster.form_cluster_matrix(i)[:,match_ind]
            m_list.append(m)
        self.m_list = m_list
        
        self.tree_cluster = tree_cluster.tree_cluster(self.g, m_list, self.num_clusters)
        self.rootnode = 0
                             
        # the following are only used in create_cut_combinations
        self.nxG = nx.from_edgelist(self.g.get_edgelist(),nx.DiGraph())
        
        self.qualified_cut_vertices = self.identify_qualified_cut_vertices()
        
    #  vertices that are good candidates for cuts
    def identify_qualified_cut_vertices(self):
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
                
        return qualified_cut_vertices
       

    # create a clustering based on random selection of vertices
    def create_random_clusters(self):
        treecluster = self.tree_cluster
        treecluster.initialize_components()
        return {'assignments':treecluster.get_assignments(),
                'residual':treecluster.compute_residual2(),
                'generated':"random"}
             
    #  form combinations of k vertices from the qualified vertices
    def create_cut_combinations(self):
 
        cut_vertices = self.qualified_cut_vertices.copy()
        cut_vertices.remove(self.rootnode)
        
        combinations = [list(combination) for combination in itertools.combinations(cut_vertices,self.num_clusters-1)]
        cutedge_combinations = [[[(u,v)] for v in comb for u in self.g.predecessors(v)] for comb in combinations]
        
        return cutedge_combinations
    
    def create_cut_combination_clusters(self, cutedge_combs):
  
        m_list = self.m_list.copy()
        p = multiprocessing.Pool(processes=self.nCPU)
        packed_results = p.starmap(computeCluster,
                                  zip(cutedge_combs,
                                      itertools.repeat(self.g),
                                      itertools.repeat(m_list),
                                      itertools.repeat(self.num_clusters)))   
        
        # need to unpack results...
        results = []
        for pr in packed_results:       
         
          residual2 = pr[0]
          assignments = pr[1]
            
          # debug
          if np.any(assignments == -1):
              sys.exit("a column has not been assigned to a cluster!")
             
          results.append({'assignments':assignments,
                          'residual':residual2,
                          'generated':"cut"})
          
        return results
    
    # generate initial clusterings    
    def generate_initial_clusters(self, nrand):
        self.results = []
        for i in range(nrand):
          self.results.append(self.create_random_clusters())
         
        cutedge_combs = self.create_cut_combinations()
        self.results = self.results \
                       + self.create_cut_combination_clusters(cutedge_combs)
       
        
      
    # load clusterings from xlsx file
    def load_results(self):
        return pd.read_csv(self.col_cluster_file,
                           sep=",")
      
  
    def save_results(self):
        assign_m = np.array([r["assignments"] for r in self.results])
        res = [r["residual"] for r in self.results]
        gen = [r["generated"] for r in self.results]
        
        tb = pd.DataFrame(assign_m,
                          columns=self.g.vs["name"])
        tb.insert(tb.shape[1], "residual", res)
        tb.insert(tb.shape[1], "generated", gen)
        
        tb.to_csv(self.col_cluster_file,
                  sep=",",
                  index=False)

class optimal_clusters:
    
    col_cluster_directory = conf.DATA_DIR + 'col_clusters/'
    
    def __init__(self, k, g, peak_cluster,
                 min_n = 100,
                 num_row_clusters=20, 
                 nCPU=4):
        
        if not os.path.isdir(self.col_cluster_directory):
            os.mkdir(self.col_cluster_directory)
        
        if not os.path.isfile(self.col_cluster_directory + str(k) + '_initial.csv'):
            sys.exit("create intiial_clusters first!")
        
        self.initial_clusters = initial_clusters(k, g, peak_cluster, num_row_clusters,nCPU)    
        self.initial_df = self.initial_clusters.load_results()
        
        self.col_cluster_file = self.col_cluster_directory + str(k) + '_optimal.csv'
        
        self.nCPU = nCPU
        self.num_clusters = k
        self.n = min_n
        
        self.g = g
        
        #form matrices for each row cluster, matching columns to vertex names
        g_names = self.g.vs["name"]
        pc_names = peak_cluster.get_cell_types()
        if not set(g_names) == set(pc_names):
            sys.exit("tree vertex names and peak cluster cell types must match")
        match_ind = utilities.match(g_names, pc_names)
        m_list = []
        for i in range(num_row_clusters+1):
            m = peak_cluster.form_cluster_matrix(i)[:,match_ind]
            m_list.append(m)
        self.m_list = m_list
        
    def optimal_col_cluster(self,different = False):    
        
        #Return indexes with minimum residual and sufficiently different cut edges, fill remainder with minimum residuals if not enough indexes with sufficiently different cut edges
        def get_min_diffindexes(df):
            def get_cut_edges(row):
                assignments = row[self.g.vs['name']].values
                cut_edges = []
                for cluster in range(1,self.num_clusters):
                    in_cluster = [vertex for vertex, a in enumerate(assignments) if a==cluster]
                    cut_edges_cluster = [(u,v) for v in in_cluster for u in self.g.predecessors(v) if u not in in_cluster]
                    cut_edges = cut_edges + cut_edges_cluster
                return cut_edges
            def is_different(cut_one,cut_two,diff_threshold):
                return len(np.intersect1d(cut_one,cut_two)) >= diff_threshold
                
            n = self.n
            k = self.num_clusters
            
            #I use the function ((x-3)^2)/14+(x-3)+2 because it has a reasonable curve for the domain [2,18]
            diff_threshold = min(1,int(2+(k-3)-((k-3)**2)/14))
            df['cut_edges'] = df.apply(lambda row: get_cut_edges(row), axis =1)
            
            index_list = [0]
            index =1
            while (index < df.shape[0]) & (len(index_list)<n):
                add=True
                index_cutedges = df.loc[index]['cut_edges']
                for index_inlist in index_list:
                    if not is_different(index_cutedges,df.loc[index_inlist]['cut_edges'],diff_threshold):
                        add=False
                if add:
                    index_list.append(index)
            
            #fill remaining indexes with min residuals that did not qualify as different
            i=1        
            while len(index_list)<n:
                if i not in index_list:
                    index_list.append(i)
                i+=1
            return index_list

        initial = self.initial_df
        initial.sort_values(by=['residual'],ascending=True,inplace=True)
        
        #subset dataframe to min residual n initial assignments
        if different:
            initial = initial.loc[get_min_diffindexes(initial)]
        else:
            initial=initial[:self.n]
            
        initial_assignments = initial[self.g.vs['name']].values
        m_list = self.m_list.copy()
        p = multiprocessing.Pool(processes=self.nCPU)
        packed_results = p.starmap(optimizeCluster,
                                  zip(initial_assignments,
                                      itertools.repeat(self.g),
                                      itertools.repeat(m_list),
                                      itertools.repeat(self.num_clusters)))  
        # need to unpack results...
        results = []
        for pr in packed_results:       
         
          residual2 = pr[0]
          assignments = pr[1]
            
          # debug
          if np.any(assignments == -1):
              sys.exit("a column has not been assigned to a cluster!")
             
          results.append({'assignments':assignments,
                          'residual':residual2})
        
        self.results = results
           
    def load_results(self):
        return pd.read_csv(self.col_cluster_file,
                           sep=",")    

    def save_results(self):
        assign_m = np.array([r["assignments"] for r in self.results])
        res = [r["residual"] for r in self.results]
        
        tb = pd.DataFrame(assign_m,
                          columns=self.g.vs["name"])
        tb.insert(tb.shape[1], "residual", res)
        
        tb.to_csv(self.col_cluster_file,
                  sep=",",
                  index=False)    
        
######################################################
# computational methods
def computeCluster(cut_edges,g,m_list,num_clusters):
    tc = tree_cluster.tree_cluster(g, m_list, num_clusters)
    tc.initialize_components(initial_cutedges = cut_edges)
    return tc.compute_residual2, tc.assignments

def optimizeCluster(assignments,g,m_list,num_clusters):
    tc = tree_cluster.tree_cluster(g, m_list, num_clusters)
    tc.initialize_components(assignments = assignments)
    tc.optimize()
    return tc.compute_residual2(), tc.assignments