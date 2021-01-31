import os, sys
import pandas as pd
import networkx as nx
import numpy as np
from pandarallel import pandarallel
import multiprocessing
import itertools
import ast
import pdb

from tree_cluster import tree_cluster
import utilities
import Yoshida
import peak_clusters
import master_peaks as mp
import configuration as conf

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
        
        self.g = Yoshida.Yoshida_tree().load_igraph()
        
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
        
        self.tree_cluster = tree_cluster(self.g, m_list, self.num_clusters)
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
                'cut_edges':treecluster.cut_edges,
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
        i = 0
        for pr in packed_results:       
         
          residual2 = pr[0]
          assignments = pr[1]
            
          # debug
          if np.any(assignments == -1):
              sys.exit("a column has not been assigned to a cluster!")
             
          results.append({'assignments':assignments,
                          'residual':residual2,
                          'cut_edges':cutedge_combs[i],
                          'generated':"cut"})
          i+=1
          
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
        cut = [r["cut_edges"] for r in self.results]
        gen = [r["generated"] for r in self.results]
        
        tb = pd.DataFrame(assign_m,
                          columns=self.g.vs["name"])
        tb.insert(tb.shape[1], "cut_edges", cut)
        tb.insert(tb.shape[1], "residual", res)
        tb.insert(tb.shape[1], "generated", gen)
        
        tb.to_csv(self.col_cluster_file,
                  sep=",",
                  index=False)

######################################################
# computational methods
def computeCluster(cut_edges,g,m_list,num_clusters):
    tc = tree_cluster(g, m_list, num_clusters)
    tc.initialize_components(cut_edges)
    return tc.compute_residual2, tc.assignments