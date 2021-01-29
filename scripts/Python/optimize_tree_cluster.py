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
    col_cluster_file = col_cluster_directory+'initial.csv'
    
    def __init__(self, k, g, peak_cluster, 
                 num_row_clusters=20, 
                 nCPU=4):
        
        if not os.path.isdir(self.col_cluster_directory):
            os.mkdir(self.col_cluster_directory)
            
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
        cluster_list = []
        for i in range(num_row_clusters+1):
            m = peak_cluster.form_cluster_matrix(i)[:,match_ind]
            m_list.append(m) 
            cluster_list.append(np.repeat(i, len(m)))
        
        self.tree_cluster = tree_cluster(self.g, m_list, self.num_clusters)
        self.rootnode = Yoshida.Yoshida_tree().get_root()
                             
        # the following are only used in create_cut_combinations
        nxG = nx.from_edgelist(self.g.get_edgelist(),nx.DiGraph())
        self.nxG_map = dict(zip(nxG.nodes(),self.g.vs['name']))
        self.nxG = nx.relabel_nodes(nxG,self.nxG_map)
        
        # create combined data matrix with cluster/row columns
        M = np.concatenate(m_list, axis=0)
        clusters = np.concatenate(cluster_list)
        # no idea why all this stuff is needed!
        mwa = pd.DataFrame(M, columns=g_names)
        mwa.insert(mwa.shape[1], "cluster", clusters)
        mwa.insert(mwa.shape[1], "row", np.arange(len(mwa)))
        mwa.set_index('row',inplace=True)
        mwa.index = mwa.index.rename('index')
        self.mwa = mwa
        
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
        combinations = [[self.rootnode] + comb for comb in combinations]
        
        return combinations
    
    def create_cut_combination_clusters(self, combinations):
  
        # debug
        mwa = self.mwa.copy()
        p = multiprocessing.Pool(processes=self.nCPU)
        packed_results = p.starmap(clusterScores,
                                  zip(combinations,
                                      itertools.repeat(self.nxG),
                                      itertools.repeat(mwa)))   
        
        # need to unpack results...
        vertex_names = self.g.vs["name"]
        assignments = np.repeat(-1, len(self.g.vs))
        results = []
        for pr in packed_results:       
         
          residual2 = pr[0]
          clusters = pr[1]
          assignments[:] = -1
 
          for i in range(len(clusters)):
            elements = clusters[i]
            assignments[utilities.match(elements, vertex_names)] = i
            
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
         
        combs = self.create_cut_combinations()
        self.results = self.results \
                       + self.create_cut_combination_clusters(combs)
       
        
      
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

######################################################
# computational methods
def optimizeColClusters(data,col_clusters):
 
    M = data.copy()
    for col_cluster in col_clusters:
        col_cluster = list(col_cluster)
        cluster_mean = pd.merge(data.reset_index()[['index','cluster']],
                                pd.DataFrame(data[col_cluster + ['cluster']].groupby('cluster').mean().mean(axis=1)),
                                left_on='cluster',right_on='cluster',how='left')
        cluster_mean.set_index('index',inplace=True)
        for col in col_cluster:
            M[col] = cluster_mean[0]
    return np.linalg.norm(M.drop('cluster',axis=1).sub(data.drop('cluster',axis=1)))**2 

def createColClusters(G,cut_vertices):
    def descendTree(vertex,G,other_cutvertices):
        if vertex in other_cutvertices:
            return []
        toreturn = [vertex]
        if len(G.out_edges(vertex))==0:
            return toreturn
        children = list(G.successors(vertex))
        for child in children:
            toreturn += descendTree(child,G,other_cutvertices)
        return toreturn
    col_clusters = []
    for vertex in cut_vertices:
        other_cutvertices = cut_vertices.copy()
        other_cutvertices.remove(vertex)
        col_clusters.append(descendTree(vertex,G,other_cutvertices))
    return col_clusters

def clusterScores(cut_vertices,G,m):
    col_clusters = createColClusters(G, cut_vertices) 
    return optimizeColClusters(m, col_clusters), col_clusters