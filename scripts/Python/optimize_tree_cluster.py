import sys, os
import pandas as pd
import networkx as nx
import numpy as np
import multiprocessing
import itertools
import pdb
import seaborn as sns
import copy

import tree_cluster as tclust
import utilities

#  Optimize tree clustering
#
# public methods
# generate_starting_points : generate random and smart starting points
# optimize : runs an optimization using each starting point
# save : saves either starting points or fits to a file
# get_starting_points
# get_fits
#
# fits and starting points are given as a list of dictionaries.  
# Dictionaries include "assignments" (1d array giving clustering),
# residual, generated (either "random", "smart" or "optimization")
#
# @ param k number of clusters in tree
# @param g igraph object, must be a tree with one root
# @param m_list list of binary matrices to be column 
# @param starting_point_file input file for starting
# points for optimization.  The file should be created
# by the save method after using generate_initial_clusters.
# If None, then generate_initial_clusters must be called
# prior to optimization
class optimize_tree_cluster:
    
    def __init__(self, k, g, m_list,
                 starting_points_file = None,
                 nCPU=1):
        
        # check that g is a tree
        np.random.seed(123)
        istree, root_node = utilities.find_root(g)
        if not istree:
            sys.exit("g does not describe a tree")
        self.rootnode = root_node
        
        self.nCPU = nCPU
        self.num_clusters = k
        self.ncol = m_list[0].shape[1]
        
        self.g = g
        self.m_list = m_list
        self.starting_points = None
        
        if not starting_points_file is None:
            self.starting_points = self.load(starting_points_file)
        else:
            self.starting_points = None
        self.fits = None
        
        self.tree_cluster = tclust.tree_cluster(self.num_clusters,
                                                      self.g, 
                                                      m_list)
                             
        # the following are only used in create_cut_combinations
        self.nxG = nx.from_edgelist(self.g.get_edgelist(),nx.DiGraph())
        self.qualified_cut_vertices = self.identify_qualified_cut_vertices()
        
        
    # constructor methods
        
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
    
    def choose_best_initial_clusters(self, results, n):
        
        def is_different(cut_one,cut_two,diff_threshold):
                return len(np.intersect1d(cut_one,cut_two)) >= diff_threshold
            
        if len(results) == 0:
            return None
        
        if n >= len(results):
            return np.arange(len(results)).tolist()
            
        assignments = np.array([r["assignments"] for r in results])
        res = np.array([r["residual"] for r in results])
        assignments = assignments[np.argsort(-1*res),:]
   
        #I use the function ((x-3)^2)/14+(x-3)+2 because it has a 
        # reasonable curve for the domain [2,18]
        k = self.num_clusters
        diff_threshold = min(1,int(2+(k-3)-((k-3)**2)/14))
        cut_edge_list = []
        for i in range(assignments.shape[0]):
          cut_edge_list.append(self.tree_cluster.assignments2cutedges(assignments[i,:]))
            
     
        index_list = [0]
        index = 1
        while (index < assignments.shape[0]) & (len(index_list)<n):
            add=True
            index_cutedges = cut_edge_list[index]
            for index_inlist in index_list:
              if not is_different(index_cutedges, 
                                  cut_edge_list[index_inlist],
                                  diff_threshold):
                  add=False
                  break
            if add:
                index_list.append(index)
            
            #fill remaining indexes with min residuals that did not qualify as different
            i=1        
            while len(index_list)<n:
                if i not in index_list:
                    index_list.append(i)
                i+=1
            
        return index_list
 
    # PUBLIC METHODS
    
    # computation methods
    def generate_starting_points(self, n_random, n_smart):
        rand_results = []
        for i in range(n_random):
          rand_results.append(self.create_random_clusters())
         
        if n_smart > 0:
         cutedge_combs = self.create_cut_combinations()
         combo_results = self.create_cut_combination_clusters(cutedge_combs)
         best_combo_results = self.choose_best_initial_clusters(combo_results,
                                                               n_smart)
         combo_results_filt = [combo_results[i] for i in best_combo_results]
        else:
         combo_results_filt = []
        
        self.starting_points = rand_results + combo_results_filt
        self.fits = None
             
    def optimize(self):    
            
        ic = self.get_starting_points()
        initial_assignments = np.array([r["assignments"] for r in ic])
        m_list = self.m_list
        
        packed_results = [] 
        if self.nCPU == 1:
          print("serial optimization...")
          for i in range(initial_assignments.shape[0]):
            packed_results.append(optimizeCluster(initial_assignments[0,:], 
                        self.g, self.m_list, self.num_clusters))
        else:
          print("parallelizing optimization...")
          p = multiprocessing.Pool(processes=self.nCPU)
          gl = [copy.deepcopy(self.g) for i in range(len(ic))]
          m_list_l = [copy.deepcopy(m_list) for i in range(len(ic))]
          packed_results = p.starmap(optimizeCluster,
                                     zip(initial_assignments,
                                         gl,
                                         m_list_l,
                                         np.repeat(self.num_clusters,len(ic))))
                                  # zip(initial_assignments,
                                  #     itertools.repeat(self.g.deepcopy),
                                  #     itertools.repeat(m_list.deepcopy),
                                  #     itertools.repeat(self.num_clusters)))  
        
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
                           'generated':'optimization'})
        
        self.fits = results
        
  
    # accessor methods
        
    def get_starting_points(self, as_dataFrame=False):
        if  not as_dataFrame:
          out =  self.starting_points
        else:
          temp_file = "temp_sp_" + str(np.random.uniform()) + ".csv"
          self.save(temp_file, starting_points=True)
          out = pd.read_csv(temp_file, sep=",")
          os.remove(temp_file)
          
        return out
    
    def get_fits(self, as_dataFrame=False):
        if  not as_dataFrame:
          out =  self.fits
        else:
          temp_file = "temp_sp_" + str(np.random.uniform()) + ".csv"
          self.save(temp_file, fits=True)
          out = pd.read_csv(temp_file, sep=",")
          os.remove(temp_file)
          
        return out
    
    def load(self, infile):
        tb = pd.read_csv(infile, sep=",")
     
        res = tb["residual"].tolist()
        gen = tb["generated"].tolist()
        
        tb = tb.drop("residual", axis=1)
        tb = tb.drop("generated", axis=1)
        
        m = tb.values
        return [{'assignments':m[i,:], 
                'residual':res[i],
                'generated':gen[i]} for i in range(len(tb))]
        
        
    # save either starting points or fits depending on flags
    def save(self, outfile, starting_points=False, fits=False):
        if not starting_points and not fits:
            sys.exit("one of starting ponts or fits flags must be set True")
            
        if starting_points:
            results = self.get_starting_points()
        else:
            results = self.get_fits()
        assign_m = np.array([r["assignments"] for r in results])
        res = [r["residual"] for r in results]
        gen = [r["generated"] for r in results]
        
        tb = pd.DataFrame(assign_m,
                          columns=self.g.vs["name"])
        tb.insert(tb.shape[1], "residual", res)
        tb.insert(tb.shape[1], "generated", gen)
        
        tb.to_csv(outfile,
                  sep=",",
                  index=False)
        
    def heatmap(self, starting_points=False, fits=False):
        if not starting_points and not fits:
            sys.exit("one of starting ponts or fits flags must be set True")
         
        if starting_points:
            results = self.get_starting_points()
        else:
            results = self.get_fits()
            
        a = np.array([r["assignments"] for r in results])
        tb = pd.DataFrame(a)
 
        sns.heatmap(tb)
      

   
######################################################
# computational methods
def computeCluster(cut_edges,g,m_list,num_clusters):
    tc = tclust.tree_cluster(num_clusters, g, m_list)
    tc.initialize_components(initial_cutedges = cut_edges)
    return tc.compute_residual2(), tc.assignments

def optimizeCluster(assignments,g,m_list,num_clusters):
    tc = tclust.tree_cluster(num_clusters, g, m_list)
    tc.optimize(assignments = assignments)
    return tc.compute_residual2(), tc.assignments