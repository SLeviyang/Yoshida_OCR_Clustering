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
import pickle
import ast
from col_cluster_suppfunctions import clusterScores
from tree_cluster import tree_cluster

class column_clusters:
    col_cluster_directory = conf.DATA_DIR + 'col_clusters/'
    
    def __init__(self,k,min_size):
        self.num_clusters = k
        self.col_cluster_file = self.col_cluster_directory +str(k) + '_colclusters.xlsx'
        self.g = Yoshida.Yoshida_tree().load_igraph()
        nxG = nx.from_edgelist(self.g.get_edgelist(),nx.DiGraph())
        self.nxG = nx.relabel_nodes(nxG,dict(zip(nxG.nodes(),self.g.vs['name'])))
        self.rootnode = self.g.vs['name'][0]
        self.clusters = peak_clusters.peak_clusters().clusters
        master_peaks = mp.master_peaks()
        self.m = pd.DataFrame(master_peaks.load_matrix().astype("float"),columns=master_peaks.get_cell_types()).loc[self.clusters.groupby('cluster').filter(lambda x: len(x) > min_size)['row']]
        if not os.path.isdir(self.col_cluster_directory):
            os.mkdir(self.col_cluster_directory)
    
    def load_column_clusters(self):
        self.get_qualified_cut_vertices()
        if os.path.isfile(self.col_cluster_file):
            df = pd.read_excel(self.col_cluster_file,index_col=0)
            self.col_clusters = [ast.literal_eval(val) for val in df['col_clusters'].values]
            self.col_cluster_scores = list(df['cluster_score'].values)
            self.cutvertex_combinations = [ast.literal_eval(val) for val in df['cut_vertices'].values]
        else:
            self.score_cut_combinations(save=True)
            
    def load_optimal_colcluster(self, n=5, plot=False):
        file = self.col_cluster_directory + str(self.num_clusters) + '_optimal_tree.pkl'
        if os.path.isfile(file):
            openfile = open(file,'rb')
            self.optimal_tree = pickle.load(openfile)
            openfile.close()
        else:
            self.optimal_col_cluster(n=n,plot=plot)
            
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
        cut_vertices = self.qualified_cut_vertices.copy()
        cut_vertices.remove(self.rootnode)
        combinations = [[self.rootnode] + list(combination) for combination in itertools.combinations(cut_vertices,self.num_clusters-1)]
        m_wclusterassign = pd.merge(self.m,self.clusters,left_index=True,right_on='row',how='left')
        m_wclusterassign.set_index('row',index=True)
        p = multiprocessing.Pool(processes=multiprocessing.cpu_count()-1)
        results = p.starmap(clusterScores,zip(combinations,itertools.repeat(self.nxG),itertools.repeat(m_wclusterassign)))   
        results = pd.DataFrame(results, columns = ['cluster_score', 'col_clusters'])
        self.col_clusters = list(results['col_clusters'].values)
        self.col_cluster_scores = list(results['cluster_score'].values)
        self.cutvertex_combinations = combinations
        
        if save:
            results['cut_vertices'] = combinations
            results.to_excel(self.col_cluster_file)
    
    def optimal_col_cluster(self, bot_pctile = .15, diff_threshold = .70, n = 1, plot=False):
        def get_min_cutvertices(scores,cutvertices, bot_pctile, diff_threshold,n):
            scores, cutvertices = [list(item) for item in zip(*sorted(zip(scores,cutvertices)))]
            threshold_value = np.percentile(scores,bot_pctile)
            cutvertices_consider = [cutvertices[i] for i in range(len(scores)) if scores[i]<=threshold_value]
            return_cutvertices = [cutvertices[0]]
            j = 1
            while (j<len(cutvertices_consider)) & (len(return_cutvertices)<n):
                cut_consider = cutvertices[j]
                add = True
                for cut in return_cutvertices:
                    #minus one in denominator to account for rootnode found in all
                    if len(np.setdiff1d(cut_consider,cut))/(len(cut_consider)-1) < diff_threshold:
                        add = False
                if add:
                    return_cutvertices.append(cut_consider)
                j+=1
            return return_cutvertices
        
        def create_rowclusterM(m_wclusterassign):
            M = []
            for cluster in m_wclusterassign.cluster.unique():
                M.append(m_wclusterassign[m_wclusterassign['cluster']==cluster].drop('cluster',axis=1).values)
            return M
        self.load_column_clusters()
        min_cutvertices = get_min_cutvertices(self.col_cluster_scores,self.cutvertex_combinations,bot_pctile,diff_threshold,n)
        m_wclusterassign = pd.merge(self.m,self.clusters,left_index=True,right_on='row',how='left')
        m_wclusterassign.set_index('row',inplace=True)
        rowclusterM = create_rowclusterM(m_wclusterassign)
        
        iterations = 5
        min_residual = np.inf
        
        for cut in min_cutvertices:
            for iteration in range(iterations):
                treecluster = tree_cluster(self.g,rowclusterM,self.num_clusters)
                treecluster.initialize_components(initial_cutvertices=[vertex for vertex in cut if vertex != self.rootnode])
                treecluster.fit()
                residual = treecluster.compute_residual2()
                if residual<min_residual:
                    min_residual = residual
                    optimal_tree = tree_cluster
                
        self.optimal_residual = min_residual
        self.optimal_tree = optimal_tree
        
        if plot:
            optimal_tree.treeplot(savepath=self.col_cluster_directory + str(self.num_clusters) + '_optimal_tree.pdf')
        
        file = open(self.col_cluster_directory + str(self.num_clusters) + '_optimal_tree.pkl','wb')
        pickle.dump(optimal_tree,file)
        file.close()    