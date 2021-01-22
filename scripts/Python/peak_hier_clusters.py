import os
import pandas as pd
import numpy as np
import igraph as ig
from sklearn.cluster import AgglomerativeClustering
import master_peaks as mp
import peak_clusters
import configuration as conf
import Yoshida

class hier_tree:
    
    def __init__(self):
        self.g = Yoshida.Yoshida_tree().load_igraph()
    
    def get_vertices(self):
        return self.g.vs['name']
        
    def treeplot(self, assignment_map, outpath, vertex_label=False, m_index=None):
        g = self.g
        vertex_labels = self.get_vertices()
        vs = {}
        
        if assignment_map is not None:
            pal = ig.drawing.colors.ClusterColoringPalette(max(assignment_map.values())+1)
            vs["vertex_color"] = pal.get_many([assignment_map.get(vertex) for vertex in vertex_labels])
                  
        vs["bbox"] = (1200, 1000)
        if m_index is None:
            vs["vertex_size"] = 20
        else:
            vs["vertex_size"] = 30*np.mean(self.m[m_index], 0) + 1
        vs["vertex_label_size"] = 20
        if vertex_label:
            vs["vertex_label"] = vertex_labels
        else:
            vs['vertex_label'] = [str(i) for i in range(g.vcount())]
        vs["vertex_label_dist"] = 1.5
           
        layout = g.layout_reingold_tilford()
     
        if outpath == '':
            pl = ig.plot(g, layout=layout, **vs)
            pl.show()
        else:
            ig.plot(g,outpath, layout=layout, **vs)
    
    
class peak_hier_clusters:
    clusters = peak_clusters.peak_clusters().clusters
    m = pd.DataFrame(mp.master_peaks().load_matrix().astype("float")).loc[clusters['row']].values
    cell_types = mp.master_peaks().get_cell_types()
    peak_hier_cluster_directory = conf.DATA_DIR + "peak_hier_clusters/"
        
    def __init__(self,n_clusters):
        self.hier_model = AgglomerativeClustering(n_clusters=n_clusters,linkage='complete')
        if not os.path.isdir(self.peak_hier_cluster_directory):
            os.mkdir(self.peak_hier_cluster_directory)

    def compute_clusters(self):
        model = self.hier_model
        clustering = model.fit(self.m.T)
        self.assignments = clustering.labels_
        
    def plot_clusters(self,save=False):
        tree = hier_tree()
        assignment_map = dict(zip(self.cell_types, self.assignments))
        tree.treeplot(assignment_map,self.peak_hier_cluster_directory + str(max(assignment_map.values())+1) + '_clusters.pdf' if save else '')