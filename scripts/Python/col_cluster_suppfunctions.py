import pandas as pd
import numpy as np

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