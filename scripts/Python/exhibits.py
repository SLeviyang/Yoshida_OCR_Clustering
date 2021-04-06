#This file contains the code for exhibits 3-7

import seaborn as sns
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import matplotlib as mpl
import configuration as conf
    
#Figure 5: Percent of variance retained by the approximation methods
def pctEnergyCharacterizedAproxx(df):
    
    #Font to match latex
    mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']},size=10)
    
    #Define figure size
    fig, ax = plt.subplots(figsize=(6,4))
    
    #Add series to the plot
    plt.plot(df.index,df['Means'],color = 'black', label = 'Means approximation')
    plt.plot(df.index,df['Binary'],color = 'red', label = 'Open/Closed approximation')
    
    #Legend location and axis labels
    plt.legend(loc = 'upper left', frameon=False)
    plt.ylabel('Percent of Across Column Variance Captured')
    plt.xlabel('Number of Cell Type Clusters')
    
    #Remove top bar from graph
    ax.tick_params(labelright=True,right=True)
    ax.spines['top'].set_visible(False)
    
    #Set x-axis ticks
    xticks = np.arange(2,13)
    locator = mpl.ticker.FixedLocator(xticks)
    ax.xaxis.set_major_locator(locator)
    
    #Set axis limits
    ax.set_xlim(1,13)
    ax.set_ylim(0,100)
    
    #Save image
    fig.savefig(conf.DATA_DIR + '/pctenergycharacterized_approx.png', dpi=400, bbox_inches='tight', pad_inches=0.01)

#Figure 6: Percent of variance across clusters for the three different clustering methods
def pctVarAcrossClusters(df):
    
    #Font to match latex
    mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']},size=10)
    
    #Define figure size
    fig, ax = plt.subplots(figsize=(6,4))
    
    #Add series to the plot
    plt.plot(df.index,df['Hier'],color = 'black', label = 'Hierarchical Clustering')
    plt.plot(df.index,df['Optimized No Tree'],color = 'blue', label = 'K-Means Clustering')
    plt.plot(df.index,df['Optimized Tree'],color = 'red', label = 'Our Clustering')
    
    #Legend location and axis labels
    plt.legend(loc='upper left',frameon=False)
    plt.ylabel('Percent of Variance Across Clusters')
    plt.xlabel('Number of Cell Type Clusters')
    
    #Remove top bar from graph
    ax.tick_params(labelright=True,right=True)
    ax.spines['top'].set_visible(False)
    
    #Set x-axis ticks
    xticks = np.arange(2,13)
    locator = mpl.ticker.FixedLocator(xticks)
    ax.xaxis.set_major_locator(locator)
    
    #Set axis limits
    ax.set_xlim(1,13)
    ax.set_ylim(0,100)
    
    #Save image
    fig.savefig(conf.DATA_DIR + '/pctvaracrossclusters.png', dpi=400, bbox_inches='tight', pad_inches=0.01)

#Figure 7: Percent of variance retained by the three different clustering methods
def pctEnergyCharacterized(df):
    
    #Font to match latex
    mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']},size=10)
    
    #Define figure size
    fig, ax = plt.subplots(figsize=(6,4))
    
    #Add series to the plot
    plt.plot(df.index,df['Hier'],color = 'black', label = 'Hierarchical Clustering')
    plt.plot(df.index,df['Optimized No Tree'],color = 'blue', label = 'K-Means Clustering')
    plt.plot(df.index,df['Optimized Tree'],color = 'red', label = 'Our Clustering')
    
    #Legend location and axis labels
    plt.legend(loc = 'upper left', frameon=False)
    plt.ylabel('Percent of Across Column Variance Captured')
    plt.xlabel('Number of Cell Type Clusters')
    
    #Remove top bar from graph
    ax.tick_params(labelright=True,right=True)
    ax.spines['top'].set_visible(False)
    
    #Set x-axis ticks
    xticks = np.arange(2,13)
    locator = mpl.ticker.FixedLocator(xticks)
    ax.xaxis.set_major_locator(locator)
    
    #Set axis labels
    ax.set_xlim(1,13)
    ax.set_ylim(0,100)
    
    #Save image
    fig.savefig(conf.DATA_DIR + '/pctenergycharacterized.png', dpi=400, bbox_inches='tight', pad_inches=0.01)

#Figure 3 - Heatmap where each bicluster is assigned its mean
def clusterMapMeans(m_list, assignments, cell_types, K):
    
    #Font to match latex
    mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']},size=10)
    
    #Convert m_list to list of means for each bicluster 
    cell_cluster_grouped = [np.where(assignments == k)[0] for k in range(K)]
    means_m_list = []
    for m in m_list:
        means_m_list.append(np.array([m[:,cell_cluster_grouped[k]].mean() for k in range(K)]))
    
    #Dataframe with mean of each bicluster
    df = pd.concat([pd.DataFrame(m).T for m in means_m_list])
    float_mlist = [m.astype('float64') for m in m_list]
    df.index = range(df.shape[0])
    
    #Group cells by cluster
    cell_cluster_grouped = [np.where(assignments == a)[0] for a in range(len(set(assignments)))]
    
    #Convert each element to the mean of the bicluster
    for i, m in enumerate(float_mlist):
        for a in list(set(assignments)):
            for col in cell_cluster_grouped[a]:
                m[:,col] = df.at[i,a]
    
    #Combine list into datafrane
    df_m_list = [pd.DataFrame(m,columns = cell_types) for m in float_mlist]
    M = pd.concat(df_m_list)
    
    #Sort cell types by cluster
    zip_tosort = zip(assignments,cell_types)
    sorted_pairs = sorted(zip_tosort)
    tuples = zip(*sorted_pairs)
    assignments_sorted, cell_types_sorted = [list(tuple) for tuple in  tuples]
    
    #Sort columns by cell type clusters
    M = M[cell_types_sorted]
    
    #Horizontal line at the end of each row cluster
    horizontal_lines = []
    
    #Y-axis label in the middle of each row cluster with a label for the size of the row cluster
    y_axis_labels = []
    y_axis_ticks = [] 
    
    #Tracks number of loci thus far
    num_loci = 0
    for m in m_list[:-1]:
        
        #Add horizontal line at the end of this row cluster
        num_loci+=m.shape[0]
        horizontal_lines.append(num_loci)
        
        #Add a label for number of loci in a cluster
        #Must have more than 700 to prevent overlapping labels at the bottom of the exhibit
        if m.shape[0]>700:
            y_axis_labels.append(m.shape[0])
            
        #y-axis ticks are in the middle of the row-cluster as opposed to at the end
        if len(horizontal_lines)==1:
            y_axis_ticks.append((0+horizontal_lines[0])/2)
        else:
            y_axis_ticks.append((horizontal_lines[-2] + horizontal_lines[-1])//2)
        
    #Vertical line at the end of each cell type cluster
    vertical_lines = []
    col = 0
    for a in list(set(assignments))[:-1]:
        col+=len(assignments[assignments==a])
        vertical_lines.append(col)

    #Plot heatmap
    ax = sns.clustermap(M.reset_index(drop=True), row_cluster=False, col_cluster=False, figsize = (18,12), linewidths = 0.0, cmap = 'bwr',cbar_pos=(0.08, 0.4, 0.05, 0.18),xticklabels=cell_types,yticklabels=y_axis_labels).ax_heatmap
    
    #Add gridlines to visually seperate clusters
    ax.hlines(horizontal_lines, linewidths = 1, *ax.get_xlim())
    ax.vlines(vertical_lines, linewidths = 1, *ax.get_ylim())
    
    #Set axis labels and location of ticks
    ax.set_ylabel('Number of Loci in Cluster', fontsize = 15)
    ax.set_xlabel('Cell Types', fontsize = 15)
    ax.set_xticks(range(len(cell_types)))
    ax.set_yticks(y_axis_ticks)
    
    #Set Labels on the colorbar
    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks([0.0,0.25,0.5,0.75,1.0])
    
    #Save image
    figure = ax.get_figure()    
    figure.savefig(conf.DATA_DIR + '/heatmap_means_' + str(K) + '.png', dpi=400, bbox_inches='tight', pad_inches=0.01)

#Figure 4 - Heatmap where each bicluster is assigned open or closed
def clusterMapBinary(m_list, characterized_df, assignments, cell_types, k):
    
    #Font to match LaTex
    mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']},size=10)
    
    #Convert row to elements of 0 or 1
    def toBinary(row):
        
        #If both means are greater than .8 then assign open to all biclusters with this row cluster
        if row.min()>.8:
            return [1]*len(row)
        
        #If both means are less than .2 then assign closed to all biclusters with this row cluster
        if row.max() <.2:
            return [0]*len(row)  
        
        #Otherwise assign open to biclusters with the higher mean and closed to lower
        max_val = row.max()
        to_return = []
        for val in row:
            if val==max_val:
                to_return.append(1)
            else:
                to_return.append(0)
        return to_return
    
    #Convert dataframe with 2 means for each row cluster to binary
    high_low_df = characterized_df.apply(lambda x: toBinary(x), axis = 1, result_type='expand')
    high_low_mlist = m_list
    
    #Group cell types by cluster
    cell_cluster_grouped = [np.where(assignments == a)[0] for a in range(len(set(assignments)))]

    #Convert elements in m_list to open or closed based on their cluster assignment
    for i, m in enumerate(high_low_mlist):
        for a in list(set(assignments)):
            for col in cell_cluster_grouped[a]:
                m[:,col] = high_low_df.at[i,a]
    
    #Combine list into single dataframe
    df_m_list = [pd.DataFrame(m,columns = cell_types) for m in high_low_mlist]
    M = pd.concat(df_m_list)
    
    #Sort cell types by cluster
    zip_tosort = zip(assignments,cell_types)
    sorted_pairs = sorted(zip_tosort)
    tuples = zip(*sorted_pairs)
    assignments_sorted, cell_types_sorted = [list(tuple) for tuple in  tuples]
    
    #Sort columns by cell type clusters
    M = M[cell_types_sorted]
    
    #Horizontal line at the end of each row cluster
    horizontal_lines = []
    
    #Y-axis label in the middle of each row cluster with a label for the size of the row cluster
    y_axis_labels = []
    y_axis_ticks = []    
    
    #Tracks number of loci thus far
    num_loci = 0
    
    for m in m_list[:-1]:
        
        #Add horizontal line at the end of each row cluster
        num_loci+=m.shape[0]
        horizontal_lines.append(num_loci)
        
        #Add a label for number of loci in a cluster
        #Must have more than 700 to prevent overlapping labels at the bottom of the exhibit
        if m.shape[0]>700:
            y_axis_labels.append(m.shape[0])
            
        #y-axis ticks are in the middle of the row-cluster as opposed to at the end
        if len(horizontal_lines)==1:
            y_axis_ticks.append((0+horizontal_lines[0])/2)
        else:
            y_axis_ticks.append((horizontal_lines[-2] + horizontal_lines[-1])//2)
        
    #Vertical line at the end of each cell type cluster
    vertical_lines = []
    col = 0
    for a in list(set(assignments))[:-1]:
        col+=len(assignments[assignments==a])
        vertical_lines.append(col)
        
    #Set colormap colors - Red for open - blue for closed
    bwr_cmap = mpl.cm.get_cmap('bwr')
    myColors = (bwr_cmap(0.0),bwr_cmap(1.0))
    cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))
    
    #Plot heatmap
    ax = sns.clustermap(M.reset_index(drop=True), row_cluster=False, col_cluster=False, figsize = (18,12), linewidths = 0.0, cmap = cmap,cbar_pos=(0.08, 0.4, 0.05, 0.18),xticklabels=cell_types,yticklabels=y_axis_labels).ax_heatmap
    
    #Add gridlines to visually seperate clusters
    ax.hlines(horizontal_lines, linewidths = 1, *ax.get_xlim())
    ax.vlines(vertical_lines, linewidths = 1, *ax.get_ylim())
    
    #Set axis labels and location of ticks
    ax.set_ylabel('Number of Loci in Cluster', fontsize = 15)
    ax.set_xlabel('Cell Types', fontsize = 15)
    ax.set_xticks(range(len(cell_types)))
    ax.set_yticks(y_axis_ticks)
    
    #Set Labels on the colorbar
    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks([0.25,0.75])
    colorbar.set_ticklabels(['Closed', 'Open'])
    
    #Save image
    figure = ax.get_figure()    
    figure.savefig(conf.DATA_DIR + '/heatmap_' + str(k) + '.png', dpi=400, bbox_inches='tight', pad_inches=0.01)
    
    
#Figure not shown - Heatmap where each bicluster is assigned one of two means from k-means clustering
def clusterMap2Means(m_list, characterized_df, assignments, cell_types, k):
    
    #Font to match LaTex
    mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']},size=10)
    
    float_mlist = [m.astype('float64') for m in m_list]
    
    #Group cell types by cluster
    cell_cluster_grouped = [np.where(assignments == a)[0] for a in range(len(set(assignments)))]
    
    #Convert elements in m_list to their respective mean based on cluster assignment
    for i, m in enumerate(float_mlist):
        for a in list(set(assignments)):
            for col in cell_cluster_grouped[a]:
                m[:,col] = characterized_df.at[i,a]
    
    #Combine list into single dataframe
    df_m_list = [pd.DataFrame(m,columns = cell_types) for m in float_mlist]
    M = pd.concat(df_m_list)
    
    #Sort cell types by cluster
    zip_tosort = zip(assignments,cell_types)
    sorted_pairs = sorted(zip_tosort)
    tuples = zip(*sorted_pairs)
    assignments_sorted, cell_types_sorted = [list(tuple) for tuple in  tuples]
    
    #Sort columns by cell type clusters
    M = M[cell_types_sorted]
    
    #Horizontal line at the end of each row cluster
    horizontal_lines = []
    
    #Y-axis label in the middle of each row cluster with a label for the size of the row cluster
    y_axis_labels = []
    y_axis_ticks = []    
    
    #Tracks number of loci thus far
    num_loci = 0
    
    for m in m_list[:-1]:
        
        #Add horizontal line at the end of each row cluster
        num_loci+=m.shape[0]
        horizontal_lines.append(num_loci)
        
        #Add a label for number of loci in a cluster
        #Must have more than 700 to prevent overlapping labels at the bottom of the exhibit           
        if m.shape[0]>700:
            y_axis_labels.append(m.shape[0])
            
        #y-axis ticks are in the middle of the row-cluster as opposed to at the end
        if len(horizontal_lines)==1:
            y_axis_ticks.append((0+horizontal_lines[0])/2)
        else:
            y_axis_ticks.append((horizontal_lines[-2] + horizontal_lines[-1])//2)
        
    #Vertical line at the end of each cell type cluster
    vertical_lines = []
    col = 0
    for a in list(set(assignments))[:-1]:
        col+=len(assignments[assignments==a])
        vertical_lines.append(col)

    #Plot heatmap
    ax = sns.clustermap(M.reset_index(drop=True), row_cluster=False, col_cluster=False, figsize = (18,12), linewidths = 0.0, cmap = 'bwr',cbar_pos=(0.08, 0.4, 0.05, 0.18),xticklabels=cell_types,yticklabels=y_axis_labels).ax_heatmap
    
    #Add gridlines to visually seperate clusters
    ax.hlines(horizontal_lines, linewidths = 1, *ax.get_xlim())
    ax.vlines(vertical_lines, linewidths = 1, *ax.get_ylim())
    
    #Set axis labels and location of ticks
    ax.set_ylabel('Number of Loci in Cluster', fontsize = 15)
    ax.set_xlabel('Cell Types', fontsize = 15)    
    ax.set_xticks(range(len(cell_types)))
    ax.set_yticks(y_axis_ticks)
    
    #Set Labels on the colorbar
    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks([0.0,0.25,0.5,0.75,1.0])
    
    #Save image
    figure = ax.get_figure()    
    figure.savefig(conf.DATA_DIR + '/heatmap_2means_' + str(k) + '.png', dpi=400, bbox_inches='tight', pad_inches=0.01)