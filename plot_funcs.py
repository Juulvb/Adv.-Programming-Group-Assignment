"""
This module contains the plot methods used in the PCA analysis
An example on how to use them can be found in the Main.py module

The methods in this module: 
    PCA_plot_2d: creating a 2d plot of the first two principal components
    PCA_plot_3d: creating a 3d plot of the first three principal components
    PCA_plot_loadings: creating loading plots of the selected principal components
    cumulative: support method for PCA_plot_cumulative_explained_variance, returning the cumulative of a list
    PCA_plot_cumulative_explained_variance: creating a bar plot of the explained variances per principal component,
                                            and a line plot of the cumulative explained variance
    PCA_plot_scree: creating a scree plot of the explained variance                    
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from AssignmentPCA import getMaxIdxs
from CellLineRMAExpressionModule import CellLineRMAExpression
    
def PCA_plot_2d(labels, targets, subspace):
    """
    Given a list of M labels, a list of target labels and a Mx2 matrix of numbers, subspace,
    Create a 2d scatter plot with M datapoints in the 2d coordinates given in the subspace columns, coloured by their corresponding labels
    from labels present in the target labels 
        
    Parameters: 
        labels, a non-empty list of strings of length M, containing the labels of the data points
        targets, a non-empty list of strings, the values to which the labels should correspond
        subspace, a Mx2 matrix of numbers, containing the coordinates of the data points in the new subspace
    """
    
    plt.figure(figsize=(12, 10))
    for target in targets: #check the datapoints for each target label
        indicesToKeep = [i for i, label in enumerate(labels) if label == target] #check which data point belongs to the specific target label
        plt.scatter(subspace[indicesToKeep,0], subspace[indicesToKeep,1], s = 50, label=target) #plot the data points belonging to the target label
        
    #make a ylabel, xlabel and title.
    plt.xlabel("Principal Component 1")
    plt.ylabel("Principal Component 2")
    plt.title("Principal Component Analysis of RMA Expression of cell lines")
    
    #make a legend with the the labels.
    plt.legend()
    plt.show()

def PCA_plot_3d(labels, targets, subspace):
    """
    Given a list of M labels, a list of target labels and a Mx3 matrix of numbers, subspace,
    Create a 3d scatter plot with M datapoints in the 3d coordinates given in the subspace columns, coloured by their corresponding labels
    from labels present in the target labels 
        
    Parameters: 
        labels, a list of strings of length M, containing the labels of the data points
        targets, a list of strings, the values to which the labels should correspond
        subspace, a Mx3 matrix of numbers, containing the coordinates of the data points in the new subspace
    """
    #initialize the 3d figure 
    fig = plt.figure(figsize=(10,10))
    ax = Axes3D(fig)
    
    for target in targets: #check the datapoints for each target label
        indicesToKeep = [i for i, label in enumerate(labels) if label == target] #check which data point belongs to the specific target label
        ax.scatter(subspace[indicesToKeep,0], subspace[indicesToKeep,1], subspace[indicesToKeep,2], s = 50, label=target) #plot the data points belonging to the target label
        
    #make a ylabel, xlabel, zlabel and title.
    ax.set_xlabel("Principal Component 1")
    ax.set_ylabel("Principal Component 2")
    ax.set_zlabel("Principal Component 3")
    ax.set_title("Principal Component Analysis of \n RMA Expression of cell lines")
    
    #make a legend with the the labels.
    plt.legend()
    plt.show()
    
def PCA_plot_loadings(loadings, nr_genes=None):
    """
    Given a non empty MxN matrix of numbers, loadings, and an integer, nr_genes,
    Create M loading plots, showing the loadings of the nr_genes maximal variables    
        
    Parameters: 
        loadings, a MxN matrix of numbers, containing the loadings of the principal components 
        nr_genes, an integer specifying how many genes should be plotted
    """
    if nr_genes is None: nr_genes = len(loadings[0]) #check how many genes to plot, if not specified plot all
    
    for i in range(len(loadings)): #loop over each principal component
        loading = abs(loadings[i])
        idxs = getMaxIdxs(loading, nr_genes) #get the indexes of the genes with the highest loading 
        ticks = [CellLineRMAExpression.allparskeys[i] for i in idxs] #set the x labels to the correct gene names
        heights = loading[idxs] #set the heights of the loadings per gene
        
        plt.figure(figsize=(10,3)) #create a new figure
        
        plt.bar(ticks, heights) #plot the results
        plt.xticks(rotation=90) #allign the x labels vertically

        #add title and ylabel
        plt.title(f"Loading plot for PC{i+1}")
        plt.ylabel("value")
        plt.show()

def cumulative(l=[1]):
    """
    Given a non empty list of numbers, l, 
    Calculate the cumulative of the list 
        
    Parameters: 
        l, a non empty list of numbers
        
    Returns: the cumulative of list l
    """
    result = [0]*(len(l)+1) #initialize an a list of zeros 
    for i in range(len(l)): result[i+1] = result[i] + l[i] #fill the list with the cumulative of list l
    return result
    
def PCA_plot_cumulative_explained_variance(explained_variance, nr_pcs):
    """
    Given a non empty list of numbers, explained_variance, and an integer, nr_pcs,
    Create a plot of the explained variance per principle component and the cumulative     
        
    Parameters: 
        explained_variance, a list of numbers, containing the explained variance of the principal components 
        nr_pcs, an integer specifying how many principal components should be plotted
    """
    cum_exp_var = cumulative(explained_variance) #calculate the cumulative explained variance
    
    end_idx = len(explained_variance) + 1
    plt.bar(range(1, end_idx), explained_variance, alpha=0.5, label='Per PC') #make a bar graph 
    
    plt.plot(range(end_idx), cum_exp_var, '-o', label='Cumulative') #make a line plot for the cumulative explained variance 
    plt.plot(range(end_idx), [0.7]*len(cum_exp_var), label='70% Threshold') #make a line plot for the 70% threshold

    #make an appropriate ylabel and xlabel
    #plt.xticks(range(end_idx+1))
    plt.xlabel(f'N of {nr_pcs} principal components')
    plt.ylabel('Variance explained (per PC & cumulative)')
    
    #make a legend
    plt.legend()
    plt.show()

def PCA_plot_scree(explained_variance):  
    """
    Given a non empty list of numbers, explained_variance, 
    Create a scree plot   
        
    Parameters: 
        explained_variance, a list of numbers, containing the explained variance of the principal components 
    """
    plt.plot(range(1, len(explained_variance)+1), explained_variance, 'o-', linewidth=2, label='PC number') #create scree plot
    
    #make an appropriate ylabel and xlabel
    plt.xlabel("Number of PC's")
    plt.ylabel("Explained Variance")
    
    #make a legend
    plt.legend()
    plt.show()