"""
This module contains the PCA of the RMA expression in cancer cells
It uses the methods from AssignmentPCA.py and plot_funcs.py, as well as the class from CellLineRMAExpressionModule.py

First the data is loaded and formatted in the right way, with the following demands: 
    * The features of the data should be in the columns, with each row being an instance
    * The indexes of metadata and rma_expr should be corresponding in attribute and format
    * Nan values, duplicates and other unclassified data should be removed

Then PCA is executed: 
    * The data is prepared 
    * Covariance, eigenvalues and eigenvectors are calculated
    * Principal components are identified and visualized in plots
    * Loadings of the principal components are visualized
    * The explained variance is calculated and visualized
"""

from AssignmentPCA import getMaxIdxs, load_RMAExp_to_CellLines, load_RMAExp_to_matrix, normalize_matrix
from plot_funcs import PCA_plot_3d, PCA_plot_2d, PCA_plot_loadings, PCA_plot_cumulative_explained_variance, PCA_plot_scree

import numpy as np
import pandas as pd

##### Data Preperation #####

# Load the data from the tsv files
metadata = pd.read_csv("Cell_lines and COSMIC_ID.tsv", sep='\t') 
rma_expr = pd.read_csv("Cell_line_RMA_proc_basalExp.tsv", sep='\t')

# Prepare the RMA_expression data
rma_expr = rma_expr.T #transform the data so that the columns represent the features and the rows the instances
rma_expr.columns = rma_expr.loc['GENE_SYMBOLS'] #set the names of the features as the column names
rma_expr.dropna(axis=1, inplace=True) #drop all genes with nan values in the data

# Prepare the metadata
metadata.index = metadata['COSMIC_ID'].astype(str) #set the cosmic_id of the cell lines as the index
metadata = metadata[~metadata.index.duplicated(keep='first')] #remove duplicated notations of the same cell line
metadata = metadata[metadata['Tissue sub-type']!='UNCLASSIFIED'] #remove cell lines which are unclassified
metadata = metadata[metadata.index.isin(rma_expr.index)] #remove cell lines of which no RMA_expression data is present in rma_expr
metadata.dropna(axis=0, inplace=True) #drop all cell lines with nan values in the data

##### PCA #####

# Preparing the data for PCA 
data = load_RMAExp_to_CellLines(metadata, rma_expr, metadata_labels=["Name", "COSMIC_ID", "Tissue sub-type"], lookup_variable="cosmic_id") #load the list of instances of class CellLineRMAExpression, and fill them with the RMA expression data
data_matrix = load_RMAExp_to_matrix(data) #load the RMA Expression data of all instances in the data list to a numpy array
data_matrix_norm = normalize_matrix(data_matrix) #normalize the RMA Expression data per gene

# Execute PCA
cov_matrix = np.cov(data_matrix_norm, rowvar=False) #covariance_matrix(data_matrix_norm) #calculate the covariance matrix for the normalized RMA expression data using the covariance_matrix method
eig_vals, eig_vecs = np.linalg.eig(cov_matrix) #calculate the eigenvalues and eigenvectors of the covariance matrix

nr_PC = 3 #determine how many principle components to investigate (normally no more then 3)
idxs = getMaxIdxs(eig_vals, nr_PC) #get the indexes of the 3 maximal eigenvalues 
loadings = [vec * val**0.5 for vec, val in zip(eig_vecs[idxs], eig_vals[idxs])] #calculate the loadings for the 3 principal components

# Read in the labels for the plots 
labels = [instance.CancerType for instance in data] #load the labels (e.g. cancer types) of the instances, by iterating over all instances and loading their CancerType attribute
targets = list(set(labels)) #extract all unique options occurring in labels, and store them in a list

# 2D plot
matrix_w_2d = np.vstack([eig_vecs[idxs[:2]]]).T #get the eigenvectors corresponding to the 2 maximal eigenvalues, and store them in a matrix 
new_2d_subspace = data_matrix_norm.dot(matrix_w_2d).real #calculate the new subspace, based on the two principal components
PCA_plot_2d(labels, targets, new_2d_subspace) #plot the cell lines in the new subspace, coloured by their label

# 3D plot
matrix_w_3d = np.vstack([eig_vecs[idxs[:3]]]).T #get the eigenvectors corresponding to the 3 maximal eigenvalues, and store them in a matrix 
new_3d_subspace = data_matrix_norm.dot(matrix_w_3d).real #calculate the new subspace, based on the two principal components
PCA_plot_3d(labels, targets, new_3d_subspace) #plot the cell lines in the new 3D subspace, coloured by their label

# Loading plot
PCA_plot_loadings(loadings, 50) #make a loading plot for each of the principal components

# Analysis of the results by an explained variance plot and a scree plot 
nr_PC = 200 #for evaluation, check more principal components
idxs = getMaxIdxs(eig_vals, nr_PC) #get the indexes of the 200 maximal eigenvalues 

explained_variance = eig_vals[idxs] / sum(eig_vals)
PCA_plot_cumulative_explained_variance(explained_variance.real, len(eig_vals))
PCA_plot_scree(explained_variance.real)


