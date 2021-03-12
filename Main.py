# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 15:02:14 2021

@author: 20164798
"""
from AssignmentPCA import getMaxIdxs, load_RMAExp_to_CellLines, load_RMAExp_to_matrix, normalize_matrix, covariance_matrix
from plot_funcs import PCA_plot_3d, PCA_plot_2d, PCA_plot_loadings, PCA_plot_cumulative_explained_variance, PCA_plot_scree

import numpy as np
import pandas as pd

import pickle

print("load data")
metadata = pd.read_csv("Cell_lines and COSMIC_ID.tsv", sep='\t') #load all metadata from the csv file to a pandas dataframe, set the first column as index
rma_expr = pd.read_csv("Cell_line_RMA_proc_basalExp.tsv", sep='\t') #load all RMA Expression data from the csv file to a pandas dataframe, set the name of the cell line as index

rma_expr = rma_expr.T
rma_expr.columns = rma_expr.loc['GENE_SYMBOLS']
rma_expr.dropna(axis=1, inplace=True)

metadata.index = metadata['COSMIC_ID'].astype(str)
metadata = metadata[~metadata.index.duplicated(keep='first')]
metadata = metadata[metadata['Tissue sub-type']!='UNCLASSIFIED']
metadata = metadata[metadata.index.isin(rma_expr.index)]
metadata.dropna(axis=0, inplace=True)

data = load_RMAExp_to_CellLines(metadata, rma_expr, metadata_labels=["Name", "COSMIC_ID", "Tissue sub-type"], lookup_variable="cosmic_id") #load the list of instances of class CellLineRMAExpression, and fill them with the RMA expression data
data_matrix = load_RMAExp_to_matrix(data) #load the RMA Expression data of all instances in the data list to a numpy array
data_matrix_norm = normalize_matrix(data_matrix) #normalize the RMA Expression data per gene

#%%
print("covariance")
cov_matrix = np.cov(data_matrix_norm, rowvar=False) #covariance_matrix(data_matrix_norm) #calculate the covariance matrix for the normalized RMA expression data using the covariance_matrix method
eig_vals, eig_vecs = np.linalg.eig(cov_matrix) #calculate the eigenvalues and eigenvectors of the covariance matrix

#%%
# with open('eig_vals.pkl', 'rb') as f: eig_vals = pickle.load(f)
# with open('eig_vecs_1.pkl', 'rb') as f: eig_vecs = pickle.load(f)
# with open('eig_vecs_2.pkl', 'rb') as f: eig_vecs2 = pickle.load(f)

# eig_vals = eig_vals[0]
# eig_vecs = np.vstack((eig_vecs[0], eig_vecs2[0]))
#%%
print("PCA")
nr_PC = 3 #determine how many principle components to investigate (normally no more then 3)
idxs = getMaxIdxs(eig_vals, nr_PC) #get the indexes of the 3 maximal eigenvalues 

labels = [instance.CancerType for instance in data] #load the labels (e.g. cancer types) of the instances, by iterating over all instances and loading their CancerType attribute
targets = list(set(labels)) #extract all unique options occurring in labels, and store them in a list

matrix_w = np.vstack([eig_vecs[idxs[:2]]]).T #get the eigenvectors corresponding to the 2 maximal eigenvalues, and store them in a matrix 
new_subspace = data_matrix_norm.dot(matrix_w).real #calculate the new subspace, based on the two principal components
PCA_plot_2d(labels, targets, new_subspace) #plot the cell lines in the new subspace, coloured by their label

matrix_w = np.vstack([eig_vecs[idxs[:3]]]).T #get the eigenvectors corresponding to the 3 maximal eigenvalues, and store them in a matrix 
new_subspace = data_matrix_norm.dot(matrix_w).real #calculate the new subspace, based on the two principal components
PCA_plot_3d(labels, targets, new_subspace) #plot the cell lines in the new 3D subspace, coloured by their label

loadings = [vec * val**0.5 for vec, val in zip(eig_vecs[idxs], eig_vals[idxs])] #calculate the loadings for the 3 principal components
PCA_plot_loadings(loadings, 50) #make a loading plot for each of the principal components

nr_PC = 250 #for evaluation, check more principal components
idxs = getMaxIdxs(eig_vals, nr_PC) #get the indexes of the 11 maximal eigenvalues 

explained_variance = eig_vals[idxs] / sum(eig_vals)
PCA_plot_cumulative_explained_variance(explained_variance.real, len(eig_vals))
PCA_plot_scree(explained_variance.real)


