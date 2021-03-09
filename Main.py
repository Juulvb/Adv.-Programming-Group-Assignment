# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 15:02:14 2021

@author: 20164798
"""
from AssignmentPCA import getMaxIdxs, load_RMAExp_to_CellLines, load_RMAExp_to_matrix, normalize_matrix, covariance_matrix
from plot_funcs import PCA_plot_3d, PCA_plot_2d, PCA_plot_loadings, PCA_plot_cumulative_explained_variance, PCA_plot_scree

import numpy as np
import pandas as pd


metadata = pd.read_csv("GDSC_metadata.csv", index_col=0) #load all metadata from the csv file to a pandas dataframe, set the first column as index
rma_expr = pd.read_csv("GDSC_RNA_expression.csv", index_col=0) #load all RMA Expression data from the csv file to a pandas dataframe, set the name of the cell line as index
#define the cell lines to inspect
ListOfCellLineNumbers= [1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 
                        17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 
                        30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 
                        43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 
                        56, 57, 58, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 
                        70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 
                        83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 
                        96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 
                        107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 
                        117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 
                        127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 
                        137, 138, 139, 140, 141, 142, 143, 144, 145, 147, 148]

data = load_RMAExp_to_CellLines(metadata, rma_expr, ListOfCellLineNumbers) #load the list of instances of class CellLineRMAExpression, and fill them with the RMA expression data
data_matrix = load_RMAExp_to_matrix(data) #load the RMA Expression data of all instances in the data list to a numpy array
data_matrix_norm = normalize_matrix(data_matrix) #normalize the RMA Expression data per gene
cov_matrix = covariance_matrix(data_matrix_norm) #calculate the covariance matrix for the normalized RMA expression data using the covariance_matrix method
cov_matrix_np = np.cov(data_matrix_norm, rowvar=False) #calculate the covariance matrix for the normalized RMA expression data using the numpy covariance method
eig_vals, eig_vecs = np.linalg.eig(cov_matrix) #calculate the eigenvalues and eigenvectors of the covariance matrix
eig_vals1, eig_vecs1 = np.linalg.eig(cov_matrix_np)


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

nr_PC = 11 #for evaluation, check more principal components
idxs = getMaxIdxs(eig_vals, nr_PC) #get the indexes of the 11 maximal eigenvalues 

explained_variance = eig_vals[idxs] / sum(eig_vals)
PCA_plot_cumulative_explained_variance(explained_variance.real, len(eig_vals))
PCA_plot_scree(explained_variance.real)


