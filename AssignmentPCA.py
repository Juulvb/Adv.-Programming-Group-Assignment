# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 11:39:51 2021

@author: 20164798
"""
from CellLineRMAExpressionModule import initclassvars, CellLineRMAExpression

import math
import numpy as np

from time import time

def load_RMAExp_to_CellLines(metadata, rma_expr, ListOfCellLineNumbers = None, metadata_labels = ["name", "COSMIC_ID", "TCGA_label"], lookup_variable = "name"):
    """
    Given: 
        a pandas dataset, metadata, containing the names, cosmic_id's and TCGA_labels of the cell lines,
        and a pandas dataset, rna_expr, containing the RMA gene expression per gene per cell line
        a non-empty list, ListOfCellLineNumbers, of indexes from metadata to use,
    Create a list of instances of class 'CellLineRMAExpression' per cell line given in ListOfCellLineNumbers, 
    and fill each instance with the correct name, cosmic_id, TCGA_label and RMAExpression values 
        
    Parameters: 
        ListOfCellLineNumbers, non-empty list of indexes corresponding to dataframe metadata
        metadata, non-empty pandas dataframe containing the name, cosmic ID and TCGA label (columns, in that order) of the cell lines (rows)
        rma_expr, non-empty pandas dataframe containing the RMAExpression per gene (columns) of each cell line (rows), indexed with the name of the cell line
        
    Returns: a list of instances of the CellLineRMAExpression class, filled with the corresponding information and RMAExpression
    """
    initclassvars(rma_expr.columns) #assign gene names to class variable 'allparskeys'

    if ListOfCellLineNumbers is None: ListOfCellLineNumbers = metadata.index #if no indexes are specified, use all indexes
    nr_instances = len(ListOfCellLineNumbers) #Get the number of instances to create
    results = [None]*nr_instances #Initialize shape size to prevent reallocation of memory
    
    for i, idx in enumerate(ListOfCellLineNumbers):  #Iterate over every cell line index, with idx the correct index and i the iteration
        name, cosmic_id, tcga_label = metadata.loc[idx, metadata_labels] #Get the name, cosmic_id and tcga_label of the cell line with index idx
        values = rma_expr.loc[str(eval(lookup_variable))].values #Get the RMA Expression values of the cell line with an index corresponding to 'name', the name of cell line i
    
        instance = CellLineRMAExpression(CellLineName=name, CosmicID=cosmic_id, CancerType=tcga_label) #load the cell line information into the instance
        instance.load_RMAExpression(values) #load the RMAExpression values into the instance
        results[i] = instance #save the instance to the list which will be used as output 

    return results

def covariance_of_two_lists(x, y):
    """
    Given two not empty lists of numbers of equal length, x and y,  
    calculate the covariance
        
    Parameters: 
        x, a non empty list of normalized numbers
        y, a non empty list of normalized numbers of equal length as x
        
    Returns: the covariance of x and y
    """
    #iterate over each corresponding element in x and y
    #for each element subtract the corresponding mean, and then multiply the two elements
    #a list, cov_list, of the multiplied variances is constructed
    
    cov_list = [a * b for (a,b) in zip(x,y)] 
    covariance = sum(cov_list) / (len(cov_list)-1) #calculate the covariance
    
    return covariance 

def load_RMAExp_to_matrix(data):
    """
    Given a non empty list of CellLineRMAExpression class instances filled with the RMA Expression values, data, 
    Create an MxN matrix containing the RMA Expression values, with M the number of instances and N the number of genes
        
    Parameters: 
        data, non-empty list of CellLineRMAExpression class instances filled with the RMA Expression values
        
    Returns: a MxN numpy array containing the RMA Expression values, with M the number of instances and N the number of genes
    """
    nr_instances = len(data) #get the number of instances
    nr_variables = len(data[0].RMAExp_dict) #get the number of genes 
    results = np.empty((nr_instances, nr_variables)) #Initialize empty numpy array to prevent reallocation of memory
    
    for i, instance in enumerate(data): #iterate over every instance, with i the iteration
        results[i] = list(instance.RMAExp_dict.values()) #extract the values of the RMAExp_dict of the instance, containing the RMAExpression values, and store them in the matrix
    
    return results

def normalize_list(variable_list):
    """
    Given a list of numbers, variable_list, create a list with the normalized values of variable_list
        
    Parameters: 
        variable_list, a list of numbers 
        
    Returns: a list containing the z-score normalized values of data
    """
    mean = sum(variable_list) / len(variable_list) #get the mean of the variable
    variance = sum(pow(x-mean,2) for x in variable_list) / len(variable_list) #get the variance of the variable
    std = math.sqrt(variance) #get the standard deviation of the variable
    
    normalized_list = [(x-mean)/std for x in variable_list] #subtract the mean and divide by the standard deviation for every element of the list
   
    return normalized_list

def normalize_matrix(data_matrix):
    """
    Given a MxN numpy array, data 
    create a MxN numpy array containing the z-score normalized values of data per column N
        
    Parameters: 
        data, a MxN numpy array of numbers
        
    Returns: a MxN numpy array containing the z-score normalized values of data
    """
    data_matrix_norm = data_matrix.copy() #copy the matrix, so the original is not overwritten
    nr_variables = data_matrix_norm.shape[1] #get the number of variables: N
    
    for i in range(nr_variables): #iterate over every variable
        data_matrix_norm[:,i] = normalize_list(data_matrix_norm[:,i]) 
   
    return data_matrix_norm


def covariance_matrix(data):
    """
    Given a MxN numpy array of numbers, with M the number of instances and N the number of variables
    Calculate the NxN covariance matrix 
        
    Parameters: 
        data, a MxN numpy array of numbers, with a variable per column
        
    Returns: a NxN numpy array containing the covariance matrix of the data
    """
    nr_variables = data.shape[1] #get the number of variables N
    covMatrix = np.zeros((nr_variables, nr_variables)) #Initialize empty numpy array to prevent reallocation of memory
    start_time = time()
    for x in range(nr_variables): #iterate over the number of variables 
        print(x, time()-start_time)
        for y in range(x, nr_variables): #iterate over the number of variables 
            covMatrix[x, y] = covariance_of_two_lists(data[:,x], data[:,y]) #Calculate the covariance of variables X and Y
            covMatrix[y, x] = covMatrix[x, y] #Assign the known covariance of Y, X to X, Y position in the covariance matrix
    return covMatrix

def calcMaxIdx(l=[1]):
    """
    Given a non-empty list l of numbers or nan's, with at least one number,
    find the index of the maximum value
    
    Parameter: l non-empty list of numbers
    
    Returns: the index of the maximum value
    """
    nx = 0 #set the index of the maximal value to 0
    mx=l[nx] #set the maximal value to the value at index nx
    
    while(math.isnan(mx)): 
        #verify that the set maximal value is a number
        #if not, iterate till a number is found, and assign it as the maximal value
        nx += 1
        mx = l[nx]
        
    #mx is the current max
    n=0
    while (n<len(l)):
        # for all 0<=i<n: l[i]<=mx and
        # there exists an i, 0<=i<n such that l[i]=mx
        # check if the element is not a nan
        if not math.isnan(l[n]) and l[n]>mx: 
            mx=l[n]
            nx=n
        n=n+1
        # for all 0<=i<n: l[i]<=mx and
        # there exists an i, 0<=i<n such that l[i]=mx
        
    # n=len(l) and for all 0<=i<n: l[i]<=mx and
    # for all 0<=i<len(l): l[i]<=mx and
    # there exists an i, 0<=i<len(l) such that l[i]=mx
    # hence, this i = nx
    
    return nx

def getMaxIdxs(l=[1], number_idxs=1):
    """
    Given a non-empty list l of numbers or nan's, with at least number_idxs numbers,
    find the indexes of the number_idxs maximum values, in order from highest to lowest
    
    Parameter: 
        l, non-empty list of numbers or nan's, with at least number_idxs numbers
        number_idxs, integer representing the amount of indexes to find
    
    Returns: a list containing the indexes of the maximum values of list l, in order of highest value
    """
    result = [None]*number_idxs #Initialize empty list to prevent reallocation of memory when saving results
    control_list = np.asarray(l.copy()).real.astype(float) #Copy the original list to prevent overwriting, and make sure the numbers are all real numbers of float type
    
    for i in range(number_idxs): #find the max index 'number_idxs' times
        idx = calcMaxIdx(control_list) #find the max index
        result[i] = idx #save the max index to the result list
        control_list[idx] = np.nan #change the value at the max index to nan so it won't be found again, thus the next max index can be found
        
    return result

