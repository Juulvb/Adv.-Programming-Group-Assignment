"""
This module contains all methods used in the PCA analysis of the RMA expression in cancer cells
An example on how to use them can be found in the Main.py module

The methods in this module: 
    load_RMAExp_to_CellLines: creating a list of CellLineRMAExpression classes and filling them with the right information 
    load_RMAExp_to_matrix: extracts the RMAExpression information and stores it in a matrix
    normalize_list: z-score normalizes a given list
    normalize_matrix: z-score normalizes a given matrix per column
    covariance_of_two_lists: calculates the covariance of two lists
    covariance_matrix: calculates the covariance of the features in a matrix 
    calcMaxIdx: gets the index of the maximum value in a list
    getMaxIdxs: gets the indexes of the N maximum values in a list                    
"""

import math
import numpy as np
from CellLineRMAExpressionModule import initclassvars, CellLineRMAExpression

def load_RMAExp_to_CellLines(metadata, rma_expr, ListOfCellLineNumbers = None, metadata_labels = ["name", "COSMIC_ID", "TCGA_label"], lookup_variable = "name"):
    """
    Create a list of instances of class 'CellLineRMAExpression' per cell line given in ListOfCellLineNumbers, 
    and fill each instance with the correct name, cosmic_id, TCGA_label from metadata and the corresponding 
    RMAExpression values from rma_expr based on the lookup_variable.  
        
    Parameters: 
        metadata, non-empty pandas dataframe containing the name, cosmic ID and TCGA label (columns) of the cell lines (rows)
        rma_expr, non-empty pandas dataframe containing the RMAExpression per gene (columns) of each cell line (rows), 
                  indexed with the lookup_variable of the cell line
        ListOfCellLineNumbers, non-empty list of indexes corresponding to dataframe metadata
        metadata_labels, list of three strings, representing the names of the columns in metadata with the name, 
                         cosmid_ID and TCGA_label information (in that order)
        lookup_variable, string containing one of: "name", "cosmic_id", "tcga_label", 
                         representing the variable to which to match metadata and rma_expr
        
    Returns: a list of instances of the CellLineRMAExpression class, filled with the corresponding information and RMAExpression
    """
    initclassvars(rma_expr.columns) #assign gene names to class variable 'allparskeys'

    if ListOfCellLineNumbers is None: ListOfCellLineNumbers = metadata.index #if no indexes are specified, use all indexes
    nr_instances = len(ListOfCellLineNumbers) #Get the number of instances to create
    List_of_cellline_classes = [None]*nr_instances #Initialize shape size to prevent reallocation of memory
    
    for i, idx in enumerate(ListOfCellLineNumbers):  #Iterate over every cell line index, with idx the correct index and i the iteration
        name, cosmic_id, tcga_label = metadata.loc[idx, metadata_labels] #Get the name, cosmic_id and tcga_label of the cell line with index idx
        values = rma_expr.loc[str(eval(lookup_variable))].values #Get the RMA Expression values of the cell line with an index corresponding to 'name', the name of cell line i
    
        instance = CellLineRMAExpression(CellLineName=name, CosmicID=cosmic_id, CancerType=tcga_label) #load the cell line information into the instance
        instance.load_RMAExpression(values) #load the RMAExpression values into the instance
        List_of_cellline_classes[i] = instance #save the instance to the list which will be used as output 

    return List_of_cellline_classes

def load_RMAExp_to_matrix(data):
    """
    Given a list of CellLineRMAExpression class instances,
    Create an MxN matrix containing the RMA Expression values, with M the number of instances and N the number of genes
        
    Parameters: 
        data, non-empty list of CellLineRMAExpression class instances filled with the RMA Expression values 
              stored in the RMAExp_dict attribute of the instance
        
    Returns: a MxN numpy array containing the RMA Expression values, with M the number of instances and N the number of genes
    """
    nr_instances = len(data) #get the number of instances
    nr_variables = len(data[0].RMAExp_dict) #get the number of genes 
    rma_exp_matrix = np.empty((nr_instances, nr_variables)) #Initialize empty numpy array to prevent reallocation of memory
    
    for i, instance in enumerate(data): #iterate over every instance, with i the iteration
        rma_exp_matrix[i] = list(instance.RMAExp_dict.values()) #extract the values of the RMAExp_dict of the instance, containing the RMAExpression values, and store them in the matrix
    
    return rma_exp_matrix

def normalize_list(variable_list):
    """
    Given a list of numbers,
    create a list with its normalized values
        
    Parameters: 
        variable_list, a non-empty list of numbers 
        
    Returns: a list containing the z-score normalized values of data
    """
    mean = sum(variable_list) / len(variable_list) #get the mean of the variable
    variance = sum(pow(x-mean,2) for x in variable_list) / len(variable_list) #get the variance of the variable
    std = math.sqrt(variance) #get the standard deviation of the variable
    
    normalized_list = [(x-mean)/std for x in variable_list] #subtract the mean and divide by the standard deviation for every element of the list
   
    return normalized_list
 
def normalize_matrix(data_matrix):
    """
    Given a MxN numpy array,
    create a MxN numpy array containing the z-score normalized values of data per column N
        
    Parameters: 
        data, a non-empty MxN numpy array of numbers
        
    Returns: a MxN numpy array containing the z-score normalized values of data
    """
    data_matrix_norm = data_matrix.copy() #copy the matrix, so the original is not overwritten
    nr_variables = data_matrix_norm.shape[1] #get the number of variables: N
    
    for i in range(nr_variables): #iterate over every variable
        data_matrix_norm[:,i] = normalize_list(data_matrix_norm[:,i]) 
   
    return data_matrix_norm

def covariance_of_two_lists(x, y):
    """
    Given two lists of numbers of equal length,
    calculate their covariance
        
    Parameters: 
        x, a non-empty list of normalized numbers
        y, a non-empty list of normalized numbers of equal length as x
        
    Returns: the covariance of x and y
    """
    
    cov_list = [a * b for (a,b) in zip(x,y)] #make an element wise multiplication
    covariance = sum(cov_list) / (len(cov_list)-1) #calculate the covariance
    
    return covariance

def covariance_matrix(data):
    """
    Given a MxN numpy array of numbers, with M the number of instances and N the number of variables
    Calculate the NxN covariance matrix 
        
    Parameters: 
        data, a non-empty MxN numpy array of normalized numbers, with M observations and N featueres
        
    Returns: a NxN numpy array containing the covariance matrix of the data
    """
    nr_variables = data.shape[1] #get the number of variables N
    covMatrix = np.zeros((nr_variables, nr_variables)) #Initialize empty numpy array to prevent reallocation of memory

    for x in range(nr_variables): #iterate over the number of variables 
        for y in range(x, nr_variables): #iterate over the number of variables 
            covMatrix[x, y] = covariance_of_two_lists(data[:,x], data[:,y]) #Calculate the covariance of variables X and Y
            covMatrix[y, x] = covMatrix[x, y] #Assign the known covariance of Y, X to X, Y position in the covariance matrix
    return covMatrix

def calcMaxIdx(l=[1]):
    """
    Given a list
    find the index of the maximum value
    
    Parameters:
        l non-empty list of numbers or nan's, with at least one number
    
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
    Given a list
    find the indexes of the number_idxs maximum values, in order from highest to lowest
    
    Parameters: 
        l, non-empty list of numbers or nan's, with at least number_idxs numbers
           number_idxs, integer representing the amount of indexes to find
    
    Returns: a list containing the indexes of the maximum values of list l, in order of highest value
    """
    max_idxs = [None]*number_idxs #Initialize empty list to prevent reallocation of memory when saving results
    control_list = np.asarray(l.copy()).real.astype(float) #Copy the original list to prevent overwriting, and make sure the numbers are all real numbers of float type
    
    for i in range(number_idxs): #find the max index 'number_idxs' times
        idx = calcMaxIdx(control_list) #find the max index
        max_idxs[i] = idx #save the max index to the result list
        control_list[idx] = np.nan #change the value at the max index to nan so it won't be found again, thus the next max index can be found
        
    return max_idxs

