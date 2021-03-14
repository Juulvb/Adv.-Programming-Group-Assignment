"""
This module holds the CellLineRMAExpression class, used to store the RMA expression information of cell lines

The methods in this module: 
    initclassvars: used to initialize the class variable allparskeys                  
"""

class CellLineRMAExpression:
    """
    Class used to store the RMA expression values of a cell line.
    """
    allparskeys = []
    
    def __init__(self, CellLineName = 'AU565', CosmicID = '910704', CancerType='BRCA'):
        """
        Initiate the values of an instance of class CellLineRMAExpression
        
        Parameters:
            CellLineName, string containing the name of this cellline 
            CosmicID, string containing the Cosmic ID of this cellline
            CancerType, string containing the Cancer type of this cellline
        """
        self.CellLineName = CellLineName #assign given name to instance
        self.CosmicID = CosmicID #assign given cosmic ID to instance
        self.CancerType = CancerType #assign given cancer type to instance 
        
    def load_RMAExpression(self, values = [0]*len(allparskeys)):
        """
        Given a list of numbers with the same length as allparskeys, values,
        assign the given values to the corresponding genes in this instance, 
        and store it in a dictionary with the gene names as keys. 
        
        Parameters:
            values, list of 244 numbers representing the RMAExpression of the genes defined in allparskeys of this instance
        """
        self.RMAExp_dict={} #create empty dictionary to store RMAExpression values in for this instance
        for i in range(len(CellLineRMAExpression.allparskeys)): #loop over all genes specified in allparskeys
            self.RMAExp_dict[CellLineRMAExpression.allparskeys[i]]=values[i] #assign corresponding values to the genes
        
def initclassvars(l):
    """
    Given a non-empty list l,
    Set the parameters names in class CellLineRMAExpression to the elements of l
    
    Parameter: l non-empty list
    """
    CellLineRMAExpression.allparskeys=l
    
    