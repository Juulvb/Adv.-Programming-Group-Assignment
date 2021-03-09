# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 10:16:14 2021

@author: 20164798
"""

default_list = ['TFPI', 'TFAP2B', 'MGST1', 'PRSS3', 'ISL1', 'SNAI2', 'SERPINB1', 'VIM', 'CD44', 'VCAN', 'CYBA', 'YBX3', 'MYLK', 'COL17A1', 'KIF26A', 'SYT1', 'GAL', 'IGF2BP2', 'CA12', 'HOXA9', 'EDN1', 'CDH17', 'CHRNA3', 'CXCL2', 'GSTP1', 'ABCB1', 'TESC', 'LYZ', 'DSP', 'LGALS1', 'PYGL', 'PLEK2', 'CHGA', 'CTSZ', 'BEX4', 'KLF5', 'ESRP2', 'QPRT', 'CA2', 'ESRP1', 'STMN2', 'TUSC3', 'TNNT1', 'CEACAM5', 'PTN', 'GSDME', 'CAV2', 'CAV1', 'SERPINE1', 'RARRES2', 'AGR2', 'ELAVL2', 'GATA3', 'ACTA2', 'DKK1', 'KRT23', 'CCL2', 'RAB34', 'PHOX2B', 'NMU', 'CPE', 'MDK', 'MGP', 'LDHB', 'TPD52L1', 'SPARC', 'C3orf14', 'GALNT3', 'EFEMP1', 'FN1', 'IGFBP2', 'RGS4', 'TSPAN1', 'CCN2', 'VAMP8', 'SPP1', 'RARRES1', 'CNRIP1', 'EPCAM', 'TGFBI', 'DUSP4', 'PLBD1', 'TNFSF10', 'TWIST1', 'SRGN', 'SLPI', 'RAB17', 'MT2A', 'BMP4', 'C3', 'IFI6', 'TSPAN8', 'GNG11', 'RPS4Y1', 'BST2', 'PXDN', 'GDF15', 'ASS1', 'AKAP12', 'MAP1B', 'PPP1R1B', 'LGALS3', 'RAMP1', 'RAB25', 'BHMT2', 'BEX2', 'RNF128', 'BEX1', 'RARRES3', 'MYCN', 'RTL8C', 'ANXA1', 'KRT7', 'SERPINE2', 'MYC', 'TUBB2B', 'IER3', 'FGFBP1', 'SYTL2', 'MMP7', 'SQOR', 'MYOF', 'ANXA3', 'ERP27', 'DUSP6', 'RTN1', 'SLC27A2', 'GCNT3', 'IFITM3', 'EMP3', 'CCN1', 'RGS5', 'SELENBP1', 'SLC10A4', 'PLK2', 'IGFBP3', 'MSN', 'MAL2', 'PLIN2', 'LCN2', 'PRSS23', 'AKR1C2', 'SPOCK1', 'PLOD2', 'UCHL1', 'MX1', 'KRTCAP3', 'TMSB15A', 'TFF3', 'TFF1', 'AZGP1', 'PDZK1IP1', 'ELAVL4', 'INHBB', 'S100A11', 'S100A9', 'IFI16', 'FABP1', 'RPL39L', 'S100P', 'HAND2', 'FOXQ1', 'ALDH1A1', 'FBP1', 'CLDN3', 'CDX2', 'IFI27', 'HTRA1', 'HDGFL3', 'CENPV', 'BEX3', 'NNMT', 'ZNF667-AS1', 'C15orf48', 'SCG5', 'PRR15L', 'IGF2', 'TUBA1A', 'TMC4', 'SPINT2', 'C19orf33', 'IGFBP6', 'SDC2', 'TM4SF4', 'TM4SF1', 'FOS', 'KRT8', 'CDH2', 'UGT2B7', 'KRT19', 'KRT20', 'LGALS4', 'SCG2', 'MUCL1', 'LRRN3', 'AGR3', 'MUC13', 'GPR160', 'HOXB2', 'SFN', 'TUBB6', 'NUPR1', 'GPX2', 'PRR15', 'BASP1', 'IRX3', 'APOBEC3B', 'RNF182', 'TSPYL5', 'MAB21L1', 'CLRN3', 'PNMA8A', 'C1S', 'TACSTD2', 'IFITM2', 'PRAME', 'IFIT1', 'AKR1C1', 'S100A16', 'CLDN4', 'S100A14', 'SYCP2', 'AKR1C3', 'S100A4', 'SULF2', 'ANXA4', 'ACSL5', 'DPP4', 'S100A10', 'S100A6', 'AKR1B10', 'EPS8L3', 'ALPK2', 'PSMB8', 'HSPA1A', 'SERPINB5', 'MLLT11', 'AKR1B10P1', 'UCA1', 'MIR205HG', 'PHGR1', 'HLA-B', 'PSMB9', 'SELENOP']

class CellLineRMAExpression:
    
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
        
def initclassvars(l=default_list):
    """
    Given a non-empty list l,
    Set the parameters names in class CellLineRMAExpression to the elements of l
    
    Parameter: l non-empty list
    """
    CellLineRMAExpression.allparskeys=l