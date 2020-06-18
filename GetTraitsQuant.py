#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is for pre-processing the raw phenotypic data of UKBB. Similar to "GetCovartiates.py" 
and "GetTraitsDichot.py" but only for quantitative traits that don't have multiple instances.
=> mainly those analysed by fastSMC.

Created on Tue Jun 11 16:32:59 2019
@author: kalantzi
18/6/20: Include both UKBB_16549 and UKBB_43206 indexing
"""
###################
# initializations #
import pandas as pd
import numpy as np
import sys

SINGLE = False # whether there are traits with single values or not
file_out = "Traits.tab"
   
#file_ukbb_data = "/well/palamara/projects/UKBB_APPLICATION_43206/DOWNLOAD_DATA/Phenotypes.csv"
#file_relabel ="/well/palamara/projects/ID_MAPPING/mapping_nonegative.csv"
    
if len(sys.argv)< 3:
    print( "Error in the number of arguments! Run:")
    print( "\tpython GetTraits.py Phenotypes.csv mapping.csv" )
    sys.exit()

####################
# Load the dataset #
print("Loading the dataset...",end='')
ukbb_data = pd.read_csv(sys.argv[1], sep=',' )
ukbb_data = ukbb_data.set_index('eid')
print("\rData loaded - Total number of variables = {0} | subjects found = {1}".format( ukbb_data.shape[1], ukbb_data.shape[0] ))

########################################################
# Select the appropriate set of subjects and load data #    
Mapping = {}
with open(sys.argv[2], 'r') as f:
    _=next(f) # skip header
    for line in f:
        tokens=line.split(sep=',')
        # map UKBB_2 to (IBD-ID, UKBB_1) and get rid of "\n":
        Mapping[ int(tokens[2][:-1]) ] = (tokens[0], tokens[1]) # Mapping has int -> (str, str)  

Samples = [x for x in Mapping if x in ukbb_data.index ] #
print("Subjects to keep = %d" % len(Samples))     
   
Traits = pd.DataFrame()
Traits['IID'] = [ Mapping[s][0] for s in Samples ]
Traits['UKBB_1'] = [ Mapping[s][1] for s in Samples ] # this and the next are used just for verification
Traits['UKBB_2'] = Samples 

#################################################################
# Select the appropriate set of columns for quantitative traits #
Quant_Dict= {
 "EosinoPerce":"30210-",
 "EosinoCount":"30150-",
 "MonocyteCount":"30130-",
 "PlateletCount":"30080-",
 "PlateletCrit":"30090-",
 "PlateletDW":	"30110-",
 "MeanPTV":"30100-",
 "MeanSCV":"30270-",
 "MeanCorpHaem":"30050-",
 "MeanCorpVol":"30040-",
 "RedBloodCellDW":"30070-",
 "RedBloodCellCount":"30010-",
}

print("Getting quantitative traits with multiple instances...")
for trait in Quant_Dict.keys(): 
    print(trait)
    temp_trait = []
    for s in Samples:
        # we have to check for NaNs and average accordingly
        subj_trait = []
        for c in range(3):
            col = Quant_Dict[trait]+str(c)+'.0'
            if  np.isnan( ukbb_data[col][s] ):
                subj_trait.append(np.nan)
            else:
                subj_trait.append( ukbb_data[col][s] )
        temp_trait.append( np.nanmean(subj_trait) ) 
    Traits[trait] = temp_trait
    # do some checking
    print("Positive values found = %d" % np.count_nonzero(temp_trait) )

if SINGLE:
    print("Getting quantitative traits with single instances...")
    singletons = {"HeelBoneMineralDensity":"3148-0.0"}
    for trait in singletons:
        print(trait)
        count=0
        subj_trait = []
        for s in Samples:
            col = singletons[trait]
            # there's only one measurement per indv
            if np.isnan( ukbb_data[col][s] ):
                subj_trait.append(np.nan)
            else:
                subj_trait.append( ukbb_data[col][s] )
                count+=1
        Traits[trait] = subj_trait
        # do some checking
        print("Positive values found = %d" % count )

Traits.to_csv(file_out, na_rep='NA', sep='\t', index=False)
# end-of-script
