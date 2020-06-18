#!/usr/bin/env python3
"""
This is for pre-processing the raw phenotypic data of UKBB. There is a bottleneck/tricky-part  
related to fields with arrays where all the entries must be parsed and aggregated.
INPUT : Original files with (1) traits and (2) PCA, and (3) file with ID mapping
OUTPUT: Tab delimited file of covariates+PCAs, with new index and sample-selection.
"""
import pandas as pd
import numpy as np
from scipy import stats


file_out="Covariates.tab"

ukbb_data = pd.read_csv("UKBB_Pheno35.csv", sep=',')
ukbb_data = ukbb_data.set_index('eid')
PCA40 = pd.read_csv("UKBB_PCA.csv")
PCA40 = PCA40.set_index('eid')

print("Data loaded - Total number of variables found = %d" % ukbb_data.shape[1])
print("              Total number of subjects  found = %d" % ukbb_data.shape[0])

Mapping = {}
with open("mapping.csv", 'r') as f:
    _=next(f) # skip header
    for line in f:
        tokens=line.split(sep=',')
        # map UKBB-ID to IBD-ID and get rid of "\n":
        Mapping[ int(tokens[2][:-1]) ] = (tokens[0], tokens[1]) # Mapping has int -> (str, str)
Samples = [x for x in Mapping if x in ukbb_data.index ] #
print("Subjects to keep = %d" % len(Samples))
print("\t Remember to check that fields with arrays are properly processed!")

Covariates = pd.DataFrame()
Covariates['IID'] = [ Mapping[s][0] for s in Samples ]
Covariates['UKBB_1'] = [ Mapping[s][1] for s in Samples ]
Covariates['UKBB_2'] = Samples # this is just for double-checking

Covariates['Sex'] = [ ukbb_data['31-0.0'][s] for s in Samples ]
Covariates['Heel_bmd'] = [ ukbb_data['3148-0.0'][s] for s in Samples ]
Covariates['Heel_fractured'] = [ ukbb_data['3082-0.0'][s] for s in Samples ]
Covariates['Foot_used'] = [ ukbb_data['3081-0.0'][s] for s in Samples ]
Covariates['Ethnic'] = [ ukbb_data['22006-0.0'][s] for s in Samples ]
Covariates['QC_failure'] = [ ukbb_data['22027-0.0'][s] for s in Samples ]
#Covariates['Age'] = [ ukbb_data['21003-0.0'][s] for s in Samples ]

print("Proceeding with array-like fields...")
Traits_Dict = {
#        "heel_bmd":"3081-",
#        "sex":"31-",
#        "testing_age":"21022-",
        "smoking_status":"20116-",
        "height":"50-",
        "weight":"21002",
        "BMI":"21001",
#        "foot_used":"3081-",
#        "heel_fracture":"3082-",
        "assessment_site":"54-",
        "menopausal":"2724",
        "age":"21003"
       }

problematic=[]
Age = []; Smoking = []; Height=[]; Weight=[]; BMI = [];  Site = []; Menopausal = [];
index = 0
for s in Samples:

    temp_age = []; temp_BMI = []; temp_smok = []; temp_site = []
    temp_height = []; temp_weight = []; temp_menop = []

    for i in range(4):
        # include age - mean(entries)
        col = '21003-'+str(i)+'.0'
        if ~np.isnan( ukbb_data[col][s]):
            temp_age.append(ukbb_data[col][s])
        # include study site - take the latest entry
        col = '54-'+str(i)+'.0'
        if ~np.isnan( ukbb_data[col][s]):
            temp_site.append(ukbb_data[col][s])
        col = '21001-'+str(i)+'.0'
        # include BMI - the average of the given values, if any!
        if ~np.isnan( ukbb_data[col][s]):
            temp_BMI.append(ukbb_data[col][s])
        # include smoking status - round(mean(entries))
        col = '20116-'+str(i)+'.0'
        if ~np.isnan( ukbb_data[col][s]):
            temp_smok.append(ukbb_data[col][s])
        # include height - mean(entries)
        col = '50-'+str(i)+'.0'
        if ~np.isnan( ukbb_data[col][s]):
            temp_height.append(ukbb_data[col][s])
        # include weight - mean(entries)
        col = '21002-'+str(i)+'.0'
        if ~np.isnan( ukbb_data[col][s]):
            temp_weight.append(ukbb_data[col][s])
        # include menopausal status - take the latest entry
        col = '2724-'+str(i)+'.0'
        if ~np.isnan( ukbb_data[col][s]):
            temp_menop.append(ukbb_data[col][s])
    # check if any of these was empty - at least one would be a problem
    if len(temp_BMI)*len(temp_site)*len(temp_smok)==0:
        problematic.append(s)#
    # the next lists all follow the same order as `Samples`
    Smoking.append( stats.mode(temp_smok)[0][0] if len(temp_smok)>0 else [] )
    BMI.append( np.mean(temp_BMI) )
    Site.append( stats.mode(temp_site)[0][0] if len(temp_site)>0 else [] )
    Height.append( np.mean(temp_height) )
    Weight.append( np.mean(temp_weight) )
    Age.append( np.nanmean(temp_age) )
    if len(temp_menop)==0:
        Menopausal.append( -3 ) # transform NaNs to "not answered" 
    else:
        Menopausal.append( temp_menop[-1] )
    
    index += 1 # this runs 0->Nsubjects (similarly to IBD)

print("%d subjects with missing values were found!" % len(problematic) )

Covariates['Age'] = Age
Covariates['BMI'] = BMI
Covariates['Smoking'] = Smoking
Covariates['Site'] = Site
Covariates['Height'] = Height
Covariates['Weight'] = Weight
Covariates['Menopausal'] = Menopausal

print("Getting PCAs...")
for i in range(4):
    Covariates['PCA'+str(i+1)] = [PCA40['22009-0.'+str(i+1)][x] for x in Samples]

print("Writing file to disk...", end='')
#Covariates = Covariates.sort_values(by=['IID'])
Covariates.to_csv(file_out, sep='\t', index=False)
print(" Bye!")
# end-of-script
