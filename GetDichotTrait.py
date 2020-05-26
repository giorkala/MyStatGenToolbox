# -*- coding: utf-8 -*-
"""
Created on Tue May 26 18:24:27 2020
@author: giork

This gets an ICD10 disease-class code, merges all the subphenotypes, and creates a new file with (ID, disease) tuples.
E.g. for I25, it will make a variable indicating every I250-I259 case.
We assume the file with all ICD10 codes as well as a "mapping" file are in the same folder.

TODO: support merging of only a few sub-phenotypes
TODO: sample selection?
"""
import sys 
import pandas as pd

ToF = sys.argv[1]
if len(sys.argv)>2:
    fout = sys.argv[2]
elif len(sys.argv)==1:
    print("\tError: run `python GetDichotTraits.py ICD10_code (output_filename)")
else:
    fout="Pheno_"+ToF+".txt"
    
print("Extracting UKBB cases for "+ToF+". Results will be saved as "+fout)

Traits = {}
with open("Pheno_ICD10.tab", 'r') as fin:
    _=next(fin)
    for line in fin:
        tokens=line.split()
        Traits[ int(tokens[0]) ]=tokens[1:]
    
print("Data for",len(Traits),"samples were loaded!")

Mapping=pd.read_csv("mapping_biased.csv", sep=',')
Samples = [x for x in list(Mapping['ukb43206']) if x in Traits ]
print("Subjects to keep = %d" % len(Samples))

total=0
with open(fout,'w') as f:
    f.write("ID\t"+ToF+"\n")    
    for x in Samples:
        T = False
        for s in Traits[x]:
            T = T or s.startswith(ToF)
        value = '1' if T else '0'
        total += int(value)
        f.write(str(x)+"\t"+value+"\n")
        
print(total,"cases were written. Bye!")
# end-of-script        
    