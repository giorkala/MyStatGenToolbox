# -*- coding: utf-8 -*-
"""
Created on Tue May 26 18:24:27 2020
@author: giork

This gets an ICD10 disease-class code, merges all the subphenotypes, and creates a new file with (ID, disease) tuples.
E.g. for I25, it will make a variable indicating every I250-I259 case.
We assume the file with all ICD10 codes as well as a "mapping" file are in the same folder.

TODO: support merging of only a few sub-phenotypes
18/6: sample selection and writing of both UKBB-1 and IBD indices.
"""
import sys 

if len(sys.argv)>3:
    ToF = sys.argv[2]
    fout = sys.argv[3]
elif len(sys.argv)<3:
    print("\tError: run `python GetDichotTraits.py mapping.csv ICD10_code (output_filename)")
    sys.exit()
else:# no filename for output was given
    ToF = sys.argv[2]
    fout="Pheno_"+ToF+".txt" 
    
print("Extracting UKBB cases for "+ToF+". Results will be saved as "+fout)

Traits = {}
with open("UKBB_ICD10.tab", 'r') as f:
    _=next(f) # skip header
    for line in f:
        tokens=line.split()
        Traits[ tokens[0] ]=tokens[1:]
    
#print("Data for",len(Traits),"samples were loaded!")

Mapping = {}
with open(sys.argv[1], 'r') as f:
    _=next(f) # skip header
    for line in f:
        tokens=line.split(sep=',')
        Mapping[ tokens[2][:-1] ] =  (tokens[0], tokens[1]) # map UKBB-2 to (IBD-ID,UKBB-1) and get rid of "\n"
Samples = [x for x in Mapping if x in Traits ]
print("Subjects to keep = %d" % len(Samples))

total=0
with open(fout,'w') as f:
    f.write("UKBB_1\tIBD_ID\t"+ToF+"\n")    
    for x in Samples:
        T = False
        for s in Traits[x]:
            T = T or s.startswith(ToF)
        value = '1' if T else '0'
        total += int(value)
        f.write(Mapping[x][1]+"\t"+Mapping[x][0]+"\t"+value+"\n")
        
print(total,"cases were written. Bye!")
# end-of-script        
