#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 14:46:18 2019 
update (5/5/2020): no need for genes-file, just variant-wise plotting
run as `python PlotManhattan.py Trait SummaryStats`
The input file should contain (at least) fields for CHR (1-22), POS (int in bp), SNPID (string) and Pvalues. Just change line-30 accordingly.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys

method= "fastSMC" 
Pheno = sys.argv[1] 
FilePvalues= sys.argv[2]
# we assume this file is ordered wrt each snp

Chr_vars={}
for i in range(1,23):
    Chr_vars[str(i)]={} # each chrom will be a dict of ID-> (pos, pvalue)

totalpvalues=0
with open(FilePvalues, 'r') as file:
    next(file) # skip the header
    for line in file:
        tokens = line.split() 
        # CHR POS SNPID Allele1 Allele2 AC_Allele2 AF_Allele2 imputationInfo N BETA SE Tstat p.value p.value.NA Is.SPA.converge varT varTstar
        Chr_vars[ tokens[0] ][ tokens[2] ] = ( int(tokens[1]), float(tokens[12]) )  # (ID, POS, Pvalue)
        totalpvalues+=1

colors = ['black','blue']

fig, ax = plt.subplots( figsize=(12, 5) )
total=0
xticks=[]
for i in range(1,23):
    c =str(i)

    names = list( Chr_vars[c].keys() )
    loci = [ Chr_vars[c][x][0] for x in names ]
    pvalues = [ Chr_vars[c][x][1] for x in names ]

    indx = np.argsort(loci)
    pvalues = [ pvalues[i] for i in indx ]
    names = [ names[i] for i in indx ]

### this is for writing tags on top of the significant variants
#    for x in [ x for x  in indx if pvalues[x]<0.05/totalpvalues]:
#        plt.text(total+indx[x]-0.1, -1*np.log10(pvalues[x])+0.15, names[x], style='italic')
        
    plt.plot([total+x for x in indx], -1*np.log10(pvalues), '.', color=colors[i%2], label=c)
    plt.xlabel( "Chromosome",fontsize=16 )
    plt.ylabel('-log10(p) for '+Pheno,fontsize=14)

    xticks.append(total + len(loci)/2)
    total += len(indx) + 0.12*len(names) # just a trick to avoid overlappings

plt.axhline(y=-np.log10(0.05/totalpvalues), color='g', linestyle='--', alpha=0.6)
plt.xticks(ticks=xticks, labels=range(1,23)) 

plt.savefig("Manhattan_SAIGE_"+method+"_"+Pheno+".png") #, transparent=True
print("\t Total number of genes tested = "+str(totalpvalues))
