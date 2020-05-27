# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 17:38:58 2020

This program creates a gene-variants index for one (given) chromosome.
It avoids an extensive O(VxG) search by using a double pointer ("leap-frogging")

"""

import gzip
import sys

if len(sys.argv)< 3:
    print( "Error in the number of arguments! Run:")
    print( "\tpython CreateGroupFiles CHR genes.txt.bed VCF" )
    sys.exit()
    
CHR=sys.argv[1]
filetosave="GROUP_files/groupFile."+CHR+".txt"
    
Genes = {}; g=0
# gene-ID maps to (name, start, stop, list-of-variants)
with open(sys.argv[2],'r') as file:
    for line in file:
        tokens = line.split()
        # chr - start - stop - name
        if (tokens[0]==CHR):
            Genes[ g ] = ( tokens[3], int(tokens[1]) , int(tokens[2]), [] )
            g +=1

Ngenes=len(Genes)
print( "The index for {0} genes in {1} was initialised!".format(Ngenes,CHR) )

g=0 # current gene index
null_genes=0
with open(filetosave, 'w') as file_out:

    with gzip.open(sys.argv[3],'r') as file:     
        for _ in range(7):
            next(file) # skip the header

        for line in file:
            tokens = line.split() # chr - pos - ID - ...
            variant=int( tokens[1] )
            # search 
            if (Genes[g][1] <= variant ) & (Genes[g][2] >= variant):
                temp="{0}:{1}_{2}/{3}".format(tokens[0].decode("utf-8"), tokens[1].decode("utf-8"), tokens[3].decode("utf-8"), tokens[4].decode("utf-8")) ## chr:pos_ref/alt
                Genes[g][3].append( temp ) 
            elif Genes[g][2] < variant :
                # variant lies after the current gene; finish processing it
                print( "Variants found for {0} = {1}".format( Genes[g][0], len(Genes[g][3]) ) )
                if len( Genes[g][3] )==0:
                    null_genes +=1
                else: # start writing
                    file_out.write( Genes[g][0] +'\t'.join(Genes[g][3])+'\n')

                # move to next gene - if any!
                while (g in Genes) & (g<=Ngenes-1):
                    # we need the first condition to avoid an error for the first variant after the last gene.  
                    # for all the other variants, this is probably just one step unless there are genes with 0 variants                 
                    if Genes[g][2] < variant:
                        g += 1 
                    else: # we found the next gene
                        break
                # now the var is a)within g b)before g c) after every gene
                if g==Ngenes:
                    print( "Variants exceeded the right-most gene!")                
                    break
                elif Genes[g][1] < variant:
                    # we know that Genes[g][2] >= variant
                    temp="{0}:{1}_{2}/{3}".format(tokens[0].decode("utf-8"), tokens[1].decode("utf-8"), tokens[3].decode("utf-8"), tokens[4].decode("utf-8")) ## chr:pos_ref/alt
                    Genes[g][3].append( temp )
#                else: continue # variant lies before the new gene
#            else: continue # variant lies before current gene
                
    if g==Ngenes-1:
        # this is when the last variant lies within the last gene
        print( "Variants found for {0} = {1}".format( Genes[g][0], len( Genes[g][3])) )
        file_out.write(Genes[g][0]+'\t'.join(Genes[g][3])+'\n')
        print( "Variants are finished!")
             
                    
print("Total number of variants found within genes = {0}".format(sum([len(Genes[g][3]) for g in Genes ])) )
print("Genes without any variants = %d" % null_genes)
# end-of-script 
