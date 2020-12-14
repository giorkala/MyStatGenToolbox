#!/bin/bash
#####
# This script gets an annotated set of variants, selects specific classes, and then creates 
# a list of carriers for each gene. It is structured for per-chromosome parallelisation.
#####

CHR=$1 # the chr to be analysed
THRES=0.01 # MAF threshold
PREFIX=deleterious
fileout=chr$CHR.gene_carriers_$PREFIX.txt # where we'll save results
rm -f $fileout

WESPATH=/well/palamara/projects/UKBB_APPLICATION_16549/EXOME_SEQUENCING_DATA

grep -w "chr$CHR" chrALL.genes.txt | while read chr GENE;
do
	# the next selects all the LoF variants -- that have "canonical=YES" -- for the current gene
	cat RESULTS_VEP/vep_$PREFIX\_chr$CHR.txt | grep $GENE | cut -f3 | while read snp;
	do 
		 grep -w $snp $WESPATH/freq/chr$CHR\_v1.FE.hg19.frq | awk -v up=$THRES '($5<up){print $2}'; # filter wrt THRES and get the "original" variant-ID
	done | sort -u > $GENE.vars 

	# use plink to create a ped file (ie, samples x variants) for the markers included in the gene:
	./plink_v1.9 --bfile $WESPATH/chr$CHR\_v1.FE.hg19 --extract $GENE.vars --recode 01 ped --output-missing-genotype . --out $GENE
	# now use awk to see who carries a mutation for this particular gene
	printf "%s" $GENE >> $fileout; # initialise row
	# for each individual (one row in the ped file) sum the genotypes and check the total
	awk '{x=$1; s=0; for(i=7; i<=NF; i++){ if($i==0){s=s+1}} if(s>0){printf "\t%s",x}}' $GENE.ped >> $fileout
	echo "" >> $fileout ; #newline
	printf "Variants in gene = "; wc -l $GENE.vars # this will save the number of variants in the log-file
	rm $GENE.* # $GENE.map $GENE.log $GENE.nosex
done

#for C in {1..22}; do qsub -N GxC.chr$C -cwd sge.gene_carriers.sh $C; done

#then do the following to create one carrier-gene index:
#cat chr*LOF* > chrALL.genes_carriers_LOF.txt
#python create_carrier_gene.py chrALL.genes_carriers_LOF.txt > get_carrier_gene.log
