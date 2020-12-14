#!/bin/bash
#####
# This script can massively submit jobs for IBD imputation using the `fastsmc-imputation` tool.
# Here we assume that IBD_SUBSET contains all the segments between (WES,nonWES) pairs. 
#####


## Remember: there's a bug with fastsmc-imputation for imputing WES samples; use '-n' to avoid that
## we need to (un) load the following for each node:
module unload GCC; module unload gcc; module load GCCcore/9.2.0

FOUT=IBD_IMPUTATION/DEL_NOWES # where to save the results

### $1,$2,$3 = CHR, START, STOP
date
mkdir -p $FOUT/$1.$2
D=$4
grep -w $1 genes.txt.bed | awk -v L=$2 -v R=$3 '!( L>$3 || R<$2 ){print $0}' | sed 's/chr//g' > $FOUT/$1.$2/genes.$1.$2.txt
if [[ ! -e $FOUT/$1.$2/matrix.dec$D.$1.$2.gz ]]
then
	# get the threshold for current chrom and decile
	THRES=$(grep $1 IBD_deciles.tab | grep $2 | awk -v d=$D '{print $(d+3)}' )
	echo "Working for $1,$2,$3,decile $D with $THRES" 
	for file in IBD_SUBSET/$1.from$2.to$3.*.bibd.gz;
	do
		./converter2.out $file 
        done | ./fastsmc-imputation -n -g $FOUT/$1.$2/genes.$1.$2.txt -c $1 -t $THRES -o $FOUT/$1.$2/matrix.dec$D.$1.$2.gz -l carriers_genes.txt; 
fi

#########################
# To resubmit broken jobs:
# grep Error LOGS/imp9*304* | cut -d: -f1 | while read x; do echo $x | sed 's/LOGS\/imp//g' | awk -F. '{print $2,$3,$4,$1}'; done > redo.txt
# and then 
# while read -r CHR START STOP D; do qsub -N imp$D.$CHR.$START.$STOP -P palamara.prjc -q short.qc -pe shmem 2 -cwd -j y -o LOGS/ -e error.$CHR.$START.txt sge.imputation.sh $CHR $START $STOP $D; done < redo.txt

# To impute, run the following
# awk 'NR>=1{print}' ChromRegions.txt | while IFS=',' read -r CHR START STOP temp; do for D in {9..0}; do qsub -N imp$D.$CHR.$START.$STOP -P palamara.prjc -q long.qc -pe shmem 2 -cwd -j y -o LOGS/ -e error.$CHR.$START.txt sge.imputation.sh $CHR $START $STOP $D; done; done
