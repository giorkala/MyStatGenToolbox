#!/bin/bash
#####
# This script can either run VEP to get annotations for sequenced variants (UKBB_WES in particular), or parse previously downloaded annotations and filter variants for next steps.
#####

# $1 is CHR 
C=$1
module load Bio-DB-HTS/2.11-foss-2018b-Perl-5.28.0
module load BioPerl/1.7.2-foss-2018b-Perl-5.28.0      

cd /well/palamara/projects/EXOME_IBD/VEP

#.vep/ensembl-vep-release-101/vep -i VCF_FILES/VEP_WES_chr$C.vcf --offline --biotype --distance 5000 --force_overwrite --sift b --polyphen b --canonical --symbol --tab --fields Uploaded_variation,SYMBOL,CANONICAL,SIFT,POLYPHEN,Consequence -o RESULTS_VEP/vep_new_chr$C.txt
#.vep/ensembl-vep-release-101/filter_vep -i RESULTS_VEP/vep_chr$C.vcf -filter "CANONICAL is YES and SIFT is deleterious and PolyPhen is probably_damaging" | grep -v "#" > RESULTS_VEP/vep_filtered_chr$C.txt
#.vep/ensembl-vep-release-101/filter_vep -i RESULTS_VEP/vep_new_chr$C.txt -filter "CANONICAL is YES and (Consequence is frameshift_variant or Consequence is splice_acceptor_variant or Consequence is splice_donor_variant or Consequence is stop_gained)" | grep -v "#" > RESULTS_VEP/vep_LOF_chr$C.txt

# the next commands prepare a set of LOF variants as those used in the Saada et al. (2020) analysis:
#cat RESULTS_VEP/vep_chr$C.vcf | egrep 'frameshift_variant|splice_acceptor_variant|splice_donor_variant|stop_gained|stop_lost|start_lost' | cut -f1-5,8 > temp.chr$C
#paste <(cut -f1-5 temp.chr$C) <(cut -f6 temp.chr$C | cut -d\| -f4 ) > RESULTS_VEP/vep_LOF_chr$C.txt
#rm temp.chr$C

# the next select a broader set of deleterious variants:
cat RESULTS_VEP/vep_chr$C.vcf | grep "deleterious" | grep "probably_damaging" | cut -f1-5,8 > temp.chr$C # this is the same as using the filter_vep command above
cut -f1-5,8 RESULTS_VEP/vep_filtered_chr$C.txt > temp.chr$C
paste <(cut -f1-5 temp.chr$C) <(cut -f6 temp.chr$C | cut -d\| -f4 ) > RESULTS_VEP/vep_newdeleterious_chr$C.txt
rm temp.chr$C

####
# then run the following
#for C in {1..22}; do qsub -N vep.chr$C -pe shmem 1 -cwd sge.vep.sh $C; done
