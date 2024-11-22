#!/bin/bash
#Removes multiallelic snps, indels, selects only those in 1000G, removes ambiguous SNPs
set -e
SUFF=$1
INPUT=$2
OUTPUT=$3
THOU_G="/scratch16/abattle4/ashton/snp_networks/reference_files/thouG_reference_snps.txt"
UK_SNPS="/scratch16/abattle4/ashton/snp_networks/reference_files/ukbb_snp_list.tsv"

function getFileName() {     FILE=$1;     f="$(basename -- $FILE)";     echo ${f%.*}; }

FN=`getFileName $INPUT`
#To remove, need the list
#echo "Selecting those in UKBB"
cut -f 1 $INPUT > partial.${FN}.${SUFF}.tmp
#awk -F ":" '(FNR == NR) {arr[$1":"$2]; next} ($1":"$2 in arr) {print $0}' partial.${FN}.${SUFF}.tmp $UK_SNPS > full_list.${FN}.${SUFF}.tmp
#we assume the snps were selected from UKBB to begin with
echo "Removing multi-allelic snps"
awk -F ":" '{print $1":"$2}' partial.${FN}.${SUFF}.tmp | uniq -d  > remove.${FN}.${SUFF}.tmp #identify multiallelic
sed -i 's/\t//g' remove.${FN}.${SUFF}.tmp
echo "FAKE:SPOTHOLDER" >> remove.${FN}.${SUFF}.tmp

awk -F ":" '(FNR == NR) {arr[$1":"$2];next} !($1":"$2 in arr) {print $0}' remove.${FN}.${SUFF}.tmp partial.${FN}.${SUFF}.tmp | cut -f 1 -d ","> filtered.${FN}.${SUFF}.tmp #no_multi-allelic. Being very lazy here and just keeping thefirst one in the list....,. best practice would be to include both.
#Get just those in 1000G
echo "Filtering to those in 1000G"
awk -F ":" '(FNR == NR) {arr[$1":"$2":"$3":"$4];next} ( $1":"$2":"$3":"$4 in arr || $1":"$2":"$4":"$3 in arr ) {print $0}' $THOU_G filtered.${FN}.${SUFF}.tmp > count_try.${FN}.${SUFF}.tmp #just in 1000G



echo "Removing ambiguous snps"
awk -F ":" '( (($3=="A"&&$4=="C") || ($3=="A"&&$4=="G") || ($3=="C"&&$4=="A") || ($3=="C"&&$4=="T") || ($3=="G"&&$4=="A") || ($3=="G"&&$4=="T") || ($3=="T"&&$4=="C") || ($3=="T"&&$4=="G"))){print $0}' count_try.${FN}.${SUFF}.tmp  > noambig.${FN}.${SUFF}.tmp
cut -f 1,2 -d ":" noambig.${FN}.${SUFF}.tmp > noambig.intermediate.${SUFF}.tmp && mv noambig.intermediate.${SUFF}.tmp noambig.${FN}.${SUFF}.tmp 

echo "Removing those in the HLA region"
#coordinates basd on Genome Reference Consortium, https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37
echo "WARNING: coordinates based only hg37, please update if using a different genome build."
awk -F ":" '!($1 == 6 && $2+0.0 > 28477796 && $2+0.0 < 33448354) {print $0}' noambig.${FN}.${SUFF}.tmp > $OUTPUT
rm *.${FN}.${SUFF}.tmp


