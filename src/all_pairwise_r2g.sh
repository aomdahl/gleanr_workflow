#!/bin/bash
set -e

while getopts ":h" option; do
   case $option in
      h) # display Help
          echo "Simple script to generate commands to estimate R2g using LDSC pairwise across many traits."
         echo "Arguments are:"
          echo -e "\t1: File listing every summary stats file you want to process pairwise"
          echo -e "\t2: the destination output directory. PLease use the FULL path."
          echo -e "\t3: Number of output files to write to."
          echo "Example command:"
          echo "bash src/all_pairwise_r2g.sh ldsr_all/pairwise.traits.txt /data/abattle4/aomdahl1/snp_networks/scratch/udler_td2/ldsr_all/ldsr_results/ 10"
         exit;;
   esac
done


 this_sucks(){
      echo `expr $1 % $2`
 }


ODIR=$2
FILE=$1
NUMO=`echo $3`
OFILE=${ODIR}ldsr.runscript.sh
#Get of list of files to run on...
N=`wc -l $FILE| awk '{print $1}'`
rm -f runames.tmp
#echo $N

for ((i=1;i<=$N;i++)); do
     #drop the current one we want
     sed -e "${i}d" $FILE | tr '\n' ',' | sed 's/,$//g' >> runames.tmp
     echo "" >> runames.tmp
done
paste $FILE runames.tmp | tr '\t' ',' > commandnames.tmp

#Get the output file names based on the reference GWAS
rm -f runids.tmp
while read p; do
     	basename $p | awk -F "." '{print $(NF-2)}' >> runids.tmp
	#basename $p | cut -f 1 -d "." >> runids.tmp
done < $FILE


#Options here- to write out to individual files each, or to a certain number of files

if [ -n "$NUMO" ] && [ "$NUMO" -eq "$NUMO" ] 2>/dev/null; then
	 echo "Splitting into a number of files"

	for ((i=0;i<$NUMO;i++)); do
	     echo "cd /data/abattle4/aomdahl1/reference_data/ldsc_ref/" > ${ODIR}_ldsc.run.${i}.sh
	     echo "set -euo pipefail" >> ${ODIR}_ldsc.run.${i}.sh
	     echo "ml anaconda; source activate py27" >> ${ODIR}_ldsc.run.${i}.sh
	done

	#Counts per file:
	 echo $N
	 echo $NUMO


	LN=1
	while read p; do
	     BLAH=`this_sucks $LN $NUMO`
	     
	     QUERY=`echo $p | sed -r 's/\s+//g'`
	     ID=`sed -n "${LN}p" runids.tmp`
	     echo "python2  ~/.bin/ldsc/ldsc.py \
	    --rg $QUERY \
	    --ref-ld-chr eur_w_ld_chr/ \
	    --w-ld-chr eur_w_ld_chr/ \
	    --out ${ODIR}${ID}" >> ${ODIR}_ldsc.run.${BLAH}.sh
	    LN=$((LN+1))
	done < commandnames.tmp
	#echo "cd -"
else
  echo "Dividing by file"
 LN=1  
  while read p; do
        ID=`sed -n "${LN}p" runids.tmp`
	echo "Prepping $ID"
	echo "cd /data/abattle4/aomdahl1/reference_data/ldsc_ref/" > ${ODIR}/${ID}_ldsc.run.sh
        echo "set -euo pipefail" >> ${ODIR}/${ID}_ldsc.run.sh
        echo "ml anaconda; source activate py27" >> ${ODIR}/${ID}_ldsc.run.sh
        QUERY=`echo $p | sed -r 's/\s+//g'`
        echo "python2  ~/.bin/ldsc/ldsc.py \
    --rg $QUERY \
    --ref-ld-chr eur_w_ld_chr/ \
    --w-ld-chr eur_w_ld_chr/ \
    --out ${ODIR}${ID}" >> ${ODIR}/${ID}_ldsc.run.sh
    LN=$((LN+1))
done < commandnames.tmp
fi

	rm runames.tmp
	rm runids.tmp
	rm commandnames.tmp
	#rm topnames.tmp

#Helper script:
