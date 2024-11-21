#!/bin/bash

buildEverything  ()
{
    mkdir -p extract_EUR_scripts
    FIN=$1
    ODIR_=$2
    RDIR=$3
    OUT="$ODIR_/${FIN%.*}.EUR_ONLY.gz"
    echo "cd $RDIR" > extract_EUR_scripts/${FIN%.*}.extract.sh
    ARG=`zcat $FIN | head -n 1 | tr '\t' '\n' | nl | grep "EUR" | cut -f 1 | tr '\n' ' '`
    AWKBUILD="1 2 3 4 arr[\$1\":\"\$2\"_\"\$3\"_\"\$4] ${ARG}"
    echo -n "awk '(FNR == NR) {arr[\$6]=\$5; next} (\$1\":\"\$2\"_\"\$3\"_\"\$4 in arr) {print " >> extract_EUR_scripts/${FIN%.*}.extract.sh
    for i in $AWKBUILD; do
    echo -n "\$$i\"\t\""  >> extract_EUR_scripts/${FIN%.*}.extract.sh
    done
    echo " }' <( zcat ../high_quality_common_variants_EUR.txt.bgz) <(zcat $FIN)| gzip > $OUT" >> extract_EUR_scripts/${FIN%.*}.extract.sh
    sed -i 's/$arr/arr/g' extract_EUR_scripts/${FIN%.*}.extract.sh
    #######3
    echo "" > extract_EUR_scripts/${FIN%.*}.extract.2.sh
    ARG=`zcat $FIN | head -n 1 | tr '\t' '\n' | nl | grep "EUR" | cut -f 1 | tr '\n' ' '`
    AWKBUILD="1 2 3 4 ${ARG}"
    echo -n "zcat $FIN | head -n 1 | awk '{print " >> extract_EUR_scripts/${FIN%.*}.extract.2.sh
    for i in $AWKBUILD; do
    echo -n "\$$i\"\t\""
    done >> extract_EUR_scripts/${FIN%.*}.extract.2.sh
    echo " }' | gzip > ${FIN%.*}.tmp.gz" >> extract_EUR_scripts/${FIN%.*}.extract.2.sh
    sed -i 's/\$4\"\\t\"/\$4\"\\tRSID\\t\"/g' extract_EUR_scripts/${FIN%.*}.extract.2.sh
    cat extract_EUR_scripts/${FIN%.*}.extract.sh extract_EUR_scripts/${FIN%.*}.extract.2.sh >> ${FIN%.*}.extract.TMP && mv ${FIN%.*}.extract.TMP extract_EUR_scripts/${FIN%.*}.extract.sh
    rm extract_EUR_scripts/${FIN%.*}.extract.2.sh
    #### join them
    echo "" >> extract_EUR_scripts/${FIN%.*}.extract.sh
    echo "cat  $OUT >> ${FIN%.*}.tmp.gz  && mv ${FIN%.*}.tmp.gz $OUT" >> extract_EUR_scripts/${FIN%.*}.extract.sh
}
  DIR=/data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/phenotype_flat_files
 cd $DIR
  ODIR="$DIR/EUR_only/"
  for i in *.bgz; do
    echo $i
    buildEverything $i $ODIR $DIR
  done

#Updated version- take a list of files to input
#while read i; do
#	echo $i
#	buildEverything $ODIR $DIR
#done < $1
