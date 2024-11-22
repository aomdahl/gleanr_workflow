#!/bin/bash
#Script to map RSIDs to chr locations
RSID=$1
OUT=$2
#awk '(FNR == NR) {arr[$1];next} ($3 in arr) {print $1":"$2":"$4":"$5"\t"$3}' $RSID /scratch16/abattle4/lab_data/hg19/variant_calls/rsid_chr_pos.txt > $OUT
awk '(FNR == NR) {arr[$1];next} ($3 in arr) {print $1":"$2":"$4":"$5"\t"$3}' $RSID /scratch16/abattle4/lab_data/hg19/variant_calls/rsid_chr_pos.SNVs_only.txt > $OUT
