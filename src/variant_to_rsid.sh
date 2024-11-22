#!/bin/bash

ID=$1
OUT=$2
awk '(FNR == NR) {arr[$1];next} ($1":"$2 in arr) {print $1":"$2"\t"$3}' $ID /scratch16/abattle4/lab_data/hg19/variant_calls/rsid_chr_pos.SNVs_only.txt > $OUT
