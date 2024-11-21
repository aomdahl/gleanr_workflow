#!/bin/bash

wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/full_variant_qc_metrics.txt.bgz

#Get just the European high quality data:
#Filter include:
#11: info score > 0.9
#34: MAF > 0.01
#9: high_quality
zcat full_variant_qc_metrics.txt.bgz |  awk '((NR == 1) || (($11 + 0.0 > 0.9) && ($34 + 0.0 > 0.01) && ($9 == "true"))) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$14"\t"$20"\t"$26"\t"$32"\t"$38}' | gzip > high_quality_common_variants_EUR.txt.bgz
