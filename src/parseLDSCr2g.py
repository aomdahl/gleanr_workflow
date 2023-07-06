


from pickle import TRUE
import re
import argparse

def getRelevantData(line, din, tn):
    if "Intercept" in line:
        din["intercept"][tn], din["intercept_se"][tn] = getLineQuery(line, "Intercept")
    elif "Lambda GC" in line:
        din["lambda_gc"][tn] = getLineQuery(line, "Lambda GC")
    elif "Mean Chi^2" in line:
        din["meanchi2"][tn] = getLineQuery(line, "Mean Chi^2")
    elif "Ratio" in line:
        din["ratio"][tn], din["ratio_se"][tn] = getLineQuery(line, "Ratio")
    else:
        return din
    return din

#major logical fallacy here you idiot

def getLineQuery(l, lookup):
    if lookup in l:
        d = l.split()
        if lookup in ["Intercept", "Ratio"]:
            if "<" in line:
                return 0, float('NaN')
            g = re.search("\w+: ([\d.NAe-]+) \(([\d.NA]+)\)",l)
            try:
                return g.group(1), g.group(2)
            except AttributeError:
                print("Unexpected..")
                print(g)
                print(l)
                print(lookup)
                return float('NaN'),float('NaN')
        elif lookup in ["Lambda GC", "Mean Chi^2"]:
            g =re.search("\w+: ([\d.NAe-]+)", l)
            return float(g.group(1))
        else:
            print('Failing here...')
            return float('NaN')
    else:
        if lookup in ["Intercept", "Ratio"]: return [float('NaN'),float('NaN')]
        else:
            return float('NaN')


def getTraitIndex(l):
    #Use regex to get the phenotpe number
    #Heritability of phenotype 1 or
    #Heritability of phenotype 2/61
    try:
        return int(re.search("phenotype (\d+)\/*", l).group(1)) -1
    except:
        print("Unable to find match...")
        return -1



parser = argparse.ArgumentParser(description = "Quick script to read in files and make ")
parser.add_argument("--input_file", help = "Path to log file output from LDSC rg")
parser.add_argument("--output", help = "output path. Name will append what is currently there and swap .log for .report.csv")
args = parser.parse_args()
import os
#query_file = "ldsr_results/CARDIoGRAM_GWAS_RESULTS.log"
query_file = args.input_file
with open(query_file, 'r') as istream:
    rel_dat = False
    for line in istream:
        line = line.strip()
        if "--rg" in line: #this has the order
            trait_order = [os.path.basename(x).replace(" ", "").replace("\\", "") for x in line[6:].split(',')]
            out_dict = {
                "study":trait_order, 
                "lambda_gc" : [""]*len(trait_order),
                "intercept":[""]*len(trait_order), 
                "intercept_se":[""]*len(trait_order),
                "meanchi2":[""]*len(trait_order),
                "ratio":[""]*len(trait_order),
                "ratio_se":[""]*len(trait_order),
                "snp_overlap":[""]*len(trait_order)
            }
        if "Computing rg for phenotype" in line:
            upcoming_trait = getTraitIndex(line)
        if "SNPs with valid" in line:
            snp_count = line.split()[0]
        if "Heritability of " in line:
            which_trait = getTraitIndex(line)
            #print("On trait", which_trait)
            rel_dat = True
            if upcoming_trait == which_trait or upcoming_trait == 1:
                out_dict["snp_overlap"][which_trait] = snp_count
            continue
        if rel_dat:
            out_dict = getRelevantData(line, out_dict, which_trait)
            if out_dict["ratio"][which_trait] != "":
                rel_dat = False

import pandas as pd
o = pd.DataFrame(out_dict)
oname = os.path.basename(args.input_file).replace(".log", ".ldsc_report.csv")
print("Writing out to", args.output + oname)
o.to_csv(args.output + oname , index= False)

"""
Helpful script to use it in

function fileName()
{
#from https://stackoverflow.com/questions/965053/extract-filename-and-extension-in-bash
    fullfile=$1
    filename=$(basename -- "$fullfile")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo $filename
}
set -eo pipefail
P=/work-zfs/abattle4/ashton/snp_networks/scratch/udler_td2/src
for i in ldsr_results/*.log; do
    echo $i
    python $P/parseLDSCr2g.py --input_file $i --output ldsr_results/
    b=`fileName $i`
    #grep -A 1000 "gcov_int_se" ldsr_results/5983_irnt.log  | grep -v "Analysis finished" | grep -v "Total time" > ldsr_results/$b.rg_report.tsv
done
"""