import argparse 
import numpy as np 
import scipy.stats as st
import glob
import gzip
parser = argparse.ArgumentParser(description = "Quickly select out a list of SNPs based on some p-value threshold, or an input list of SNPs.") 
parser.add_argument("--gwas_list", help = "list of gwas studies to extract from") 
parser.add_argument("--gwas_id", help = "Single gwas study name to extract from.", default = None)
parser.add_argument("--gwas_dir", help = "gwas directory to look in.") 
parser.add_argument("--output", help = "output path") 
parser.add_argument("--type", help = "if its ldsc or not", default = "NOT")
parser.add_argument("--pval", help = "pvalue threshold", type = float, default = 1e-5) 
parser.add_argument("--maf", help = "MAF threshold (greater than)", type = float, default = 0.01) 
parser.add_argument("--snp_list", help = "Optional: list of SNPs to include in the search", default = "") 
parser.add_argument("--extension", help = "File extension to grab, default is none.", default = "") 
args = parser.parse_args() 

if args.gwas_id is None:
    with open(args.gwas_list, 'r') as istream: 
        file_list = [x.strip().split()[0] for x in istream] 
else:
    file_list = [args.gwas_id]
matchlist = False
if args.snp_list != "":
    print("We are looking for a SNP list")
    matchlist = True
    with open(args.snp_list, 'r') as istream:
        snp_set = set([x.strip() for x in istream])
    print("SNP list has length", str(len(snp_set)))

print("Read in files to read") 
nstudies = len(file_list) 
var_set = set()
print("Designated threshold is associated with pvalue", args.pval)
for f in file_list:
    fname = args.gwas_dir + "/" + f + "*" + args.extension
    try:
        open_file = glob.glob(fname)[0]
        print(open_file)

    except:
        print(fname)
        print(glob.glob(fname))
        input()
    if args.type == "ldsc":
        #2 sided distribution, so look at half the pvalue
        adjp = args.pval/2.0
        thresh = abs(st.norm.ppf(adjp))
        #print("Designated threshold is a abs z score of", thresh)
    if (".gz" in args.extension) or (".gz" in f):
      file_obj = gzip.open(open_file, 'rt')
    else:
      file_obj = open(open_file, 'r')
    with file_obj as istream: 
        for line in istream:
            line_ = line.strip().split('\t')
            if line_[0]== "SNP" or line_[0] == "variant": 
              if args.type == "ldsc_custom": #The header names vary here
                pval_i=line_.index("P")
                if "FRQ" in line_:
                    maf_i=line_.index("FRQ")
                else:
                    maf_i = None
              continue #first line
            if args.type == "std":
                if abs(float(line_[-1])) < args.pval:#looking at a pvalue
                    if matchlist and line_[0] in snp_set: #we want to filter and there is the filter
                        var_set.add(line_[0])
                    if not matchlist:
                        var_set.add(line_[0])
            elif args.type == "ldsc":
                if(len(line_) < 2):
                  continue
                if abs(float(line_[-1])) > thresh: #looking at a zscore
                    if matchlist and line_[0] in snp_set: #we want to filter and there is the filter
                        var_set.add(line_[0])
                    if not matchlist:
                        var_set.add(line_[0])
            elif args.type == "ldsc_custom": #we modified the ldsc thing...
                if(len(line_) < 4):
                  continue #missing entry
                if float(line_[pval_i]) < args.pval: #looking at a zscore
                    if (maf_i is not None) and (float(line_[maf_i]) < args.maf):
                        continue
                    else:
                        if matchlist and line_[0] in snp_set: #we want to filter and there is the filter
                            var_set.add(line_[0])
                        if not matchlist:
                            var_set.add(line_[0])
            else:
                print("havent implemented yet")
                exit
    file_obj.close()
    #print(len(var_set))
    print("Finished with " + open_file)
with open(args.output, 'w') as ostream: 
    for var in var_set: 
        ostream.write(var + '\n')                      




        #python src/unionVariants.py --gwas_list ./trait_selections/seed1_thresh0.9_h2-0.05.studies.tsv  --output test.txt --type ldsc --pval 1e-5 --extension ".both_sexes.tsv" --gwas_dir /work-zfs/abattle4/lab_data/UKBB/GWAS_Neale/highly_heritable_traits_2/ldsr_format/unzipped/
        #python src/unionVariants.py --gwas_list ./trait_selections/seed1_thresh0.9_h2-0.05.studies.tsv  --output test.txt --type ldsc --pval 1e-5 --extension ".both_sexes.tsv" --gwas_dir /work-zfs/abattle4/lab_data/UKBB/GWAS_Neale/highly_heritable_traits_2/ldsr_format/unzipped/
