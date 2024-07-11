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
parser.add_argument("--output_counts", help = "output path", default = None) 
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
count_set = dict()
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
        #May not contain p-value but does contain z-score, so use that.
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
            if "SNP" in line.upper()  or "VARIANT" in line.upper():
              snp_i = 0
              if "SNP" in line.upper(): snp_i = line_.index("SNP")
              if "VARIANT" in line.upper(): snp_i = line_.index("variant")
              #print("Found the snp index", snp_i)
              if args.type == "ldsc_custom": #The header names vary here
                pval_i=line_.index("P")
                #print("found the pvalue index,",pval_i)
                if "FRQ" in line_:
                    maf_i=line_.index("FRQ")
                else:
                    maf_i = None
              continue #first line
            if args.type == "std":
                if abs(float(line_[-1])) < args.pval:#looking at a pvalue
                    if line_[0] not in count_set: #Update- track the number of times it hits.
                        count_set[line_[snp_i]] = 0
                    if matchlist and line_[snp_i] in snp_set: #we want to filter and there is the filter
                        var_set.add(line_[snp_i])
                        count_set[line_[snp_i]]+=1 #Update- track the number of times it hits.
                    if not matchlist:
                        var_set.add(line_[snp_i])
                        count_set[line_[snp_i]]+=1 #Update- track the number of times it hits.
            elif args.type == "ldsc":
                if(len(line_) < 2):
                  continue
                if abs(float(line_[-1])) > thresh: #looking at a zscore
                    if line_[snp_i] not in count_set: #Update- track the number of times it hits.
                        count_set[line_[snp_i]] = 0
                    if matchlist and line_[snp_i] in snp_set: #we want to filter and there is the filter
                        var_set.add(line_[snp_i])
                        count_set[line_[snp_i]]+=1 #Update- track the number of times it hits.
                    if not matchlist:
                        var_set.add(line_[snp_i])
                        count_set[line_[snp_i]]+=1 #Update- track the number of times it hits.
            elif args.type == "ldsc_custom": #we modified the ldsc thing...
                if(len(line_) < 4):
                  continue #missing entry
                if float(line_[pval_i]) < args.pval: #looking at a zscore
                    if line_[snp_i] not in count_set: #Update- track the number of times it hits.
                        count_set[line_[snp_i]] = 0
                    if "rs28463616" in line_[snp_i]:
                        print(line_)
                        print("Found the error")
                        input()
                    if (maf_i is not None) and (float(line_[maf_i]) < args.maf): #Omit by MAF
                        continue
                    else:
                        if matchlist and line_[snp_i] in snp_set: #we want to filter and there is the filter
                            var_set.add(line_[snp_i])
                            count_set[line_[snp_i]]+=1 #Update- track the number of times it hits.
                        if not matchlist:
                            var_set.add(line_[snp_i])
                            count_set[line_[snp_i]]+=1 #Update- track the number of times it hits.
            else:
                print("havent implemented yet")
                exit
    file_obj.close()
    #print(len(var_set))
    print("Finished with " + open_file)
with open(args.output, 'w') as ostream: 
    for var in var_set: 
        ostream.write(var + '\n')                      

#The count data, if requested
if args.output_counts is not None:
  print("Count data will be represented as 1-(frequency/[sum(frequencies)]), to put the smallest p-value on the most pleiotropic SNP")
  with open(args.output_counts, 'w') as ostream: 
      for var in var_set: 
          ostream.write(str(var) + "\t" + str((1-(count_set[var]/(nstudies+1)))) + "\n")

#python src/unionVariants.py --gwas_list ./trait_selections/seed1_thresh0.9_h2-0.05.studies.tsv  --output test.txt --type ldsc --pval 1e-5 --extension ".both_sexes.tsv" --gwas_dir /work-zfs/abattle4/lab_data/UKBB/GWAS_Neale/highly_heritable_traits_2/ldsr_format/unzipped/
#python src/unionVariants.py --gwas_list ./trait_selections/seed1_thresh0.9_h2-0.05.studies.tsv  --output test.txt --type ldsc --pval 1e-5 --extension ".both_sexes.tsv" --gwas_dir /work-zfs/abattle4/lab_data/UKBB/GWAS_Neale/highly_heritable_traits_2/ldsr_format/unzipped/
