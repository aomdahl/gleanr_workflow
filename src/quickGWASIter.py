import sys
import argparse
import numpy as np
import glob
import pickle
def readInSNPs(snp_list, ind):
    snps = list()
    with open(snp_list, 'r') as istream:
        for l in istream:
            dat = l.strip().split()
            snps.append(dat[ind])
    return(snps)

#quick helper to make a dictionary of the first column
def colNameKeyBuilder(f):
    with open(f, 'r') as istream:
        first = istream.readline().strip().split()
    return dict(zip(first, range(0, len(first))))

parser = argparse.ArgumentParser(description = "Quickly extract relevant data from Neale lab PAN-GWAS data for certain snps. Default extracts maf, pval, beta, standard error, z and N for full sum stats. For LDSC preformatted ones, only does Z-scores and N")
parser.add_argument("--snp_list", help = "list of snps to extract")
parser.add_argument("--gwas_list", help = "list of gwas studies to extract from")
parser.add_argument("--ancestries", help = "list of ancestries to extract", default="EUR,AFR,EAS")
parser.add_argument("--pickle", default = "", help = "Path to a pickle file.")
parser.add_argument("--output", help = "output path")
parser.add_argument("--gwas_dir", help = "gwas directory to look in.")
parser.add_argument("--type", help = "if its ldsc or not", default = "NOT")
parser.add_argument("--extension", help = "File extension:q to grab")
args = parser.parse_args()

snps =list()
if args.type == "ldsc":
    IND=1
else:
    IND=0
snps = readInSNPs(args.snp_list, IND)
#Quick lookup
lookup_snps = dict(zip(snps, list(range(0, len(snps)))))
if "rs" in list(lookup_snps.keys())[0] and args.type != "ldsc":
    print("We have rsid list but not the LDSC format. Converting now:")
    import subprocess
    check_call("bash src/variant_lookup.sh " + args.snp_list + " " + args.snp_list + ".convert", shell = True)
    #read them in again
    snps = readInSNPs(args.snp_list + ".convert", IND)
    input("Please confirm the conversion went correction- see", args.snp_list + ".convert")
ancestries = args.ancestries.split(",")
ancestries.append("meta")
print(ancestries)
#Get the list of gwas files to examine.
with open(args.gwas_list, 'r') as istream:
        file_list = [x.strip() for x in istream]
print("read in files to read")
nsnps = len(snps)
nstudies = len(file_list)
file_names = [x.split("/")[-1].split(".")[0] for x in file_list]
#TODO: fill the list with NAna, not zeros.

if args.pickle == "" :
    if args.type != "ldsc":
        ret_dat = dict()
        for a in ancestries: #specify these
            ret_dat["af_" + a] = dict(zip(file_names, [np.zeros(nsnps)]*nstudies))
            ret_dat["beta_" + a] = dict(zip(file_names, [np.zeros(nsnps)]*nstudies))
            ret_dat["se_" + a] = dict(zip(file_names, [np.zeros(nsnps)]*nstudies))
            ret_dat["pval_" + a] = dict(zip(file_names, [np.zeros(nsnps)]*nstudies))
        #ret_dat = {"se":dict(zip(file_names, [np.zeros(nsnps)]*nstudies)), "p":dict(zip(file_names, [np.zeros(nsnps)]*nstudies)), "z":dict(zip(file_names, [np.zeros(nsnps)]*nstudies)), "beta":dict(zip(file_names, [np.zeros(nsnps)]*nstudies)), "maf": dict(zip(file_names, [np.zeros(nsnps)]*nstudies))}

    else:
        print("not implemented here...")
        ret_dat = {"z":dict(zip(file_names, [np.zeros(nsnps)]*nstudies)), "n":dict(zip(file_names,[np.zeros(nsnps)]*nstudies))}
else:
    dbfile = open(args.pickle, 'rb')
    ret_dat = pickle.load(dbfile)
    files_so_far = set(ret_dat['pval_meta'].keys())


snp_ids = list()
first = True
print("Reading through files")
for f in file_list:
    print(f)
    fname = f.split("/")[-1].split(".")[0]
    fname = args.gwas_dir + "/" + fname + "*" + args.extension
    if args.pickle != "" and ret_dat["pval_meta"][f][1] != 0:
        #check to see if all 0
        print("we have", f)
        continue
    try:
        open_file = glob.glob(fname)[0]
    except: #cannot find the specified file.
        print('cannot find specified file')
        print(fname)
        print(glob.glob(fname))
        input()
    key = colNameKeyBuilder(open_file) #to get the columns you need
    print(key)
    with open(open_file, 'r') as istream:
        line_num = 0
        print(open_file)

        #fill all entries for that one with nas.
        for df in ret_dat:
            new_list = np.zeros(nsnps)
            new_list = [float('nan') for x in new_list]
            ret_dat[df][f] = new_list

        for line in istream:
            tab = line.strip().split()
            k = key
            if tab[0] == "variant" or tab[0] == "SNP" or tab[0] == "chr":
                continue
            if args.type == "ldsc":
                fid = tab[0]
            else:
                fid = tab[0] + ":" + tab[1]
            if fid in lookup_snps:
                snp_ref_index = lookup_snps[fid]
                for dt in ret_dat:
                    try:
                        #here there are cases where a certain ancestry isn't included in the trait.
                        if dt in k:
                            ret_dat[dt][f][snp_ref_index] = tab[k[dt]]
                        else:
                            ret_dat[dt][f][snp_ref_index] = float('nan')
                    except:
                        print("dt", dt)
                        print("f",f)
                        print("snp index", snp_ref_index)
                        print("data table", tab)
                        print(len(ret_dat[dt][f]))
                        print(nsnps)
                        print(len(ret_dat[dt][f][snp_ref_index]))
                        input()
                line_num += 1
    #Store a local backup so you aren't working everything over and over and over
    #print(ret_dat)
    #input()
    dbfile = open(args.output + 'tmp_results.dump', 'wb')
    pickle.dump(ret_dat, dbfile)
    first = False
import pandas as pd
for dt in ret_dat:
    d = ret_dat[dt]
    d["ids"] = snps
    #put Ids first:
    out = pd.DataFrame.from_dict(d)
    cols = out.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    out = out[cols]
    out.to_csv(args.output + "." + dt + ".tsv" ,index = False, sep = '\t')



