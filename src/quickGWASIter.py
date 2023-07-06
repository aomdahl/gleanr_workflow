import sys
import argparse
import numpy as np
import glob
import pickle
import gzip
import os
def readInSNPs(snp_list, ind):
    snps = list()
    with open(snp_list, 'r') as istream:
        for l in istream:
            dat = l.strip().split()
            snps.append(dat[ind])
    return(snps)

def learnHeaders(tl):
    ret_dict = dict()
    header_map = {"SIGNED_SUMSTAT":"beta","SE":"se","P":"p", "N":"n","Z":"z", "FRQ":"maf"}
    for i,entry in enumerate(tl):
        if(entry in header_map): #just the headers we care about
            ret_dict[header_map[entry]] = i
    return(ret_dict)

parser = argparse.ArgumentParser(description = "Quickly extract relevant data from Neale lab GWAS data for certain snps. Default extracts maf, pval, beta, standard error, z and N for full sum stats. For LDSC preformatted ones, only does Z-scores and N")
parser.add_argument("--snp_list", help = "list of snps to extract")
parser.add_argument("--gwas_list", help = "list of gwas studies to extract from. give the full path to the GWAS file.")
parser.add_argument("--pickle", default = "", help = "Path to a pickle file.")
parser.add_argument("--output", help = "output path")
parser.add_argument("--gwas_dir", help = "gwas directory to look in.")
parser.add_argument("--type", help = "if its ldsc or not", default = "NOT")
parser.add_argument("--keep_dump", help = "Should we keep a backup dump file as we go?", default =False, action = "store_true")
parser.add_argument("--extension", help = "File extension to grab")
args = parser.parse_args()

snps =list()
if args.type == "ldsc" or args.type == "ldsc_custom":
    IND=1
else:
    IND=0
snps = readInSNPs(args.snp_list, IND)
#Quick lookup
lookup_snps = dict(zip(snps, list(range(0, len(snps)))))
if "rs" in list(lookup_snps.keys())[0] and ("ldsc" not in args.type):
    print("We have rsid list but not the LDSC format. Converting now:")
    import subprocess
    check_call("bash src/variant_lookup.sh " + args.snp_list + " " + args.snp_list + ".convert", shell = True)
    #read them in again
    snps = readInSNPs(args.snp_list + ".convert", IND)
    input("Please confirm the conversion went correction- see", args.snp_list + ".convert")

#Get the list of gwas files to examine.
with open(args.gwas_list, 'r') as istream:
    file_list = [x.strip() for x in istream]
print("read in files to read")
nsnps = len(snps)
nstudies = len(file_list)
#file_names = [x.split("/")[-1].split(".")[0] for x in file_list]
file_names = [os.path.basename(x.split()[0]) for x in file_list]
#TODO: fill the list with NAna, not zeros.

if args.pickle == "" : 
    if  "ldsc" not in args.type:
        ret_dat = {"se":dict(zip(file_names, [np.zeros(nsnps)]*nstudies)), "p":dict(zip(file_names, [np.zeros(nsnps)]*nstudies)), "z":dict(zip(file_names, [np.zeros(nsnps)]*nstudies)), "beta":dict(zip(file_names, [np.zeros(nsnps)]*nstudies)), "maf": dict(zip(file_names, [np.zeros(nsnps)]*nstudies)), "n":dict(zip(file_names,[np.zeros(nsnps)]*nstudies))}
    if args.type == "ldsc_custom":    
        ret_dat = {"se":dict(zip(file_names, [np.zeros(nsnps)]*nstudies)), "p":dict(zip(file_names, [np.zeros(nsnps)]*nstudies)), "z":dict(zip(file_names, [np.zeros(nsnps)]*nstudies)), "beta":dict(zip(file_names, [np.zeros(nsnps)]*nstudies)), "n":dict(zip(file_names,[np.zeros(nsnps)]*nstudies))}

    else:
        ret_dat = {"z":dict(zip(file_names, [np.zeros(nsnps)]*nstudies)), "n":dict(zip(file_names,[np.zeros(nsnps)]*nstudies))}
else:
    dbfile = open(args.pickle, 'rb')     
    ret_dat = pickle.load(dbfile)
    files_so_far = set(ret_dat['z'].keys())
#keys for each table type
cont = {"maf":2, "p":10 , "beta":7 , "se":8  , "z":9, "n":4}
disc = {"maf":2, "p":11 , "beta":8 , "se":9  , "z":11, "n":5}
ldsc  = {"z":4, "n":3}
#ldsc custom can vary, depending on if MAF data is available or not
ldsc_custom = {"z":7, "n":3, "se":5, "p":6,"beta":4}
#keep this here, but adjust it leater
key = {11:cont, 12:disc, 5:ldsc, 8:ldsc_custom}
snp_ids = list()
first = True
print("Reading through files")
print(file_list)
for f_line in file_list:
    #its possible its a list, so clean up
    #f = f_line.split()[0]
    f= os.path.basename(f_line.split()[0])
    print(f)
    fname = f.split("/")[-1].split(".")[0]
    fname = args.gwas_dir + "/" + fname + "*" + args.extension
    if args.pickle != "" and ret_dat["z"][f][1] != 0:
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
    #account for zipped files.
    if ".gz" in args.extension:
      file_obj = gzip.open(open_file, 'rt')
    else:
      file_obj = open(open_file, 'r')
      
    with file_obj as istream:
        line_num = 0
        #print(open_file)
        #fill all entries for that one with nas.
        for df in ret_dat:
            new_list = np.zeros(nsnps)
            new_list = [float('nan') for x in new_list]
            ret_dat[df][f] = new_list
        for line in istream:
            tab = line.strip().split()
            if len(tab) == 1: #only have rsid dad
                fid=tab[0]
            else:
              k = key[len(tab)] #if its ldsc, continuous, or binary
            if tab[0] == "variant" or tab[0] == "SNP":
                if "ldsc" in args.type:
                    key[len(tab)] = learnHeaders(tab) #extend this to all types?
                    k = key[len(tab)]
                continue 
            if args.type == "ldsc" or args.type == "ldsc_custom":
                fid = tab[0]
            else:
                id = tab[0].split(":")
                try:
                    fid = id[0] + ":" + id[1]
                except IndexError:
                    print(id)
                    print("Error in either 0 or 1")
                    input("Skipping...")
            if fid in lookup_snps:
               snp_ref_index = lookup_snps[fid]
               for dt in ret_dat:
                    if len(tab) == 1:
                      ret_dat[dt][f][snp_ref_index]=np.NaN
                    else:
                      try:
                          ret_dat[dt][f][snp_ref_index] = tab[k[dt]]
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
    if args.keep_dump:
        dbfile = open(args.output + 'tmp_results.dump', 'wb')
        pickle.dump(ret_dat, dbfile)
        dbfile.close()
    first = False
#print(ret_dat["beta"]["GWAS_abdominal_hernia_SAIGE_550.txt.vcf.abdominal_hernia.sumstats.gz"][lookup_snps["rs10489588"]] )
#print(ret_dat["beta"].keys())
import pandas as pd 
for dt in ret_dat:
    d = ret_dat[dt]
    d["ids"] = snps
    #if dt == "beta":
    #    #matching_index = [i if x is "rs10489588" for x,i in enumerate(d["ids"])]
    #    print(lookup_snps["rs10489588"])
    #    print(d["ids"][lookup_snps["rs10489588"]])
    #    print(d["GWAS_abdominal_hernia_SAIGE_550.txt.vcf.abdominal_hernia.sumstats.gz"][lookup_snps["rs10489588"]])
    #put Ids first:
    out = pd.DataFrame.from_dict(d)
    cols = out.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    out = out[cols]
    out.to_csv(args.output + "." + dt + ".tsv" ,index = False, sep = '\t')


