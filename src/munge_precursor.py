""" 
for f in files
    read in first few lines
        identify headers, and each assignment
        if is chr:00 and hg37:
            read in, do conversion lookup
        if effect/alt is lowercase, make uppercase
        call to mung 
"""
#NOTE that using effectivee sample size isn't recommended for Ldsc, but whatever.
n_options =["N_SAMPLES", "EFFSAMPLESIZE", "N_COMPLETE_SAMPLES", "NEFF", "TOTALSAMPLESIZE"]
n_options_LDSC=['N','NCASE','CASES_N','N_CASE', 'N_CASES','N_CONTROLS','N_CAS','N_CON','N_CASE','NCONTROL','CONTROLS_N','N_CONTROL', "NSTUDY", "N_STUDY", "NSTUDIES", "N_STUDIES"] #These are ones in LDSC
snp_options = ["MARKER",  "VARIANT_ID", "SNP-ID", "RSIDS"] 
a1_options = ["RISK_ALLELE", "ALT"] #effect allele
a2_options = ["OTHER_ALLELE", "REF"] #non-effect allele
sum_stats = ["MAINEFFECTS", "ESTIMATE", "ODDS_RATIO"]
maf_options = ["MINOR_AF", "REF_ALLELE_FREQUENCY", "FREQ1"]
se_options = ["MAINSE", "SDE"]
pval_options = ["MAINP", "P_GC"]
PVAL_FULL=["MAINP", "P_GC", "P", "PVALUE", "P_VALUE", "PVAL", "P_VAL", "GC_PVALUE"]
IDM={"chr_opts" : ["CHR", "CHROMOSOME"], "pos_opts" :["POS", "base_pair_location".upper(), "BP", "POSITION"], "snpid_opts":["VARIANT_ID", "VARIANT", "SNP", "NAME", "MARKERNAME"]}
#chromosome	base_pair_location
labels = {"n" : n_options, "snp" : snp_options, "maf" : maf_options, "a1" : a1_options, "a2" :a2_options, "ss" : sum_stats, "se" :se_options, "p" : pval_options, "n_existing" : n_options_LDSC}
import sys
import argparse
import gzip
import os
from subprocess import check_call
import subprocess
import shlex
import csv
import numpy as np
import re
import pdb

which = lambda lst:list(np.where(lst)[0])
def any_in(x,y):
     return np.isin(np.array(x), np.array(y)).any()

def isVariantID(text):
    """
    Helper function to detect if this has addresses, not rsids, which we want.
    """
    split = text.split(":")
    if (split[0] in [str(x) for x in list(range(0,24))]) or ("CHR" in split[0].upper()): 
        return True
    else:
        return False


def mergeAlleles(argin):
    if argin == "":
        return ""
    else:
        return " --merge-alleles " + argin

def openFileContext(path):
    if path[-3:] == ".gz" or path[-4:] == ".bgz":
        return gzip.open(path, 'rb')
    else:
        return open(path, 'r')

def outFileContext(path):
    if path[-3:] == ".gz" or path[-4:] == ".bgz":
        return gzip.open(path, 'wt')
    else:
        return open(path, 'w')
    
def detectDelimiter(line):
    #from https://stackoverflow.com/questions/3952132/how-do-you-dynamically-identify-unknown-delimiters-in-a-data-file
    sniffer = csv.Sniffer()
    return sniffer.sniff(line).delimiter

#File handler- if gz or not.
def fh(st, p):
    try:
        if ".bgz" in  p or ".gz" in p:
            return st.decode(errors='replace')
        else:
            return st
    except TypeError:
        print(p)
        print(str(p))
        input()
        return p

def processN(intext):
    if "/" in intext: #this means a case/control one was given:
        fl = intext.split("/")
        return (" --N-cas " + str(fl[0]) + " --N-con " + str(fl[1]))  
    else:
        return " --N " + str(intext)


def labelSpecify(header, checkfor, labels):
    """
    Specifies which column maps to which thing in cases when munge won't pick it up
    @header- the LIST of header items
    @checkfor: which item we are looking for
    @labels: the map of labels to lookup options.
    """
    str_map = {
        "snp" : " --snp ",
        "n" : " --N-col ",
        "a1": " --a1 ",
        "a2": " --a2 ",
        "ss": " --signed-sumstats ",
        "p" : " --p ",
        "maf":" --frq ",
    }
    #Check- is the header one of the aberrant ones?
    relevant_col = which([x.upper() in labels[checkfor] for x in header])
    if(len(relevant_col) > 0):
        header_label = str(header[relevant_col[0]])
        ret_snp = str_map[checkfor] + header_label
        #Special case for providing sum stats, need to know what the null value is
        if checkfor == "ss":
            size_options = {"OR":"1", "ODDS_RATIO": "1","BETA":"0", "MAINEFFECTS":"0", "EFFECTS":"0", "B":"0", "ZSCORE":"0", "Z":"0", "ESTIMATE":"0"}
            ret_snp = ret_snp + "," + size_options[header_label.upper()]
        ret_snp = ret_snp + " "
    else:
        ret_snp = ""
        if checkfor == "n":
            #Check if its in the existing LDSC list.
            #If not, return the signal to specify it
            relevant_col = which([x.upper() in labels["n_existing"] for x in header])
            if len(relevant_col) == 0:
                ret_snp =  "USE_PROVIDED"
            else:
                ret_snp = ""
        else:
            ret_snp = ""
    return ret_snp

def filePeek(readin, argv):
    """
    This looks at the summary stats file and determines if munge_stats will be compatible with it.
    If not, it makes necessary changes or notifies the user. Yup.
    @return cleanup_protocol: instructions for cleaning up the file, if any.
    @return ret_n the string to append to the munge command to deal with sample sizes
    TODO: if chr and pos are included, maybe favor those. Some of the other stuff is janky wanky.
    """
    fpath = readin[0]
    pheno = readin[1]
    cleanup_protocol = ""
    ret_n = ""
    try:
        with openFileContext(fpath) as istream:
            header = fh(istream.readline(), fpath)
            header = header.strip()
            first = fh(istream.readline(),fpath)
            first_tab = first.strip().split()
            delim  = detectDelimiter(header)

            if delim == ",":
                cleanup_protocol = "CSV"
            
            if ("PaxHeader" in header) or ("GIANT" in fpath):
                #This is a GIANT file, need to do cleanup
                cleanup_protocol = cleanup_protocol + "GIANT"
                return cleanup_protocol, ret_n

            header_dat = (header.upper()).split(delim)
            #if "RSID" in header.upper() and "RS" in first.upper():
            #    cleanup_protocol = cleanup_protocol + "NONE"
            if "hm_variant_id".upper() in header_dat[0]:
                cleanup_protocol = cleanup_protocol + "HEADER_HM"

            if header_dat[0] == "variant".upper() and header_dat[-1] == "pval".upper(): #we suspect this is UKBB full format file.
                if isVariantID(first_tab[0]):
                    #print("The current file at", fpath, "appears to be a full Neale Lab UKBB file. We don't recommend using this at this stage, use the LDSC pre-formatted one available onlinel")
                    cleanup_protocol = cleanup_protocol + "UKBB"
                    #return cleanup_protocol, ret_n
                    #need to deal with UKBB cases of both continuous and case-control
            #the split column case
            if any_in(header_dat, IDM["chr_opts"]) and any_in(header_dat, IDM["pos_opts"]) and "RSID" not in header.upper(): #changed this.
                if "RS" not in first.upper(): #sometimes its under another header label we can recognize
                    cleanup_protocol = cleanup_protocol + "TO_RSID_2"
            else:
                if ":" in header_dat[0] or ("MARKERNAME" in header_dat[0].upper()):
                    if "MARKERNAME" in header_dat[0].upper():
                        #weird format with rs3131972:752721:A:G
                        cleanup_protocol = cleanup_protocol + "FIRST_RSID"
                    else:
                        cleanup_protocol = cleanup_protocol + "TO_RSID_1"
                if ":" in first and "RS" not in first.upper(): #No rsid but likely something else
                    cleanup_protocol = cleanup_protocol + "TO_RSID_1"
            #Get the sample size if its there...
            ret_n = labelSpecify(header_dat, "n", labels)
            if ret_n == "USE_PROVIDED": 
                ret_n =  processN(readin[3])
            #Check for the marker data
            ret_a1 = labelSpecify(header_dat, "a1", labels)
            ret_a2 = labelSpecify(header_dat, "a2", labels)
            ret_ss = labelSpecify(header_dat, "ss", labels)
            ret_snp = labelSpecify(header_dat, "snp", labels)
            ret_p = labelSpecify(header_dat, "p", labels)
            if(argv.keep_maf):
                ret_maf = labelSpecify(header_dat, "maf", labels)
            else:
                ret_maf = ""
    except FileNotFoundError:
        print("Specified path does not exist")
        
        cleanup_protocol = "ERROR"
        return cleanup_protocol, [fpath]
    
    if cleanup_protocol == "": #nothing to change.
        cleanup_protocol = "NONE"

    return cleanup_protocol, [ret_a1, ret_a2, ret_ss, ret_snp, ret_n, ret_p, ret_maf] 

def makeInterFileName(fpath, protocol):
    """
    If no change to be made, don't modify the file.
    """
    remakes = ["TO_RSID_2", "GIANT", "CSV", "UKBB","HEADER_HM","FIRST_RSID", "RSID_1" ]
    #os.path.basename()
    if any([x in protocol for x in remakes]):
        r = os.path.splitext(fpath)[0] + ".INTERMEDIATE" + os.path.splitext(fpath)[1]
        if ".csv" in fpath:
            return r.replace("csv", "tsv")
        else:
            return r
    else:
        return fpath

def buildOutName(outpath, inpath, inpheno):
    """
    Specifies the name of the output file based on the output directory, the phenotype name, etc.
    @param outpath: the provided output directory
    @param inpath: the input file info
    """
    fn = os.path.basename(inpath)
    with_ext = os.path.splitext(fn)

    return outpath + with_ext[0] + "." + inpheno


def validVersionExists(fn):
    #file exists
    import os.path
    if not os.path.isfile(fn):
        return False
    
    #file has some reasonable number of lines (i.e. over 500,000)
    try:
        num_lines = sum(1 for line in openFileContext(fn))
        if num_lines < 850000:
            return False
    except gzip.BadGzipFile:
        return False
    #if its a gzip file, make sure its valid
    #code snapped from https://stackoverflow.com/questions/41998226/python-checking-integrity-of-gzip-archive
    if fn[-3:] == ".gz" or fn[-4:] == ".bgz":    
        import os
        r=os.popen('gunzip -t ' + fn).read()
        if r != "":
            return False
    return True
    #also do some kind of RSID-based sanity check? Or look at the last line in the file?

def doCleanup(readin, protocol, rsids, correct_files, force_update = False):
    """
    Performs the cleanup step on files that need it
    UKBB, we should recommend ldsc
    Giant: remove teh first few lines, and get the header, keep the rest of the file the same
    TO_RSID: convert the chr/things to rsids; make sure to check using the right genome build (!)
    """
    print("Correction protocol is", protocol)
    fpath = readin[0]
    intermediate_file_name = makeInterFileName(fpath, protocol)
    print("curr file", fpath)
    print("new file", intermediate_file_name)
    if correct_files:
        print("Not correcting files, just writing out.")
        return intermediate_file_name
    if validVersionExists(intermediate_file_name) and not force_update and protocol != "NONE": #if no changes to be made, don't bother checking..
        print("Valid version of file already exists. Will not force update")
        return intermediate_file_name
    else:
        print("Proceeding with update...")
    delim = ""
    T='\t'
    if protocol == "ERROR":
        print("Print unable to process current file")
        return ""
    if protocol != "NONE":
        print("Changes being made to", fpath)
        print("Updated file name is ", intermediate_file_name)
        print("File change protocol code is ", protocol)
    if "GIANT" in protocol :
        #remove the first few lines of the file
        #done_header="\"MarkerName\\tChr\\tPos\\tAllele1\\tAllele2\\tFreqAllele1HapMapCEU\\tb\\tse\\tp\\tN\""
        done_header="MarkerName\tChr\tPos\tAllele1\tAllele2\tFreqAllele1HapMapCEU\tb\tse\tp\tN\n"
        print_switch = False
        with openFileContext(fpath) as istream:
            with outFileContext(intermediate_file_name) as ostream:
                for i, line in enumerate(istream):
                    line = fh(line.strip(), fpath)
                    
                    if print_switch:
                        ostream.write(line + '\n')
                    if "MarkerName" in line:
                        ostream.write(done_header)
                        print_switch = True
        fpath = intermediate_file_name #changes updated
    if "UKBB" in protocol:
        adjacent=False
        with openFileContext(fpath) as istream:
            with outFileContext(intermediate_file_name) as ostream:
                for i, line in enumerate(istream):
                    line = fh(line.strip(), fpath)
                    if i == 0:
                        header = line
                        if("effect_allele" in header):
                            print("This one is UKBB-like, but not the same")
                            adjacent = True 
                            #no change to header
                        else:
                            #header = header.replace("variant\t", "SNP\teffect_allele\tother_allele")
                            header = re.sub('^variant\s', "SNP\teffect_allele\tother_allele\t", header)
                        ostream.write(header + '\n')
                    else:
                        dat = line.split('\t')
                        var_dat = dat[0].split(":")
                        varid = var_dat[0] + ":" + var_dat[1]
                        if varid in rsids:                     
                            snp = rsids[varid]
                            effect = var_dat[3]
                            other = var_dat[2]
                            if adjacent:
                                ostream.write(snp + T + T.join(dat[1:]) + '\n')
                            else:
                                ostream.write(snp + T + effect + T + other + T + T.join(dat[1:]) + '\n')
                        #otherwise we skip it
        fpath = intermediate_file_name #changes updated

           
    if "TO_RSID" in protocol and "UKBB" not in protocol:
        #If its TO_RSID, its a 2 column thing. if its TO_RSID_1 its a 1 column thing.
        dropped_snps = 0
        T='\t'
        print("Converting SNP IDs to RSIDs, constraining to only those in RSID list (default is Hapmap3 non-HLA SNPs)")
        with openFileContext(fpath) as istream:
            with outFileContext(intermediate_file_name) as ostream:
                for i, line in enumerate(istream):
                    line = fh(line.strip(), fpath)                     
                    if i == 0: 
                        delim = detectDelimiter(line)
                        #ostream.write("SNP" + T + T.join(line.split(delim)[1:]) + '\n')
                        ostream.write("SNP" + T + T.join(line.split(delim)) + '\n')
                        if ("TO_RSID_2" in protocol) and ("TO_RSID_1" not in protocol):
                            which_chr = which([x.upper() in IDM["chr_opts"] for x in line.split(delim)])[0]
                            which_snp = which([x.upper() in IDM["pos_opts"] for x in line.split(delim)])[0]
                        else: #protocol == "TO+RSID1"
                            which_id = which([x.upper() in IDM["snpid_opts"] for x in line.split(delim)])[0]
                    else:
                        line = line.split(delim)
                        if ("TO_RSID_2" in protocol) and ("TO_RSID_1" not in protocol):
                            address = (line[which_chr] + ":" + line[which_snp]).replace("chr", "")
                        else:
                            address = line[which_id].replace("chr", "")
                            if ":" in line[which_id]:
                                address = ":".join(line[which_id].split(":")[0:2]).replace("chr", "")
                        if address in rsids:
                            snp = rsids[address]
                            ostream.write(snp + T + T.join(line) +  '\n')
                            #ostream.write(snp + T + T.join(line[1:]) +  '\n')
                        else:
                            dropped_snps += 1
                            continue
        print(dropped_snps, "SNPs omitted.")
        fpath = intermediate_file_name #changes updated
        if(protocol == 'UKBBTO_RSID_1'):
            return fpath
    #For all remaining protocols, we don't do anything.                    
        #determine which columnes have the info we need, look it up, and convert it.
    
    if "CSV" in protocol:
        if "TO_RSID" not in protocol: #if it is, we've already made the adjustments to delimiter
            print("Changing CSV to TSV...")
            with openFileContext(fpath) as istream:
                with outFileContext(intermediate_file_name) as ostream:
                    for i, line in enumerate(istream):
                        line = fh(line.strip(), fpath)
                        ostream.write(line.replace(",", T) + '\n') 
            fpath = intermediate_file_name #changes updated                   
        
    if "HEADER_HM" in protocol: #just need to change the first line
        with openFileContext(fpath) as istream:
            with outFileContext(intermediate_file_name) as ostream:
                for i, line in enumerate(istream):
                    line = fh(line.strip(), fpath)
                    if i == 0:
                        header = line.replace("hm_", "")
                        ostream.write(header + '\n')
                    else:
                        ostream.write(line + '\n')
        fpath = intermediate_file_name #changes updated

    #This is an unusual type of header where they have RSID:CHR:POS or something like that.
    if "FIRST_RSID" in protocol and "UKBB" not in protocol and ("TO_RSID_1" not in protocol): #this last option should have been done above
        with openFileContext(fpath) as istream:
            with outFileContext(intermediate_file_name) as ostream:
                for i, line in enumerate(istream):
                    line = fh(line.strip(), fpath)
                    if i == 0:
                        header = line.upper().replace("MARKERNAME", "SNP")
                        ostream.write(header + '\n')
                    else:
                        dat = line.split()
                        if "RS" in line.upper(): #the wierd one
                            rsid = dat[0].split(":")[0]
                        else: #another common setup.
                            #May be a 4 parter:
                            snpdat = ":".join(dat[0].split(":")[0:2])
                            if(snpdat in rsids):
                                rsid = rsids[snpdat]
                            else:
                                #print("Missed", snpdat)
                                continue
                        ostream.write(rsid + '\t' + '\t'.join(dat[1:]) + '\n')
        fpath = intermediate_file_name #changes updated
    return intermediate_file_name


def importRSDict(p, invert = False):
    """
    By default uses a list of no hLA region hm3 snps
    """
    ret_dict = dict()
    with open(p, 'r') as istream:
        for line in istream:
            line = line.strip().split()
            #ret_dict[line[0] + ":" line[1]] = line[2]
            if not invert:
                ret_dict[line[0]] = line[1]
            else:
                ret_dict[line[1]] = line[0]
    return ret_dict

#Helper function to upldate the rference list if a sublist is submitted.
def buildReferenceList(argsin):
    if "," not in args.rsid_ref:
        build_list = {"hg37":importRSDict(args.rsid_ref)}
        
        if argsin.rsid_ref != "/data/abattle4/aomdahl1/reference_data/hm3_nohla.snpids":
            print("non-default list is entered; going to do a quick liftover to get hg38 too...")
            temp_new = importRSDict("/data/abattle4/aomdahl1/reference_data/hm3_nohla.hg38.snpids", invert = True)
            build_list["hg38"] = dict()
            for i in build_list["hg37"].values():
                if i in temp_new:
                    build_list["hg38"][temp_new[i]] = i 
            print("Converted list has", len(build_list["hg38"]), "entries")
            print("Original list has", len(build_list["hg37"]), "entries")
            #breakpoint()
        else:
            build_list["hg38"] = importRSDict("/data/abattle4/aomdahl1/reference_data/hm3_nohla.hg38.snpids")
            build_list["hg36"] = None
    else:
        print("Multiple input references recieved. Assuming build is specified on each (i.e. HG37:PATH,HG38:PATH")
        paths = args.rsid_ref.split(",")
        build_list = dict()
        for p in paths:
            g = p.split(":")
            build_list[g[0].lower()] = importRSDict(g[1])
    #breakpoint()
    return build_list

parser = argparse.ArgumentParser(description = "Quick script to get many GWAS summary stats converted into an LDSC-friendly format. A nice wrapper for ldsc's munge script.")
parser.add_argument("--study_list", help = "list of studies to work on. Columns are study path, phenotype name, genome build, N samples. Last 2 columns not required, in which case we assume hg37 and that N is in the data. Note that if N is in the data, that is prioritized above a specified N.")
parser.add_argument("--rsid_ref", help = "path to the RSID reference", default = "/data/abattle4/aomdahl1/reference_data/hm3_nohla.snpids")
parser.add_argument("--ldsc_path", help = "path to LDSC code", default="~/.bin/ldsc/")
parser.add_argument("--merge_alleles", help = "Path to set up a predefined allele reference; if is empty string, then none used", default = "/data/abattle4/aomdahl1/reference_data/hm3_snps_ldsc_ukbb.tsv")
parser.add_argument("--keep_maf", help = "Use this to get the MAF from files.", default = False, action = "store_true")
parser.add_argument("--force_update", help = "Specify this if you want to force all intermediate files to be re-created. If intermediate files exist already and are a reasonable length, they aren't recreated.", default = False, action = "store_true")
parser.add_argument("--mod", help = "Specify this if you want to run the modified version of LDSC that lets you extract B, SE, etc.", default = False, action = "store_true")
parser.add_argument("--calls_only", help = "Use this if you want only to generate the commands, not run the actual correction.", default = False, action = "store_true")

parser.add_argument("--output", help = "output path")
args = parser.parse_args()

#get the hg37 build dict in memory....
build_list = buildReferenceList(args)

snp_list_pval_opt = dict()
with open(args.study_list ,'r') as istream:
    with open(args.output + "munge_sumstats.all.sh", 'w') as ostream:
        for line in istream:
            CONVERT_GENOME_BUILD=False
            dat = line.strip().split()
            print("Currently processing file ", dat[0], "(", dat[1], ")")
            """
            if args.pval_filter == "":
                #peek and determine which column actually has the pvalue
                #filter by the corresponding threshold, return just the list of snps you desire as RSIDs for passage into the next phase.
                add_snps = pvalPeekAndClean(dat) #get the RSIDs of all snps that pass threshold thank you very much.
                snp_list_pval_opt = updatePvalSNPs(add_snps,snp_list_pval_opt)
                #Or run fastAsset instead...
                #write this out to output.

                sys.exit()
            """
            #Check the first lines of the file, see if its UKBB template, if it has N, etc.
            #Check to see if its got RSIDs or now...
            fix_instructions, fixed_labels = filePeek(dat, args)
            #detect the genome build if needed....
            if (len(dat) > 2) and (dat[2] != "" )and ("TO_RSID" in fix_instructions):
                print("Input specifies", dat[2], "build.")
                if(dat[2] == "hg36"):
                    print("not yet implemented")
                    continue
                else:
                    build_dict = build_list[dat[2]]
                    if dat[2] not in ["hg37", "hg19"]:
                            CONVERT_GENOME_BUILD = True
            else:
                print("Assuming hg37 build, using hapmap3 SNPS only...")
                build_dict = build_list["hg37"]
            intermediate_file = doCleanup(dat, fix_instructions,build_dict, args.calls_only, force_update = args.force_update)
            #build commands for it
            if fix_instructions and fix_instructions != "ERROR":
                call = " munge_sumstats.py "
                if args.mod:
                    call = "munge_sumstats_MOD.py "
                lc = "python2 " + args.ldsc_path + call + "--sumstats " + intermediate_file + " --out " + buildOutName(args.output, dat[0], dat[1]) + "".join(fixed_labels) + mergeAlleles(args.merge_alleles)
                if args.keep_maf:
                    lc = lc + " --keep-maf "
                ostream.write(lc + '\n')
            else:
                print("Skipping for incorrect file path,", fixed_labels[0])
                with open(args.output +"_errorlog.txt",'w') as ostream:
                    ostream.write("File could not be found: " + fixed_labels[0] + '\n')
            print()
        

