
#NOTE that using effectivee sample size isn't recommended for Ldsc, but whatever.
case_maf_options = ["AF_CASES_EUR","AF_CASES"]
control_maf_options = ["AF_CONTROLS_EUR","AF_CONTROLS"]
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

def makeNum(str):
    if str.upper() == "NA" or str == "NaN":
        return float("nan")
    else:
        return float(str)

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


def getAFIndices(dat_tab,case_names, control_names):
    case_i = -1
    control_i=-1
    for i, term in enumerate(dat_tab):
        if term.upper() in case_names:
            case_i=i
        if term.upper() in control_names:
            control_i = i
        if case_i != -1 and control_i != -1:
            print("Found the c/c af columns")
            break
    return case_i, control_i

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
        "maf":" --frq "
    }
    #Check- is the header one of the aberrant ones?
    relevant_col = which([x.upper() in labels[checkfor] for x in header])
    if(len(relevant_col) > 0):
        header_label = str(header[relevant_col[0]])
        ret_snp = str_map[checkfor] + header_label
        #Special case for providing sum stats, need to know what the null value is
        if checkfor == "ss":
            size_options = {"OR":"1", "ODDS_RATIO": "1","BETA":"0", "MAINEFFECTS":"0", "beta_EUR".upper():"0","EFFECTS":"0", "B":"0", "ZSCORE":"0", "Z":"0", "ESTIMATE":"0"}
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

def missingAnyAlleleInformation(header, IDM):
    """
    Ensures the header has some kind of allele information. Can customize by adjusting the IDM lists....
    """
    if any_in([h.upper() for h in header], IDM["a1_opts"]):
        if any_in([h.upper() for h in header], IDM["a2_opts"]): #both conditions satisfied, NP
            return False
    return True

def whichMissingAlleleInformation(header, IDM):
    """
    Ensures the header has some kind of allele information. Can customize by adjusting the IDM lists....
    """
    if any_in([h.upper() for h in header], IDM["a1_opts"]):
        aoi = which([h.upper() in IDM["a1_opts"]  for h in header])[0]
        if any_in(header, IDM["a2_opts"]): #both conditions satisfied, NP
            return ""
        else:
            
            return "ALLELE_REPAIR_2:" + str(aoi)  #you have 1, so fix 2. Index of 1 is aoi
    else: # you don't have 1
        if any_in([h.upper() for h in header], IDM["a2_opts"]): #both conditions satisfied, NP
            ati = which([h.upper() in IDM["a2_opts"]  for h in header])[0]
            return "ALLELE_REPAIR_1:" + str(ati) #have 2, need 1
        else:
            return "ALLELE_REPAIR_B"

def readInRefDict(path):
    ret_d = dict()
    if path == "":
        return None
    with open(path, 'r') as istream:
        for l in istream:
            d = l.strip().split()
            ret_d[d[0]] = [d[1], d[2]]
    return ret_d


# filePeek(dat, args)
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
            upper_header = header.upper()
            if delim == ",":
                cleanup_protocol = "CSV"
            
            if ("PaxHeader" in header) or ("GIANT" in fpath):
                #This is a GIANT file, need to do cleanup
                cleanup_protocol = cleanup_protocol + "GIANT"
                return cleanup_protocol, ret_n

            header_dat = (upper_header).split(delim)
            #if "RSID" in header.upper() and "RS" in first.upper():
            #    cleanup_protocol = cleanup_protocol + "NONE"
            if "hm_variant_id".upper() in header_dat[0]:
                cleanup_protocol = cleanup_protocol + "HEADER_HM"

            if "neglog10_pval".upper() in upper_header:
                cleanup_protocol = cleanup_protocol + "CONVLOGPTOSTD"
            if case_control_afs[1] in upper_header and case_control_afs[2] in upper_header:
                cleanup_protocol = cleanup_protocol + "CC_JOINALLELEFREQ"
            if header_dat[0] == "variant".upper() and header_dat[-1] == "pval".upper(): #we suspect this is UKBB full format file.
                if isVariantID(first_tab[0]):
                    #print("The current file at", fpath, "appears to be a full Neale Lab UKBB file. We don't recommend using this at this stage, use the LDSC pre-formatted one available onlinel")
                    cleanup_protocol = cleanup_protocol + "UKBB"
                    #return cleanup_protocol, ret_n
                    #need to deal with UKBB cases of both continuous and case-control
            #the split column case
            if any_in(header_dat, IDM["chr_opts"]) and any_in(header_dat, IDM["pos_opts"]) and "RSID" not in upper_header: #changed this.
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
            if missingAnyAlleleInformation(header_dat, IDM):
                    print("Missing allele information....")
                    cleanup_protocol = cleanup_protocol + whichMissingAlleleInformation(header_dat, IDM)

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
    remakes = ["TO_RSID_2", "GIANT", "CSV", "UKBB","HEADER_HM","FIRST_RSID", "RSID_1", "ALLELE_REPAIR", "CONVLOGPTOSTD", "CC_JOINALLELEFREQ"]
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
    #Alternative option
    size=os.path.getsize(fn)
    if size < 5000000:
        return False

    """
    try:
        num_lines = sum(1 for line in openFileContext(fn))
        if num_lines < 850000:
            return False
    except gzip.BadGzipFile:
        return False
    """
    #if its a gzip file, make sure its valid
    #code snapped from https://stackoverflow.com/questions/41998226/python-checking-integrity-of-gzip-archive
    if fn[-3:] == ".gz" or fn[-4:] == ".bgz":    
        import os
        r=os.popen('gunzip -t ' + fn).read()
        if r != "":
            return False
    return True
    #also do some kind of RSID-based sanity check? Or look at the last line in the file?

def alleleRepair(fix_allele, lookup_dict, rs_index, allele_index, dat):
        """
        Give the missing allele information based on the lookup_dict reference
        fix_allele: which allele is missing- either 1, 2, or both
        lookup_dict: the reference to which we are aligning (snp:[a1,a2])
        rs_index: which index has the rsid we use as a reference
        allele_index: which index in the table has the current allele information. if its -1, it means neither.
        dat: the current line of the file.
        """
        
        fix_code  = {"1":"A1", "2": "A2", "B":"A1\tA2"}
        if dat == "HEADER":
            return fix_code[fix_allele]
        else:
            if dat[rs_index] in lookup_dict:
                alleles = lookup_dict[dat[rs_index]] #a tuple, [A1, A2]
                if fix_allele == "B":
                    return "\t".join(alleles)
                else:
                    if dat[allele_index] == alleles[0]:
                        return alleles[1]
                    elif dat[allele_index] == alleles[1]:
                        return alleles[0]
                    else: #They don't match
                        return "NA"
            else:
                #print("No allele data for SNP:", dat[rs_index])
                if fix_allele == "b":
                    return "NA\tNA"
                return "NA"


def getRepairAllele(protocol):
    f = re.search("ALLELE_REPAIR_([12B]):*(\d*)", protocol)
    if f.group(1) == "B":
        return "B", None
    else:
        return f.group(1), int(f.group(2))

def doCleanup(readin, protocol, rsids, correct_files,allele_reference, force_update = False):
    """
    Performs the cleanup step on files that need it
    UKBB, we should recommend ldsc
    Giant: remove teh first few lines, and get the header, keep the rest of the file the same
    TO_RSID: convert the chr/things to rsids; make sure to check using the right genome build (!)
    ALLELE_REPAIR: Add in both the effect and other allele, based on the provided reference.... 
        Note that this requires that the File already have rsids...
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
                            header = re.sub('^variant\s', "SNP_\teffect_allele\tother_allele\t", header)
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
                        ostream.write("SNP_" + T + T.join(line.split(delim)[1:]) + '\n')
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
                            ostream.write(snp + T + T.join(line[1:]) +  '\n')
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

    if "CONVLOGPTOSTD" in protocol  and "UKBB" not in protocol and ("TO_RSID_1" not in protocol) and "CSV" not in protocol: 
        with openFileContext(fpath) as istream:
            print(intermediate_file_name)
            with outFileContext(intermediate_file_name) as ostream:
                for i, line in enumerate(istream):
                    line = fh(line.strip(), fpath)
                    if i == 0:
                        header = line.replace("hm_", "")
                        pval_index= ["PVAL" in x.upper() for x in header.split()].index(True)
                        header = header.split()
                        header[pval_index] = "PVAL"
                        ostream.write(T.join(header) +  '\n')
                        #ostream.write(header + '\n')
                    else:
                        line=line.split()
                        line[pval_index] = str(10**(-1 * makeNum(line[pval_index])))
                        ostream.write(T.join(line) +  '\n')

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
                        header = line.upper().replace("MARKERNAME", "SNP_")
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
    
    if "ALLELE_REPAIR" in protocol and ("UKBB" not in protocol):
        print("Repairing allele data based on provided reference...")
        #At this point, we assume the output file has the rsid information in it.
        repair_a, allele_i = getRepairAllele(protocol) #a quick regtex, formate ALLELE_REPAIR:2:I
        print(repair_a, allele_i)
        #set the file name
        infile = intermediate_file_name if os.path.isfile(intermediate_file_name) else fpath
        ofile = os.path.splitext(intermediate_file_name)[0] + ".REPAIR" + os.path.splitext(intermediate_file_name)[1]
        print(ofile)
        input()
        rs_i = None
        with openFileContext(infile) as istream:
            with outFileContext(ofile) as ostream:
                for i, line in enumerate(istream):
                    lined = line.strip().split()
                    if i == 0:
                        ostream.write("\t".join(lined) +'\t' + alleleRepair(repair_a,allele_reference,0, allele_i,"HEADER" ) + '\n')
                        if("SNP" in line):
                            rs_i = which(["SNP_" in d for d in lined])[0]
                    else:
                        if not rs_i:
                            rs_i = which(["rs" in d for d in line])[0] #
                        ostream.write("\t".join(lined) +'\t'+ alleleRepair(repair_a,allele_reference, rs_i,allele_i, line ) + '\n')
        
        #final step- move the modified file into the correct file...
        #pdb.set_trace()
        d = os.popen("mv " + ofile + " " + intermediate_file_name)
        fpath =  intermediate_file_name
    
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
parser.add_argument("--output", help = "output path")
args = parser.parse_args()

#Pick out all the files to edit
files_to_update = list()
ncase=None
ncontrol = None

with open(args.study_list ,'r') as istream:
    for line in istream:
        CONVERT_GENOME_BUILD=False
        dat = line.strip().split()
        if "/" in dat[3]:
            print("Will process file: ", dat[0], "(", dat[1], ")", "with case/control of sample size ", dat[3])
            files_to_update.append(dat)
            ncase = float(dat[3].split("/")[0])
            ncontrol = float(dat[3].split("/")[1])
            print(ncase, ncontrol)
print("Processing a total of ", len(files_to_update), " files")
#Edit them if possible
for file_dat in files_to_update:
    tmp_file_name = args.output + file_dat[0] + ".tmp.gz"
    fpath = file_dat[0]
    existing_intermediate_file = os.path.splitext(fpath)[0] + ".INTERMEDIATE" + os.path.splitext(fpath)[1]
    if os.path.exists(existing_intermediate_file):
        fpath = existing_intermediate_file
    with openFileContext(fpath) as istream:
        with outFileContext(tmp_file_name) as ostream:
            for i, line in enumerate(istream):
                line = fh(line.strip(), fpath)
                if i ==0:
                    #Get the indices, if its there
                    case_i, control_i = getAFIndices(line.split(), case_maf_options, control_maf_options)
                    if case_i != -1 and control_i != -1:
                        print("Found both, proceeding with edit")
                        ostream.write(line + '\t' + "FRQ" + "\n")
                    else:
                        print("AF data is unavailable, not going to edit file")
                        break
                else:
                    dat=line.split()
                    if dat[case_i] == "NA" or dat[control_i] == "NA":
                        frq="NA"
                    else:
                        frq= round(((ncase*float(dat[case_i])) + (ncontrol*float(dat[control_i])))/(ncase + ncontrol), 5)
                    #print("frq", str(frq))
                    #print("Ncase", str(ncase))
                    #print("Ncontrol", str(ncontrol))
                    #print("AF_case", str(dat[case_i]))
                    #print("AF_control", str(dat[control_i]))
                    #print(line + '\t' + str(frq) + "\n")
                    #input()
                    ostream.write(line + '\t' + str(frq) + "\n")
    print("Finished processing", fpath)
    #mov intermediate_file_name into fpath
    os.rename(tmp_file_name, fpath)
    print("Moved into current INTERMEDIATE file")



