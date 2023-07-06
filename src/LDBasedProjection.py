import numpy as np
import argparse
import pdb
def adjustVarExplained(rin, thresh): #change to 0.7
    """
    basically, just threshold these at some value
    """
    if rin > thresh:
        return rin
    return 0


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--ld_ref', type=str, 
                    help='path to the ld_reference file onto which we wish to project.')
parser.add_argument('--snp_map', type = str, 
                    help='Conversion between CHR:ID and RSID, as relevant')
parser.add_argument('--r_thresh', type = float, help = "Specify the cutoff for R2 inclusion.", default = 1e-5)
parser.add_argument('--loading', type = str, 
                    help='Loadings from which to projection from.')
parser.add_argument('--output', type = str, 
                    help='Where to save output .csv file.')

args = parser.parse_args()


#ld_ref = "hm3_infertility_p0.001.ld"
ld_ref = args.ld_ref
#ilist="snp.order.factors.txt"
ilist = args.snp_map
#loading_path="K10_top_run_AUG_A0.036854_L32.894377_B_SE.1.loadings.txt"
loading_path = args.loading
target_snp_list=""
r_thresh = args.r_thresh

#map each snp to its index
snp_ref = dict()
snp_rside = dict()

#get the input snps
with open(ilist, 'r') as ostream:
    i=0
    for line in ostream:
        line = line.strip().split()
        if line[0] not in snp_ref:
            snp_ref[line[0]]=i
        else:
            print("Duplicate uh oh")
        i = i+1
print("Input list read in!")
#Map each projectee snp to those its associafted with 
indices = dict()
correlations = dict()
with open(ld_ref, 'r') as ostream:
    for line in ostream:
        if "SNP" in line:
            continue
        line = line.strip().split()
        projectee = line[5]
        r = float(line[6])
        projector = line[2]
        #5 is the sNP
        #6 is the R2
        #being lazy here, just appending as we go
        if projectee not in indices:
            indices[projectee] = list()
            correlations[projectee] = list()
        indices[projectee].append(snp_ref[projector])
        correlations[projectee].append(adjustVarExplained(r, r_thresh))
print(len(indices))
print(len(correlations))
print("Data loaded in...")
#so far so good....
import pandas as pd
loadings = pd.read_table(loading_path,delimiter=",")
projections = dict()
projections["SNP"] = list()
nfactors = loadings.shape[1]
print(loadings.shape)
for colnum in range(1, nfactors): #first index is the snps, so don't bother with that
    projections[colnum] = list()
    # TODO put in a sanity check, that the orders of the SNPs are the same between loadings and ilist...
    for projectee in indices:
        curr_lookups = indices[projectee]
        if projectee in projections:
            print("Error- duplicate id encountered. shouldn't be possible")
        projections[colnum].append(np.dot(loadings.iloc[curr_lookups, colnum].to_numpy(), correlations[projectee]).round(decimals = 6))
        #update the column with snp ides
        if colnum == 1:
            projections["SNP"].append(projectee) #attach the rsid


out_dict = pd.DataFrame.from_dict(projections)
#out_dict.to_csv("surya_projection_full.csv", index = False)
out_dict.to_csv(args.output, index = False)
#print(out_dict)
#Question now about these scores- what is their distribution?
#Our factor effect sizes come from a Laplace distribution, we assume R is constant
#So our projected scores will be linear sums of laplace distributions
#We should plot what they look like, would be helpful
    #We are hoping the distribution 2 is chi-squared, that's what they rely on...
#Need to think about the scaling here... what is the right scaling?
#Is R2 even the right metric? as in should it be R instead?

