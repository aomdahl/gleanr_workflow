import numpy as np
def adjustVarExplained(rin, thresh = 1e-5):
    """
    basically, just threshold these at some value
    """
    if rin > thresh:
        return rin
    return 0

ld_ref = "hm3_infertility_p0.001.ld"
#ilist = "toy.input.snp.list.txt"
ilist="snp.order.factors.txt"
#loading_path="toy.loadings.txt"
loading_path="K10_top_run_AUG_A0.036854_L32.894377_B_SE.1.loadings.txt"
target_snp_list=""
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
        correlations[projectee].append(adjustVarExplained(r))

#so far so good....
import pandas as pd
loadings = pd.read_table(loading_path,delimiter=",")
projections = dict()
projections["SNP"] = list()
nfactors = loadings.shape[1]
for colnum in range(1, nfactors): #first index is the snps, so don't bother with that
    print("on factor ", colnum)
    projections[colnum] = list()
    # TODO put in a sanity check, that the orders of the SNPs are the same between loadings and ilist...
    for projectee in indices:
        curr_lookups = indices[projectee]
        if projectee in projections:
            print("Error- duplicacte id encountered. shouldn't be possible")
        projections[colnum].append(np.dot(loadings.iloc[curr_lookups, colnum].to_numpy(), correlations[projectee]).round(decimals = 6))
        #update the column with snp ides
        if colnum == 1:
            projections["SNP"].append(projectee) #attach the rsid

out_dict = pd.DataFrame.from_dict(projections)
out_dict.to_csv("surya_projection_full.csv", index = False)
print(out_dict)
#Question now about these scores- what is their distribution?
#Our factor effect sizes come from a Laplace distribution, we assume R is constant
#So our projected scores will be linear sums of laplace distributions
#We should plot what they look like, would be helpful
    #We are hoping the distribution 2 is chi-squared, that's what they rely on...
#Need to think about the scaling here... what is the right scaling?
#Is R2 even the right metric? as in should it be R instead?

