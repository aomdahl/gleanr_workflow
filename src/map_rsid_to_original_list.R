#!/usr/bin/env Rscript
argv = commandArgs(trailingOnly=TRUE)
pacman::p_load(magrittr, dplyr, data.table)
setwd()
argv<-c("gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.1.pruned_rsids.intermediate.txt",
       "gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.1.union.txt",
       "gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.1.pruned_rsids.txt")
t#packages
############################################################################
# Fix case where there are duplicate RSID-snp mappings- as in a single SNP ID from 1KG maps onto multiple RSIDs. Need to make these
# Ashton Omdahl April 2024
#
# Arguments are: 1) the output map of SNP ids to RSIDs, no header, 2 columns
#                2) The original input list of RSIDs eligible for pruning, should be longer than argument 1
#                3) The desired output file
#
############################################################################
#procedure:
# If we can map duplicates to a unique entry in the original list, keep that one
# Otherwise, just drop them out.

list.to.extract <- fread(argv[1], header=FALSE)
original.query.list <- fread(argv[2], header=FALSE)

dups.ids <- list.to.extract$V1[duplicated(list.to.extract$V1)]
dups.rsid <- (list.to.extract %>% filter(V1 %in% dups.ids))$V2

#side check -are any of the rsids duplicated?
if(length(which(duplicated(list.to.extract$V2))) > 1)
{
  message("Some of the RSIDs are also redundant, beware!")
}

fout.list <- list.to.extract %>% filter(!(V1 %in% dups.ids)) #default behavior is to remove them

#Check if any of the remainng ones already map to what's left....
if(any(fout.list$V2 %in% dups.rsid))
{
  message("Concern that may be mapping to existing RSID, recommend a manual review.")
}

offending.members.rsids <- (filter(original.query.list, V1 %in% dups.rsid))$V1 #these are the rsids that exist in the original data
if(length(offending.members.rsids) > 0)
{
  #paste the original RSID as the option for the final ouutput.
  if(length(offending.members.rsids) > length(dups.ids))
  {
    message("Unexpected scenario- we have redundant SNP-IDs but multiple RSIDs have been mapped onto them from the original data.")
    message("Strongly recommend a manual review")
    quit()
  }
  add.in <- list.to.extract %>% filter(V2 %in% offending.members.rsids)
  if(nrow(add.in) > length(offending.members.rsids))
  {
    message("Redundant RSIDs have arisen. Please perform manual review")
    quit()
  }
  if(length(offending.members.rsids) < length(dups.ids))
  {
    message("We are missing some of the mapped ids. Please perform manual review")
    quit()
  }
  fout.list <- rbind(fout.list,add.in ) %>% arrange()
  
}
if(any(!(fout.list$V2) %in% original.query.list$V1))
{
  
  message("There are RSIDs NOT in the original query list. This will not work for SNP extraction.")
  message("Please provide a mapping of RSIDs to sNP IDs that works for this data")
  quit()
}
message("Writing out final RSID list to ", argv[3])
write.table(fout.list, file=argv[3],col.names = FALSE, row.names = FALSE, quote = FALSE )
