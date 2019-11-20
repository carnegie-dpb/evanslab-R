##
## load in the VCFSegregation data for ig1 calls
##
B73 = read.table("B73v4.tsv", header=TRUE)
B73.sig = B73[B73$p<0.01,]
