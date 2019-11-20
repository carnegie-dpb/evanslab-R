##
## load in the VCFSegregation data for ig1 calls
##
Mo17 = read.table("Mo17.tsv", header=TRUE)
Mo17.sig = Mo17[Mo17$p<0.01,]
