##
## load the segregation data
##
segFile = "gsnap-Zm-B73-REFERENCE-GRAMENE-4.0-combined.seg.txt.gz"
seg = read.table(file=segFile, header=TRUE)
seg$log10OR = log10((seg$d/seg$c)/(seg$b/seg$a))
