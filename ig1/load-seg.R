##
## load the segregation data
##

seg = read.table(file="gsnap-B73v4-merged.seg.txt", header=TRUE)

seg$log10OR = log10((seg$d/seg$c)/(seg$b/seg$a))



