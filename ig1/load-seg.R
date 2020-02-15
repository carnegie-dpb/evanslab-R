##
## load the segregation data from a seg file
##
## a=NON-TWIN REF
## b=NON-TWIN ALT
## c=TWIN REF
## d=TWIN ALT
##
## NOTE: there may be consecutive calls at the same position, e.g.
## contig     start                                       REF ALT  a  b  c  d size          p   mlog10p signif
##      1 102591307                                         T   C 11 50  5 34  122 0.17886000 0.7474868  false
##      1 102591307 TAGCCAGGGCTTGGTCTCGGGTGAGTCGTGACTATGTCTCC   T 93 44 75 29  274 0.08815873 1.0547347  false

segFile = "gsnap-Zm-B73-REFERENCE-GRAMENE-4.0-combined.seg.txt.gz"

seg = read.table(file=segFile, header=TRUE)

seg$log10OR = log10((seg$a*seg$d)/(seg$b*seg$c))
