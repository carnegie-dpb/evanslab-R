## load the zygosity data

zyg=read.table(file="zygosity.tsv", header=TRUE)

A = 45
B = 48

## redo odds ratios
zyg$lnAOR = log( (zyg$Ahom+0.5)/(A-zyg$Ahom+0.5) / ( (zyg$Ahet+0.5)/(A-zyg$Ahet+0.5) ) )
zyg$lnBOR = log( (zyg$Bhom+0.5)/(B-zyg$Bhom+0.5) / ( (zyg$Bhet+0.5)/(B-zyg$Bhet+0.5) ) )
