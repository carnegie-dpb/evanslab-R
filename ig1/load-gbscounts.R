## load the vcfcounts GBS data file

##         contig   pos ref alt TWRF TWRR TWAF TWAR NTRF NTRR NTAF NTAR TWRef TWHet TWHom NTRef NTHet NTHom        p   mlog10p  OR log10OR
## 1_51732      1 51732  T*   C    0    0    0    0    2    2    2    2    48     0     0    44     1     0 0.483871 0.3152704   0    -Inf
## 1_51750      1 51750  A*   T    0    0    0    0    2    2    2    2    48     0     0    44     1     0 0.483871 0.3152704   0    -Inf
## 1_52772      1 52772  G*   A    2    2    0    7    0    0    0    0    47     1     0    45     0     0 0.516129 0.2872417 Inf     Inf

gbscounts = read.table(file="gsnap-Zm-B73-REFERENCE-GRAMENE-4.0.GBS.vars.counts.txt.gz", header=TRUE)
gbscounts$log10OR = log10(gbscounts$OR)
rownames(gbscounts) = paste(gbscounts$contig,"_",gbscounts$pos,sep="")
