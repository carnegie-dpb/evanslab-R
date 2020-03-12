## load the GBS counts file
##
## contig  pos     ref     alt     TWRF    TWRR    TWAF    TWAR    NTRF    NTRR    NTAF    NTAR    TWRef   TWHet   TWHom   NTRef   NTHet   NTHom   p       mlog10p OR

gbscountsFile = readline(prompt="GBS counts file: ")
gbscounts = read.table(file=gbscountsFile, header=TRUE)
gbscounts = gbscounts[!is.nan(gbscounts$p),]

gbscounts$log10OR = log10(gbscounts$OR)
##rownames(gbscounts) = paste(gbscounts$chr,"_",gbscounts$pos,sep="")
