## load the GBS counts file
##
## chr     pos     ref     alt     twRef   twHet   twHom   ntRef   ntHet   ntHom   p                       mlog10p                 OR
## 1       33879   C*      T       44      1       3       42      0       3       0.5172413793103448      0.2863067388432749      Infinity
## 1       51732   T*      C       48      0       0       44      1       0       0.48387096774189314     0.3152704347786294      0.0
## 1       51750   A*      T       48      0       0       44      1       0       0.48387096774189314     0.3152704347786294      0.0

gbscountsFile = readline(prompt="GBS counts file: ")
gbscounts = read.table(file=gbscountsFile, header=TRUE)
gbscounts = gbscounts[!is.nan(gbscounts$p),]

gbscounts$log10OR = log10(gbscounts$OR)
rownames(gbscounts) = paste(gbscounts$chr,"_",gbscounts$pos,sep="")
