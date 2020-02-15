## load the vcfcounts GBS data file
## chr	pos	ref	alt	twRef	twHet	twHom	ntRef	ntHet	ntHom	p	                mlog10p	                OR
## 1	33879	C*	T	44	1	3	42	0	3	0.5172413793103448	0.2863067388432749	Infinity
## 1	52665	A*	C	29	0	19	27	0	18	1.0	                -0.0	                NaN
## 1	52671	T*	C	29	2	17	26	2	17	0.38620074440922747	0.4131868934450027	0.896551724137931

vcfcounts = read.table(file="vcfcounts.txt", header=TRUE)
vcfcounts$log10OR = log10(vcfcounts$OR)

rownames(vcfcounts) = paste(vcfcounts$chr,"_",vcfcounts$pos,sep="")
