## plot y vs x TPM for genes from htseq-counts

## x = htseq.B73$P.B73.B73.tpm
## y = htseq.B73$S.W22.P.B73.B73.tpm

x = htseq.B73$S.B73.B73.tpm
y = htseq.B73$S.B73.P.W22.B73.tpm

plot(x, y,
     main="B73 alignment: B73 silk+W22 pollen vs. B73 silk",
     xlab="(TPM)",
     ylab="(TPM)",
     cex=0.5,
     log="xy")

lines(c(1e-4,1e5), c(1e-4,1e5)*1e0, col="gray", lwd=3)
lines(c(1e-4,1e5), c(1e-4,1e5)*1e1, col="green", lwd=1)
lines(c(1e-4,1e5), c(1e-4,1e5)*1e2, col="green", lwd=1)
lines(c(1e-4,1e5), c(1e-4,1e5)*1e3, col="green", lwd=1)
lines(c(1e-4,1e5), c(1e-4,1e5)*1e4, col="green", lwd=1)
lines(c(1e-4,1e5), c(1e-4,1e5)*1e-1, col="red", lwd=1)
lines(c(1e-4,1e5), c(1e-4,1e5)*1e-2, col="red", lwd=1)
lines(c(1e-4,1e5), c(1e-4,1e5)*1e-3, col="red", lwd=1)

text(par()$xaxp[1], 10^(par()$usr[4]-0.5), "Lines at factors of ten", pos=4, offset=0)

