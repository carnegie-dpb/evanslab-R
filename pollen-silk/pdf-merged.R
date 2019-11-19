## save the merged plots to appropriately-named PDFs
source("plot.merged.R")
pdf(file=paste(reference,mixSampleName,refSampleName,"pdf",sep="."), width=8, height=8)
plot.merged(FALSE)
dev.off()
pdf(file=paste(reference,mixSampleName,altSampleName,"pdf",sep="."), width=8, height=8)
plot.merged(TRUE)
dev.off()

