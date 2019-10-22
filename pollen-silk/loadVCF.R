## Loads VCF data from a text file that has been created with  edu.carnegiescience.dpb.evanslab.VCFLoader.
## This is VASTLY faster than parsing the DP4 field in R!
loadVCF = function(vcfTxtFile) {

    vcf = read.table(vcfTxtFile, header=FALSE, sep="\t")

    ##CHROM  POS REF ALT REFfor REFrev ALTfor ALTrev
    colnames(vcf) = c("CHR","POS","REF","ALT","REFfor","REFrev","ALTfor","ALTrev")

    ## total forward+reverse
    vcf$REFtot = vcf$REFfor + vcf$REFrev
    vcf$ALTtot = vcf$ALTfor + vcf$ALTrev

    return(vcf)
}
