##
## Write out a wiggle file
##

write.wig = function(file="wiggle.wig", chroms, positions, values) {

    chrom = chroms[1]
    write(file=file, paste("variableStep","\t","chrom=",chrom, sep=""), append=FALSE)

    for (i in 1:length(positions)) {
        if (chroms[i]!=chrom) {
            ## output new chrom line
            chrom = chroms[i];
            write(file=file, paste("variableStep","\t","chrom=",chrom, sep=""), append=TRUE)
        }
        write(file=file, paste(positions[i],values[i], sep="\t"), append=TRUE)
    }

}
