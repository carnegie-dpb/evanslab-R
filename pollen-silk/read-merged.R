## Read the merged DF in from a text dump.

sampleTissues = c(
    "SS354+355"="B73 silk", 
    "YX24"="B73 pollen", 
    "SS364"="W22 silk", 
    "PS422"="W22 pollen",
    "S364_354"="B73 pollen x W22 silk",
    "S35464"="W22 pollen x B73 silk"
)

## reference = "B73"
## refMixRatio = 1.81769
## altMixRatio = 0.04887
## mixSampleName = "S35464"    # B73 silk <- W22 pollen
## refSampleName = "SS354+355" # B73 silk REF
## altSampleName = "PS422"     # W22 pollen ALT

## reference = "W22"
## refMixRatio = 0.57338
## altMixRatio = 0.22939
## mixSampleName = "S364_354"  # W22 silk <- B73 pollen
## refSampleName = "SS364"     # W22 silk REF
## altSampleName = "YX24"      # B73 pollen ALT

reference = "B73"
refMixRatio = 0.23941
altMixRatio = 0.55649
mixSampleName = "S364_354"  # B73 pollen -> W22 silk
refSampleName = "YX24"      # B73 pollen REF
altSampleName = "SS364"     # W22 silk ALT

## reference = "W22"
## refMixRatio = 0.04977
## altMixRatio = 1.81575
## mixSampleName = "S35464"    # W22 pollen -> B73 silk
## refSampleName = "PS422"     # W22 pollen REF
## altSampleName = "SS354+355" # B73 silk ALT

merged = read.table(file=paste(reference,mixSampleName,refSampleName,"merged.txt.gz",sep="."))

## create the filter arrays
source("filter-merged.R")

