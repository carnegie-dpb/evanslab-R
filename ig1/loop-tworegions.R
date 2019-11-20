##
## convenience to spin through the two-region plots from a dataframe of synteny
##


for (i in 1:length(nine.on.three$length)) {
    plot.logp.tworegions(df=nine.on.three, index=i)
    readline(prompt="Press [enter] to continue")
}
    
