## create various filtered index lists

## "clean" filter:
##  >5 ALT reads in alt sample (so we know it's an ALT)
##  >5 REF reads in ref sample (so we know it's a REF)
##  >9 total reads in mix sample (decent stats)
##  >0 REF reads in mix sample (for finite ratios)
##  >0 ALT reads in mix sample (for finite ratios)
##  <1% REF reads in alt sample
##  <1% ALT reads in ref sample
clean = 
    merged$ALTtot.a>5 &
    merged$REFtot.r>5 &
    merged$TOT.m>9 &
    merged$REFtot.m>0 &
    merged$ALTtot.m>0 &
    merged$REFtot.a/(merged$REFtot.a+merged$ALTtot.a)<0.01 &
    merged$ALTtot.r/(merged$REFtot.r+merged$ALTtot.r)<0.01

## "bigref" filter:
##  clean
##  >100 REF reads for ref sample
bigref =
    clean &
    merged$REFtot.r>100

## "bigalt" filter:
##  clean
##  >100 ALT reads for alt sample
bigalt =
    clean &
    merged$ALTtot.a>100

## "refdeup" filter:
##  clean
##  median-adjusted REF ratio > 16
refdeup =
    clean &
    merged$REFratio/refMixRatio > 16

## "refdedn" filter:
##  clean
##  median-adjusted REF ratio < 1/16
refdedn =
    clean &
    merged$REFratio/refMixRatio < (1/16)

## "altdeup" filter:
##  clean
##  median-adjusted ALT ratio > 16
altdeup =
    clean &
    merged$ALTratio/altMixRatio > 16

## "altdedn" filter:
##  clean
##  median-adjusted ALT ratio < 1/16
altdedn =
    clean &
    merged$ALTratio/altMixRatio < (1/16)
