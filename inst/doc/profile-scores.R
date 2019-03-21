## ----Set up, echo=FALSE--------------------------------------------------
require(phangorn)
require(TreeSearch)
data(referenceTree)
data(congreveLamsdellMatrices)
dataset <- congreveLamsdellMatrices[[1]]
suppressWarnings(RNGversion("3.5.0")) # Until we can require R3.6.0
set.seed(0)

## ----Precision testing, cache=TRUE---------------------------------------
preci1 <- PrepareDataProfile(dataset, precision=2e+05) # Quick, imprecise
preci2 <- PrepareDataProfile(dataset, precision=4e+05)
preci3 <- PrepareDataProfile(dataset, precision=8e+05)
info1 <- attr(preci1, 'info.amounts')
info2 <- attr(preci2, 'info.amounts')
info3 <- attr(preci3, 'info.amounts')
diff32 <- as.double(info3 - info2)
hist (diff32, breaks=seq(min(diff32) - 0.002, max(diff32) + 0.005, by=0.002))

if (all_the_time_in_the_world <- FALSE) {
preci4 <- PrepareDataProfile(dataset, precision=1.6e+06)
preci5 <- PrepareDataProfile(dataset, precision=3.2e+06) # Slow, more precise

info4 <- attr(preci4, 'info.amounts')
info5 <- attr(preci5, 'info.amounts')

diff42 <- as.double(info4 - info2)
diff43 <- as.double(info4 - info3)
diff54 <- as.double(info5 - info4)
nonzero <- info4 > 0.00001

hist (diff43)
hist (thisDiff <- diff54); quantile(thisDiff, probs=c(0, 5, 10, 50, 90, 95, 100)/100)
hist (diff42)
hist(100*(diff32 / info4)[nonzero])
hist(100*(diff42 / info4)[nonzero])
hist(100*(diff43 / info4)[nonzero])
}

## ----More histograms-----------------------------------------------------
diff12 <- info1[1:10, ] - info2

hist(diff12, breaks=seq(min(diff12)-0.01, max(diff12)+0.01, by=0.01))

hist(info3 - info2)
hist(info3 - info1[1:10, ])
if (all_the_time_in_the_world) {
hist(info4 - info2)
}

