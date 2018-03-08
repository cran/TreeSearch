## ----Load library--------------------------------------------------------
library('TreeSearch')

## ----Load data-----------------------------------------------------------
my.data <- TreeSearch::inapplicable.datasets[['Vinther2008']]
my.phyDat <- phangorn::phyDat(my.data, type='USER', levels=c(0:9, '-'))

## ----Random tree---------------------------------------------------------
# Set random seed so that random functions will generate consistent output in this document
set.seed(0)
random.tree <- RandomTree(my.phyDat)
par(mar=rep(0.25, 4), cex=0.75) # make plot easier to read
plot(random.tree)
Fitch(random.tree, my.phyDat)

## ----Neighbour-joining---------------------------------------------------
nj.tree <- NJTree(my.phyDat)
par(mar=rep(0.25, 4), cex=0.75) # make plot easier to read
plot(nj.tree)
Fitch(nj.tree, my.phyDat)

## ----Outgroup------------------------------------------------------------
outgroup <- c('Nemertean', 'Lingula', 'Terebratulina')
rooted.tree <- EnforceOutgroup(nj.tree, outgroup)
par(mar=rep(0.25, 4), cex=0.75) # make plot easier to read
plot(rooted.tree)

## ----Basic NNI search----------------------------------------------------
better.tree <- TreeSearch(tree=rooted.tree, dataset=my.phyDat, 
                           EdgeSwapper=RootedNNISwap, verbosity=3)

## ----SPR and TBR---------------------------------------------------------
better.tree <- TreeSearch(better.tree, my.phyDat, maxHits=8, maxIter=10000,
                           EdgeSwapper=RootedSPRSwap, verbosity=2)
better.tree <- TreeSearch(better.tree, my.phyDat, maxHits=20, maxIter=40000, 
                           EdgeSwapper=RootedTBRSwap, verbosity=2)

## ----Ratchet search------------------------------------------------------
best.tree <- Ratchet(better.tree, my.phyDat, verbosity=0, ratchIter=5,
                     swappers=list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap))
attr(best.tree, 'score') # Each tree is labelled with its score during tree search

## ----Plot best tree------------------------------------------------------
par(mar=rep(0.25, 4), cex=0.75) # make plot easier to read
plot(best.tree)

## ----Consensus-----------------------------------------------------------
my.consensus <- RatchetConsensus(best.tree, my.phyDat, swappers=
                                 list(RootedTBRSwap, RootedNNISwap))

## ----Plot consensus------------------------------------------------------
par(mar=rep(0.25, 4), cex=0.75) # make plot easier to read
plot(ape::consensus(my.consensus))

