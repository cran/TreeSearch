## ----Load-library-------------------------------------------------------------
library('TreeSearch')

## ----Load-data----------------------------------------------------------------
my.data <- TreeSearch::inapplicable.datasets[['Vinther2008']]
my.phyDat <- phangorn::phyDat(my.data, type = 'USER', levels = c(0:9, '-'))

## ----RNG-version--------------------------------------------------------------
# Set a random seed so that random functions in this document are reproducible
suppressWarnings(RNGversion("3.5.0")) # Until we can require R3.6.0
set.seed(0)

## ----Random-tree--------------------------------------------------------------
random.tree <- TreeTools::RandomTree(my.phyDat, root = TRUE)
par(mar = rep(0.25, 4), cex = 0.75) # make plot easier to read
plot(random.tree)
Fitch(random.tree, my.phyDat)

## ----Neighbour-joining--------------------------------------------------------
nj.tree <- TreeTools::NJTree(my.phyDat)
par(mar = rep(0.25, 4), cex = 0.75) # make plot easier to read
plot(nj.tree)
Fitch(nj.tree, my.phyDat)

## ----Outgroup-----------------------------------------------------------------
outgroup <- c('Nemertean', 'Lingula', 'Terebratulina')
rooted.tree <- TreeTools::EnforceOutgroup(nj.tree, outgroup)
par(mar = rep(0.25, 4), cex = 0.75) # make plot easier to read
plot(rooted.tree)

## ----Basic-NNI-search---------------------------------------------------------
better.tree <- TreeSearch(tree = rooted.tree, dataset = my.phyDat, 
                           EdgeSwapper = RootedNNISwap, verbosity = 3)

## ----SPR-and-TBR--------------------------------------------------------------
better.tree <- TreeSearch(better.tree, my.phyDat, maxHits = 8, maxIter = 10000,
                           EdgeSwapper=RootedSPRSwap, verbosity = 2)
better.tree <- TreeSearch(better.tree, my.phyDat, maxHits = 20, maxIter = 40000, 
                           EdgeSwapper=RootedTBRSwap, verbosity = 2)

## ----Ratchet-search-----------------------------------------------------------
best.tree <- Ratchet(better.tree, my.phyDat, verbosity=0, ratchIter=5,
                     swappers=list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap))
attr(best.tree, 'score') # Each tree is labelled with its score during tree search

## ----Plot-best-tree-----------------------------------------------------------
par(mar = rep(0.25, 4), cex = 0.75) # make plot easier to read
plot(best.tree)

## ----Consensus----------------------------------------------------------------
some.optimal.trees <- MultiRatchet(best.tree, nSearch = 20, my.phyDat,
                                   swappers = 
                                   list(RootedTBRSwap, RootedNNISwap))

## ----Plot-strict-consensus----------------------------------------------------
par(mar = rep(0.25, 4), cex = 0.75) # make plot easier to read
strict.consensus <- ape::consensus(some.optimal.trees)
plot(strict.consensus)

## ----Plot-majority-rule-consensus---------------------------------------------
par(mar=rep(0.25, 4), cex=0.75) # make plot easier to read
majority.consensus <- ape::consensus(some.optimal.trees, p = 0.5)
plot(majority.consensus)

## ----Jackknife-annotations----------------------------------------------------
jack.trees <- Jackknife(best.tree, my.phyDat, EdgeSwapper = RootedTBRSwap,
                        jackIter = 20, verbosity = 0)
par(mar = rep(0.25, 4), cex = 0.75) # make plot easier to read

JackLabels(strict.consensus, jack.trees) -> XX

