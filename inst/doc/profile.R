## ----Load-library, message=FALSE----------------------------------------------
library('TreeSearch')

## ----RNG-version--------------------------------------------------------------
# Set a random seed so that random functions in this document are reproducible
suppressWarnings(RNGversion('3.5.0')) # Until we can require R3.6.0
set.seed(888)

## ----Load Longrich data-------------------------------------------------------
data(congreveLamsdellMatrices)
myMatrix <- congreveLamsdellMatrices[[10]]

## ----Addition tree------------------------------------------------------------
additionTree <- AdditionTree(myMatrix, concavity = 'profile')
TreeLength(additionTree, myMatrix, 'profile')
par(mar = rep(0.25, 4), cex = 0.75) # make plot easier to read
plot(additionTree)

## ----Random tree--------------------------------------------------------------
randomTree <- TreeTools::RandomTree(myMatrix, root = TRUE)
TreeLength(randomTree, myMatrix, 'profile')
njTree <- TreeTools::NJTree(myMatrix)
TreeLength(njTree, myMatrix, 'profile')

## ----Starting score, message = FALSE------------------------------------------
betterTrees <- MaximizeParsimony(myMatrix, additionTree, concavity = 'profile',
                                 ratchIter = 3, tbrIter = 3, maxHits = 8)

## ----Ratchet-search-results---------------------------------------------------
TreeLength(betterTrees[[1]], myMatrix, 'profile')
par(mar = rep(0.25, 4), cex = 0.75) # make plot easier to read
plot(ape::consensus(betterTrees))

## ----Suboptimal sampling, message = FALSE-------------------------------------
suboptimals <- MaximizeParsimony(myMatrix, betterTrees, tolerance = 3,
                                 ratchIter = 2, tbrIter = 3,
                                 maxHits = 25,
                                 concavity = 'profile')

## ----Plot suboptimal consensus------------------------------------------------
par(mar = rep(0.25, 4), cex = 0.75)
table(signif(TreeLength(suboptimals, myMatrix, 'profile')))
plot(ape::consensus(suboptimals))

