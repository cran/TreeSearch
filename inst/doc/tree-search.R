## ----Load-library-------------------------------------------------------------
library("TreeSearch")

## ----Load-data----------------------------------------------------------------
rawData <- TreeSearch::inapplicable.datasets[['Vinther2008']]
vinther <- phangorn::phyDat(rawData, type = 'USER', levels = c(0:9, '-'))

## ----RNG-version--------------------------------------------------------------
# Set a random seed so that random functions in this document are reproducible
suppressWarnings(RNGversion("3.5.0")) # Until we can require R3.6.0
set.seed(0)

## ----first-pass, message = FALSE----------------------------------------------
bestTrees <- MaximizeParsimony(vinther)

## ----inspect-progress---------------------------------------------------------
firstHit <- attr(bestTrees, 'firstHit')
firstHit

## ----map-search, fig.asp = 1--------------------------------------------------
distances <- TreeDist::ClusteringInfoDistance(bestTrees)
searchStages <- length(firstHit)
map <- cmdscale(distances, k = 3)
cols <- hcl.colors(searchStages, alpha = 0.8)
par(mar = rep(0, 4))
TreeDist::Plot3(map,
                col = cols[rep(seq_along(firstHit), firstHit)],
                pch = 16, cex = 2,
                axes = FALSE, xlab = '', ylab = '', asp = 1)
TreeTools::MSTEdges(distances, plot = TRUE, map[, 1], map[, 2],
                    col = '#00000030', lty = 2)
legend('topright', names(firstHit), col = cols, pch = 16, bty = 'n')

## ----second-pass, message = FALSE---------------------------------------------
bestTrees <- MaximizeParsimony(vinther, tree = bestTrees,
                               ratchIter = 6L,
                               tbrIter = 4L, 
                               finalIter = 3L,
                               maxHits = 80L)

## ----plot-tree----------------------------------------------------------------
par(mar = rep(0.25, 4), cex = 0.75) # make plot easier to read
plot(ape::consensus(bestTrees))
TreeLength(bestTrees[[1]], vinther)

## ----plot-label-nodes---------------------------------------------------------
par(mar = rep(0.25, 4), cex = 0.75) # make plot easier to read
cons <- ape::consensus(bestTrees, p = 0.5)
splitFreqs <- TreeTools::SplitFrequency(cons, bestTrees) / length(bestTrees)
plot(cons)
TreeTools::LabelSplits(cons, round(splitFreqs * 100), unit = '%',
                       col = TreeTools::SupportColor(splitFreqs),
                       frame = 'none', pos = 3L)

## ----Jackknife-annotations----------------------------------------------------
nReplicates <- 10
jackTrees <- lapply(logical(nReplicates), function (x)
  Resample(vinther, bestTrees, ratchIter = 0, tbrIter = 1, maxHits = 4,
          verbosity = 0)
)

par(mar = rep(0, 4), cex = 0.8)
# Take the strict consensus of all trees for each replicate
JackLabels(cons, lapply(jackTrees, ape::consensus)) -> XX

## ----concordance--------------------------------------------------------------
concordance <- QuartetConcordance(cons, vinther)

# Alternative measures:
# concordance <- ClusteringConcordance(cons, vinther)
# concordance <- PhylogeneticConcordance(cons, vinther)

par(mar = rep(0, 4), cex = 0.8)
plot(cons)
TreeTools::LabelSplits(cons, signif(concordance, 3),
                       col = TreeTools::SupportColor(concordance / max(concordance)),
                       frame = 'none', pos = 3L)

## ----iw-search, message = FALSE-----------------------------------------------
iwTrees <- MaximizeParsimony(vinther, concavity = 10)
par(mar = rep(0.25, 4), cex = 0.75) # make plot easier to read
plot(ape::consensus(iwTrees))

## ----simple-constraints, message = FALSE--------------------------------------
library("TreeTools", quietly = TRUE, warn.conflicts = FALSE)
constraint <- MatrixToPhyDat(c(a = 1, b = 1, c = 0, d = 0, e = 0, f = 0))
characters <- MatrixToPhyDat(matrix(
  c(0, 1, 1, 1, 0, 0,
    1, 1, 1, 0, 0, 0), ncol = 2,
  dimnames = list(letters[1:6], NULL)))
plot(MaximizeParsimony(characters, constraint = constraint,
                       verbosity = -1)[[1]])

## ----complex-constraints, message = FALSE-------------------------------------
characters <- MatrixToPhyDat(matrix(
  c(0, 0, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 0, 0, 0), ncol = 2,
  dimnames = list(letters[1:7], NULL)))
constraint <- MatrixToPhyDat(matrix(
  c(0, 0, 1, '?', 1, 1,
    1, 1, 1,   1, 0, 0), ncol = 2,
  dimnames = list(letters[1:6], NULL)))
plot(MaximizeParsimony(characters, constraint = constraint,
                       verbosity = -1)[[1]])

