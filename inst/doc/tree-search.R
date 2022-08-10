## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(fig.width = 7.2, fig.asp = 0.7) 

## ----Load-library-------------------------------------------------------------
library("TreeSearch")

## ----Load-data----------------------------------------------------------------
vinther <- TreeSearch::inapplicable.phyData[["Vinther2008"]]

## ----RNG-version--------------------------------------------------------------
# Set a random seed so that random functions in this document are reproducible
suppressWarnings(RNGversion("3.5.0")) # Until we require R v3.6.0
set.seed(0)

## ----first-pass, message = FALSE----------------------------------------------
bestTrees <- MaximizeParsimony(vinther)

## ----inspect-progress---------------------------------------------------------
firstHit <- attr(bestTrees, "firstHit")
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
                axes = FALSE, xlab = "", ylab = "", asp = 1)
TreeTools::MSTEdges(distances, plot = TRUE, map[, 1], map[, 2],
                    col = "#00000030", lty = 2)
legend("topright", names(firstHit), col = cols, pch = 16, bty = "n")

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
majCons <- ape::consensus(bestTrees, p = 0.5)
splitFreqs <- TreeTools::SplitFrequency(majCons, bestTrees) / length(bestTrees)
plot(majCons)
TreeTools::LabelSplits(majCons, round(splitFreqs * 100), unit = "%",
                       col = TreeTools::SupportColor(splitFreqs),
                       frame = "none", pos = 3L)

## ----Jackknife-annotations----------------------------------------------------
nReplicates <- 10
jackTrees <- lapply(logical(nReplicates), function (x)
  Resample(vinther, bestTrees, ratchIter = 0, tbrIter = 1, maxHits = 4,
           verbosity = 0)
)

strict <- ape::consensus(bestTrees, p = 1)

par(mar = rep(0, 4), cex = 0.8)
# Take the strict consensus of all trees for each replicate
JackLabels(strict, lapply(jackTrees, ape::consensus)) -> XX

## ----concordance--------------------------------------------------------------
concordance <- QuartetConcordance(strict, vinther)

# Alternative measures:
# concordance <- ClusteringConcordance(strict, vinther)
# concordance <- PhylogeneticConcordance(strict, vinther)

par(mar = rep(0, 4), cex = 0.8)
plot(strict)
TreeTools::LabelSplits(strict, signif(concordance, 3),
                       col = TreeTools::SupportColor(concordance / max(concordance)),
                       frame = "none", pos = 3L)

## ----stability----------------------------------------------------------------
par(mar = rep(0, 4), cex = 0.8)

plot(strict, tip.color = Rogue::ColByStability(bestTrees))

## ----find-rogues--------------------------------------------------------------
Rogue::QuickRogue(bestTrees, p = 1)

## ----cons-without-halk--------------------------------------------------------
par(mar = rep(0, 4), cex = 0.8)
noWiwaxia <- lapply(bestTrees, TreeTools::DropTip, "Wiwaxia")
plot(ape::consensus(noWiwaxia), tip.color = Rogue::ColByStability(noWiwaxia))

## ----restore-wiwaxia----------------------------------------------------------
par(mar = rep(0, 4), cex = 0.8)
xx <- TreeTools::RoguePlot(bestTrees, "Wiwaxia", p = 1)

## ----iw-search, message = FALSE-----------------------------------------------
iwTrees <- MaximizeParsimony(vinther, concavity = 10)
par(mar = rep(0.25, 4), cex = 0.75) # make plot easier to read
plot(ape::consensus(iwTrees))

## ----simple-constraints, fig.width = 4, fig.align = "center", message = FALSE----
library("TreeTools", quietly = TRUE)
constraint <- MatrixToPhyDat(c(a = 1, b = 1, c = 0, d = 0, e = 0, f = 0))
characters <- MatrixToPhyDat(matrix(
  c(0, 1, 1, 1, 0, 0,
    1, 1, 1, 0, 0, 0), ncol = 2,
  dimnames = list(letters[1:6], NULL)))
plot(MaximizeParsimony(characters, constraint = constraint,
                       verbosity = -1)[[1]])

## ----complex-constraints, fig.width = 4, fig.align = "center", message = FALSE----
characters <- MatrixToPhyDat(matrix(
  c(0, 0, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 0, 0, 0), ncol = 2,
  dimnames = list(letters[1:7], NULL)))
constraint <- MatrixToPhyDat(matrix(
  c(0, 0, 1, "?", 1, 1,
    1, 1, 1,   1, 0, 0), ncol = 2,
  dimnames = list(letters[1:6], NULL)))
plot(MaximizeParsimony(characters, constraint = constraint,
                       verbosity = -1)[[1]])

