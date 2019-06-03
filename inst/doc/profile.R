## ----Load-library, message=FALSE-----------------------------------------
library('TreeSearch')

## ----RNG-version, warning=FALSE------------------------------------------
# Set a random seed so that random functions in this document are reproducible
RNGversion("3.5.0") # Until we can require R3.6.0
set.seed(888)

## ----Load Longrich data--------------------------------------------------
data(congreveLamsdellMatrices)
my.data <- congreveLamsdellMatrices[[10]]
my.phyDat <- phangorn::phyDat(my.data, type='USER', levels=c(1, 2))

## ----Prepare the data for analysis, warning=FALSE------------------------
my.prepdata <- PrepareDataProfile(my.phyDat, precision=4e+04)

## ----Random tree, eval=FALSE---------------------------------------------
#  tree <- RandomTree(my.phyDat)

## ----NJ Tree-------------------------------------------------------------
tree <- NJTree(my.phyDat)
par(mar=rep(0.25, 4), cex=0.75) # make plot easier to read
plot(tree)

## ----Starting score------------------------------------------------------
ProfileScore(tree, my.prepdata)

better.tree <- ProfileTreeSearch(tree, my.prepdata, EdgeSwapper=RootedTBRSwap)

## ----longwinded, eval=FALSE----------------------------------------------
#  # Longwinded approach:
#  better.tree <- Ratchet(better.tree, my.prepdata, searchHits=10, searchIter=100, ratchIter=5,
#                         swappers=list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap),
#                         InitializeData=ProfileInitMorphy, CleanUpData=ProfileDestroyMorphy,
#                         TreeScorer=ProfileScoreMorphy, Bootstrapper=ProfileBootstrap)

## ----Cached results of ratchet search, echo=FALSE------------------------
# The ratchet search takes a little while to run,
# so I've cached its results
better.tree <- structure(list(edge = structure(c(23L, 23L, 24L, 25L, 26L, 26L, 
27L, 28L, 28L, 29L, 29L, 30L, 30L, 27L, 31L, 31L, 25L, 32L, 32L, 
33L, 33L, 34L, 34L, 35L, 35L, 36L, 36L, 37L, 37L, 38L, 38L, 39L, 
39L, 40L, 40L, 41L, 41L, 24L, 42L, 42L, 43L, 43L, 1L, 24L, 25L, 
26L, 13L, 27L, 28L, 14L, 29L, 15L, 30L, 17L, 16L, 31L, 18L, 19L, 
32L, 12L, 33L, 11L, 34L, 10L, 35L, 9L, 36L, 8L, 37L, 7L, 38L, 
6L, 39L, 5L, 40L, 4L, 41L, 3L, 2L, 42L, 20L, 43L, 22L, 21L), .Dim = c(42L, 
2L)), tip.label = c("1", "2", "3", "4", "5", "6", "7", "8", "9", 
"10", "11", "12", "19", "18", "17", "15", "16", "13", "14", "20", 
"21", "22"), Nnode = 21L), class = "phylo", order = "cladewise", score = -309.16150317298, hits = 2)

## ----Ratchet search, eval=FALSE------------------------------------------
#  # Equivalent, but less typing!
#  RootedSwappers <- list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap)
#  better.tree <- ProfileRatchet(better.tree, my.prepdata,
#                                swappers=RootedSwappers,
#                                searchHits=10, searchIter=100, ratchIter=5)
#  

## ----Ratchet-search-results----------------------------------------------
attr(better.tree, 'score')
par(mar=rep(0.25, 4), cex=0.75) # make plot easier to read
plot(better.tree)

## ----Suboptimal sampling, eval=FALSE-------------------------------------
#  suboptimals <- ProfileRatchet(better.tree, my.prepdata,
#                                swappers=list(RootedTBRSwap),
#                                returnAll=TRUE, suboptimal=5,
#                                ratchHits=25, ratchIter=500,
#                                bootstrapHits=15, bootstrapIter=450,
#                                searchHits=10, searchIter=100)

## ----cached-SoS, echo=FALSE----------------------------------------------
# Use cached results of last block to reduce compilation time.
# A slow compile means that CRAN can't check incoming vignettes.
# 
tip.labels <- c("1", "2", "3", "4", 
"5", "6", "7", "8", "9", "10", "11", "12", "19", "18", "17", 
"15", "16", "13", "14", "20", "21", "22")
suboptimals <- 
structure(list(structure(list(edge = structure(c(23L, 23L, 24L, 
25L, 26L, 26L, 27L, 28L, 28L, 29L, 29L, 30L, 30L, 27L, 31L, 31L, 
25L, 32L, 32L, 33L, 33L, 34L, 34L, 35L, 35L, 36L, 36L, 37L, 37L, 
38L, 38L, 39L, 39L, 40L, 40L, 41L, 41L, 24L, 42L, 42L, 43L, 43L, 
1L, 24L, 25L, 26L, 13L, 27L, 28L, 14L, 29L, 15L, 30L, 16L, 17L, 
31L, 19L, 18L, 32L, 12L, 33L, 11L, 34L, 10L, 35L, 9L, 36L, 8L, 
37L, 7L, 38L, 6L, 39L, 5L, 40L, 4L, 41L, 2L, 3L, 42L, 20L, 43L, 
21L, 22L), .Dim = c(42L, 2L)), tip.label = tip.labels, Nnode = 21L), class = "phylo", order = "cladewise", score = -307.586176470711, hits = 2), 
    structure(list(edge = structure(c(23L, 23L, 24L, 25L, 26L, 
    26L, 27L, 28L, 28L, 29L, 29L, 30L, 30L, 27L, 31L, 31L, 25L, 
    32L, 32L, 33L, 33L, 34L, 34L, 35L, 35L, 36L, 36L, 37L, 37L, 
    38L, 39L, 39L, 38L, 40L, 40L, 41L, 41L, 24L, 42L, 42L, 43L, 
    43L, 1L, 24L, 25L, 26L, 13L, 27L, 28L, 14L, 29L, 15L, 30L, 
    17L, 16L, 31L, 18L, 19L, 32L, 12L, 33L, 11L, 34L, 10L, 35L, 
    9L, 36L, 8L, 37L, 7L, 38L, 39L, 5L, 6L, 40L, 4L, 41L, 3L, 
    2L, 42L, 20L, 43L, 22L, 21L), .Dim = c(42L, 2L)), tip.label = tip.labels,
    Nnode = 21L), class = "phylo", order = "cladewise", score = -304.566939215937, hits = 2), 
    structure(list(edge = structure(c(23L, 23L, 24L, 25L, 26L, 
    26L, 27L, 27L, 28L, 28L, 29L, 29L, 30L, 30L, 31L, 31L, 32L, 
    32L, 33L, 33L, 34L, 34L, 35L, 35L, 36L, 36L, 25L, 37L, 38L, 
    38L, 39L, 39L, 40L, 40L, 37L, 41L, 41L, 24L, 42L, 42L, 43L, 
    43L, 1L, 24L, 25L, 26L, 13L, 27L, 12L, 28L, 11L, 29L, 10L, 
    30L, 9L, 31L, 8L, 32L, 7L, 33L, 6L, 34L, 5L, 35L, 4L, 36L, 
    3L, 2L, 37L, 38L, 14L, 39L, 15L, 40L, 17L, 16L, 41L, 18L, 
    19L, 42L, 20L, 43L, 22L, 21L), .Dim = c(42L, 2L)), tip.label = tip.labels,
    Nnode = 21L), class = "phylo", order = "cladewise", score = -305.74605645655, hits = 2), 
    structure(list(edge = structure(c(23L, 23L, 24L, 25L, 26L, 
    26L, 27L, 28L, 28L, 29L, 29L, 30L, 30L, 27L, 31L, 31L, 25L, 
    32L, 32L, 33L, 33L, 34L, 34L, 35L, 35L, 36L, 36L, 37L, 37L, 
    38L, 38L, 39L, 39L, 40L, 40L, 41L, 41L, 24L, 42L, 42L, 43L, 
    43L, 1L, 24L, 25L, 26L, 13L, 27L, 28L, 14L, 29L, 15L, 30L, 
    17L, 16L, 31L, 18L, 19L, 32L, 12L, 33L, 11L, 34L, 10L, 35L, 
    9L, 36L, 7L, 37L, 8L, 38L, 6L, 39L, 5L, 40L, 4L, 41L, 3L, 
    2L, 42L, 20L, 43L, 22L, 21L), .Dim = c(42L, 2L)), tip.label = tip.labels,
    Nnode = 21L), class = "phylo", order = "cladewise",
    score = -304.733211428379, hits = 2), 
    structure(list(edge = structure(c(23L, 23L, 24L, 25L, 25L, 
    26L, 27L, 28L, 28L, 27L, 29L, 29L, 30L, 30L, 31L, 31L, 26L, 
    32L, 32L, 33L, 33L, 34L, 34L, 35L, 35L, 36L, 36L, 37L, 37L, 
    38L, 38L, 39L, 39L, 40L, 40L, 41L, 41L, 24L, 42L, 42L, 43L, 
    43L, 1L, 24L, 25L, 13L, 26L, 27L, 28L, 19L, 18L, 29L, 14L, 
    30L, 15L, 31L, 17L, 16L, 32L, 12L, 33L, 11L, 34L, 10L, 35L, 
    9L, 36L, 8L, 37L, 7L, 38L, 6L, 39L, 5L, 40L, 4L, 41L, 3L, 
    2L, 42L, 20L, 43L, 22L, 21L), .Dim = c(42L, 2L)), tip.label = tip.labels,
    Nnode = 21L), class = "phylo", order = "cladewise", 
    score = -305.74605645655, hits = 2), 
    structure(list(edge = structure(c(23L, 23L, 24L, 25L, 26L, 
    27L, 27L, 28L, 28L, 29L, 29L, 30L, 30L, 31L, 31L, 32L, 32L, 
    33L, 33L, 34L, 34L, 35L, 35L, 36L, 36L, 37L, 37L, 26L, 38L, 
    38L, 25L, 39L, 39L, 40L, 40L, 41L, 41L, 24L, 42L, 42L, 43L, 
    43L, 1L, 24L, 25L, 26L, 27L, 13L, 28L, 12L, 29L, 11L, 30L, 
    10L, 31L, 9L, 32L, 8L, 33L, 7L, 34L, 6L, 35L, 5L, 36L, 4L, 
    37L, 3L, 2L, 38L, 19L, 18L, 39L, 14L, 40L, 15L, 41L, 16L, 
    17L, 42L, 20L, 43L, 21L, 22L), .Dim = c(42L, 2L)), tip.label = tip.labels,
    Nnode = 21L), class = "phylo", order = "cladewise", 
    score = -306.478924175555, hits = 2), 
    structure(list(edge = structure(c(23L, 23L, 24L, 25L, 26L, 
    26L, 27L, 28L, 28L, 29L, 29L, 30L, 30L, 27L, 31L, 31L, 25L, 
    32L, 32L, 33L, 33L, 34L, 34L, 35L, 35L, 36L, 36L, 37L, 37L, 
    38L, 38L, 39L, 39L, 40L, 40L, 41L, 41L, 24L, 42L, 42L, 43L, 
    43L, 1L, 24L, 25L, 26L, 13L, 27L, 28L, 14L, 29L, 15L, 30L, 
    17L, 16L, 31L, 18L, 19L, 32L, 11L, 33L, 12L, 34L, 10L, 35L, 
    9L, 36L, 8L, 37L, 7L, 38L, 6L, 39L, 5L, 40L, 4L, 41L, 3L, 
    2L, 42L, 20L, 43L, 22L, 21L), .Dim = c(42L, 2L)), tip.label = tip.labels,
    Nnode = 21L), class = "phylo", order = "cladewise",
    score = -305.69089947461, hits = 2), 
    structure(list(edge = structure(c(23L, 23L, 24L, 25L, 26L, 
    27L, 27L, 28L, 28L, 29L, 29L, 30L, 30L, 31L, 31L, 32L, 32L, 
    33L, 33L, 34L, 34L, 35L, 35L, 36L, 36L, 37L, 37L, 26L, 38L, 
    38L, 39L, 39L, 40L, 40L, 25L, 41L, 41L, 24L, 42L, 42L, 43L, 
    43L, 1L, 24L, 25L, 26L, 27L, 13L, 28L, 12L, 29L, 11L, 30L, 
    10L, 31L, 9L, 32L, 8L, 33L, 7L, 34L, 6L, 35L, 5L, 36L, 4L, 
    37L, 2L, 3L, 38L, 14L, 39L, 15L, 40L, 16L, 17L, 41L, 19L, 
    18L, 42L, 20L, 43L, 21L, 22L), .Dim = c(42L, 2L)), tip.label = tip.labels,
    Nnode = 21L), class = "phylo", order = "cladewise", 
    score = -305.74605645655, hits = 2), 
    structure(list(edge = structure(c(23L, 23L, 24L, 25L, 25L, 
    26L, 26L, 27L, 28L, 28L, 29L, 29L, 30L, 30L, 31L, 31L, 32L, 
    32L, 33L, 33L, 34L, 34L, 35L, 35L, 36L, 36L, 37L, 37L, 38L, 
    38L, 27L, 39L, 39L, 40L, 40L, 41L, 41L, 24L, 42L, 42L, 43L, 
    43L, 1L, 24L, 25L, 19L, 26L, 18L, 27L, 28L, 13L, 29L, 12L, 
    30L, 11L, 31L, 10L, 32L, 9L, 33L, 8L, 34L, 7L, 35L, 6L, 36L, 
    5L, 37L, 4L, 38L, 3L, 2L, 39L, 14L, 40L, 15L, 41L, 17L, 16L, 
    42L, 20L, 43L, 22L, 21L), .Dim = c(42L, 2L)), tip.label = tip.labels,
    Nnode = 21L), class = "phylo", order = "cladewise",
    score = -304.193173481336, hits = 2)), class = "multiPhylo",
    old.index = c(1L, 
    2L, 3L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 20L, 21L, 1L, 1L, 1L, 20L, 2L, 1L, 1L, 29L, 1L, 1L, 21L, 
    1L, 3L, 1L, 1L, 37L, 1L, 20L, 21L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 21L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 20L, 1L, 88L, 20L, 1L, 1L, 1L, 2L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 29L, 1L, 1L, 1L, 20L, 20L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 20L, 
    1L, 1L, 88L, 1L, 1L, 2L, 1L, 1L, 3L, 1L, 1L, 20L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 21L, 1L, 1L, 1L, 1L, 1L, 
    20L, 1L, 1L, 1L, 1L, 1L, 1L, 178L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    29L, 1L, 1L, 1L, 37L, 1L, 37L, 3L, 1L, 1L, 1L, 20L, 1L, 2L, 1L, 
    1L, 1L, 1L, 1L, 1L, 29L, 1L, 1L, 1L, 20L, 1L, 1L, 1L, 1L, 21L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 37L, 1L, 1L, 1L, 1L, 1L, 1L, 20L, 37L, 1L, 1L, 1L, 1L, 
    20L, 1L, 21L, 1L, 88L, 1L, 1L, 1L, 1L, 1L, 1L, 3L, 37L, 1L, 3L, 
    1L, 1L, 1L))

## ----Plot suboptimal consensus-------------------------------------------
par(mar=rep(0.25, 4), cex=0.75)
plot(my.consensus <- ape::consensus(suboptimals))

