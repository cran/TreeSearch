
########Sectorial <- function (tree, dataset, TreeScorer = TODO, sectRearrangements=list(RootedNNI), 
########                             searchRearrangements=list(RootedNNI, RootedTBR, RootedNNI), 
########                             maxHits=c(30, 40, 60), maxIter=2000, verbosity=3, ...) {
########  best.score <- attr(tree, "score")
########  if (is.null(treeOrder <- attr(tree, "order")) || treeOrder != "preorder") tree <- Preorder(tree)
######## 
########  tree <- RenumberTips(tree, names(dataset))
########  if (length(best.score) == 0) best.score <- TreeScorer(tree, dataset, ...)
########  sect <- DoSectorial(tree, dataset, verbosity=verbosity-1, maxit=30, 
########    maxIter=max(maxIter), maxHits=15, smallestSector=6, 
########    largestSector=dim(tree[["edge"]])[1]*0.25, Rearrangements=sectRearrangements)
########  for (i in seq_along(sectRearrangements)) {
########    iters <- if (length(maxIter) <= i) maxIter[[i]] else min(maxIter)
########    hits  <- if (length(maxHits) <= i) maxHits[[i]] else min(maxHits)
########    sect <- TreeSearch(sect, dataset, TreeScorer, sectRearrangements[[i]], maxIter=iters,
########                         maxHits=hits, verbosity=verbosity-1)
########  }
########  if (attr(sect, "score") <= best.score) {
########    return (sect)
########  } else return (tree)
########}
########
#########' Sector Data
#########' Check that chosen sector contains parsimony-informative data
#########'
#########' The function simply checks that some characters have more than one state.
#########" It"s crude, but the cost of a false positive is low.
#########'
#########' @param dataset the dataset to subsample
#########' @param tips character vector listing tips that exist in the sector 
#########'
#########' @keywords internal
#########' @export
########SectorHasData <- function (dataset, tips) {
########  if (class(dataset) =="phyDat") {
########    levs <- attr(dataset, "levels")
########    contrast <- attr(dataset, "contrast")
########    index <- as.integer(contrast %*% 2L ^ (seq_along(attr(dataset, "levels")) - 1))
########    characters <- vapply(dataset, function (X) index[X], integer(attr(dataset, "nr")))
########  } else if (names(dataset)) {
########    characters <- vapply(dataset, c, integer(length(dataset[[1]])))
########  } else if(rownames(dataset)) {
########    characters <- t(dataset)
########  } else characters <- dataset
########  tokens <- apply(characters[, tips], 2, function (x) length(unique(x)))
########  return (any(tokens > 1))
########}
########
#######' Sectorial Search with inapplicable data
#######'
#######' \code{DoSectorial} is called by Sectorial
#######'
#######' @template preorderTreeParam
#######' @param dataset a dataset in the format expected by \code{TreeScorer}
#######' @param TreeScorer a function that will score a tree topology
#######' @param maxSectIter maximum number of sectorial iterations to perform
#######' @param maxIter maximum number of iterations to perform in tree rearrangement functions
#######' @param maxImprovements maximum number of times to find an optimal score before ending sectorial search
#######' @param Rearrangements a list of tree rearrangement functions that retain the root of the tree
#######'                       (e.g. \code{list(RootedSPR)})
#######' @param smallestSector integer giving size of smallest sector to rearrange
#######' @param largestSector integer giving size of largest sector to rearrange (rounded down if non-integral)
#######' @template verbosityParam
#######' @template treeScorerDots
#######'  
#######' @return a tree of class \code{phylo} with a \code{TreeScorer} score as good or better than that of \code{tree}
#######' 
#######' @author Martin R. Smith
#######' @importFrom ape root
#######' @export
######MorphySectorial <- function (parent, child, dataset, TreeScorer = MorphyLength, maxSectIter=100, 
######                         maxIter=500, maxImprovements=5, smallestSector=4, largestSector=1e+06, 
######                         Rearrangements=list(RootedNNI), verbosity=0, ...) {
######  if (verbosity >= 0) message(" - Sectorial search: optimizing sectors of", smallestSector, "to", floor(largestSector), "tips")
######  nEdge <- length(parent)
######  nTip <- (nEdge / 2) + 1
######  nonRootNodes <- (nTip + 2):(nEdge + 1)
######  
######  epsilon <- 1e-08
######  improvements <- 1
######  
######  NodeChildren <- function (parent, child) {
######    result <- integer(nEdge + 1L)
######    result[1:nTip] <- 1
######    edgeCounted <- child <= nTip
######    while (any(!counted)) {
######      
######      result[!counted] <- result[!counted] + 
######      
######    }
######  }
######  for (i in seq_len(maxSectIter)) {
######    
######    nodeLengths <- CladeSizes(tree, nonRootNodes)
######    candidateNodes <- nonRootNodes[nodeLengths >= smallestSector & nodeLengths <= largestSector]
######    if (verbosity >= 0) message("\n - Iteration ", i, "- attempting sectorial search on node ")
######    repeat {
######      sector <- sample(candidateNodes, 1)
######      candidate <- Subtree(tree, sector)
######      crownTips <- candidate[["tip.label"]]
######      sectorSize <- length(crownTips)
######      message(sector, "size", sectorSize, "...")
######      
######      if (SectorHasData(dataset, crownTips)) break else message("unsuitable (no dataset); trying")
######      
######      candidateNodes <- candidateNodes[-match(sector, candidateNodes)]
######      if (length(candidateNodes == 0)) stop("No selectable sectors contain parsimony information! Either "largestSector" is close to "smallestSector" or your dataset is short of parsimony information.")
######    }
######    if (verbosity >= 0) message(" Sector OK.")
######    
######    crown <- root(AddTip(crown, 0, "SECTOR_ROOT"), length(crown[["tip.label"]]) + 1, resolve.root=TRUE)
######    initialScore <- TreeScorer(candidate, dataset, ...)
######    attr(candidate, "score") <- initialScore
######    
######    if (verbosity >= 0) message("\n - Rearranging sector", sector)
######    for (Rearrange in Rearrangements) {
######      candidate <- TreeSearch(candidate, dataset, TreeScorer, Rearrange,
######                                verbosity=verbosity-1, maxIter=maxIter, ...) 
######    }
######    candidateScore <- attr(candidate, "score")
######    
######    if((candidateScore + epsilon) < initialScore) {
######      improvements <- improvements + 1
######      
######      subtree.labels <- crownTips
######      subtree.nTips <- sectorSize
######      subtree.edge <- candidate[["edge"]]
######      subtree.parent <- subtree.edge[, 1]
######      subtree.child  <- subtree.edge[, 2]
######           
######      isTip <- subtree.child <= subtree.nTips
######      subtree.child[isTip] <- match(crownTips[subtree.child[isTip]], tipOrder)
######      nodeAdjust <- sector - (subtree.nTips + 1)
######
######      subtree.child[!isTip] <- subtree.child[!isTip] + nodeAdjust
######      edges <- which(tree[["edge"]][, 2] == sector) + seq_along(subtree.parent)
######      tree[["edge"]][edges, 1] <- subtree.parent + nodeAdjust
######      tree[["edge"]][edges, 2] <- subtree.child
######      
######      if (verbosity > 0) message(" : improved local pscore, updated tree")
######    } else if (verbosity > 0) message (" : no improvement to local pscore")
######    if (improvements == maxImprovements) break()
######  } # for
######  if (verbosity >= 0)
######    message ("\nCompleted sectorial rearrangements.\n")
######  attr(tree, "score") <- NULL
######  attr(tree, "hits") <- NULL
######  # Return:
######  tree
######}  # DoSectorial
######
#######' Sectorial Search
#######'
#######' \code{SectorialSearch} performs a sectorial search on a tree.
#######' 
#######' A sectorial search detaches a random part of the tree, performs rearrangments
#######'  on this subtree, then reattaches it to the main tree (Goloboff, 1999).
#######' The improvement to local \var{score} hopefully (but not necessarily) improves 
#######' the overall \var{score}.
#######' 
#######' \code{SectorialSearch} then uses this tree as a starting point for a series
#######'  of tree rearrangements, by default using NNI, SPR and TBR swappers.
#######' It returns the new tree, unless the starting tree had a better score, 
#######' in which case the starting tree is returned.
#######'
#######''
#######' @inheritParams TreeTools::Renumber
#######' @param dataset a dataset in the format required by `TreeScorer()`.
#######' @importParams TreeSearch
#########' @param Bootstrapper Function to perform bootstrapped rearrangements of tree. 
#########'                     First arguments will be an edgeList and a dataset, initialized using \code{InitializeData}
#########'                     Should return a rearranged edgeList.
#######' @param SectorialSwapper Function such as \code{\link{RootedNNISwap}} to use 
#######'                         to rearrange sector.
#######' @param sectIter stop sectorial search when this many rearrangements have been performed.
#######' @param sectHits stop sectorial search when this many rearrangements have 
#######'                 found the same best score.
#######' @param searchIter maximum rearrangements for subsequent search on whole tree.
#######' @param searchHits maximum times to hit best score before terminating whole-tree search.
#######' @template verbosityParam
#######' @template treeScorerDots
#######' 
#######' @return a rooted tree of class \code{phylo}.
#######' 
#######' @references \insertRef{Goloboff1999}{TreeSearch}
#######' 
#######' @author Martin R. Smith
#######' 
#######' @seealso \code{\link{TreeSearch}}
#######' @seealso \code{\link{MorphyRatchet}}
#######' 
#######' @examples
#######" data("Lobo", package="TreeTools')
#######' njtree <- TreeTools::NJTree(Lobo.phy)
#######'
#######' \dontrun{
#######' SectorialSearch(njtree, Lobo.phy, maxIter=20, EdgeSwapper=NNISwap,
#######'  ratchIter=1, maxIter=50, largest.sector=7)} # Will be time-consuming }
#######' 
#######' @keywords  tree 
#######' @export
######SectorialSearch <- function (tree, dataset, 
######                             InitializeData = PhyDat2Morphy,
######                             CleanUpData    = UnloadMorphy,
######                             TreeScorer     = MorphyLength,
######                             Bootstrapper   = MorphyBootstrap,
######                             swappers = list(TBRSwap, SPRSwap, NNISwap),
######                             SectorialSwapper = swappers[[length(swappers)]],
######                             returnAll=FALSE, stopAtScore=NULL,
######                             sectIter=400L, sectHits=20L,
######                             searchIter=sectIter * 5L, searchHits=sectHits * 2L,
######                             verbosity=1L, ...) {
######  epsilon <- 1e-08
######  hits <- 0L
######  # initialize tree and data
######  if (dim(tree[["edge"]])[1] != 2 * tree[["Nnode"]]) stop("tree must be bifurcating; try rooting with ape::root")
######  tree <- RenumberTips(tree, names(dataset))
######  edgeList <- tree[["edge"]]
######  edgeList <- RenumberEdges(edgeList[, 1], edgeList[, 2])
######  
######  initializedData <- InitializeData(dataset)
######  on.exit(initializedData <- CleanUpData(initializedData))
######  
######  bestScore <- if (is.null(attr(tree, "score"))) {
######    TreeScorer(edgeList[[1]], edgeList[[2]], initializedData, ...)
######  } else {
######    attr(tree, "score")
######  }
######  
######  if (verbosity > 0L) message("\n* Beginning Sectorial Search, with initial score", bestScore)
######  if (!is.null(stopAtScore) && bestScore < stopAtScore + epsilon) return(tree)
######  
######  
######  # Rearrange a sector of the tree:
######  sect <- MorphySectorial(edgeList[[1]], edgelist[[2]], morphyObj, verbosity=verbosity-1, ratchIter=30, 
######    maxIter=maxIter, maxHits=15, smallest.sector=6, 
######    largest.sector=length(edgeList[[1]]) * 0.25, rearrangements=rearrangements)
######  
######  if (class(swappers) == "function") swappers <- list(swappers)
######  # Now use sectorially rearranged tree as starting point for conventional search
######  edgeList <- EdgeListSearch(edgeList, initializedData, TreeScorer=TreeScorer,
######                             EdgeSwapper=swappers, maxIter = maxIter, 
######                             maxHits = maxHits, verbosity = verbosity - 1L)
######  
######  if (edgeList[[3]] <= bestScore) {
######    sect[["edge"]] <- cbind(edgeList[[1]], edgeList[[2]])
######    attr(sect, "score") <- edgeList[[3]]
######    attr(sect, "hits") <- edgeList[[4]]
######    # Return:
######    sect
######  } else {
######    # Return:
######    tree
######  } 
######}
######
