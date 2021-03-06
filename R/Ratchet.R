#' Parsimony Ratchet
#'
#' \code{Ratchet} uses the parsimony ratchet (Nixon 1999) to search for a more parsimonious tree.
#'
#' @template treeParam 
#' @param dataset a dataset in the format required by `TreeScorer()`.
#' @template InitializeDataParam
#' @template CleanUpDataParam
#' @template treeScorerParam
#' @param Bootstrapper Function to perform bootstrapped rearrangements of tree.
#' First arguments will be an `edgeList` and a dataset, initialized using 
#' `InitializeData()`. Should return a rearranged `edgeList`.
#' @template swappersParam
#' @param BootstrapSwapper Function such as \code{\link{RootedNNISwap}} to use 
#' to rearrange trees within `Bootstrapper()`.
#' @param returnAll Set to \code{TRUE} to report all MPTs encountered during the
#'  search, perhaps to analyse consensus.
#' @param ratchIter Stop when this many ratchet iterations have been performed.
#' @param ratchHits Stop when this many ratchet iterations have found the same
#' best score.
#' @param searchIter Integer specifying maximum rearrangements to perform on each bootstrap or 
#' ratchet iteration.  
#' To override this value for a single swapper function, set e.g. 
#' `attr(SwapperFunction, 'searchIter') <- 99`
#' @param searchHits Integer specifying maximum times to hit best score before terminating a tree
#' search within a ratchet iteration.
#' To override this value for a single swapper function, set e.g. 
#' `attr(SwapperFunction, 'searchHits') <- 99`
#' @param bootstrapIter Integer specifying maximum rearrangements to perform on each bootstrap
#'  iteration (default: `searchIter`).
#' @param bootstrapHits Integer specifying maximum times to hit best score on each bootstrap 
#' iteration (default: `searchHits`).
#' @template stopAtScoreParam
#' @template stopAtPeakParam
#' @template stopAtPlateauParam
#' @template verbosityParam
#' @param suboptimal retain trees that are suboptimal by this score. Defaults to 1e-08 to counter rounding errors.
#' @template treeScorerDots
#' 
#' @return This function returns a tree modified by parsimony ratchet iterations.
#'
#' @references 
#' 
#' - \insertRef{Nixon1999}{TreeSearch}
#' 
#' - \insertRef{SmithTern}{TreeSearch}
#'
#' @examples
#' data('Lobo', package='TreeTools')
#' njtree <- TreeTools::NJTree(Lobo.phy)
#' # Increase value of ratchIter and searchHits to do a proper search
#' quickResult <- Ratchet(njtree, Lobo.phy, ratchIter = 2, searchHits = 3)
#' 
#' # Plot result (legibly)
#' oldPar <- par(mar = rep(0, 4), cex = 0.75)
#' plot(quickResult)
#' par(oldPar)
#' 
#' # IW search is currently much slower:
#' quickIWResult <- IWRatchet(quickResult, Lobo.phy, concavity = 2.5,
#'                            ratchIter = 1, searchIter = 25, searchHits = 2,
#'                            swappers = RootedTBRSwap, verbosity = 5)
#'  
#' @author Martin R. Smith
#' 
#' @seealso
#' - [`TreeSearch()`]
#' 
#' - Adapted from \code{\link[phangorn:parsimony]{pratchet}} in the 
#' \pkg{phangorn} package.
#' 
#' @keywords tree
#' @importFrom TreeTools RenumberEdges RenumberTips
#' @export
Ratchet <- function (tree, dataset, 
                     InitializeData = PhyDat2Morphy,
                     CleanUpData    = UnloadMorphy,
                     TreeScorer     = MorphyLength,
                     Bootstrapper   = MorphyBootstrap,
                     swappers = list(TBRSwap, SPRSwap, NNISwap),
                     BootstrapSwapper = if (class(swappers) == 'list')
                      swappers[[length(swappers)]] else swappers,
                     returnAll=FALSE, stopAtScore=NULL,
                     stopAtPeak=FALSE, stopAtPlateau=0L, 
                     ratchIter=100, ratchHits=10,
                     searchIter=4000, searchHits=42,
                     bootstrapIter=searchIter, bootstrapHits=searchHits,
                     verbosity=1L, 
                     suboptimal=1e-08, ...) {
  epsilon <- 1e-08
  hits <- 0L
  # initialize tree and data
  if (dim(tree$edge)[1] != 2 * tree$Nnode) stop("tree must be bifurcating; try rooting with ape::root")
  tree <- RenumberTips(tree, names(dataset))
  edgeList <- tree$edge
  edgeList <- RenumberEdges(edgeList[, 1], edgeList[, 2])

  initializedData <- InitializeData(dataset)
  on.exit(initializedData <- CleanUpData(initializedData))
  
  bestScore <- if (is.null(attr(tree, 'score'))) {
    TreeScorer(edgeList[[1]], edgeList[[2]], initializedData, ...)
  } else {
    attr(tree, 'score')
  }
  if (verbosity > 0L) {
    message("* Beginning Parsimony Ratchet, with initial score ", bestScore,
            if (!is.null(stopAtScore)) "; will stop at score ", stopAtScore)
  }
  if (!is.null(stopAtScore) && bestScore < stopAtScore + epsilon) {
    if (verbosity > 1L) {
      message("*** Target score of ", stopAtScore, " met.")
    }
    return(tree)
  }
  if (class(swappers) == 'function') swappers <- list(swappers)
  
  if (returnAll) {
    nullForest <- vector('list', ratchIter)
    forest <- nullForest
    forestScores <- rep(NA, ratchIter)
  }

  iterationsWithBestScore <- 0
  BREAK <- FALSE
  for (i in 1:ratchIter) {
    if (verbosity > 1L) {
      message("\n* Ratchet iteration ", i, '.')
      if (verbosity > 2L) {
        message(" - Generating new candidate tree by bootstrapping dataset.")
      }
    }
    candidate <- Bootstrapper(edgeList, initializedData, maxIter=bootstrapIter,
                              maxHits=bootstrapHits, verbosity=verbosity-2L,
                              EdgeSwapper=BootstrapSwapper, stopAtPeak=stopAtPeak,
                              stopAtPlateau=stopAtPlateau, ...)
    candScore <- 1e+08
    
    if (verbosity > 2L) message(" - Rearranging from new candidate tree:")
    for (EdgeSwapper in swappers) {
      at <- attributes(EdgeSwapper)
      Argument <- function (arg) if (!is.null(at[[arg]])) at[[arg]] else get(arg)
      candidate <- EdgeListSearch(candidate, dataset = initializedData,
                                  TreeScorer = TreeScorer, 
                                  EdgeSwapper = EdgeSwapper, 
                                  maxIter = Argument("searchIter"),
                                  stopAtScore = Argument("stopAtScore"), 
                                  stopAtPeak = Argument("stopAtPeak"),
                                  stopAtPlateau = Argument("stopAtPlateau"),
                                  maxHits = Argument("searchHits"),
                                  verbosity = verbosity - 2L, ...)
      candScore <- candidate[[3]]
      if (!is.null(stopAtScore) && candScore < stopAtScore + epsilon) {
        BREAK <- TRUE
        if (verbosity > 1L) {
          message("  * Target score ", stopAtScore, 
                  " met; terminating tree search.")
        }
        bestScore <- candScore
        break
      }
    }
    if (BREAK) {
      break
    }
    
    if (verbosity > 2L) message(" - Rearranged candidate tree scored ", candScore)
    if (returnAll && candScore < (bestScore + suboptimal)) { # Worth saving this tree in forest
      forest[[i]] <- candidate
      forestScores[i] <- candScore
    }
    if ((candScore + epsilon) < bestScore) {
      # New 'best' tree
      edgeList <- candidate
      bestScore <- candScore
      iterationsWithBestScore <- 1L
    } else if (bestScore + epsilon > candScore) { # i.e. best == cand, allowing for floating point error
      iterationsWithBestScore <- iterationsWithBestScore + 1L
      edgeList <- candidate
    }
    if (verbosity > 1L) {
      message("* Best score after ", i, "/", ratchIter, 
              " ratchet iterations: ", signif(bestScore), " (hit ", 
              iterationsWithBestScore, "/", ratchHits, ")\n")
    }
    if ((!is.null(stopAtScore) && bestScore < stopAtScore + epsilon) 
    || (iterationsWithBestScore >= ratchHits)) {
      break
    }
  } # end for

  if (verbosity > 0L) {
    message("Completed parsimony ratchet after ", i, " iterations with score ",
            bestScore, "\n")
  }
  
  if (returnAll) {
    keepers <- !is.na(forestScores) & forestScores < bestScore + suboptimal
    forestScores <- forestScores[keepers]
    forest <- forest[keepers]
    if (verbosity > 1L) {
      message("\n - Keeping ", sum(keepers), 
              " trees from iterations numbered:\n   ", which(keepers))
    }
    if (length(forest) > 1) {
      forest <- structure(lapply(forest, function (phy) {
        x <- tree
        x$edge <- cbind(phy[[1]], phy[[2]])
        attr(x, 'score') <- phy[[3]]
        # Return to lapply: 
        x}), class = 'multiPhylo')
      ret <- unique(forest)
      if (verbosity > 1L) {
        message(" - Removing duplicates leaves ", length(ret), " unique trees")
      }
      uniqueScores <- vapply(ret, attr, double(1), 'score')
    } else if (length(forest) == 1) {
      ret <- tree
      newEdge <- forest[[1]]
      ret$edge <- cbind(newEdge[[1]], newEdge[[2]])
      uniqueScores <- newEdge[[3]]
    } else {
      stop("\nNo trees!? Is suboptimal set to a sensible (positive) value?")
    }
    if (verbosity > 0L) {
      message('\nFound ', sum(uniqueScores == min(uniqueScores)),
              ' unique MPTs and ',
              length(ret) - sum(uniqueScores == min(uniqueScores)), 
              ' suboptimal trees.\n')
    }
    # Return:
    ret
  } else {
    tree$edge <- cbind(edgeList[[1]], edgeList[[2]])
    attr(tree, 'score') <- bestScore
    # Return:
    tree  
  }
}

#' @describeIn Ratchet Shortcut for Ratchet search under Profile Parsimony
#' @export
ProfileRatchet <- function (tree, dataset,
                            swappers = list(TBRSwap, SPRSwap, NNISwap),
                            BootstrapSwapper = if (class(swappers) == 'list')
                              swappers[[length(swappers)]] else swappers,
                            returnAll=FALSE, stopAtScore=NULL,
                            stopAtPeak=FALSE, stopAtPlateau=0L, 
                            ratchIter=100, ratchHits=10, 
                            searchIter=2000, searchHits=40,
                            bootstrapIter=searchIter, bootstrapHits=searchHits, verbosity=1L, 
                            suboptimal=1e-08, ...) {
  Ratchet(tree=tree, dataset=dataset,
          InitializeData=ProfileInitMorphy, CleanUpData=ProfileDestroyMorphy,
          TreeScorer=ProfileScoreMorphy, Bootstrapper=ProfileBootstrap,
          swappers=swappers, BootstrapSwapper=BootstrapSwapper,
          returnAll=returnAll, suboptimal=suboptimal, stopAtScore=stopAtScore,
          ratchIter=ratchIter, ratchHits=ratchHits,
          searchIter=searchIter, searchHits=searchHits,
          stopAtPeak=stopAtPeak, stopAtPlateau=stopAtPlateau, 
          bootstrapIter=searchIter, bootstrapHits=bootstrapHits, 
          verbosity=verbosity,  ...)
}

#' @describeIn Ratchet Shortcut for Ratchet search using implied weights
#' @template concavityParam
#' @export
IWRatchet <- function (tree, dataset, concavity = 10,
                            swappers = list(TBRSwap, SPRSwap, NNISwap),
                            BootstrapSwapper = if (class(swappers) == 'list')
                              swappers[[length(swappers)]] else swappers,
                            returnAll=FALSE, stopAtScore=NULL,
                            stopAtPeak=FALSE, stopAtPlateau=0L, 
                            ratchIter=100, ratchHits=10, searchIter=2000, searchHits=40,
                            bootstrapIter=searchIter, bootstrapHits=searchHits, verbosity=1L, 
                            suboptimal=1e-08, ...) {
  dataset <- PrepareDataIW(dataset)
  if (verbosity > 1L) {
    message("* Using implied weighting with concavity constant k = ", concavity)
  }
  
  Ratchet(tree=tree, dataset=dataset, 
          concavity=concavity, minLength=attr(dataset, 'min.length'), 
          InitializeData=IWInitMorphy, CleanUpData=IWDestroyMorphy,
          TreeScorer=IWScoreMorphy, Bootstrapper=IWBootstrap,
          swappers=swappers, BootstrapSwapper=BootstrapSwapper,
          returnAll=returnAll, suboptimal=suboptimal, stopAtScore=stopAtScore,
          ratchIter=ratchIter, ratchHits=ratchHits,
          stopAtPeak=stopAtPeak, stopAtPlateau=stopAtPlateau, 
          searchIter=searchIter, searchHits=searchHits,
          bootstrapIter=searchIter, bootstrapHits=bootstrapHits, 
          verbosity=verbosity, ...)
}

#' Unique trees (ignoring 'hits' attribute)
#' @author Martin R. Smith
#' @keywords internal
#' @export
UniqueExceptHits <- function (trees) {
  unique(lapply(trees, function(tree) {
    attr(tree, 'hits') <- NULL
    tree
  }))
}

#' @rdname Ratchet 
#' @return `MultiRatchet()` returns a list of optimal trees 
#' produced by `nSearch` 
#' ratchet searches, from which a consensus tree can be generated using 
#'  [`ape::consensus()`] or [`TreeTools::ConsensusWithout()`].
#' @param nSearch Number of Ratchet searches to conduct
#' (for `RatchetConsensus()`)
#' @export
MultiRatchet <- function (tree, dataset, ratchHits=10, 
                              searchIter=500, searchHits=20, verbosity=0L, 
                              swappers=list(RootedNNISwap), nSearch=10, 
                              stopAtScore=NULL, ...) {
  trees <- lapply(seq_len(nSearch), function (i) {
    if (verbosity > 1L) message("\nRatchet search ", i, '/', nSearch, ':')
    Ratchet(tree, dataset, ratchIter = 1, ratchHits = 0L, 
            searchIter = searchIter, searchHits = searchHits, 
            verbosity = verbosity, swappers = swappers, 
            stopAtScore = stopAtScore, ...)
  })
  scores <- vapply(trees, function (x) attr(x, 'score'), double(1))
  trees <- UniqueExceptHits(trees[scores == min(scores)])
  message("Found ", length(trees), ' unique trees from ', nSearch, ' searches.')
  
  # Return:
  structure(trees, class = 'multiPhylo')
}

#' @rdname Ratchet 
#' @return `IWMultiRatchet` returns a list of optimal trees produced by `nSearch` 
#'                     Ratchet searches, using implied weighting.
#' @export
IWMultiRatchet <- function (tree, dataset, ratchHits=10, concavity=4,
                              searchIter=500, searchHits=20, verbosity=0L, 
                              swappers=list(RootedNNISwap), nSearch=10, 
                              suboptimal=suboptimal,
                              stopAtScore=NULL, ...) {
  trees <- lapply(seq_len(nSearch), function (i) {
    if (verbosity > 1L) message("\nRatchet search ", i, '/', nSearch, ':')
    IWRatchet(tree, dataset, ratchIter = 1L, ratchHits = 0L, 
              concavity = concavity, 
              searchIter = searchIter, searchHits = searchHits, 
              verbosity = verbosity, swappers = swappers,
              stopAtScore = stopAtScore, ...)
  })
  scores <- vapply(trees, function (x) attr(x, 'score'), double(1))
  trees <- UniqueExceptHits(trees[scores == min(scores)])
  message("Found ", length(trees), ' unique trees from ', nSearch, ' searches.')
  
  # Return:
  structure(trees, class = 'multiPhylo')
}

#' @describeIn Ratchet deprecated alias for `MultiRatchet()`
#' @export
RatchetConsensus <- MultiRatchet
#' @describeIn Ratchet deprecated alias for `MultiRatchet()`
#' @export
IWRatchetConsensus <- IWMultiRatchet