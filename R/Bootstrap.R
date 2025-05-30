#' @inheritParams TreeSearch
#' @inheritParams EdgeListSearch
#' @inheritParams MorphyTreeLength
#' @param maxIter Numeric specifying maximum number of iterations to perform in
#' tree search.
#' @param maxHits Numeric specifying maximum number of hits to accomplish in
#' tree search.
#' @param stopAtScore stop search as soon as this score is hit or beaten.
#' @param \dots further parameters to send to `TreeScorer()`
#'
#' @return `MorphyBootstrap()` returns a tree that is optimal under a random
#' sampling of the original characters.
#' 
#' @rdname Ratchet
#' @export
MorphyBootstrap <- function (edgeList, morphyObj, EdgeSwapper = NNISwap, 
                             maxIter, maxHits, verbosity = 1L, 
                             stopAtPeak = FALSE, stopAtPlateau=0L, ...) {
  startWeights <- MorphyWeights(morphyObj)["exact", ]
  eachChar <- seq_along(startWeights)
  deindexedChars <- rep.int(eachChar, startWeights)
  resampling <- tabulate(sample(deindexedChars, replace = TRUE),
                         length(startWeights))
  errors <- vapply(eachChar, function (i) 
            mpl_set_charac_weight(i, resampling[i], morphyObj), integer(1))
  
  if (any(errors)) { # nocov start
    stop ("Error resampling morphy object: ",
          mpl_translate_error(unique(errors[errors < 0L])))
  }
  if (mpl_apply_tipdata(morphyObj) -> error) {
    stop("Error applying tip data: ", mpl_translate_error(error))
  } # nocov end
  
  res <- EdgeListSearch(edgeList[1:2], morphyObj, EdgeSwapper = EdgeSwapper,
                        maxIter = maxIter, maxHits = maxHits,
                        stopAtPeak = stopAtPeak, stopAtPlateau = stopAtPlateau,
                        verbosity = verbosity - 1L, ...)
  errors <- vapply(eachChar, function (i) 
         mpl_set_charac_weight(i, startWeights[i], morphyObj), integer(1))
  if (any(errors)) { # nocov start
    stop ("Error resampling morphy object: ",
          mpl_translate_error(unique(errors[errors < 0L])))
  }
  if (mpl_apply_tipdata(morphyObj) -> error) {
    stop("Error applying tip data: ", mpl_translate_error(error))
  } # nocov end
  
  # Return:
  res[1:2]
}
