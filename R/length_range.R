#' Minimum length
#' 
#' The smallest length that a character can obtain on any tree.
#' 
#' 
#' @param x An object of class `phyDat`;
#' or a string to be coerced to a `phyDat` object via 
#' [`TreeTools::StringToPhyDat()`];
#' or an integer vector listing the tokens that may be present at each 
#' tip along a single character, with each token represented as a binary digit;
#' e.g. a value of 11 ( = 2^0 + 2^1 + 2^3) means that
#' the tip may have tokens 0, 1 or 3.
#' 
#' Inapplicable tokens should be denoted with the integer `0` (not 2^0).
#' 
#' @template compressParam
#' 
#' @return `MinimumLength()` returns a vector of integers specifying the 
#' minimum number of steps that each character must contain.
#'
#' @examples
#' data("inapplicable.datasets")
#' myPhyDat <- inapplicable.phyData[[4]]
#' 
#' # load your own data with
#' # my.PhyDat <- as.phyDat(read.nexus.data("filepath"))
#' # or Windows users can select a file interactively using:
#' # my.PhyDat <- as.phyDat(read.nexus.data(choose.files()))
#' 
#' class(myPhyDat) # phyDat object
#' 
#' # Minimum length of each character in turn
#' MinimumLength(myPhyDat)
#' 
#' # Collapse duplicate characters, per phyDat compression
#' MinimumLength(myPhyDat, compress = TRUE)
#' 
#' # Calculate length of a single character from its textual representation
#' MinimumLength("-{-1}{-2}{-3}2233")
#' @template MRS
#' @family tree scoring
#' @export
MinimumLength <- function (x, compress = FALSE) UseMethod("MinimumLength")

#' @rdname MinimumLength
#' @export
MinimumLength.phyDat <- function (x, compress = FALSE) {
  
  at <- attributes(x)
  nLevel <- length(at$level)
  nChar <- at$nr
  nTip <- length(x)
  cont <- at$contrast
  if (is.null(colnames(cont))) {
    colnames(cont) <- as.character(at$levels)
  }
  
  inappLevel <- at$levels == "-"
  powersOf2 <- 2L ^ (seq_len(nLevel - sum(inappLevel)) - 1L)
  
  # Treat {-, 1} as {1}
  unlisted <- unlist(x, use.names = FALSE)
  tmp <- as.integer(cont[, colnames(cont) != "-"] %*% powersOf2)
  ambigIsApp <- matrix(tmp[unlisted], nChar, nTip)
  
  if (any(inappLevel)) {
    # Treat {-, 1} as {-}
    tmp[cont[, "-"] == 1] <- 0
    ambigIsInapp <- matrix(tmp[unlisted], nChar, nTip)
    
    inappCount <- rowSums(matrix(unlisted %in% which(at$allLevels == "-"),
                                 nChar, nTip))
    binaryMatrix <- ambigIsApp
    binaryMatrix[inappCount > 1, ] <- ambigIsInapp[inappCount > 1, ]
  } else {
    binaryMatrix <- ambigIsApp
  }
  
  ret <- apply(binaryMatrix, 1, MinimumLength)
  
  # Return:
  if (compress) {
    ret
  } else {
    ret[attr(x, "index")]
  }
}

#' @rdname MinimumLength
#' @export
MinimumLength.numeric <- function (x, compress = NA) {
  
  uniqueStates <- unique(x[x > 0])
  if (length(uniqueStates) < 2) return (0)
  tokens <- vapply(uniqueStates, intToBits, raw(32)) != 00
  tokens <- tokens[apply(tokens, 1, any), ]
  
  lastDim <- dim(tokens)
  tokensUsed <- 0
  
  repeat {
    tokens <- tokens[!duplicated(tokens), , drop = FALSE]
    unambiguous <- colSums(tokens) == 1
    tokenNecessary <- apply(tokens[, unambiguous, drop = FALSE], 1, any)
    statesRemaining <- !unambiguous
    statesRemaining[statesRemaining] <- colSums(
      tokens[tokenNecessary, statesRemaining, drop = FALSE]) == 0
    tokensUsed <- tokensUsed + sum(tokenNecessary)
    
    if (!any(statesRemaining)) {
      # Return:
      return (tokensUsed - 1L)
    }
    
    tokens <- tokens[!tokenNecessary, statesRemaining, drop = FALSE]
    if (identical(dim(tokens), lastDim)) {
      occurrences <- rowSums(tokens)
      unnecessary <- occurrences == 1
      if (any(unnecessary)) {
        tokens <- tokens[!unnecessary, , drop = FALSE]
      } else {
        squish <- which.max(occurrences)
        tokensUsed <- tokensUsed + 1L
        tokens <- tokens[, tokens[!squish], drop = FALSE]
      }
    }
    lastDim <- dim(tokens)
  }
}

#' @rdname MinimumLength
#' @importFrom TreeTools NexusTokens StringToPhyDat
#' @export
MinimumLength.character <- function (x, compress = TRUE) {
  nTip <- length(NexusTokens(x[1]))
  vapply(x, function (x) MinimumLength(StringToPhyDat(x, nTip)),
         1, USE.NAMES = FALSE)
}


#' @rdname MinimumLength
MinimumSteps <- function(x) {
  .Deprecated("MinimumLength",
              msg = "Renamed to `MinimumLength()` and recoded to better support inapplicable tokens")
  MinimumLength(x, compress = TRUE)
}
