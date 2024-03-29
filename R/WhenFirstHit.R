#' When was a tree topology first hit?
#' 
#' Reports when each tree in a list was first found by tree search.
#' This information is read from the `firstHit` attribute if present.
#' If not, trees are taken to be listed in the order in which they were found,
#' and named according to the search iteration in which they were first hit - 
#' the situation when trees found by [`MaximizeParsimony()`] are saved to file.
#' 
#' @param trees A list of trees, or a `multiPhylo` object.
#' @return `trees`, with a `firstHit` attribute listing the number of trees hit
#' for the first time in each search iteration.
#' @template MRS
#' 
#' @examples
#' library("TreeTools", quietly = TRUE)
#' trees <- list(
#'    seed_00 = as.phylo(1, 8),
#'    ratch1_01 = as.phylo(2, 8),
#'    ratch1_02 = as.phylo(3, 8),
#'    ratch4_44 = as.phylo(4, 8),
#'    final_99 = as.phylo(5, 8)
#' )
#' attr(WhenFirstHit(trees), "firstHit")
#' @family utility functions
#' @seealso
#' - [`MaximizeParsimony()`]
#' @export
WhenFirstHit <- function(trees) {
  if (is.null(attr(trees, "firstHit"))) {
    treeNames <- names(trees)
    pattern <- "(seed|start|ratch\\d+|final)_\\d+"
    if (length(grep(pattern, treeNames, perl = TRUE)) == length(trees)) {
      whenHit <- gsub(pattern, "\\1", treeNames, perl = TRUE)
      
      attr(trees, "firstHit") <- table(whenHit)[unique(whenHit)]
    }
  }
  
  # Return:
  trees
}
