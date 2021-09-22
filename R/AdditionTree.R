#' Addition tree
#' 
#' Generates a starting tree by adding each taxon in turn to the most
#' parsimonious location.
#' 
#' @inheritParams MaximizeParsimony
#' @param sequence Character or numeric vector listing sequence in which to add
#' taxa. Randomized if not provided.
#' @examples 
#' data('Lobo', package = 'TreeTools')
#' AdditionTree(Lobo.phy, concavity = 10)
#' @template MRS
#' @return `AdditionTree()` returns a tree of class `phylo`, rooted on
#' `sequence[1]`.
#' @importFrom TreeTools AddUnconstrained AddTipEverywhere PectinateTree
#' @importFrom cli cli_progress_bar cli_progress_update
#' @family tree generation functions
#' @export
AdditionTree <- function (dataset, concavity = Inf, constraint, sequence) {
  # Initialize missing parameters
  taxa <- names(dataset)
  if (missing(sequence)) {
    sequence <- taxa[1]
  } else if (is.numeric(sequence)) {
    sequence <- taxa[sequence]
  }
  nTaxa <- length(taxa)
  if (length(taxa) < 4) {
    return(PectinateTree(taxa))
  }
  unlisted <- setdiff(taxa, sequence)
  if (length(unlisted) > 0) {
    sequence <- c(sequence, sample(unlisted))
  }
  if (!missing(constraint)) {
    constraint <- AddUnconstrained(constraint, taxa)
  }
  # PrepareDataXXX attributes only valid for full dataset
  attr(dataset, 'info.amounts') <- NULL
  attr(dataset, 'min.length') <- NULL
  attr(dataset, 'informative') <- NULL
  
  # Starting tree, rooted on first element in sequence
  tree <- PectinateTree(sequence[1:3])
  
  cli_progress_bar('Addition tree', total = sum(2 * (4:nTaxa) - 5))
  for (addition in sequence[4:nTaxa]) {
    candidates <- AddTipEverywhere(tree, addition)
    nCands <- length(candidates)
    
    theseTaxa <- candidates[[1]]$tip.label
    theseData <- .Recompress(dataset[theseTaxa])
    if (is.finite(concavity)) {
      theseData <- PrepareDataIW(theseData)
    } else if (is.character(concavity)) {
      theseData <- suppressMessages(PrepareDataProfile(theseData))
    }
    
    if (!missing(constraint)) {
      thisConstr <- constraint[theseTaxa]
      morphyConstr <- PhyDat2Morphy(thisConstr)
      # Calculate constraint minimum score
      constraintLength <- sum(MinimumLength(thisConstr, compress = TRUE) *
                              attr(thisConstr, 'weight'))
      
      .Forbidden <- function (edges) {
        preorder_morphy(edges, morphyConstr) != constraintLength
      }
      
    
      candidates <- candidates[!vapply(lapply(candidates, `[[`, 'edge'),
                                       .Forbidden, logical(1))]
      UnloadMorphy(morphyConstr)
    }
    
    # Score remaining candidates
    if (length(theseData)) {
      scores <- TreeLength(candidates, theseData, concavity)
      minScore <- which.min(scores)
      nMin <- length(minScore)
      if (nMin > 1) {
        minScore <- minScore[sample.int(nMin, 1)]
      }
      tree <- candidates[[minScore]]
    } else {
      tree <- sample(candidates, 1)[[1]]
    }
    cli_progress_update(nCands)
  }
  tree
}

.Recompress <- function (dataset) {
  MatrixToPhyDat(PhyDatToMatrix(dataset))
}