% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_length.R
\name{CharacterLength}
\alias{CharacterLength}
\alias{FitchSteps}
\alias{FastCharacterLength}
\title{Character length}
\usage{
CharacterLength(tree, dataset, compress = FALSE)

FitchSteps(tree, dataset)

FastCharacterLength(tree, dataset)
}
\arguments{
\item{tree}{A tree of class \code{\link[ape:read.tree]{phylo}}.}

\item{dataset}{A phylogenetic data matrix of \pkg{phangorn} class
\code{phyDat}, whose names correspond to the labels of any accompanying tree.
Perhaps load into R using \code{\link[TreeTools]{ReadAsPhyDat}}.
Additive (ordered) characters can be handled using
\code{\link[TreeTools]{Decompose}}.}

\item{compress}{Logical specifying whether to retain the compression of a
\code{phyDat} object or to return a vector specifying to each individual
character, decompressed using the dataset's \code{index} attribute.}
}
\value{
\code{CharacterLength()} returns a vector listing the contribution of each
character to tree score, according to the algorithm of
\insertCite{Brazeau2018;textual}{TreeTools}.
}
\description{
Homoplasy length of each character in a dataset on a specified tree.
}
\section{Functions}{
\itemize{
\item \code{FastCharacterLength()}: Do not perform checks.  Use with care: may cause
erroneous results or software crash if variables are in the incorrect format.

}}
\examples{
data("inapplicable.datasets")
dataset <- inapplicable.phyData[[12]]
tree <- TreeTools::NJTree(dataset)
CharacterLength(tree, dataset)
CharacterLength(tree, dataset, compress = TRUE)
}
\references{
\insertAllCited{}
}
\seealso{
Other tree scoring: 
\code{\link{IWScore}()},
\code{\link{LengthAdded}()},
\code{\link{MinimumLength}()},
\code{\link{MorphyTreeLength}()},
\code{\link{TaxonInfluence}()}
}
\author{
\href{https://smithlabdurham.github.io/}{Martin R. Smith}
(\href{mailto:martin.smith@durham.ac.uk}{martin.smith@durham.ac.uk})
}
\concept{tree scoring}
