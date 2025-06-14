% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TBR.R
\name{TBR}
\alias{TBR}
\alias{TBRMoves}
\alias{TBRMoves.phylo}
\alias{TBRMoves.matrix}
\alias{TBRSwap}
\alias{RootedTBR}
\alias{RootedTBRSwap}
\title{Tree bisection and reconnection (TBR)}
\usage{
TBR(tree, edgeToBreak = NULL, mergeEdges = NULL)

TBRMoves(tree, edgeToBreak = integer(0))

\method{TBRMoves}{phylo}(tree, edgeToBreak = integer(0))

\method{TBRMoves}{matrix}(tree, edgeToBreak = integer(0))

TBRSwap(
  parent,
  child,
  nEdge = length(parent),
  edgeToBreak = NULL,
  mergeEdges = NULL
)

RootedTBR(tree, edgeToBreak = NULL, mergeEdges = NULL)

RootedTBRSwap(
  parent,
  child,
  nEdge = length(parent),
  edgeToBreak = NULL,
  mergeEdges = NULL
)
}
\arguments{
\item{tree}{A bifurcating tree of class \code{\link[ape]{phylo}}, with all nodes resolved;}

\item{edgeToBreak}{(optional) integer specifying the index of an edge to bisect/prune,
generated randomly if not specified.
Alternatively, set to \code{-1} to return a complete list
of all trees one step from the input tree.}

\item{mergeEdges}{(optional) vector of length 1 or 2, listing edge(s) to be joined:
In SPR, this is where the pruned subtree will be reconnected.
In TBR, these edges will be reconnected (so must be on opposite
sides of \code{edgeToBreak}); if only a single edge is specified,
the second will be chosen at random}

\item{parent}{Integer vector corresponding to the first column of the edge
matrix of a tree of class \code{\link[ape]{phylo}}, i.e.
\code{tree[["edge"]][, 1]}}

\item{child}{Integer vector corresponding to the second column of the edge
matrix of a tree of class \code{\link[ape]{phylo}}, i.e.
\code{tree[["edge"]][, 2]}.}

\item{nEdge}{(optional) Number of edges.}
}
\value{
\code{TBR()} returns a tree in \code{phyDat} format that has undergone one
\acronym{TBR} iteration.

\code{TBRMoves()} returns a \code{multiPhylo} object listing all trees one
\acronym{TBR} move away from \code{tree}, with edges and nodes in preorder,
rooted on the first-labelled tip.

\code{TBRSwap()} returns a list containing two elements corresponding
to the rearranged \code{parent} and \code{child} parameters.
}
\description{
\code{TBR} performs a single random \acronym{TBR} iteration.
}
\details{
Branch lengths are not (yet) supported.

All nodes in a tree must be bifurcating; \link[ape:collapse.singles]{ape::collapse.singles} and
\link[ape:multi2di]{ape::multi2di} may help.
}
\section{Functions}{
\itemize{
\item \code{TBRSwap()}: faster version that takes and returns parent and child
parameters

\item \code{RootedTBR()}: Perform \acronym{TBR} rearrangement, retaining position of root

\item \code{RootedTBRSwap()}: faster version that takes and returns parent and child parameters

}}
\examples{
library("ape")
tree <- rtree(20, br=NULL)
TBR(tree)
}
\references{
The \acronym{TBR} algorithm is summarized in
\insertRef{Felsenstein2004}{TreeSearch}
}
\seealso{
\code{\link[=RootedTBR]{RootedTBR()}}: useful when the position of the root node should be retained.

Other tree rearrangement functions: 
\code{\link{NNI}()},
\code{\link{SPR}()}
}
\author{
\href{https://smithlabdurham.github.io/}{Martin R. Smith}
(\href{mailto:martin.smith@durham.ac.uk}{martin.smith@durham.ac.uk})
}
\concept{tree rearrangement functions}
