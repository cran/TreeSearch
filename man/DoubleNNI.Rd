% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NNI.R
\name{DoubleNNI}
\alias{DoubleNNI}
\title{Double NNI}
\usage{
DoubleNNI(parent, child, edgeToBreak)
}
\arguments{
\item{parent}{Integer vector corresponding to the first column of the edge
matrix of a tree of class \code{\link[ape]{phylo}}, i.e.
\code{tree[["edge"]][, 1]}}

\item{child}{Integer vector corresponding to the second column of the edge
matrix of a tree of class \code{\link[ape]{phylo}}, i.e.
\code{tree[["edge"]][, 2]}.}

\item{edgeToBreak}{In (\code{Rooted}/\code{Double})\code{NNI()}, an optional integer
specifying the index of an edge to rearrange, generated randomly if not
specified.
If \code{-1}, a complete list of all trees one step from the input tree
will be returned.
In \code{cNNI()}, an integer from zero to \code{nEdge(tree) - nTip(tree) - 2},
specifying which internal edge to break.}
}
\value{
the \code{tree[["edge"]]} parameter of the two trees consistent with the specified rearrangement
}
\description{
Returns the edge parameter of the two trees consistent with the speficied \acronym{NNI} rearrangement
}
\author{
Martin R. Smith
}
\keyword{internal}
