% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SPR.R
\name{cSPR}
\alias{cSPR}
\title{\code{cSPR()} expects a tree rooted on a single tip.}
\usage{
cSPR(tree, whichMove = NULL)
}
\arguments{
\item{tree}{A tree of class \code{\link[ape:read.tree]{phylo}}.}

\item{whichMove}{Integer specifying which SPR move index to perform.}
}
\description{
\code{cSPR()} expects a tree rooted on a single tip.
}
\examples{
tree <- TreeTools::BalancedTree(8)

# Tree must be rooted on leaf
tree <- TreeTools::RootTree(tree, 1)

# Random rearrangement
cSPR(tree)

# Specific rearrangement
cSPR(tree, 9)
}
\author{
\href{https://smithlabdurham.github.io/}{Martin R. Smith}
(\href{mailto:martin.smith@durham.ac.uk}{martin.smith@durham.ac.uk})
}
