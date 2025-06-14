% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Jackknife.R
\name{Jackknife}
\alias{Jackknife}
\title{Jackknife resampling}
\usage{
Jackknife(
  tree,
  dataset,
  resampleFreq = 2/3,
  InitializeData = PhyDat2Morphy,
  CleanUpData = UnloadMorphy,
  TreeScorer = MorphyLength,
  EdgeSwapper = TBRSwap,
  jackIter = 5000L,
  searchIter = 4000L,
  searchHits = 42L,
  verbosity = 1L,
  ...
)
}
\arguments{
\item{tree}{A fully-resolved starting tree in \code{\link[ape]{phylo}} format,
with the desired outgroup.
Edge lengths are not supported and will be removed.}

\item{dataset}{a dataset in the format required by \code{TreeScorer()}.}

\item{resampleFreq}{Double between 0 and 1 stating proportion of characters
to resample.}

\item{InitializeData}{Function that sets up data object to prepare for tree search.
The function will be passed the \code{dataset} parameter.
Its return value will be passed to \code{TreeScorer()} and \code{CleanUpData()}.}

\item{CleanUpData}{Function to destroy data object on function exit.
The function will be passed the value returned by \code{InitializeData()}.}

\item{TreeScorer}{function to score a given tree.
The function will be passed three parameters, corresponding to the
\code{parent} and \code{child} entries of a tree's edge list, and a dataset.}

\item{EdgeSwapper}{a function that rearranges a parent and child vector,
and returns a list with modified vectors; for example \code{\link[=SPRSwap]{SPRSwap()}}.}

\item{jackIter}{Integer specifying number of jackknife iterations to conduct.}

\item{searchIter}{Integer specifying maximum rearrangements to perform on
each bootstrap or ratchet iteration.
To override this value for a single swapper function, set e.g.
\code{attr(SwapperFunction, "searchIter") <- 99}}

\item{searchHits}{Integer specifying maximum times to hit best score before
terminating a tree search within a ratchet iteration.
To override this value for a single swapper function, set e.g.
\code{attr(SwapperFunction, "searchHits") <- 99}}

\item{verbosity}{Numeric specifying level of detail to display in console:
larger numbers provide more verbose feedback to the user.}

\item{...}{further parameters to send to \code{TreeScorer()}}
}
\value{
\code{Jackknife()} returns a list of trees recovered after jackknife
iterations.
}
\description{
Resample trees using Jackknife resampling, i.e. removing a subset of
characters.
}
\details{
The function assumes
that \code{InitializeData()} will return a morphy object; if this doesn't hold
for you, post a \href{https://github.com/ms609/TreeSearch/issues/new/}{GitHub issue}
or e-mail the maintainer.
}
\seealso{
\itemize{
\item \code{\link[=Resample]{Resample()}}: Jackknife resampling for non-custom searches performed
using \code{MaximizeParsimony()}.
\item \code{\link[=JackLabels]{JackLabels()}}: Label nodes of a tree with jackknife supports.
}

Other split support functions: 
\code{\link{JackLabels}()},
\code{\link{MaximizeParsimony}()},
\code{\link{SiteConcordance}}

Other custom search functions: 
\code{\link{EdgeListSearch}()},
\code{\link{MorphyBootstrap}()},
\code{\link{SuccessiveApproximations}()}
}
\author{
\href{https://smithlabdurham.github.io/}{Martin R. Smith}
(\href{mailto:martin.smith@durham.ac.uk}{martin.smith@durham.ac.uk})
}
\concept{custom search functions}
\concept{split support functions}
