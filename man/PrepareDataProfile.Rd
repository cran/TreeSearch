% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_manipulation.R
\encoding{UTF-8}
\name{PrepareDataProfile}
\alias{PrepareDataProfile}
\alias{PrepareDataIW}
\title{Prepare data for Profile Parsimony}
\usage{
PrepareDataProfile(dataset)

PrepareDataIW(dataset)
}
\arguments{
\item{dataset}{dataset of class \code{phyDat}}
}
\value{
An object of class \code{phyDat}, with additional attributes.
\code{PrepareDataProfile} adds the attributes:
\itemize{
\item \code{info.amounts}: details the information represented by each
character when subject to N additional steps.
\item \code{informative}: logical specifying which characters contain any
phylogenetic information.
\item \code{bootstrap}: The character vector
\code{c("info.amounts", "split.sizes")}, indicating attributes to sample
when bootstrapping the dataset (e.g. in Ratchet searches).
}

\code{PrepareDataIW} adds the attribute:
\itemize{
\item \code{min.length}: The minimum number of steps that must be present in each
transformation series.
}
}
\description{
Calculates profiles for each character in a dataset.  Will also simplify
characters, with a warning, where they are too complex for the present
implementation of profile parsimony:
\itemize{
\item inapplicable tokens will be replaced with the ambiguous token
(i.e. \code{-} \ifelse{html}{\out{&rarr;}}{\eqn{\rightarrow}{-->}} \verb{?});
\item Ambiguous tokens will be treated as fully ambiguous
(i.e. \code{{02}} \ifelse{html}{\out{&rarr;}}{\eqn{\rightarrow}{-->}} \verb{?})
\item Where more than two states are informative (i.e. unambiguously present in
more than one taxon), states beyond the two most informative will be
ignored.
}
}
\section{Functions}{
\itemize{
\item \code{PrepareDataIW()}: Prepare data for implied weighting

}}
\examples{
data("congreveLamsdellMatrices")
dataset <- congreveLamsdellMatrices[[42]]
PrepareDataProfile(dataset)
}
\seealso{
Other profile parsimony functions: 
\code{\link{Carter1}()},
\code{\link{StepInformation}()},
\code{\link{WithOneExtraStep}()},
\code{\link{profiles}}
}
\author{
Martin R. Smith; written with reference to
\code{phangorn:::prepareDataFitch()}
}
\concept{profile parsimony functions}
