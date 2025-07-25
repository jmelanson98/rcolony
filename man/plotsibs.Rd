% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotsibs.R
\name{plotsibs}
\alias{plotsibs}
\title{A function to plot sibships (full and half-sibs) in an x-y style plot.}
\usage{
plotsibs(colony.object, maintitle = "sibships", prob = 0.8, pairwise = FALSE)
}
\arguments{
\item{colony.object}{A list object derived from the colony data
(\code{\link{get.colony.data}}).}

\item{maintitle}{The plot's title.}

\item{prob}{Assignment probability threshold. The probability of sibship
must be greater than this value to be included in the plot.}

\item{pairwise}{TRUE/FALSE: should pairwise sibship estimates be used.
(TRUE= pairwise, FALSE=full likelihood)}
}
\value{
A plot showing full- or half-sibships.
}
\description{
This function plots the estimated sibships from Colony data.  Both the x and
y axis represent individuals in the population and points indicate a sibship
assignment between individuals x and y.
}
\details{
Full sibships are shown in the upper left part of the plot, and
half-sibships are shown in the lower right part. The data can be derived
from pairwise methods, or from the full likelihood method.
}
\examples{

data(testdata)

plotsibs(testdata,maintitle="Sibships (pairwise)",pairwise=TRUE)
plotsibs(testdata,maintitle="Sibships (full)",pairwise=FALSE)

}
\references{
Wang, J. (2004) Sibship reconstruction from genetic data with
typing errors. Genetics 166: 1963-1979.
}
\seealso{
\code{\link{get.colony.data}},\code{\link{run.colony}}
}
\author{
Owen R. Jones
}
\keyword{aplot}
