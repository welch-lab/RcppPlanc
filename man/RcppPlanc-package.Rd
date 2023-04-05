\name{RcppPlanc-package}
\alias{RcppPlanc-package}
\alias{RcppPlanc}
\docType{package}
\title{
\packageTitle{RcppPlanc}
}
\description{
\packageDescription{RcppPlanc}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{RcppPlanc}
\packageIndices{RcppPlanc}
~~ An overview of how to use the package, including the most ~~
~~ important functions ~~
}
\author{
\packageAuthor{RcppPlanc}

Maintainer: \packageMaintainer{RcppPlanc}
}
\references{
~~ Literature or other references for background information ~~
}

\keyword{ package }

\examples{
  if (require("Matrix")) {
  x <- Matrix(runif(n=20, min=0, max=2), 4,5, sparse = TRUE)
  rcppplanc_bppnmf(t(x), 10, 10)
  }
}
