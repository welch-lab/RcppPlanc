\name{RcppPlanc-Functions}
\alias{rcppplanc_bppnmf}
\title{Set of functions to call NMF routines from planc.}
\description{
  These functions return factorized W and H matrices
  for a given input matrix X.
}
\usage{
rcppplanc_bppnmf(x, k, niter)
}
\arguments{
  \item{x}{Input matrix for factorization}
  \item{k}{Factor matrix rank}
  \item{niter}{Maximum number of NMF iterations}
}
\value{
  \code{rcppplanc_bppnmf()} calls planc NMF with the ANLS-BPP algorithm.
}
\details{
  This is a wrapper of planc, originally written in C++. It can be extended
  as any other Armadillo library.
}
\references{
  See the documentation and the c++ repository for more information.
}
\examples{
  if (require("Matrix")) {
  x <- Matrix(runif(n=20, min=0, max=2), 4,5, sparse = TRUE)
  rcppplanc_bppnmf(t(x), 10, 10)
  }
}
\author{Andrew Robbins}

