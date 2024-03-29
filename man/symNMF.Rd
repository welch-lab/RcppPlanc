% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{symNMF}
\alias{symNMF}
\title{Perform Symmetric Non-negative Matrix Factorization}
\usage{
symNMF(
  x,
  k,
  niter = 30L,
  lambda = 0,
  algo = "gnsym",
  nCores = 2L,
  Hinit = NULL
)
}
\arguments{
\item{x}{Input matrix for factorization. Must be symmetric. Can be either
dense or sparse.}

\item{k}{Integer. Factor matrix rank.}

\item{niter}{Integer. Maximum number of symNMF interations.
Default \code{30}}

\item{lambda}{Symmetric regularization parameter. Must be
non-negative. Default \code{0.0} uses the square of the maximum value in
\code{x}.}

\item{algo}{Algorithm to perform the factorization, choose from "gnsym" or
"anlsbpp". Default \code{"gnsym"}}

\item{nCores}{The number of parallel tasks that will be spawned. Only applies to anlsbpp.
Default \code{2}}

\item{Hinit}{Initial right-hand factor matrix, must be of size n x k.
Default \code{NULL}.}
}
\value{
A list with the following elements:
\itemize{
 \item{\code{W} - the result left-hand factor matrix, non-empty when using
 \code{"anlsbpp"}}
 \item{\code{H} - the result right hand matrix.}
 \item{\code{objErr} - the objective error of the factorization.}
}
}
\description{
Symmetric input matrix \eqn{X} of size \eqn{n \times n} is required. Two
approaches are provided. Alternating Non-negative Least Squares Block
Principal Pivoting algorithm (ANLSBPP) with symmetric regularization, where
the objective function is set to be \eqn{\arg\min_{H\ge0,W\ge0}||X-WH||_F^2+
\lambda||W-H||_F^2}, can be run with \code{algo = "anlsbpp"}.
Gaussian-Newton algorithm, where the objective function is set to be
\eqn{\arg\min_{H\ge0}||X-H^\mathsf{T}H||_F^2}, can be run with \code{algo =
"gnsym"}. In the objectives, \eqn{W} is of size \eqn{n \times k} and \eqn{H}
is of size \eqn{k \times n}. The returned results will all be
\eqn{n \times k}.
}
\references{
Srinivas Eswar and et al., Distributed-Memory Parallel Symmetric Nonnegative
Matrix Factorization, SC '20, 2020, 10.5555/3433701.3433799
}
