% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create.R
\name{mkn_create_gaussian}
\alias{mkn_create_gaussian}
\title{Generator of Multiple Gaussian Model-X Knockoffs}
\usage{
mkn_create_gaussian(X, k, mu, Sigma, method = c("sdp", "equi"),
  diag_s = NULL, ...)
}
\arguments{
\item{X}{covariate matrix}

\item{k}{positive integer. The number of knockoffs}

\item{mu}{vector. Mean vector of each row of X}

\item{Sigma}{matrix. Covariance matrix of each row of X}

\item{method}{string. Should be "asdp" (approximate SDP) or "sdp" or "equi". See \code{\link[knockoff]{create.solve_asdp}}, \code{\link[knockoff]{create.solve_sdp}} and \code{\link[knockoff]{create.solve_equi}} for details}

\item{diag_s}{an optional vector of the diagonal elements of S}

\item{...}{other arguments passed into \code{\link[knockoff]{create.solve_asdp}}, \code{\link[knockoff]{create.solve_sdp}} or \code{\link[knockoff]{create.solve_equi}}}
}
\value{
a matrix of size (n, pk) consisting all knockoff variables
}
\description{
\code{mkn_create_gaussian} generates \code{k} gaussian model-X knockoffs Xk such that the joint distribution of each row of (X, Xk) is multivariate-gaussian with mean \code{rep(mu, k)} and variance being a matrix with diagonal blocks \code{Sigma} and off-diagonal blocks Sigma-S where S is a diagonal matrix that satisfies S <= (k+1) / k Sigma.
}
\examples{
\donttest{
## Generate X from an AR(1) process (with Toeplitz covariance matrix)
n <- 100
p <- 50
rho <- 0.5
Sigma <- rho^stats::toeplitz(0:(p - 1))
mu <- rep(0, p)
k <- 10
X <- matrix(rnorm(n * p), n, p) \%*\% Matrix::chol(Sigma)

## Generate knockoffs via SDP
Xk <- mkn_create_gaussian(X, k, mu, Sigma, method = "sdp")
}
}
