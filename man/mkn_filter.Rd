% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{mkn_filter}
\alias{mkn_filter}
\title{The Multiple-Knockoffs Filter}
\usage{
mkn_filter(X, y, k, nref = 1, knockoffs_fun = mkn_create_gaussian,
  scores_fun = mkn_scores_glmnet_coef, fstats_fun = mkn_fstats_max,
  pofstats_fun = NULL, knockoffs_args = list(), scores_args = list(),
  fstats_args = list(), pofstats_args = list(), randomized = FALSE,
  winnow = TRUE, inter = TRUE, offset = 1, verbose = TRUE,
  return_data = TRUE)
}
\arguments{
\item{X}{matrix or data.frame. Should be compatible to \code{knockoffs_fun} and \code{scores_fun}}

\item{y}{vector. Response/Outcome}

\item{k}{positive integer. The number of knockoffs}

\item{nref}{positive integer. The number of reference points. See Details}

\item{knockoffs_fun}{function to generate knockoffs. See Details}

\item{scores_fun}{function to calculate scores. See Details}

\item{fstats_fun}{function to calculate filtering statistics. See Details}

\item{pofstats_fun}{function to polish filtering statistics. See Details}

\item{knockoffs_args}{list. Extra arguments passed into knockoffs_fun}

\item{scores_args}{list. Extra arguments passed into scores_fun}

\item{fstats_args}{list. Extra arguments passed into fstats_fun}

\item{pofstats_args}{list. Extra arguments passed into pofstats_fun}

\item{randomized}{logical. Indicate whether the p-values are randomized when scores have ties.}

\item{winnow}{logical. Indicate whether the knockoffs are winnowed after initial score calculation}

\item{inter}{logical. Indicate whether interactive re-fitting is performed}

\item{offset}{0/1. The correction factor in the formula of estimated FDP. The values 1 yields Knockoffs+ procedure that is guaranteed to control FDR. The value 0 yields a more liberal procedure that is proved to control modified-FDR = E[V / (R + 1 / alpha)] where alpha is the target FDR level}

\item{verbose}{logical. Indicate whether the progress is printed to the console}

\item{return_data}{logical. Indicate whether to output the data in the form of \code{list(X = , Xk = , y = )} where Xk is the knockoff variables before winnowing.}
}
\value{
\item{qvals}{vector. The i-th element gives the minimum target FDR level for which the i-th hypothesis can be rejected. The value Inf refers to the case that the i-th hypothesis is never rejected}
\item{data}{list. NULL if return_date = FALSE. }
}
\description{
\code{mkn_filter} is a framework for model selection via multiple knockoff variables.
}
\details{
\code{mkn_filter} consists of five main steps:
\itemize{
 \item{Step 1: } generate \code{k} knockoffs for each variable in \code{X} by \code{knockoffs_fun};
 \item{Step 2: } calculate the initial scores for all variables and their knockoffs by \code{scores_fun} and then calculate the p-values, as the normalized rank;
 \item{Step 3: } winnow the original knockoff variables if \code{winnow = TRUE} and maintain one knockoff for each variable;
 \item{Step 4: } at each time when a certain number of hypotheses are revealed, re-calculate the scores for each variable by \code{scores_fun}, calculate the filtering statistics by \code{fstats_fun} and (optionally) polished the filtering statistics using the masked p-values;
 \item{Step 5: } reveal the hypotheses based on the filtering statistics, and the masked p-values, i.e. min\{p, 1 - p\}.
}

\code{knockoffs_fun}, \code{scores_fun}, \code{fstats_fun} should be function objects, and \code{pofstats_fun} could be either NULL or a function object:
\itemize{
 \item{\code{knockoffs_fun}} {should take at least \code{X} as the input and output an object of "mkn_knockoffs" class, which includes all necessary inputs for \code{scores_fun}. By default, \code{knockoffs_fun = mkn_create_gaussian}, which will generate model-X gaussian knockoffs in which case \code{mu} and \code{Sigma} should appear in \code{knockoffs_args}. See \code{\link{mkn_create_gaussian}} for other arguments and outputs.}
 \item{\code{scores_fun}} {should take at least four inputs: \code{X}, \code{y}, \code{knockoffs} as an object of class "mkn_knockoffs", and \code{subset} as a logical vector of size p which indicates the set of masked hypotheses. Given the input, \code{scores_fun} will calculate the scores based on all columns in \code{X} and the columns that correspond to \code{subset} in the knockoff variables. The output should be of class "mkn_scores" and in the form of \code{list(mask = , unmask = )} where \code{mask} is a matrix of size (s, k + 1) (\code{s = sum(subset)}) that gives the score for each variable and its knockoff variables in \code{subset}, and \code{unmask} is a vector of size (p - s) that gives the score for each variable in the complement of \code{subset}. By default, \code{scores_fun = mkn_scores_glmnet_coef}. See \code{\link{mkn_scores_glmnet_coef}} for details.}
\item{\code{fstats_fun}} {should take at least \code{scores}, of class "mkn_scores", as input. It should output an object of class "mkn_fstats" as a vector of scores that summarize the columns of \code{scores} in a symmetric way. By default, \code{fstats_fun = mkn_fstats_range} which will output the difference of the maximum and the minimum scores. If \code{fstats_fun = mkn_fstats_max}, it will output the maximum scores}
}
}
\examples{
\donttest{
## Generate data
n <- 1000
p <- 500
X <- matrix(rnorm(n * p), n, p)
beta <- c(rep(0.1, 50), rep(0, p - 50))
y <- X \%*\% beta + rnorm(n)
Sigma <- diag(p)

knockoffs_args <- list(mu = rep(0, p), Sigma = Sigma,
                       method = "sdp")

## Run mkn_filter
k <- 10
set.seed(1)
res <- mkn_filter(X, y, 10, knockoffs_args = knockoffs_args)

## Get the set of rejections when alpha = 0.1
which(res$qvals <= 0.1)
}
}
