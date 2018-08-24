fdp_hat <- function(A, R, offset){
    (offset + A) / pmax(1, R)
}

get_scores_info <- function(scores){
    m <- ncol(scores)
    if (m %% 2 == 0){
        info <- apply(scores, 1, function(x){
            rk <- rank(-x, ties.method = "last")
            pval <- (rk[1] - 1) / (m - 1)
            reflectid <- which(rk == m + 1 - rk[1]) - 1
            c(pval, reflectid)
        })
    } else {
        info <- apply(scores, 1, function(x){
            rk <- rank(-x, ties.method = "last")
            pval <- rk[1] / m
            if (rk[1] == m){
                reflectid <- 0
            } else {
                reflectid <- which(rk == m - rk[1]) - 1
            }
            c(pval, reflectid)
        })
    }
    t(info)
}

get_nreveals <- function(nmasks, ninter){
    fit_stamps <- c(seq(nmasks, 0, -max(floor(nmasks / ninter), 1))[1:ninter], 0)
    diff_stamps <- -diff(fit_stamps)
    diff_stamps <- diff_stamps[!is.na(diff_stamps)]
    diff_stamps[diff_stamps > 0]
}

#' The Multiple-Knockoffs Filter
#'
#' \code{mkn_filter} is a framework for model selection via multiple knockoff variables.
#'
#' \code{mkn_filter} consists of five main steps:
#' \itemize{
#'  \item{Step 1: } generate \code{k} knockoffs for each variable in \code{X} by \code{knockoffs_fun};
#'  \item{Step 2: } calculate the initial scores for all variables and their knockoffs by \code{scores_fun} and then calculate the p-values, as the normalized rank;
#'  \item{Step 3: } winnow the original knockoff variables if \code{winnow = TRUE} and maintain one knockoff for each variable;
#'  \item{Step 4: } at each time when a certain number of hypotheses are revealed, re-calculate the scores for each variable by \code{scores_fun} and calculate the filtering statistics by \code{fstats_fun}.
#'  \item{Step 5: } reveal the hypotheses based on the filtering statistics, and the masked p-values, i.e. min\{p, 1 - p\}, if \code{use_masked_pvals = TRUE}.
#' }
#' 
#' \code{knockoffs_fun}, \code{scores_fun}, \code{fstats_fun} should be function objects:
#' \itemize{
#'  \item{\code{knockoffs_fun}} {should take at least \code{X} as the input and output a matrix of size (n, pk) where (n, p) = dim(X). By default, \code{knockoffs_fun = mkn_create_gaussian}, which will generate model-X gaussian knockoffs in which case \code{mu} and \code{Sigma} should appear in \code{knockoffs_args}. See \code{\link{mkn_create_gaussian}} for other arguments.}
#'  \item{\code{scores_fun}} {should take at least four inputs: \code{X}, \code{Xk} (for the output of \code{knockoffs_fun}), \code{y} and \code{subset} as a logical vector of size p which indicates the set of masked hypotheses. Given the input, \code{scores_fun} will calculate the scores based on all columns in \code{X} and the columns that correspond to \code{subset} in \code{Xk}. The output should be in the form of \code{list(mask = , unmask = )} where \code{mask} is a matrix of size (s, k + 1) (\code{s = sum(subset)}) that gives the score for each variable and its knockoff variables in \code{subset}, and \code{unmask} is a vector of size (p - s) that gives the score for each variable in the complement of \code{subset}. By default, \code{scores_fun = mkn_scores_glmnet_coef}. See \code{\link{mkn_scores_glmnet_coef}} for details.}
#' \item{\code{fstats_fun}} {should take at least three inputs: \code{score}, a matrix that corresponds to the original variables and their knockoff variables, \code{mask}, a logical vector that indicates whether each hypothesis is masked, and \code{masked_pvals}, a numeric vector that gives the masked p-values (min\{p, 1 - p\}). It should output a vector of scores by combining the columns of \code{score} in a symmetric way. By default, \code{fstats_fun = mkn_fstats_range} which will output the difference of the maximum and the minimum scores. If \code{fstats_fun = mkn_fstats_max}, it will output the maximum scores}
#' }
#'
#' If \code{use_masked_pvals = TRUE}, by default, the hypotheses will be first sorted by masked p-values from the largest (least promising) to the smallest (least promising) and then sorted by the filtering statistics from the smallest. If \code{use_masked_pvals = TRUE}, the hypotheses will only be sorted by the filter statistics.
#' 
#' @param X matrix or data.frame. Should be compatible to \code{knockoffs_fun} and \code{scores_fun}
#' @param y vector. Response/Outcome
#' @param k positive integer. The number of knockoffs
#' @param knockoffs_fun function to generate knockoffs. See Details
#' @param scores_fun function to calculate scores. See Details
#' @param fstats_fun function to calculate filtering statistics. See Details
#' @param knockoffs_args list. Extra arguments passed into knockoffs_fun
#' @param scores_args list. Extra arguments passed into scores_fun
#' @param fstats_args list. Extra arguments passed into fstats_fun
#' @param use_masked_pvals logical. Indicate whether the masked p-values are incorporated into the ranking. See Details
#' @param winnow logical. Indicate whether the knockoffs are winnowed after initial score calculation
#' @param inter logical. Indicate whether interactive re-fitting is performed 
#' @param ninter positive integer. Number of rounds of interactive re-fitting. See Details
#' @param offset 0/1. The correction factor in the formula of estimated FDP. The values 1 yields Knockoffs+ procedure that is guaranteed to control FDR. The value 0 yields a more liberal procedure that is proved to control modified-FDR = E[V / (R + 1 / alpha)] where alpha is the target FDR level
#' @param verbose logical. Indicate whether the progress is printed to the console
#' @param return_data logical. Indicate whether to output the data in the form of \code{list(X = , Xk = , y = )} where Xk is the knockoff variables before winnowing.
#' 
#' @return
#' \item{qvals}{vector. The i-th element gives the minimum target FDR level for which the i-th hypothesis can be rejected. The value Inf refers to the case that the i-th hypothesis is never rejected}
#' \item{data}{list. NULL if return_date = FALSE. }
#' 
#' @examples
#' \donttest{
#' ## Generate data
#' n <- 1000
#' p <- 500
#' X <- matrix(rnorm(n * p), n, p)
#' beta <- c(rep(0.1, 50), rep(0, p - 50))
#' y <- X %*% beta + rnorm(n)
#' Sigma <- diag(p)
#'
#' knockoffs_args <- list(mu = rep(0, p), Sigma = Sigma,
#'                        method = "sdp")
#'
#' ## Run mkn_filter
#' k <- 10
#' set.seed(1)
#' res <- mkn_filter(X, y, 10, knockoffs_args = knockoffs_args)
#'
#' ## Get the set of rejections when alpha = 0.1
#' which(res$qvals <= 0.1)
#' }
#' @export
mkn_filter <- function(X, y, k,
                       knockoffs_fun = mkn_create_gaussian,
                       scores_fun = mkn_scores_glmnet_coef,
                       fstats_fun = mkn_fstats_range,
                       knockoffs_args = list(),
                       scores_args = list(),
                       fstats_args = list(),
                       use_masked_pvals = TRUE,
                       winnow = TRUE,
                       inter = TRUE,
                       ninter = 10,
                       offset = 1,
                       verbose = TRUE,
                       return_data = TRUE){
    if (!is.function(knockoffs_fun) ||
        !"X" %in% methods::formalArgs(knockoffs_fun)){
        stop("knockoffs_fun should be a valid function with at least X as input")
    }
    if (!is.function(scores_fun) ||
        any(!c("X", "Xk", "y", "subset") %in% methods::formalArgs(scores_fun))){
        stop("scores_fun should be a valid function with at least X, Xk, y and subset as inputs")
    }
    if (!is.function(fstats_fun) ||
        any(!c("scores", "mask", "masked_pvals") %in% methods::formalArgs(fstats_fun))){
        stop("fstats_fun should be a valid function")
    }
    
    if (!offset %in% c(0, 1)){
        stop("offset must be either 1 or 0")
    }
    if (k != floor(k) || k <= 0){
        stop("k must be a positive integer")
    }
    if (k == 1){
        winnow <- FALSE
    }
    n <- nrow(X)
    if (length(y) != n){
        stop("y should have the same length as nrow(X)")
    }
    if (!inter){
        ninter <- 1
    }
        
    p <- ncol(X)
    knockoffs_args <- c(list(X = X, k = k), knockoffs_args)
    if (verbose){
        cat("Generating knockoff variables\n")
    }
    Xk <- do.call(knockoffs_fun, knockoffs_args)
    if (verbose){
        cat("Knockoff variables generated\n")
        cat("\n")
    }

    mask <- rep(TRUE, p)

    scores_args_root <- c(list(X = X, Xk = Xk, y = y), scores_args)
    scores_args <- c(list(subset = mask), scores_args_root)
    if (verbose){
        cat("Fitting the initial scores\n")
    }
    scores <- do.call(scores_fun, scores_args)
    if (verbose){
        cat("Initial scores obtained\n")
        cat("\n")
    }

    info <- get_scores_info(scores$mask)
    pvals <- info[, 1]
    winnow_inds <- info[, 2]
    if (winnow){
        knockoff_inds <- (1:p) + pmax(winnow_inds - 1, 0) * p
        scores_args_root$Xk <- Xk[, knockoff_inds]
    }

    A <- sum(pvals > 0.5)
    R <- sum(pvals < 0.5)
    minfdp <- fdp_hat(A, R, offset)
    fdp_list <- minfdp
    reveal_order <- vector("integer")

    if (any(winnow_inds == 0)){
        init_inds <- which(winnow_inds == 0)
        init_scores <- scores$mask[init_inds, 1]
        init_inds <- init_inds[order(init_scores)]
        fdp_list <- c(fdp_list,
                      fdp_hat(A - 1:length(init_inds), R, offset))
        A <- A - length(init_inds)
        reveal_order <- c(reveal_order, init_inds)
        mask[init_inds] <- FALSE
    }

    nmasks <- sum(mask)
    if (nmasks > 0){
        nreveal_list <- get_nreveals(nmasks, ninter)
    }
    if (use_masked_pvals){
        masked_pvals <- pmin(pvals, 1 - pvals)
    }

    if (verbose && nmasks > 0){
        cat("Interactive re-fitting starts:\n")
        pb <- utils::txtProgressBar(min = 0, max = ninter,
                                    style = 3, width = 50)
    }

    if (use_masked_pvals){
        fstats_args_root <- c(list(masked_pvals = masked_pvals), fstats_args)
    } else {
        fstats_args_root <- c(list(masked_pvals = NULL), fstats_args)
    }

    for (i in 1:ninter){
        nmasks <- sum(mask)
        if (nmasks == 0){
            break
        }
        
        if (verbose){
            utils::setTxtProgressBar(pb, i)
        }

        nreveal <- nreveal_list[i]
        scores_args <- c(list(subset = mask), scores_args_root)
        scores <- do.call(scores_fun, scores_args)
        ## if (!winnow && k > 1){
        ##     knockoff_inds <- cbind(1:nmasks, winnow_inds[mask])
        ##     scores$mask <- cbind(scores$mask[, 1], scores$mask[knockoff_inds])
        ## }
        fstats_args <- c(list(scores = scores, mask = mask), fstats_args_root)
        fstats <- do.call(fstats_fun, fstats_args)

        rm_inds <- which(mask)[order(fstats)[1:nreveal]]

        Adecre <- cumsum(pvals[rm_inds] > 0.5)
        Rdecre <- cumsum(pvals[rm_inds] < 0.5)
        fdp_list <- c(fdp_list,
                      fdp_hat(A - Adecre, R - Rdecre, offset))
        reveal_order <- c(reveal_order, rm_inds)
        A <- A - max(Adecre)
        R <- R - max(Rdecre)

        mask[rm_inds] <- FALSE
    }
    if (verbose){
        cat("\n")
        cat("Interactive re-fitting ends\n")
    }

    qvals <- rep(1, p)
    qvals[reveal_order] <- cummin(fdp_list[1:p])
    qvals[pvals > 0.5] <- Inf
    if (return_data){
        data <- list(X = X, Xk = Xk, y = y)
    } else {
        data <- NULL
    }

    structure(list(call = match.call(),
                   qvals = qvals, data = data),
              class = "mkn_result")
}
