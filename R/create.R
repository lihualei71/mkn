block_t <- function(Z, n, p, k){
    Z <- matrix(Z, nrow = p)
    N <- n * k
    ind <- (0:(N-1)) %% k * n + floor((0:(N-1)) / k) + 1
    Z[, ind]
}

#' Generator of Multiple Gaussian Model-X Knockoffs
#'
#' \code{mkn_create_gaussian} generates \code{k} gaussian model-X knockoffs Xk such that the joint distribution of each row of (X, Xk) is multivariate-gaussian with mean \code{rep(mu, k)} and variance being a matrix with diagonal blocks \code{Sigma} and off-diagonal blocks Sigma-S where S is a diagonal matrix that satisfies S <= (k+1) / k Sigma.
#'
#' \code{mkn_create_gaussian} first calculates the diagonal matrix S that is smaller than \code{Sigma}, without a correction factor. Then it modifies S by multiplying it by \code{s_const}, which should be smaller than (k + 1) / k by definition. 
#' 
#' @param X covariate matrix
#' @param k positive integer. The number of knockoffs
#' @param mu vector. Mean vector of each row of X
#' @param Sigma matrix. Covariance matrix of each row of X
#' @param method string. Should be "asdp" (approximate SDP) or "sdp" or "equi". See \code{\link[knockoff]{create.solve_asdp}}, \code{\link[knockoff]{create.solve_sdp}} and \code{\link[knockoff]{create.solve_equi}} for details
#' @param s_const positive real. Should be smaller than (k + 1) / k. See Details
#' @param diag_s an optional vector of the diagonal elements of S
#' @param ... other arguments passed into \code{\link[knockoff]{create.solve_asdp}}, \code{\link[knockoff]{create.solve_sdp}} or \code{\link[knockoff]{create.solve_equi}}
#'
#' @return a matrix of size (n, pk) consisting all knockoff variables
#' 
#' @examples
#' \donttest{
#' ## Generate X from an AR(1) process (with Toeplitz covariance matrix)
#' n <- 100
#' p <- 50
#' rho <- 0.5
#' Sigma <- rho^stats::toeplitz(0:(p - 1))
#' mu <- rep(0, p)
#' k <- 10
#' X <- matrix(rnorm(n * p), n, p) %*% Matrix::chol(Sigma)
#' 
#' ## Generate knockoffs via SDP
#' Xk <- mkn_create_gaussian(X, k, mu, Sigma, method = "sdp")
#' }
#' @export
mkn_create_gaussian <- function(X, k, mu, Sigma,
                                method = c("asdp", "sdp", "equi"),
                                s_const = (k + 1) / k,
                                diag_s = NULL,
                                ...){
    if (!is.matrix(X)){
        X <- as.matrix(X)
    }
    
    if (s_const > (k + 1)){
        warning("The multiplier of S cannot exceed (k + 1) / k. Replace it by (k + 1) / k.")
    }
    s_const <- min(s_const, (k + 1) / k - 1e-6)

    method <- match.arg(method)[1]
    if ((nrow(Sigma) <= 500) && method == "asdp") {
        method <- "sdp"
    }
    if (is.null(diag_s)) {
        diag_s = switch(match.arg(method),
            equi = mkn_create_solve_equi(Sigma, ...),
            sdp = mkn_create_solve_sdp(Sigma, ...),
            asdp = mkn_create_solve_asdp(Sigma, ...))
    }
    diag_s <- diag_s * s_const

    n <- nrow(X)
    p <- nrow(Sigma)
    Sigma_inv_S <- solve(Sigma, diag(diag_s))
    mu_cond <- X - sweep(X, 2, mu, "-") %*% Sigma_inv_S

    #### Naive approach
    ## time.start <- Sys.time()
    ## SigmaInv_S <- solve(Sigma, diag(diag_s))
    ## Lambda <- diag(diag_s) - diag_s * SigmaInv_S
    ## cov_cond <- Matrix::kronecker(matrix(1, k, k), Lambda) + Matrix::kronecker(diag(k), diag(diag_s))
    ## Z <- matrix(stats::rnorm(n * p * k), nrow = n) %*% chol(cov_cond)
    ## time.end <- Sys.time()
    ## time.end - time.start

    sqrtS <- sqrt(diag_s)
    if (s_const <= 1){
        ## ~ 28x faster than the naive approach with (n, p, k) = (10p00, 500, 10)
        off_diag <- diag(diag_s) - diag_s * Sigma_inv_S
        addon <- matrix(stats::rnorm(n * p), nrow = n) %*%
            Matrix::chol(off_diag) + mu_cond
        Xk <- rep(sqrtS, k) * matrix(stats::rnorm(n * p * k), ncol = n)
        Xk <- t(Xk) + as.numeric(addon)
    } else {
        ## ~ 14x faster than the naive approach with (n, p, k) = (1000, 500, 10)
        p <- length(diag_s)
        res <- base::eigen(sqrtS * solve(Sigma, diag(sqrtS)), symmetric = TRUE)
        D <- res$values
        U <- sqrtS * res$vectors
        ## Step 1
        Xk <- matrix(stats::rnorm(n * p * k), ncol = k)
        ## Step 2
        Xk[, 1] <- Xk[, 1] * sqrt(rep(1 + k * (1 - D), n))
        ## Step 3
        V <- cbind(rep(1, k) / sqrt(k), MASS::Null(rep(1, k)))
        Xk <- Xk %*% t(V)
        ## Step 4
        Xk <- block_t(Xk, n, p, k)
        ## Step 5
        Xk <- U %*% Xk
        ## Step 6
        Xk <- matrix(Xk, nrow = p * k)
        ## Step 7
        Xk <- t(Xk) + as.numeric(mu_cond)
    }

    structure(list(Xk = Xk, k = k,
                   Sigma = Sigma, diag_s = diag_s),
              class = "mkn_knockoffs")
}
