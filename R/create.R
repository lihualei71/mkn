block_t <- function(Z, n, p, k){
    Z <- matrix(Z, nrow = p)
    N <- n * k
    ind <- (0:(N-1)) %% k * n + floor((0:(N-1)) / k) + 1
    Z[, ind]
}

mkn_create.gaussian <- function(X, k, mu, Sigma,
                                method = c("asdp", "sdp", "equi"),
                                s_const = 1,
                                diag_s = NULL,
                                ...){
    if (s_const > (k + 1) / k - 1e-9){
        stop("The multiplier of S cannot exceed (k + 1) / k.")
    }

    method <- match.arg(method)[1]
    if ((nrow(Sigma) <= 500) && method == "asdp") {
        method <- "sdp"
    }
    if (is.null(diag_s)) {
        diag_s = switch(match.arg(method),
            equi = mkn_create.solve_equi(Sigma, ...), 
            sdp = mkn_create.solve_sdp(Sigma, ...),
            asdp = mkn_create.solve_asdp(Sigma, ...))
        diag_s <- diag_s * s_const
    }

    n <- nrow(X)
    p <- nrow(Sigma)
    mu_cond <- t(
        solve(k * Sigma - (k - 1) * diag(diag_s)) %*%
        (Sigma - diag(diag_s)) %*% (t(X) - mu) + mu
        )

    #### Naive approach
    ## time.start <- Sys.time()
    ## SigmaInv_S <- solve(Sigma, diag(diag_s))
    ## Lambda <- diag(diag_s) - diag_s * SigmaInv_S
    ## cov_cond <- Matrix::kronecker(matrix(1, k, k), Lambda) + Matrix::kronecker(diag(k), diag(diag_s))
    ## Z <- matrix(rnorm(n * p * k), nrow = n) %*% chol(cov_cond)
    ## time.end <- Sys.time()
    ## time.end - time.start

    sqrtS <- sqrt(diag_s)    
    if (s_const <= 1){
        ## ~ 28x faster than the naive approach with (n, p, k) = (1000, 500, 10)
        off_diag <- diag(diag_s) - diag_s * (solve(Sigma, diag(diag_s)))
        addon <- matrix(rnorm(n * p), nrow = n) %*%
            chol(off_diag) + mu_cond
        Xk <- rep(sqrtS, k) * matrix(rnorm(n * p * k), ncol = n)
        return(t(Xk) + as.numeric(addon))
    } else {
        ## ~ 14x faster than the naive approach with (n, p, k) = (1000, 500, 10)       
        p <- length(diag_s)
        res <- base::eigen(sqrtS * solve(Sigma, diag(sqrtS)), symmetric = TRUE)
        D <- res$values
        U <- sqrtS * res$vectors
        ## Step 1
        Xk <- matrix(rnorm(n * p * k), ncol = k)
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
        return(t(Xk) + as.numeric(mu_cond))
    }
}
