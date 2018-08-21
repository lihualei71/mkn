mkn_stat.glmnet_coef <- function(X, Xk, y,
                                 subset = rep(TRUE, ncol(X)),
                                 family = "gaussian",
                                 cores = 2,
                                 ...){
    if (!requireNamespace("glmnet", quietly = TRUE)) 
        stop("glmnet is not installed", call. = FALSE)
    parallel <- TRUE
    if (!requireNamespace("doMC", quietly = TRUE)) {
        warning("doMC is not installed. Without parallelization, the statistics will be slower to compute", 
                call. = FALSE, immediate. = TRUE)
        parallel <- FALSE
    }
    if (!requireNamespace("parallel", quietly = TRUE)) {
        warning("parallel is not installed. Without parallelization, the statistics will be slower to compute.", 
                call. = FALSE, immediate. = TRUE)
        parallel <- FALSE
    }
    if (parallel) {
        ncores <- parallel::detectCores(all.tests = TRUE, logical = TRUE)
        if (cores == 2) {
            cores <- min(2, ncores)
        }
        else {
            if (cores > ncores) {
                warning(paste("The requested number of cores is not available. Using instead", 
                              ncores, "cores"), immediate. = TRUE)
                cores <- ncores
            }
        }
        if (cores > 1) {
            doMC::registerDoMC(cores = cores)
            parallel <- TRUE
        }
        else {
            parallel <- FALSE
        }
    }
    ## Z <- knockoff:::cv_coeffs_glmnet(cbind(X, Xk), y,
    ##                                  family = family,
    ##                                  parallel = parallel, ...)
    p <- ncol(X)    
    if (any(!is.logical(subset))){
        tmp <- rep(FALSE, p)
        tmp[subset] <- TRUE
        subset <- tmp
    }
    k <- floor(ncol(Xk) / p)
    if (ncol(Xk) > p * k){
        stop("The number of knockoffs is not an integer.")
    }
    inds <- rep(subset, k)
    Xfull <- cbind(X, Xk[, inds])
    fit <- glmnet::cv.glmnet(Xfull, y,
                             standardize = FALSE,
                             family = family,
                             parallel = parallel, ...)

    lambda <- fit$lambda.min
    Z <- abs(coef(fit, s = "lambda.min")[-1])
    resids <- y - predict(fit, Xfull, s = "lambda.min")    
    LR <- abs(Matrix::crossprod(Xfull, resids) / n)
    LR[Z != 0] <- lambda
    score <- Z + LR

    if (sum(!subset) > 0){
        mask <- matrix(score[-which(!subset)], nrow = sum(subset))
    } else {
        mask <- matrix(score, nrow = sum(subset))
    }
    unmask <- score[1:p][!subset]
    
    structure(list(mask = mask, unmask = unmask), class = "mkn_score")
}
