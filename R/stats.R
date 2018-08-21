#' @export
mkn_stat.glmnet_coef <- function(X, Xk, y,
                                 subset = rep(TRUE, ncol(X)),
                                 family = "gaussian",
                                 cores = 2,
                                 ...){
    if (!requireNamespace("glmnet", quietly = TRUE))
        stop("glmnet is not installed", call. = FALSE)

    parallel <- (cores > 1)
    n <- nrow(X)
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
    Z <- abs(glmnet::coef.cv.glmnet(fit, s = "lambda.min")[-1])
    resids <- y - glmnet::predict.cv.glmnet(fit, Xfull, s = "lambda.min")
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
