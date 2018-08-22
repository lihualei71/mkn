mkn_scores_glmnet_coef <- function(X, Xk, y,
                                   subset = rep(TRUE, ncol(X)),
                                   use_LR = FALSE,
                                   nlambda = 100,
                                   lambda.min.ratio = 0.0001,
                                   cvtype = "lambda.min",
                                   family = "gaussian",
                                   cores = 1,
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
    if (family == "gaussian") {
        lambda_max <- max(abs(t(Xfull) %*% y))/n
        lambda_min <- lambda_max * lambda.min.ratio
        k <- (0:(nlambda - 1)) / nlambda
        lambda <- lambda_max * (lambda_min / lambda_max)^k
    }
    else {
        lambda <- NULL
    }
    fit <- glmnet::cv.glmnet(Xfull, y,
                             family = family,
                             parallel = parallel, ...)

    lambda <- fit$lambda.min
    scores <- abs(glmnet::coef.cv.glmnet(fit, s = cvtype)[-1])
    resids <- y - glmnet::predict.cv.glmnet(fit, Xfull, s = cvtype)
    if (use_LR){
        LR <- abs(Matrix::crossprod(Xfull, resids) / n)
        LR[scores != 0] <- lambda
        scores <- scores + LR
    }

    if (sum(!subset) > 0){
        mask <- matrix(scores[-which(!subset)], nrow = sum(subset))
    } else {
        mask <- matrix(scores, nrow = sum(subset))
    }
    unmask <- scores[1:p][!subset]

    structure(list(mask = mask, unmask = unmask), class = "mkn_scores")
}
