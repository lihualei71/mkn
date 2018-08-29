mineig <- function(A){
    min(eigen(A, symmetric = TRUE)$values)
}

mkn_solve_sdp <- function(Sigma, k, gaptol = 1e-06,
                          maxit = 1000, psdtol = 1e-09){
    stopifnot(isSymmetric(Sigma))
    ratio <- (k + 1) / k
    G <- stats::cov2cor(Sigma)
    p <- dim(G)[1]
    if (mineig(G) < psdtol) {
        stop("The covariance matrix is not positive-definite: cannot solve SDP", 
            immediate. = T)
    }
    Cl1 <- rep(0, p)
    Al1 <- -Matrix::Diagonal(p)
    Cl2 <- rep(1, p)
    Al2 <- Matrix::Diagonal(p)
    d_As <- c(diag(p))
    As <- Matrix::Diagonal(length(d_As), x = d_As)
    As <- As[which(Matrix::rowSums(As) > 0), ]
    Cs <- c(ratio * G)
    A <- cbind(Al1, Al2, As)
    C <- matrix(c(Cl1, Cl2, Cs), 1)
    K <- NULL
    K$s <- p
    K$l <- 2 * p
    b <- rep(1, p)
    OPTIONS <- NULL
    OPTIONS$gaptol <- gaptol
    OPTIONS$maxit <- maxit
    OPTIONS$logsummary <- 0
    OPTIONS$outputstats <- 0
    OPTIONS$print <- 0
    sol <- Rdsdp::dsdp(A, b, C, K, OPTIONS)
    if (!identical(sol$STATS$stype, "PDFeasible")) {
        warning("The SDP solver returned a non-feasible solution")
    }
    s <- sol$y
    s[s < 0] <- 0
    s[s > 1] <- 1
    psd <- 0
    s_eps <- 1e-08
    while (mineig(ratio * G - diag(s * (1 - s_eps), length(s))) < psdtol) {
            s_eps <- s_eps * 10
    }
    s <- s * (1 - s_eps)
    if (max(s) == 0) {
        warning("In creation of SDP knockoffs, procedure failed. Knockoffs will have no power.", 
            immediate. = T)
    }
    return(s * diag(Sigma))
}

mkn_solve_equi <- function(Sigma, k, psdtol = 1e-9){
    ratio <- (k + 1) / k
    G <- stats::cov2cor(Sigma)
    lambda_min <- mineig(G)
    Smat <- rep(lambda_min, nrow(Sigma))
    eps <- psdtol / 10
    while (mineig(G - diag(Smat * (1 - eps))) < psdtol){
        eps <- eps * 10
        if (eps > 1){
            stop("Solution not found!")
        }
    }
    return(Smat * (1 - eps) * diag(Sigma) * ratio)
}
