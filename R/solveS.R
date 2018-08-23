mineig <- function(A){
    min(eigen(A, symmetric = TRUE)$values)
}

mkn_create_solve_sdp <- function(Sigma,
                                 rho = 0.9, gaptol = 1e-06,
                                 maxit = 1000, psdtol = 1e-9){
    Smat <- knockoff::create.solve_sdp(Sigma, gaptol, maxit)
    while (TRUE){
        if (mineig(Sigma - diag(Smat)) > psdtol){
            break
        }
        Smat <- Smat * rho
    }
    return(Smat)
}

mkn_create_solve_asdp <- function(Sigma,
                                  rho = 0.9,
                                  nBlocks = 10, cores = 1,
                                  gaptol = 1e-06,
                                  maxit = 1000, psdtol = 1e-9){
    Smat <- knockoff::create.solve_asdp(Sigma, nBlocks, cores, gaptol, maxit)
    while (TRUE){
        if (mineig(Sigma - diag(Smat)) > psdtol){
            break
        }
        Smat <- Smat * rho
    }
    return(Smat)
}

mkn_create_solve_equi <- function(Sigma,
                                  tol = 1e-9, psdtol = 1e-9){
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
    return(Smat * (1 - eps) * diag(Sigma))
}
