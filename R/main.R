fdp_hat <- function(A, R, offset){
    (offset + A) / pmax(1, R)
}

get_score_info <- function(score){
    m <- ncol(score)
    if (m %% 2 == 0){
        info <- apply(score, 1, function(x){
            rk <- rank(-x, ties.method = "last")
            pval <- (rk[1] - 1) / (m - 1)
            reflectid <- which(rk == m + 1 - rk[1]) - 1
            c(pval, reflectid)
        })
    } else {
        info <- apply(score, 1, function(x){
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
    fit_stamps <- c(seq(nmasks, 0, -floor(nmasks / ninter))[1:ninter], 0)
    -diff(fit_stamps)
}

mkn.filter <- function(X, y, k,
                       knockoffs_fun = mkn_create.gaussian,
                       stats_fun= mkn_stat.glmnet_coef,
                       summarize_fun = mkn_summarize.max,
                       knockoffs_args = list(),
                       stats_args = list(),
                       summarize_args = list(),
                       purify = TRUE,
                       inter = TRUE,
                       use_masked_pvals = TRUE,
                       ninter = 20,
                       offset = 1,
                       verbose = TRUE){
    stopifnot(is.function(knockoffs_fun))
    stopifnot(is.function(stats_fun))
    stopifnot(offset %in% c(0, 1))
    stopifnot(k == floor(k) && k > 0)

    n <- nrow(X)
    p <- ncol(X)
    knockoffs_args <- c(list(X = X, k = k), knockoffs_args)
    Xk <- do.call(knockoffs_fun, knockoffs_args)    
    mask <- rep(TRUE, p)
    
    stats_args_root <- c(list(X = X, Xk = Xk, y = y), stats_args)
    stats_args <- c(list(subset = mask), stats_args_root)
    score <- do.call(stats_fun, stats_args)

    info <- get_score_info(score$mask)
    pvals <- info[, 1]
    inds <- info[, 2]
    A <- sum(pvals > 0.5)
    R <- sum(pvals < 0.5)
    minfdp <- fdp_hat(A, R, offset)
    fdp_list <- minfdp
    reveal_order <- vector("integer")    
    
    if (any(inds == 0)){
        init_inds <- which(inds == 0)
        init_score <- score$mask[init_inds, 1]
        init_inds <- init_inds[order(init_score)]
        fdp_list <- c(fdp_list,
                      fdp_hat(A - 1:length(init_inds), R, offset))
        A <- A - length(init_inds)
        reveal_order <- c(reveal_order, init_inds)
        mask[init_inds] <- FALSE
    }

    if (purify){
        knockoff_inds <- (1:p) + pmax(info[, 2] - 1, 0) * p
        Xk <- Xk[, knockoff_inds]
        stats_args_root$Xk <- Xk
    }

    nmasks <- sum(mask)
    if (!inter){
        ninter <- 1
    }
    nreveal_list <- get_nreveals(nmasks, ninter)
    if (use_masked_pvals){
        masked_pvals <- pmin(pvals, 1 - pvals)
    }

    if (verbose){
        pb <- utils::txtProgressBar(min = 0, max = ninter,
                                    style = 3, width = 50)
    }
    
    for (i in 1:ninter){
        if (verbose){
            utils::setTxtProgressBar(pb, i)
        }
        
        nreveal <- nreveal_list[i]
        stats_args <- c(list(subset = mask), stats_args_root)
        score <- do.call(stats_fun, stats_args)
        stats <- summarize_fun(score$mask)

        if (use_masked_pvals){
            stats <- stats - masked_pvals[mask] * max(stats) * (k + 1)
        }
        rm_inds <- which(mask)[order(stats)[1:nreveal]]

        Adecre <- cumsum(pvals[rm_inds] > 0.5)
        Rdecre <- cumsum(pvals[rm_inds] < 0.5)
        fdp_list <- c(fdp_list,
                      fdp_hat(A - Adecre, R - Rdecre, offset))
        reveal_order <- c(reveal_order, rm_inds)

        mask[rm_inds] <- FALSE
    }

    print(sum(mask))
    qvals <- rep(1, p)
    qvals[reveal_order] <- cummin(fdp_list[1:p])
    return(qvals)
}
