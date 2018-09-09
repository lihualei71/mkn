get_scores_info <- function(scores, randomized, nref,
                             winnow = TRUE){
    p <- nrow(scores)
    m <- ncol(scores)
    ties.method <- ifelse(randomized, "random", "last")
    rkmat <- t(apply(scores, 1, function(x){
        rank(-x, ties.method = ties.method)
    }))
    flag <- as.numeric(m %% (nref + 1) == 0)
    pvals <- (rkmat[, 1] - flag) / (m - flag)
    if (!winnow){
        return(list(pvals = pvals))
    }
    if (m %% (nref + 1) == 0){
        q <- floor(m / (nref + 1))
        refid <- apply(rkmat, 1, function(rk){
            tmp <- rk[1] - 1
            anchor <- ifelse(tmp < q, tmp, q - 1 - floor((tmp - q) / nref))
            reflist <- c(anchor, q + nref * (q - anchor) - 1:nref) + 1
            reflist <- reflist[!reflist == rk[1]]
            which(rk %in% reflist) - 1
        })
    } else {
        q <- floor(m / (nref + 1))
        r <- m - (nref + 1) * q
        refid <- apply(rkmat, 1, function(rk){
            tmp <- rk[1] - 1
            if (tmp >= (nref + 1) * q){
                return(rep(0, nref))
            }
            anchor <- ifelse(tmp < q, tmp, q - 1 - floor((tmp - q) / nref))
            reflist <- c(anchor, q + nref * (q - anchor) - 1:nref) + 1
            reflist <- reflist[!reflist == rk[1]]
            which(rk %in% reflist) - 1
        })
    }
    if (nref > 1){
        refid <- t(refid)
    }
    return(list(pvals = pvals, refid = as.matrix(refid)))
}

mkn_winnow <- function(knockoffs, scores,
                       winnow = TRUE,
                       nref = 1,
                       randomized = FALSE){
    p <- floor(ncol(knockoffs$Xk) / knockoffs$k)
    info <- get_scores_info(scores$mask, randomized, nref, winnow)
    mask <- rep(TRUE, p)
    if (winnow){
        winnow_knockoffs <- lapply(1:nref, function(i){
            id <- info$refid[, i]
            knockoffs_inds <- (1:p) + pmax(id - 1, 0) * p
            knockoffs$Xk[, knockoffs_inds]
        })
        winnow_knockoffs <- do.call(cbind, winnow_knockoffs)
        knockoffs$Xk <- winnow_knockoffs
        knockoffs$k <- nref
        mask <- (rowSums(info$refid) > 0)
    }
    return(list(pvals = info$pvals, mask = mask, knockoffs = knockoffs))
}
