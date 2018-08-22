mkn_fstats_max <- function(scores, mask, masked_pvals){
    fstats <- apply(abs(scores$mask), 1, max)
    if (!is.null(masked_pvals) && any(masked_pvals > 0)){
        fstats <- fstats - masked_pvals[mask] * max(fstats) /
            max(min(masked_pvals[masked_pvals > 0]), 1e-4)
    }
    return(fstats)
}

mkn_fstats_diff <- function(scores, mask, masked_pvals){
    fstats <- abs(scores$mask[, 1] - scores$mask[, 2])
    if (!is.null(masked_pvals) && any(masked_pvals > 0)){
        fstats <- fstats - masked_pvals[mask] * max(fstats) /
            max(min(masked_pvals[masked_pvals > 0]), 1e-4)
    }
    return(fstats)    
}
