mkn_pofstats_pvals_last <- function(fstats, mask, pvals){
    if (class(fstats) != "mkn_fstats"){
        stop("The input fstats must be of class mkn_fstats")
    }
    
    masked_pvals <- ifelse(mask, pmin(pvals, 1 - pvals), pvals)
    incr <- diff(sort(fstats))
    incr <- max(min(incr[incr > 0]), 1e-15)
    fstats <- fstats - masked_pvals[mask] * incr

    structure(fstats, class = "mkn_fstats")    
}

mkn_pofstats_pvals_first <- function(fstats, mask, pvals){
    if (class(fstats) != "mkn_fstats"){
        stop("The input fstats must be of class mkn_fstats")
    }
    
    masked_pvals <- ifelse(mask, pmin(pvals, 1 - pvals), pvals)    
    incr <- max(fstats) / max(min(masked <- pvals[masked <- pvals > 0]), 1e-4)
    fstats <- fstats - masked_pvals[mask] * incr

    structure(fstats, class = "mkn_fstats")
}

mkn_pofstats_gam <- function(fstats, mask, pvals,
                             thresh = min(pvals) + 1e-9){
    if (class(fstats) != "mkn_fstats"){
        stop("The input fstats must be of class mkn_fstats")
    }
    
    masked_pvals <- ifelse(mask, pmin(pvals, 1 - pvals), pvals)    
    if (!requireNamespace("mgcv")){
        stop("mgcv not installed")
    }

    x <- rank(fstats, ties.method = "average")
    y <- as.numeric(pvals <= thresh)
    df <- data.frame(x = x, y = y)
    mod <- mgcv::gam(y ~ s(x), family = "binomial", data = df)
    fstats <- as.numeric(mgcv::predict.gam(mod, df, type = "response"))

    structure(fstats, class = "mkn_fstats")                         
}
