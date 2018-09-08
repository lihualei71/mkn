#' Calculate Filtering Statistics
#'
#' Calculate filtering statistics based on the scores, obtained for the original variables and their winnowed knockoff variables, and masked p-values.
#'
#' Given a sequence of filtering statistics, one could either sort the hypotheses by them directly or first sort them by masked p-values and then sort the hypotheses with the same masked p-values by the filtering statistics. The latter takes advantage of higher accuracy from multiple knockoffs. For \code{mkn_fstats_range} and \code{mkn_fstats_max}, the filtering statistics are simply modified by substracting \code{masked_pvals * M} from them for some large enough M. 
#' 
#' The following a list of available functions:
#' \itemize{
#' \item{\code{mkn_fstats_range}} {Calculate the difference between the maximum and the minimum absolute scores, with the above modification if \code{masked_pvals} is not \code{NULL}}
#' \item{\code{mkn_fstats_max}} {Calculate the maximum score, with the above modification if \code{masked_pvals} is not \code{NULL}}
#' }
#'
#' @param scores matrix with 2 columns. \code{nrow(scores)} should equal to \code{sum(mask)}
#' @name fstats
NULL

#' @rdname fstats
#'
#' @export
mkn_fstats_max <- function(scores){
    if (class(scores) != "mkn_scores"){
        stop("The input scores must be of class mkn_scores")
    }
    
    fstats <- apply(scores$mask, 1, max)

    structure(fstats, class = "mkn_fstats")
}

#' @rdname fstats
#'
#' @export
mkn_fstats_range <- function(scores){
    if (class(scores) != "mkn_scores"){
        stop("The input scores must be of class mkn_scores")
    }
    
    fstats <- apply(scores$mask, 1, function(x){
        max(x) - min(x)
    })

    structure(fstats, class = "mkn_fstats")    
}
