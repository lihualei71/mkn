mkn_summarize.max <- function(score){
    apply(score, 1, max)
}

mkn_summarize.diff <- function(score){
    abs(score[, 1] - score[, 2])
}
