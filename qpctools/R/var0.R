#' Calculate variance where the mean is set to zero
#'
#' This function takes a string of numbers and calculates the variance of these numbers, assuming that the mean is 0.
#' @param x a string of vectors
#' @export


var0 <- function(x){  #variance where mean is set to 0
return(sum(x^2)/length(x))
}





