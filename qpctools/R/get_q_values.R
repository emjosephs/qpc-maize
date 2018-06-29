#' Calculate qvalues for a table of values
#'
#' Calculates qvalues for a table of values
#' @param ptable table of values
#' @export

get_q_values <- function(ptable){
qobj = qvalue(p = c(ptable))
myqvals = matrix(qobj$qvalues, nrow=dim(ptable)[1])
return(myqvals)
}

