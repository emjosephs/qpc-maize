#' Kinship matrix function for complete data using the estimated variance across all loci
#'
#' This function makes a kinship matrix using the cov function. Unlike make_k_complete, it standardizes by the estimated genic variance across all loci, not each locus individually
#' @param myG matrix where the rows are individuals/populations and the columns are loci and the values are the allele frequency (not the # of copies present in an individual!!!).
#' @export


make_k_E <- function (myG) 
{
  scaleFactor = sqrt(mean(colMeans(myG) * (1 - colMeans(myG))))
    myM = dim(myG)[1]
    myT = matrix(data = -1/myM, nrow = myM - 1, ncol = myM)
    diag(myT) = (myM - 1)/myM
    myGstand = (myT %*% myG)/scaleFactor
    myK = cov(t(myGstand))
    return(myK)
}



