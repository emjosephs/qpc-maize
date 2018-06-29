#' Calculate Va
#'
#' This function calculates additive genetic variation (Va) from a set of loci effect sizes and mean allele frequencies
#' @param afs vector of mean allele frequencies (mafs or afs are ok)
#' @param betas vector of effect sizes

#' @export

calcVa <-function(afs, betas){
  return(sum(2*afs*(1-afs)*(betas^2)))
}

 
