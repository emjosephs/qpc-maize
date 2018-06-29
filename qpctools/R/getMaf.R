#' Calculate minor allele frequency
#'
#' This function takes an allele frequency (range 0 to 1) and calculates the minor allele frequency (0 to 0.5)
#' @param af an allele frequency
#' @export


getMaf = function(af){if (af <= 0.5){maf=af}
  else{ maf = 1-af}
  return(maf)
}






