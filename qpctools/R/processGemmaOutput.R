#' Process Gemma Output
#'
#' This function takes an output file from gemma and turns it into a useable data.frame
#' @param x the path to the gemma file

#' @export

processGemmaOutput <- function(x){
  gemmao = read.table(x, header=T, stringsAsFactors=F)
  gemmao$chr = sapply(gemmao$rs, function(x){as.numeric(substr(strsplit(x, '_')[[1]][1],2,100))})
  gemmao$ps = sapply(gemmao$rs, function(x){as.numeric(strsplit(x, '_')[[1]][2])})
  return(gemmao)
}


 
