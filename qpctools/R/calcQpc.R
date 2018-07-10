#' Calculate Qpc
#'
#' This function calculates Qpc given data about the relatedness matrix, and a set of trait values
#' @param myZ vector of traits. Not normalized yet.
#' @param myU matrix of eigenvectors of the kinship matrix (each column is an eigenvector)
#' @param myLambdas vector of eigenvalues of the kinship matrix 
#' @param myPCcutoff a value that determines how many PCs you want to look at. For example, 0.5 would mean that you would look at the first set of PCs that explain 0.5 of the variation. The default here is 0.3
#' @export
#' @examples
#' calcQpc()

calcQpc <- function(myZ, myU, myLambdas, myPCcutoff = 0.3){
  tailCutoff = round(.9*length(myLambdas))
  pcm = which(sapply(1:length(myLambdas), function(x){sum(myLambdas[1:x])/sum(myLambdas)}) > myPCcutoff)[1]
  
  myZ = myZ[-length(myZ)] - mean(myZ)
  myCm = (myZ %*% myU)/sqrt(myLambdas)
  myQm = sapply(1:pcm, function(n){
    var0(myCm[n])/var0(myCm[(tailCutoff-50):tailCutoff])
  })
  myPs = sapply(1:pcm, function(x){pf(myQm[x], 1, 50, lower.tail=F)})
  retdf = data.frame(qm = myQm, pvals = myPs)
  return(retdf)
  }



 
