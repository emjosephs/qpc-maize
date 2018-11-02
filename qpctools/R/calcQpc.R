#' Calculate Qpc
#'
#' This function calculates Qpc given data about the relatedness matrix, and a set of trait values
#' @param myZ vector of traits. Not normalized yet.
#' @param myU matrix of eigenvectors of the kinship matrix (each column is an eigenvector)
#' @param myLambdas vector of eigenvalues of the kinship matrix 
#' @param myPCcutoff a value that determines how many PCs you want to look at. For example, 0.5 would mean that you would look at the first set of PCs that explain 0.5 of the variation. The default here is 0.3
#' @param tailCutoff is there if you don't want to use the last PCs to estimate Va because of excess noise. The default value is 0.9, which means that you're not using the last 10% of your PCs. Set to 1 if you want to use all PCs
#' @param vapcs is the number of pcs used to estimate Va. Default is 50.
#' @export
#' @examples
#' calcQpc()

calcQpc <- function(myZ, myU, myLambdas, myPCcutoff = 0.3, tailCutoff = 0.9, vapcs = 50){
  myTailCutoff = round(tailCutoff*length(myLambdas)) #picks the end of the set of pcs used to calculate va
  pcm = which(sapply(1:length(myLambdas), function(x){sum(myLambdas[1:x])/sum(myLambdas)}) > myPCcutoff)[1] #the number of pcs tested
  
  myZ = myZ[1:dim(myU)[1]] - mean(myZ) #mean center phenotypes
  myCm = (myZ %*% myU)/sqrt(myLambdas) #project + standardize by the eigenvalues
  myQm = sapply(1:pcm, function(n){
    var0(myCm[n])/var0(myCm[(myTailCutoff-vapcs):myTailCutoff])
  })  #test for selection
  myPs = sapply(1:pcm, function(x){pf(myQm[x], 1, vapcs, lower.tail=F)}) #get a pvalue
  retdf = list(cm = myCm, qm = myQm, pvals = myPs)
  return(retdf)
  }



 
