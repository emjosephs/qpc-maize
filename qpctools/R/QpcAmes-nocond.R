#' Calculate Qpc on Ames panel
#'
#' This function calculates Qpc on Ames panel.
#' @param myI the trait number that we're looking at
#' @param myM the number of lines in the genotyping panel
#' @param gwasPrefix path prefix for the GWAS results
#' @param sigPrefix path prefix for the genotypes of GWAS loci in the genotyping panel
#' @param cutoff p value cutoff for including SNPs. Default is 1.
#' @param mysigma the kinship matrix for the genotyping panel.
#' @param myLambda a list of eigenvalues of the genotyping panel kinship matrix
#' @param myU a matrix of eigenvectors of the genotyping panel kinship matrix
#' @export

Qpcames_nocond <- function(myI, myM = 2704, gwasPrefix = 'data/281-gwas-results/ldfiltered.', sigPrefix='data/281-gwas-results/sigSnps.',
                           mypcmax = pcmax, myU = amesEig$vectors, myLambda = amesEig$values)
  {
  tailCutoff = round(0.9 * myM)
  
#read in data
gwasHits = read.table(paste(gwasPrefix,myI,sep=""), stringsAsFactors=F)
names(gwasHits) = c('x','y',strsplit('chr     rs      ps      n_miss  allele1 allele0 af      beta    se      l_remle l_mle   p_wald  p_lrt   p_score', split=' +')[[1]])
sigGenos = read.table(paste(sigPrefix,myI, sep=""), header=T, stringsAsFactors=F)

#combine table of GWAS results with genotypes in the GWAS set
combData = dplyr::left_join(sigGenos, gwasHits, by = c('locus'='rs'))
myBetas = as.matrix(combData$beta)
myG = t(as.matrix(sigGenos[,4:(myM+3)]))

#center genotype matrix
myT = matrix(data = -1/myM, nrow = myM - 1, ncol = myM)
diag(myT) = (myM - 1)/myM
myGcent = myT %*% myG

#calculate breeding values
allZ = myGcent %*% myBetas

#project breeding values onto PCs and standardize
myBm = t(allZ) %*% myU

myCmprime = sapply(1:(myM-1), function(x){t(myBm[,x]/sqrt(myLambda[x]))})
myQm = sapply(1:mypcmax, function(n){
    var0(myCmprime[n])/var0(myCmprime[(tailCutoff-50):tailCutoff])
  })
myPsprime = sapply(1:mypcmax, function(x){pf(myQm[x], 1, 50, lower.tail=F)})

outList = list(muprime = allZ, cmprime = myCmprime, pprime = myPsprime)
return(outList)
}





