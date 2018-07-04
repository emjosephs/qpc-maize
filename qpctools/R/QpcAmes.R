#' Calculate Qpc on Ames panel
#'
#' This function calculates Qpc the Ames panel (referred to as the genotyping panel).
#' @param myI the trait number that we're looking at
#' @param myM the number of lines in the genotyping panel
#' @param gwasPrefix path prefix for the GWAS results
#' @param sigPrefix path prefix for the genotypes of GWAS loci in the genotyping panel
#' @param cutoff p value cutoff for including SNPs. Default is 1.
#' @param mysigma the kinship matrix for the combined genotyping and GWAS panels.
#' @param myLambda a list of eigenvalues of the conditional kinship matrix
#' @param myU a matrix of eigenvectors of the conditional kinship matrix
#' @export


Qpcames <- function(myI, myM = 2704, gwasPrefix = 'data/281-gwas-results/ldfiltered.', sigPrefix='data/281-gwas-results/sigSnps.',
                    mysigma=myF, mypcmax = pcmax, myLambda = cEigValues, myU = cEigVectors){

#calculate matrices
sigma22 = as.matrix(mysigma[(myM + 1):dim(mysigma)[1], (myM + 1):dim(mysigma)[1]])
sigma12 = as.matrix(myF[1:myM, (myM + 1):dim(mysigma)[1]])
tailCutoff = round(myM*0.9)

#read in data
gwasHits = read.table(paste(gwasPrefix,myI,sep=""), stringsAsFactors=F)
names(gwasHits) = c('x','y',strsplit('chr     rs      ps      n_miss  allele1 allele0 af      beta    se      l_remle l_mle   p_wald  p_lrt   p_score', split=' +')[[1]])
sigGenos = read.table(paste(sigPrefix,myI, sep=""), header=T, stringsAsFactors=F)

#combine table of GWAS results with genotypes in the GWAS set
combData = dplyr::left_join(sigGenos, gwasHits, by = c('locus'='rs'))
myBetas = as.matrix(combData$beta)
myG = t(as.matrix(sigGenos[,4:ncol(sigGenos)]))

#center genotype matrix
m = nrow(myG)
myT = matrix(data = -1/m, nrow = m - 1, ncol = m)
diag(myT) = (m - 1)/m
myGcent = myT %*% myG

#calculate breeding values
allZ = myGcent %*% myBetas
z1 = allZ[1:myM]
z2 = allZ[(myM+1):length(allZ)]
zcond = mean(allZ) + sigma12 %*% solve(sigma22) %*%  z2 #calculating the conditional prediction for Z (j)

#project breeding values onto PCs and standardize
myBm = t(z1 - zcond) %*% as.matrix(myU) #z1 - zcond is the observed - expected under conditional

myCmprime = sapply(1:(myM-1), function(x){t(myBm[,x]/sqrt(myLambda[x]))})
myQm = sapply(1:mypcmax, function(n){
    var0(myCmprime[n])/var0(myCmprime[(tailCutoff-50):tailCutoff])
  })
myPsprime = sapply(1:mypcmax, function(x){pf(myQm[x], 1, 50, lower.tail=F)})

outList = list(muprime = zcond, cmprime = myCmprime, pprime = myPsprime)
#return the data in a reasonable way
}





