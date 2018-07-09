#' Calculate Qpc on European landraces
#'
#' This function calculates Qpc on European landraces.
#' @param myI the trait number that we're looking at
#' @param myM the number of lines in the genotyping panel
#' @param gwasPrefix path prefix for the GWAS results
#' @param sigPrefix path prefix for the genotypes of GWAS loci in the genotyping panel
#' @param cutoff p value cutoff for including SNPs. Default is 1.
#' @param mysigma the kinship matrix for the genotyping panel.
#' @param myLambda a list of eigenvalues of the genotyping panel kinship matrix
#' @param myU a matrix of eigenvectors of the genotyping panel kinship matrix
#' @export



Qpceuro_nocond <- function(myI, myM = 906, gwasPrefix = "data/263-gwas-results/ldfiltered.assoc.", 
    sigPrefix = "data/263-gwas-results/sigSnpsEuro.", mysigma = euroOnlyF, 
    mypcmax = pcmax, myU = euroOnlyeigen$vectors, myLambdas = euroOnlyeigen$values)
  {
  
tailCutoff = round(0.9 * myM)

#read in data
gwasHits = read.table(paste(gwasPrefix,myI,sep=""), stringsAsFactors=F) #gwas results
names(gwasHits) = c('x','y',strsplit('chr     rs      ps      n_miss  allele1 allele0 af      beta    se      l_remle l_mle   p_wald  p_lrt   p_score', split=' +')[[1]])
gwasHits$locus =  sapply(gwasHits$rs, function(x){paste('s',gsub(":","_",x),sep="")})
sigGenos = read.table(paste(sigPrefix,myI, sep=""), header=T, stringsAsFactors=F) #genotypes of gwas snps in the European landraces

#combine table of GWAS results with genotypes in the GWAS set
combInfo = dplyr::left_join(sigGenos, gwasHits, by = 'locus')
combInfo$mybetas = ifelse(combInfo$allele1 == combInfo$ALT, combInfo$beta, -combInfo$beta)
myBetas = as.matrix(combInfo$mybetas)
myG = t(as.matrix(sigGenos[,6:(myM+5)]))

#center genotype matrix 
m = nrow(myG)
myT = matrix(data = -1/m, nrow = m - 1, ncol = m)
diag(myT) = (m - 1)/m
myGcent = myT %*% myG

#calculate breeding values
allZ = myGcent %*% myBetas

#project breeding values onto PCs and standardize by eigenvalue
myBm = t(allZ) %*% myU

#do Qpc
myCmprime = sapply(1:(myM-1), function(x){t(myBm[,x]/sqrt(myLambdas[x]))})
myQm = sapply(1:mypcmax, function(n){
    var0(myCmprime[n])/var0(myCmprime[(tailCutoff-50):tailCutoff])
  })
myPsprime = sapply(1:mypcmax, function(x){pf(myQm[x], 1, 50, lower.tail=F)})

outList = list(cmprime = myCmprime, pprime = myPsprime, n.sites = nrow(combInfo))
return(outList)

}






