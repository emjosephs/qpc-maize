#' Calculate Qpc on European landraces
#'
#' This function calculates Qpc on European landraces.
#' @param myI the trait number that we're looking at
#' @param myM the number of lines in the genotyping panel
#' @param gwasPrefix path prefix for the GWAS results
#' @param sigPrefix path prefix for the genotypes of GWAS loci in the genotyping panel
#' @param cutoff p value cutoff for including SNPs. Default is 1.
#' @param mysigma the kinship matrix for the combined genotyping and GWAS panels.
#' @param myLambda a list of eigenvalues of the conditional kinship matrix
#' @param myU a matrix of eigenvectors of the conditional kinship matrix
#' @export


Qpceuro <- function(myI, myM = 906, cutoff=1, gwasPrefix = 'data/263-gwas-results/ldfiltered.assoc.', 
                    sigPrefix = 'data/263-gwas-results/sigSnpsEuro.', mysigma = myF, mypcmax = pcmax,
                    myLambda = cEigValues, myU = cEigVectors){ 
  
#remove the last end of PCs 
tailCutoff = round(.9*myM)
  
#generate sigmas
sigma22 = as.matrix(mysigma[(myM + 1):dim(mysigma)[1],(myM + 1):dim(mysigma)[1]])
sigma12 = as.matrix(myF[1:myM,(myM+1):dim(mysigma)[1]])


#read in data
gwasHits = read.table(paste(gwasPrefix,myI,sep=""), stringsAsFactors=F)
names(gwasHits) = c('x','y',strsplit('chr     rs      ps      n_miss  allele1 allele0 af      beta    se      l_remle l_mle   p_wald  p_lrt   p_score scaf', split=' +')[[1]])
gwasHits$locus =  sapply(gwasHits$rs, function(x){paste('s',gsub(":","_",x),sep="")})
sigGenos = read.table(paste(sigPrefix,myI, sep=""), header=T, stringsAsFactors=F)

##filter based on p cutoff
gwasHits = dplyr::filter(gwasHits, p_lrt < cutoff)

#combine table of GWAS results with genotypes in the GWAS set
combInfo = dplyr::inner_join(sigGenos, gwasHits, by = c('locus'))
combInfo$mybetas = ifelse(combInfo$allele1 == combInfo$ALT, combInfo$beta, -combInfo$beta)
myBetas = as.matrix(combInfo$mybetas)

#center genotype matrix
myG = t(as.matrix(combInfo[,6:1174]))
#myG = t(as.matrix(sigGenos[,6:ncol(sigGenos)]))
m = nrow(myG)
myT = matrix(data = -1/m, nrow = m - 1, ncol = m)
diag(myT) = (m - 1)/m
myGcent = myT %*% myG

#calculate breeding values
allZ = myGcent %*% myBetas
z1 = allZ[1:myM]
z2 = allZ[(myM+1):length(allZ)]
#z2cent = z2 - mean(z2)
zcond = mean(allZ) + sigma12 %*% solve(sigma22) %*%  z2 #calculating the conditional prediction for Z
#zcond = zcond - mean(zcond)#center zcond
#z1 = z1 - mean(z1)

#project breeding values onto PCs and standardize
myBm = t(z1 - zcond) %*% as.matrix(myU) #z1 - zcond is the observed - expected under conditional

#do PC specific test -- here still using Va from the loci effect sizes and frequency
myCmprime = sapply(1:(myM-1), function(x){t(myBm[,x]/sqrt(myLambda[x]))})
myQm = sapply(1:pcmax, function(n){
    var0(myCmprime[n])/var0(myCmprime[(tailCutoff-50):tailCutoff])
  })
myPsprime = sapply(1:mypcmax, function(x){pf(myQm[x], 1, 50, lower.tail=F)})

outList = list(muprime = zcond, bv = z1, cmprime = myCmprime, pprime = myPsprime, n.sites = nrow(combInfo))
return(outList)

}





