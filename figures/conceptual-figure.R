
library(MASS)
library(quaint)
library(LaCroixColoR)

myBound = function(x){
if (x>1) {x=1}
if (x<0) {x=0}
  return(x)
}

calcVa = function (afs, betas) 
{
    return(sum(2 * afs * (1 - afs) * (betas^2)))
}

getPopGenos <- function(p, pops, popSize){ #turns allele frequencies into individual genotypes, p is row, pops is the matrix of frequencies, popSize is the # of individuals per pop
  indFreqs = t(matrix(rep(pops[p,],popSize),nrow=nloci,ncol=popSize))
  popGenos = apply(indFreqs,c(1,2),function(x){sum(sample(c(1,0),2,replace=TRUE, prob = c(x,1-x)))})
  return(popGenos)
}

getPopPhenos <- function(pg,b){
  pp = lapply(pg, function(x){x %*% b}) 
  return(pp)}

getBetas <- function(res, colu, loci){  #get betas including 0s that weren't tested
  labs = data.frame(1:loci)
  names(labs) = 'lab'
  res$lab = as.numeric(substr(res$rs, 3, 100))
  allB = merge(labs,res, by="lab", all=TRUE)
  betas = allB[,colu]
  betas[is.na(betas)] <- 0
  return(betas)
}

#i=1
#x1 <- runif(1) #getting a random number so there's a seed
#save(".Random.seed", file=paste("data/figure_sims/randomSeed.",i, sep=""))
load('../data/figure_sims/randomSeed.1')

npops=3
Faa = 0.15

sigma = matrix(0,nrow=3, ncol=3)
sigma[1:2,1:2] = matrix(Faa/2, nrow=2, ncol=2)
diag(sigma) = Faa
sigma2 = matrix(0, nrow=9, ncol=9)

for (x in 0:2){
 sigma2[(1:3)+x*3,(1:3)+x*3] = sigma + sigma[x+1,x+1]
  }
sigma2[4:6,1:3] = 0.075
sigma2[1:3,4:6] = 0.075




#simulate allele freqs in pops at these loci
nloci = 500
ancPop = runif(nloci, min=0, max=1)
presentPops1 = sapply(ancPop, function(x){mvrnorm(n=1, mu = rep(x,npops*3), x*(1-x)*sigma2)})
presentPops = apply(presentPops1, c(1,2), myBound) #deal with numbers greater or less than 0 (the outer bounds are sticky)

#get the population genotypes
npop = 30
popGenos = lapply(1:(npops*3), function(x) getPopGenos(x, presentPops, npop)) #a list of elements, each is a population

#make a kinship matrix with the last 400 sites
myG = do.call(rbind, lapply(popGenos, function(x){x[,101:500]}))
myK = make_k(myG/2)

myEig = eigen(myK)

beetas = matrix(c(rnorm(100), rep(0, 400)), ncol=1, nrow=nloci) 
popPhenos = getPopPhenos(popGenos, beetas)

nind = 269


library(dendextend)
nloci = 500
mycol = lacroix_palette('Lime')[c(2,5,6)]
ancPopd = runif(nloci, min=0, max=1)
presentPops3 = sapply(ancPopd, function(x){mvrnorm(n=1, mu = rep(x,npops), x*(1-x)*sigma)})
dd = dist(presentPops3, method='euclidean')
hc = hclust(dd)
hcd = as.dendrogram(hc)

myPhenos = unlist(popPhenos) - mean(unlist(popPhenos))

postscript("conceptfigure.eps",height=5,width=10,paper="special",horizontal=FALSE,colormodel="cymk")
#png('conceptfigure.png', height=250, width=500)
par(mfrow=c(1,4), mar = c(15,5,5,2), cex.lab = 1.5, cex.axis=1.5, xpd=T)

plot(rotate(hcd, c(3,2,1)), type="rectangle", leaflab = 'none', yaxt = "n", edgePar = c(lwd=3))
points(c(1,2,3), c(0,0,0), pch = 1, lwd=8, col = mycol, cex=5)

par(mar=c(5,5,5,2))

#heatmap(sigma, col = "white", rowv = NULL, labCol = "", labRow = "", Rowv= NA)



plot(-myEig$vectors[,1], myPhenos[1:nind], col = c(rep(mycol[1], 90),rep(mycol[2], 90), rep(mycol[3], 90)), ylim = c(-20,30),bty="n", xlab = "PC 1", ylab ="Trait", lwd=2)
#myl = lm(myPhenos[1:nind]~ myEig$vectors[,1])
#legend('topleft', c('Pop 1','Pop 2','Pop 3'), fill = mycol, border="white", bty="n", cex=1.5)
mtext('A', side=3, adj=-0.1, cex=2, line=1)

#abline(myl, lwd=2, col = mycol[5])


plot(myEig$vectors[,2], myPhenos[1:nind],  col = c(rep(mycol[1], 90),rep(mycol[2], 90), rep(mycol[3], 90)), bty="n", xlab = "PC 2", ylim = c(-20,30),bty="n", ylab ="Trait", lwd=2)
#myl = lm(myPhenos[1:nind]~ myEig$vectors[,2])
mtext('B', side=3, adj=-0.1, cex=2, line=1)

#abline(myl, lwd=2, col = mycol[5])


plot(myEig$vectors[,3], myPhenos[1:nind],  col = c(rep(mycol[1], 90),rep(mycol[2], 90), rep(mycol[3], 90)),  bty="n", xlab = "PC 3", ylab ="Trait",lwd=2,
ylim = c(-20,30),bty="n")
mtext('C', side=3, adj=-0.1, cex=2, line=1)

#abline(myl, lwd=2, col = mycol[5])


dev.off()


