---
title: "conceptual_figure"
author: "em"
date: "September 24, 2018"
output:
  html_document:
    keep_md: yes
---







3 populations


```r
#i=1
#x1 <- runif(1) #getting a random number so there's a seed
#save(".Random.seed", file=paste("data/figure_sims/randomSeed.",i, sep=""))
load('data/figure_sims/randomSeed.1')

npops=3
Faa = 0.15

sigma = matrix(0,nrow=3, ncol=3)
sigma[1:2,1:2] = matrix(Faa/2, nrow=2, ncol=2)
diag(sigma) = Faa
library(viridis)
```

```
## Loading required package: viridisLite
```

```r
heatmap(sigma, col = viridis(4))
```

![](conceptual-figure_files/figure-html/3pops-1.png)<!-- -->

```r
#simulate allele freqs in two pops at these loci
nloci = 500
ancPop = runif(nloci, min=0, max=1)
presentPops1 = sapply(ancPop, function(x){mvrnorm(n=1, mu = rep(x,npops), x*(1-x)*sigma)})
presentPops = apply(presentPops1, c(1,2), myBound) #deal with numbers greater or less than 0 (the outer bounds are sticky)

#plot(ancPop, bty='n', xlab = "locus", ylab = "p")
#sapply(1:500, function(x){lines(c(x,x), c(presentPops[1,x], presentPops[2,x]))})

#get the population genotypes
npop = 50
popGenos = lapply(1:npops, function(x) getPopGenos(x, presentPops, npop)) #a list of elements, each is a population

#make a kinship matrix
myG = rbind(popGenos[[1]], popGenos[[2]], popGenos[[3]])/2
myK = make_k(myG)
heatmap(myK)
```

![](conceptual-figure_files/figure-html/3pops-2.png)<!-- -->

```r
myEig = eigen(myK)

mycol = lacroix_palette('Mango')
plot(myEig$values, lwd=2, bty="n", col = mycol[6], xlab = "PC", ylab = "Eigenvalue")
```

![](conceptual-figure_files/figure-html/3pops-3.png)<!-- -->

```r
#plot(myEig$vectors[,1], myEig$vectors[,2], col = c(rep(mycol[1], 50),rep(mycol[2], 50), rep(mycol[3], 50), rep(mycol[4], 50)), bty='n', xlab= 'PC1', ylab = 'PC2', lwd=2, xlim = c(-.2, .2))
#legend('topright', c('pop1','pop2', 'pop3','pop4'), bty="n", pch=1, pt.lwd=2, col = mycol)

#plot(myEig$vectors[,3], myEig$vectors[,4], col = c(rep(mycol[1], 50),rep(mycol[2], 50), rep(mycol[3], 50), rep(mycol[4], 50)), bty='n', xlab= 'PC3', ylab = 'PC4', lwd=2, xlim = c(-.2, .2))
#legend('topright', c('pop1','pop2', 'pop3','pop4'), bty="n", pch=1, pt.lwd=2, col = mycol)

beetas = matrix(c(rnorm(nloci)), ncol=1, nrow=nloci) 
popPhenos = getPopPhenos(popGenos, beetas)
#individual noise shifts
#popPhenosNoise = lapply(popPhenos, function(x){x + rnorm(length(x), mean=0, sd=sd(x)/2 )})

#plot(popPhenos[[1]], popPhenosNoise[[1]])

nind = npops*npop - 1
#myPhenos = unlist(popPhenosNoise)
myPhenos = unlist(popPhenos) - mean(unlist(popPhenos))
plot(myEig$vectors[,1], myPhenos[1:nind], col = c(rep(mycol[1], 50),rep(mycol[2], 50), rep(mycol[3], 50)), bty="n", xlab = "PC1", ylab ="Sim phenos", lwd=2)
myl = lm(myPhenos[1:nind]~ myEig$vectors[,1])
abline(myl)
```

![](conceptual-figure_files/figure-html/3pops-4.png)<!-- -->

```r
plot(myEig$vectors[,2], myPhenos[1:nind], col = c(rep(mycol[1], 50),rep(mycol[2], 50), rep(mycol[3], 50), rep(mycol[4], 50)), bty="n", xlab = "PC1", ylab ="Sim phenos", lwd=2)
myl = lm(myPhenos[1:nind]~ myEig$vectors[,2])
abline(myl)
```

![](conceptual-figure_files/figure-html/3pops-5.png)<!-- -->

```r
plot(myEig$vectors[,3], myPhenos[1:nind], col = c(rep(mycol[1], 50),rep(mycol[2], 50), rep(mycol[3], 50), rep(mycol[4], 50)), bty="n", xlab = "PC1", ylab ="Sim phenos",lwd=2)
myl = lm(myPhenos[1:nind]~ myEig$vectors[,3])
abline(myl)
```

![](conceptual-figure_files/figure-html/3pops-6.png)<!-- -->


Now to show how we're estimating Va with the lower pcs


```r
## estimate Va with the loci and allele freqs
myPs =  colMeans(presentPops) #get allele freq across all pops (could do with anc pop too)
myVaGeno = calcVa(myPs, beetas)
myVaGeno
```

```
## [1] 158.4759
```

```r
myVaAnc = calcVa(ancPop, beetas)
myVaAnc
```

```
## [1] 170.628
```

```r
#estimate Va with all PCs
ninds = npops*npop - 1
myCs = sapply(1:ninds, function(x){
  (myPhenos[1:ninds] %*% myEig$vectors[,x])/sqrt(myEig$values[x])
}) #project traits onto PCs


myVaAll = var0(myCs[1:(nind-1)])/2
myVaAll 
```

```
## [1] 169.7287
```

```r
#plot Va estimates
myVaPC = sapply(myCs[1:(nind-1)], var0)/2
hist(myVaPC, col = mycol[4], border="white", main = "") ##this is the right distribution, right?
```

![](conceptual-figure_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

```r
plot(myVaPC, bty="n", lwd=2, col = mycol[4], xlab = 'PC', ylab = "Va from one PC")
abline(h=myVaAll, col = mycol[3], lwd=2)
```

![](conceptual-figure_files/figure-html/unnamed-chunk-1-2.png)<!-- -->

```r
summary(myVaPC)
```

```
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##    0.0016   12.7532  105.6428  169.7287  223.4823 1334.2860
```


Do simulations

```r
i=2

doSim <- function(i){
x1 <- runif(i) #getting a random number so there's a seed
save(".Random.seed", file=paste("data/figure_sims/randomSeed.",i, sep=""))

npops=3
Faa = 0.15

sigma = matrix(0,nrow=3, ncol=3)
sigma[1:2,1:2] = matrix(Faa/2, nrow=2, ncol=2)
diag(sigma) = Faa

#simulate allele freqs in two pops at these loci
nloci = 500
ancPop = runif(nloci, min=0, max=1)
presentPops1 = sapply(ancPop, function(x){mvrnorm(n=1, mu = rep(x,npops), x*(1-x)*sigma)})
presentPops = apply(presentPops1, c(1,2), myBound) #deal with numbers greater or less than 0 (the outer bounds are sticky)

#get the population genotypes
npop = 50
popGenos = lapply(1:npops, function(x) getPopGenos(x, presentPops, npop)) #a list of elements, each is a population

#make a kinship matrix
myG = rbind(popGenos[[1]], popGenos[[2]], popGenos[[3]])/2
myK = make_k(myG)

myEig = eigen(myK)

beetas = matrix(c(rnorm(nloci)), ncol=1, nrow=nloci) 
popPhenos = getPopPhenos(popGenos, beetas)

nind = npops*npop - 1
myPhenos = unlist(popPhenos) - mean(unlist(popPhenos))

myPs =  colMeans(presentPops) #get allele freq across all pops (could do with anc 
myVaGeno = calcVa(myPs, beetas)
myVaAnc = calcVa(ancPop, beetas)

#estimate Va with all PCs
ninds = npops*npop - 1
myCs = sapply(1:ninds, function(x){
  (myPhenos[1:ninds] %*% myEig$vectors[,x])/sqrt(myEig$values[x])
}) #project traits onto PCs


myVaAll = var0(myCs[1:(nind-1)])/2
myVaPC = sapply(myCs[1:(nind-1)], var0)/2

return(list(vapc = myVaPC, vaanc = myVaAnc, vageno = myVaGeno, vaall = myVaAll))
}

vaSims = sapply(2:202, doSim)
save(vaSims, file = "data/figure_sims/vaSimresults")
```

Look at sims

```r
load('data/figure_sims/vaSimresults')

#plot out the mean Va for each PC across all PCs, divided by what you'd get for total
allVaPCs = rowMeans(sapply(1:200, function(x){vaSims[,x]$vapc/vaSims[,x]$vaall}))

plot(allVaPCs)
```

![](conceptual-figure_files/figure-html/simslook-1.png)<!-- -->

So what does the final figure look like??


```r
postscript("conceptfigure.eps",height=5,width=10,paper="special",horizontal=FALSE,colormodel="cymk")
par(mfrow=c(1,4), mar = c(5,5,2,2), cex.lab = 1.5, cex.axis=1.5)
myPhenos = unlist(popPhenos) - mean(unlist(popPhenos))
plot(myEig$vectors[,1], myPhenos[1:nind], col = c(rep(mycol[1], 50),rep(mycol[2], 50), rep(mycol[3], 50)), bty="n", xlab = "PC1", ylab ="Trait", lwd=2)
myl = lm(myPhenos[1:nind]~ myEig$vectors[,1])
abline(myl, lwd=2, col = mycol[5])


plot(myEig$vectors[,2], myPhenos[1:nind], col = c(rep(mycol[1], 50),rep(mycol[2], 50), rep(mycol[3], 50), rep(mycol[4], 50)), bty="n", xlab = "PC1", ylab ="Trait", lwd=2)
myl = lm(myPhenos[1:nind]~ myEig$vectors[,2])
abline(myl, lwd=2, col = mycol[5])


plot(myEig$vectors[,3], myPhenos[1:nind], col = c(rep(mycol[1], 50),rep(mycol[2], 50), rep(mycol[3], 50), rep(mycol[4], 50)), bty="n", xlab = "PC1", ylab ="Trait",lwd=2)
myl = lm(myPhenos[1:nind]~ myEig$vectors[,3])
abline(myl, lwd=2, col = mycol[5])


plot(allVaPCs, bty="n", lwd=2, col = mycol[4], xlab = "PC", ylab = "Estimated Va")

dev.off()
```

TODO add CIs
