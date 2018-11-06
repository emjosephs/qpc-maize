---
title: "Simulations-polygenicqpc.Rmd"
author: "em"
date: "July 5 2017"
output:
  html_document:
    keep_md: yes
---




## Doing GWAS
### Ames

```r
#as before, only run this once
#the bimbam file is too big to share on github, but I will share the GWAS results
cutoff=0.1
doGwas <- function(myI, cutoff=0.1){
  set.seed(myI)
  
  #read in file and pull out 'significant hits'
  gemma.out = processGemmaOutput(paste('data/simFiles/gemmaout.',myI,'.assoc.txt',sep=""))
  sigs = dplyr::filter(gemma.out, p_lrt < cutoff)
  
  #LD prune
  linkm = read.table('data/FileS3.csv', header=T) 
  
  #make ranges basked on LD distance
  myranges = sapply(1:10, function(j){
    mychr = dplyr::filter(linkm, chromosome==j)
    windowStarts <- mychr[seq(1, nrow(mychr), 5),]
    windowEnds <- mychr[c(seq(1, nrow(mychr)-5, 5)+5,nrow(mychr) ),]
    mywin = IRanges(start=windowStarts$position, end=windowEnds$position)
    return(mywin)})
  ldwindows = GRanges()
  for (k in 1:10){
    chrld = GRanges(k, strand = "+",myranges[[k]])
    suppressWarnings(ldwindows <- append(ldwindows,chrld))
  }
  sig.ranges = GRanges(seqname = sigs$chr, ranges = IRanges(start=sigs$ps, width=1))
  myOverlaps = as.matrix(findOverlaps(ldwindows,sig.ranges))
  sigs$myIndex = as.numeric(row.names(sigs))
  myTop = dplyr::inner_join(as.data.frame(myOverlaps), sigs, by = c("subjectHits" = "myIndex")) %>% group_by(queryHits) %>% filter(p_lrt == min(p_lrt)) %>% sample_n(1) #pull out ones with lowest pvalues
  top.nodup = myTop[!duplicated(myTop$rs),]#remove duplicates
  write.table(top.nodup$rs, file=paste('data/simFiles/ldfiltered.',myI, sep=""),quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(top.nodup, file=paste('data/simFiles/ldfiltered.assoc.',myI, sep=""),quote=FALSE, row.names=FALSE, col.names=FALSE)
  return(top.nodup)
}

mytest = sapply(1:200, doGwas)
```


Then pull out the hits from the larger dataset

```r
getSigSnps <- function(myI, sitePrefix = 'data/simFiles/ldfiltered.'){
  gwasSites = read.table(paste(sitePrefix,myI,sep=""), stringsAsFactors=F)
  system(paste('head -n 1 data/ames.281.june14.geno95.maf01.randomimpute > data/simFiles/sigSnps.ames.',myI,sep=""))
  #sapply(gwasSites$V1, function(x) {system(paste('grep ',x,' data/ames.281.june14.geno95.maf01.randomimpute >> data/simFiles/sigSnps.ames.',myI, sep=""))})
 sapply(gwasSites$V1, function(x) {system(paste('awk \'$3 == ', '"', x, '"' ,'\' data/ames.281.june14.geno95.maf01.randomimpute >> data/simFiles/sigSnps.ames.',myI, sep=""))})

  }

sapply(1:200, getSigSnps)
```



### European lines
And on the European sims


```r
cutoff=0.1
doGwas <- function(myI){
  set.seed(myI)
  
  system(paste("~/Apps/GEMMA/bin/gemma -bfile data/chip263 -k data/chip263.sXX.txt -n ",
                myI,
                "  -lmm 4 -miss 0.1 -k -o gemmaout.euro.",
                myI,
                sep=""))
}

system("mv output/* data/simFiles/")
  

pruneGwas <- function(myI, cutoff=0.1){
   
  #read in file and pull out 'significant hits'
  #filter out significant hits with awk???
  gemma.out = read.table(paste('data/simFiles/gemmaout.euro.',myI,'.assoc.txt',sep=""), header=T)
  sigs = dplyr::filter(gemma.out, p_lrt < cutoff)
  
  #LD prune
  linkm = read.table('data/FileS3.csv', header=T) 
  
  #make ranges basked on LD distance
  myranges = sapply(1:10, function(j){
    mychr = dplyr::filter(linkm, chromosome==j)
    windowStarts <- mychr[seq(1, nrow(mychr), 5),]
    windowEnds <- mychr[c(seq(1, nrow(mychr)-5, 5)+5,nrow(mychr) ),]
    mywin = IRanges(start=windowStarts$position, end=windowEnds$position)
    return(mywin)})
  ldwindows = GRanges()
  for (k in 1:10){
    chrld = GRanges(k, strand = "+",myranges[[k]])
    suppressWarnings(ldwindows <- append(ldwindows,chrld))
  }
  sig.ranges = GRanges(seqname = sigs$chr, ranges = IRanges(start=sigs$ps, width=1))
  myOverlaps = as.matrix(findOverlaps(ldwindows,sig.ranges))
  sigs$myIndex = as.numeric(row.names(sigs))
  myTop = dplyr::inner_join(as.data.frame(myOverlaps), sigs, by = c("subjectHits" = "myIndex")) %>% group_by(queryHits) %>% filter(p_lrt == min(p_lrt)) %>% sample_n(1) #pull out ones with lowest pvalues
  top.nodup = myTop[!duplicated(myTop$rs),]#remove duplicates
  #write.table(top.nodup$rs, file=paste('data/simFiles/ldfiltered.euro.',myI, sep=""),quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(top.nodup, file=paste('data/simFiles/ldfiltered.euro.assoc.',myI, sep=""),quote=FALSE, row.names=FALSE, col.names=FALSE)
  #return(top.nodup)
}

#sapply(1:1000, doGwas)

sapply(1:200, pruneGwas)
sapply(201:1000, pruneGwas)



#pull out significant snps -- do this on Farm actually
getSigSnpsEuro <- function(myI, sitePrefix = 'data/simFiles/ldfiltered.euro.assoc.'){
  gwasSites = read.table(paste(sitePrefix,myI,sep=""), stringsAsFactors=F)
  system(paste('head -n 1 data/merged263Landraces.alleles.gmri.randomimpute > data/simFiles/sigSnps.euro.',myI,sep=""))
  gwasSiteNames = sapply(gwasSites$V4, function(x){paste('s',gsub(":","_",x),sep="")})
  #sapply(gwasSites$V1, function(x) {system(paste('grep ',x,' data/ames.281.june14.geno95.maf01.randomimpute >> data/simFiles/sigSnps.ames.',myI, sep=""))})
 sapply(gwasSiteNames, function(x) {system(paste('awk \'$3 == ', '"', x, '"' ,'\' data/merged263Landraces.alleles.gmri.randomimpute >> data/simFiles/sigSnps.euro.',myI, sep=""))})

  }

sapply(1:200, getSigSnpsEuro)
sapply(201:1000, getSigSnpsEuro)

#system("rsync -avz -e 'ssh -p 2022' farm:/home/emjo/euro-maize/data/simFiles/sigSnps* simFiles/")
```


## Non-conditional test 
### Ames


```r
#do the nonconditional test using the same function used in real data (in Qpc-ames.md)
load('data/amesOnly.eig.rda')

ncamesOut = lapply(1:200, function(x){
  Qpcames_nocond(myI = x, gwasPrefix = 'data/simFiles/ldfiltered.assoc.', sigPrefix = 'data/simFiles/sigSnps.ames.', mypcmax = 100,
    myU = amesEig$vectors, myLambda = amesEig$values
  )})
save(ncamesOut, file='data/simFiles/qxpc_nonconditional_ames_200')  
```

### European lines


```r
#do noncoditional test using the same function used for the test on real data (in Qpc-euro.md)

load('data/euroOnlyK.rda')

nceuroOut <- lapply(1:200, function(x){
  Qpceuro_nocond(myI=x, gwasPrefix = 'data/simFiles/ldfiltered.euro.assoc.', sigPrefix = 'data/simFiles/sigSnps.euro.', mypcmax=100,
                 myU = euroOnlyeigen$vectors, myLambdas = euroOnlyeigen$values
                 )})

save(nceuroOut, file='data/simFiles/qxpc_nonconditional_euro_200') 
```

## Conditional test
### Ames


```r
## load data
load("data/ames.281E.K.rda")
load('data/ames.281E.condeig.rda')
ames281=myF
sigma11 = as.matrix(ames281[1:2704,1:2704])
sigma12 = as.matrix(ames281[1:2704,2705:2984])
sigma21 = as.matrix(ames281[2705:2984,1:2704])
sigma22 = as.matrix(ames281[2705:2984,2705:2984]) #we are dropping the last row
sigma.cond = sigma11 - sigma12 %*% solve(sigma22) %*% sigma21 
pcmax = 182

#run the test with same function used for real data (in Qpc-ames.md)
camesOut = lapply(1:200,function(x) {Qpcames(myI = x, gwasPrefix = 'data/simFiles/ldfiltered.assoc.', sigPrefix = 'data/simFiles/sigSnps.ames.', mypcmax = 100,
        myLambda = condEig$values, myU = condEig$vectors)})
save(camesOut, file='data/simFiles/qpc_ames_200')
```



```r
###load data
load('data/euro.282.E.rda')

myM=906
euro282 = myF
sigma11 = as.matrix(euro282[1:myM,1:myM])
sigma12 = as.matrix(euro282[1:myM,(myM+1):ncol(euro282)])
sigma21 = as.matrix(euro282[(myM+1):ncol(euro282),1:myM])
sigma22 = as.matrix(euro282[(myM+1):ncol(euro282),(myM+1):ncol(euro282)]) #we are dropping the last row
sigma.cond = sigma11 - sigma12 %*% solve(sigma22) %*% sigma21 
condEig = eigen(sigma.cond)


ceuroOut = lapply(1:200, function(x){
  Qpceuro(myI = x, gwasPrefix = 'data/simFiles/ldfiltered.euro.assoc.', sigPrefix = 'data/simFiles/sigSnps.euro.', mysigma = euro282, mypcmax = 100,
          myLambda = condEig$values, myU = condEig$vectors)
  })
save(ceuroOut, file='data/simFiles/qpc_euro_200') 
```



## Look at the results


```r
### 500 SNPs
load('data/simFiles/qxpc_nonconditional_ames_200') #ncamesOut
load('data/simFiles/qxpc_nonconditional_euro_200') #nceuroOut
load('data/simFiles/qpc_euro_200') #ceuroOut
load('data/simFiles/qpc_ames_200') #camesOut

####compare nc and c for ames and euro (with inflation factor)
nap = sapply(ncamesOut, function(x) {x$pprime}) #matrix, rows are pvals, columns are traits
cap = sapply(camesOut, function(x) {x$pprime}) #matrix, rows are pvals, columns are traits
nep = sapply(nceuroOut, function(x) {x$pprime}) #matrix, rows are pvals, columns are traits
cep = sapply(ceuroOut, function(x) {x$pprime}) #matrix, rows are pvals, columns are traits

###heatmaps
par(mar=c(4,4,2,2), xpd=TRUE, mfrow=c(1,1), mfrow=c(1,2))
mysig2 =  cut((1:1000/1000), c(0,0.001,0.01,0.05,0.1,1)) #for legend
mycol = c(viridis(6, direction=1)[1:4], "white")
image(nap, col=mycol, xaxt="n", yaxt="n", bty="l", breaks=c(0,0.001,0.01,0.05,0.1,1), ylab = "simulated trait", xlab="PC", main="Non-conditional")
legend(0,-0.15, levels(mysig2), fill=mycol, bty="n", ncol=3)
axis(1, at = c(0,0.2,0.4,0.6,0.8,1), labels=round(c(0,0.2,0.4,0.6,0.8,1)*nrow(nap)))
image(cap, col=mycol, xaxt="n", yaxt="n", bty="l", breaks=c(0,0.001,0.01,0.05,0.1,1), ylab = "simulated trait", xlab = "PC", main="Conditional")
axis(1, at = c(0,0.2,0.4,0.6,0.8,1), labels=round(c(0,0.2,0.4,0.6,0.8,1)*nrow(cap)))
```

![](Simulations-polygenicqpc_files/figure-html/gwasnocond-1.png)<!-- -->

```r
image(nep, col=mycol, xaxt="n", yaxt="n", bty="l", breaks=c(0,0.001,0.01,0.05,0.1,1), ylab = "simulated trait", xlab="PC", main="Non-conditional")
legend(0,-0.15, levels(mysig2), fill=mycol, bty="n", ncol=3)
axis(1, at = c(0,0.2,0.4,0.6,0.8,1), labels=round(c(0,0.2,0.4,0.6,0.8,1)*99+1))
image(cep[,1:200], col=mycol, xaxt="n", yaxt="n", bty="l", breaks=c(0,0.001,0.01,0.05,0.1,1), ylab = "simulated trait", xlab = "PC", main="Conditional")
axis(1, at = c(0,0.2,0.4,0.6,0.8,1), labels=round(c(0,0.2,0.4,0.6,0.8,1)*99+1))
```

![](Simulations-polygenicqpc_files/figure-html/gwasnocond-2.png)<!-- -->

```r
### bar plots of the proportion of sims that have p < 0.05
par(xpd=F, mfrow = c(2,1), mar=c(5,5,2,2))
mycol = lacroix_palette('Lime')

#get proportion of ps that are below (just first 10 PCs)
prop05 <- function(pvals){apply(pvals, 1, function(x){sum(x< 0.05)/length(x)})}

plot(-1,-1, ylim = c(0,1), xlim = c(1,30), bty="n", xlab = "PC", ylab = "Proportion significant tests", xaxt="n", yaxt = "n")
abline(h=0.05, col = mycol[6], lwd=2)
test = barplot(rbind(prop05(nap)[1:10], prop05(cap)[1:10]), beside=T, border=NA, col = mycol[c(2,5)], ylim=c(0,1), add=T)
axis(1, at = test[1,]+ 0.5, lab = 1:10, cex=1.5)
legend('topleft', c('Ames Non-conditional test', 'Ames Conditional test'), fill = mycol[c(2,5)], border="white", bty="n")
## show affect of inflation factor (compared to genic Va test), for supps

plot(-1,-1, ylim = c(0,1), xlim = c(1,30), bty="n", xlab = "PC", ylab = "Proportion significant tests", xaxt="n", yaxt = "n")
abline(h=0.05, col = mycol[6], lwd=2)
test = barplot(rbind(prop05(nep)[1:10], prop05(cep)[1:10]), beside=T, border=NA, col = mycol[c(4,6)], ylim=c(0,1), add=T)
axis(1, at = test[1,]+ 0.5, lab = 1:10, cex=1.5)
legend('topleft', c('Europe Non-conditional test', 'Europe Conditional test'), fill = mycol[c(4,6)], border="white", bty="n")
```

![](Simulations-polygenicqpc_files/figure-html/gwasnocond-3.png)<!-- -->

```r
### 50 SNPs
load('data/simFiles/qxpc_nc_ames50_200.rda') #ncamesOut
load('data/simFiles/qxpc_nc_euro50_200.rda') #nceuroOut
load('data/simFiles/qxpc_euro50_200.rda') #ceuroOut
load('data/simFiles/qxpc_ames50_200.rda') #camesOut

nap50 = sapply(ncamesOut50, function(x) {x$pprime}) #matrix, rows are pvals, columns are traits
cap50 = sapply(qxpcames_sims50, function(x) {x$pprime}) #matrix, rows are pvals, columns are traits
nep50 = sapply(nceuroOut50, function(x) {x$pprime}) #matrix, rows are pvals, columns are traits
cep50 = sapply(euroOut50, function(x) {x$pprime}) #matrix, rows are pvals, columns are traits

plot(-1,-1, ylim = c(0,1), xlim = c(1,30), bty="n", xlab = "PC", ylab = "Proportion significant tests", xaxt="n", yaxt = "n")
abline(h=0.05, col = mycol[6], lwd=2)
test = barplot(rbind(prop05(nap50)[1:10], prop05(cap50)[1:10]), beside=T, border=NA, col = mycol[c(2,5)], ylim=c(0,1), add=T)
axis(1, at = test[1,]+ 0.5, lab = 1:10, cex=1.5)
legend('topleft', c('Ames Non-conditional test 50', 'Ames Conditional test 50'), fill = mycol[c(2,5)], border="white", bty="n")

plot(-1,-1, ylim = c(0,1), xlim = c(1,30), bty="n", xlab = "PC", ylab = "Proportion significant tests", xaxt="n", yaxt = "n")
abline(h=0.05, col = mycol[6], lwd=2)
test = barplot(rbind(prop05(nep50)[1:10], prop05(cep50)[1:10]), beside=T, border=NA, col = mycol[c(4,6)], ylim=c(0,1), add=T)
axis(1, at = test[1,]+ 0.5, lab = 1:10, cex=1.5)
legend('topleft', c('Europe Non-conditional test', 'Europe Conditional test'), fill = mycol[c(4,6)], border="white", bty="n")
```

![](Simulations-polygenicqpc_files/figure-html/gwasnocond-4.png)<!-- -->

