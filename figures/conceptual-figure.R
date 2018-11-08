
library(gap)
library(qpctools)
library(dendextend)
library(LaCroixColoR)

##below is the command for the ms simulations
#system('~/Apps/msdir/ms 180 1 -t 500 -r 500 10000 -I 6 30 30 30 30 30 30 -ej 0.06 1 2 -ej 0.06 3 4 -ej 0.06 5 6 -ej 0.1 4 6 -ej 0.15 2 6  > pop3')

msout <- read.ms.output('pop3')
myn = dim(msout$gametes[[1]])[1]
causalLociIndex = sample(1:myn, 50) #pick the 50 to be causal for the phenotype
myG = msout2$gametes[[1]][-causalLociIndex,] #get the genotypes for the kinship matrix
myK = make_k_E(t(myG)) #make the kinship matrix
myEig = eigen(myK) #eigendecompose K

#make traits
beetas = as.matrix(rnorm(50), nrow=1)
causalLoci = t(msout2$gametes[[1]][causalLociIndex,])
myZ = causalLoci %*% beetas



#make the dendrogram
myPops = data.frame(pop1 = rowMeans(myG[1:1000,1:70]), pop2 = rowMeans(myG[1:1000,61:120]), pop3 = rowMeans(myG[1:1000,121:180]))

dd = dist(t(myPops), method='euclidean')
hc = hclust(dd)
hcd = as.dendrogram(hc)



mycol = lacroix_palette('Lime')
popcol =  c(rep(mycol[2], 60),rep(mycol[5], 60), rep(mycol[6], 60))


postscript("conceptfigure.eps",height=10,width=10,paper="special",horizontal=FALSE,colormodel="cymk")
par(mfrow=c(2,2), mar = c(15,5,5,2), cex.lab = 1.5, cex.axis=1.5, xpd=T)
plot(rotate(hcd, c(3,2,1)), type="rectangle", leaflab = 'none', yaxt = "n", edgePar = c(lwd=3))
points(c(1,2,3), c(0,0,0), pch = 1, lwd=8, col = mycol[c(6,5,2)], cex=5)
par(mar=c(5,5,5,2))

plot(-myEig$vectors[,1], myZ[-180], col =popcol,bty="n", xlab = "PC 1", ylab ="Trait", lwd=2, cex=2, xaxt="n", yaxt="n")
axis(1, at = c(-.1,0,.1))
axis(2, las=2)
mtext('A', side=3, adj=-0.1, cex=2, line=1)


plot(myEig$vectors[,2], myZ[-180],  col = popcol, bty="n", xlab = "PC 2",bty="n", ylab ="Trait", lwd=2, cex=2, xaxt="n", yaxt="n")
mtext('B', side=3, adj=-0.1, cex=2, line=1)
axis(1, at = c(-.15,-.05,0.05,.15))
axis(2, las=2)

plot(myEig$vectors[,3], myZ[-180],  col = popcol,  bty="n", xlab = "PC 3", ylab ="Trait",lwd=2,bty="n", cex=2, xaxt="n", yaxt="n")
mtext('C', side=3, adj=-0.1, cex=2, line=1)
axis(1, at = c(-.1,0,.1))
axis(2, las=2)
dev.off()
