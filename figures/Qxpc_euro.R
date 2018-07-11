
library(LaCroixColoR)
library(viridis)
library(qpctools)
library(qvalue)


##read in data from Qpc-maize.Rmd
load('../data/euro_qpc_data.rda')
load('../data/qxpc_euro_output.rda')
#cEigVectors = as.matrix(read.table('../data/euro.cond.eigenvectors', header=F, stringsAsFactors = F))
#cEigValues = as.matrix(read.table('../data/euro.cond.eigenvalues', header=F, stringsAsFactors = F))
load('../data/euro.282.condeig.rda')
pcpvalsprime = sapply(qxpceuroOut, function(x) {x$pprime})[1:17,] #matrix, rows are pvals, columns are traits
qvalsprime = get_q_values(pcpvalsprime)
myCm = sapply(qxpceuroOut, function(x){x$cmprime})
mycols = viridis(5)[2:5]
tailCutoff = round(.9*906)


##load in power sim info
prop05 <- function(pvals){apply(pvals, 1, function(x){sum(x< 0.05)/length(x)})}
load('../data/simFiles/ceuroOut_200_power01')
cpvalsprime = sapply(ceuroOut_power01, function(x) {x$pprime}) #matrix, rows are pvals, columns are traits
load('../data/simFiles/ceuroOut_200_power005')
cpvalsprime005 = sapply(ceuroOut_power005, function(x) {x$pprime}) #matrix, rows are pvals, columns are traits
load('../data/simFiles/ceuroOut_200_power05')
cpvalsprime05 = sapply(ceuroOut_power05, function(x) {x$pprime}) #matrix, rows are pvals, columns are traits



load('../data/simFiles/powerSimLatCors.rda')
mycol = lacroix_palette('Lime')
latpcols = ifelse(latps<0.05, mycol[2], mycol[5])

ycords = seq(0,1, length=dim(qvalsprime)[2])
xcords = seq(0,1,length=dim(qvalsprime)[1])


#plot stuff

postscript("Qxpc_euro.eps",height=8,width=16,paper="special",horizontal=FALSE,colormodel="cymk")

par(mar=c(10,14,5,3), xpd=TRUE, mfrow=c(1,1), cex.axis=1.5, cex.lab=2)
layout(
  matrix(c(1,1,1,2,2,2,3,4,4),3,3,byrow=F)
  )

### Plot A is a heatmap
mysig2 =  cut((1:1000/1000), c(0,0.001,0.01,0.05,0.1,1)) #for legend
image(pcpvalsprime, col=c(mycols,'white'), xaxt="n", yaxt="n", bty="n", breaks=c(0,0.001,0.01,0.05,0.1,1), xlab = "PC")
#sigPCs = c(1:18)[apply(pcpvalsprime, 1, function(x){min(x)<0.05*18})]
axis(1, at=seq(0,1, length=nrow(pcpvalsprime)), label=1:nrow(pcpvalsprime))
#axis(1, cex.axis=1.5, label=sigPCs, las=2)
axis(2, at=(0:21)/21, labels = niceTraitnames, las=2)
legend(-0.6,-0.13, levels(mysig2), fill=c(mycols,'white'), bty="n", horiz=F, cex=1.8, ncol=3, text.width = 0.4)
mtext('A', side=3, adj=0, cex=2, line=1)

for (x in 1:17){ #read through PCs
  for (y in 1:22){ # read through traits
  if (qvalsprime[x,y] < 0.05){
    points(xcords[x],ycords[y], pch=16, col="white")
  }
  }
}


#mycol = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", '#D55E00', '#CC79A7')
#palette(mycol[c(2,6)])
palette(mycol[c(1,2)])

### B 


par(mar=c(10,5,5,3), xpd=FALSE)
myI=22
myVa = var0(myCm[(tailCutoff-50):tailCutoff, myI])
myCIs = sapply(1:21, function(x){1.96*sqrt(myVa*cEigValues[x])})
bonfP = 0.05/(18*22)
bonfZ = qnorm(bonfP/2, lower.tail=F)
bonfCIs = sapply(1:21, function(x){bonfZ*sqrt(myVa*cEigValues[x])})

plot(cEigVectors[,1], mymerge[,niceTraitnames[myI]], bty="n", xlab = "Conditional PC 1", ylab = expression(paste('Brace root X - ',mu,',', sep="")),
  col=mycol[2], lwd=2, main="", yaxt="n")
axis(2, las=2)
abline(lm(mymerge[,niceTraitnames[myI]]~cEigVectors[,1]), lwd=3, col= mycol[6])
abline(a=mean(mymerge[,niceTraitnames[myI]]), b = myCIs[1], lty=2, col=mycol[4], lwd=3)
abline(a=mean(mymerge[,niceTraitnames[myI]]), b = -myCIs[1], lty=2, col=mycol[4], lwd=3)
abline(a=mean(mymerge[,niceTraitnames[myI]]), b = bonfCIs[1], lty=2, col=mycol[5], lwd=3)
abline(a=mean(mymerge[,niceTraitnames[myI]]), b = -bonfCIs[1], lty=2, col=mycol[5], lwd=3)
par(xpd=T)
#legend(-0.055, -6.1, levels(as.factor(mymerge$Type)), col = palette()[1:2], pch=1, pt.lwd=2, bty="n", cex=1.8, horiz=T, text.width = .04)
#legend(-0.06, -6.5, c('slope', '95% CI'), col = mycol[c(6,5)], lwd=3, lty = c(1,2), bty="n", cex=1.8, horiz=T, text.width = .03)
legend(-0.06, -6, c('slope', '95% CI', '99.99% CI'), col = mycol[c(6,4, 5)], lwd=3, lty = c(1,2,2), bty="n", cex=1.8, ncol=2, text.width = .02)

mtext('B', side=3, adj=0, cex=2, line=1)
### add in neutral expectation


### C 

par(mar = c(4,5,5,1))

plot(abs(latcors)[1:5], bty="n", xlab = "PC", ylab = "R", xaxt="n", yaxt = "n", col = mycol[5], cex=2, lwd=4, ylim = c(0,0.8), xlim = c(0.6,5.4))
axis(1, lab = 1:5, at=1:5)
axis(2, las=2)
mtext('C', side=3, adj=0, cex=2, line=1)

par(mar=c(10,5,4,1))
plot(-1,-1, ylim = c(0,1.1), xlim = c(1,20), bty="n", xlab = "PC", ylab = "Proportion significant tests", xaxt="n", yaxt = "n")
test = barplot(rbind(prop05(cpvalsprime005)[1:5],prop05(cpvalsprime)[1:5], prop05(cpvalsprime05)[1:5]), beside=T, col = mycol[c(2,4,5)], border = mycol[c(2,4,5)], add=T, yaxt="n")
axis(1, at = test[1,]+ 1, lab = 1:5)
axis(2, at = (0:5)/5,las=2)
mtext('D', side=3, adj=0, cex=2, line=-2)
legend(0.5, -0.28, c(expression(paste(alpha, '= 0.005')), expression(paste(alpha, '= 0.01')), expression(paste(alpha, '= 0.05'))),fill = mycol[c(2,4,5)], 
       border=NA,cex=2, bty="n", horiz=T, text.width=3.8)


dev.off()



