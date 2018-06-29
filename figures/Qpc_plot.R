
##read in data from Qpc-maize.Rmd
load('../data/qpc-maize_results.rda')
library(viridis)
mycols = viridis(5)[2:5]

library(LaCroixColoR)
#palette(c('darkgray',lacroix_palette('Lime')[c(3,1,4,5,6)]))
palette(c('darkgray',lacroix_palette('Pamplemousse')[c(3,1,4,5,6)]))

#plot stuff
#mycol = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", '#D55E00', '#CC79A7')
#palette(mycol[c(1,4,5,3,7,8)])
nicepops = c('mixed','non-stiff stalk','popcorn','stiff-stalk','sweet','tropical') #readable names
postscript("Qpc_results.eps",height=6,width=12,paper="special",horizontal=FALSE,colormodel="cymk")

#par(mar=c(8,15,5,3), xpd=TRUE, mfrow=c(1,3), cex.axis=1.5, cex.lab=1.5)
par(mar=c(11,15,7,3), xpd=TRUE,  cex.axis=1.5, cex.lab=1.5)
#layout(matrix(c(1,2,1,3),2,2,byrow=T))
layout(matrix(c(1,1,1,2,2,3,3, 1,1, 1, 4, 4, 4, 4), 2,7,byrow=T), heights = c(0.85,0.15))

### Plot A is a heatmap
mysig2 =  cut((1:1000/1000), c(0,0.001,0.01,0.05,0.1,1)) #for legend
pcnum = dim(allpvals)[2]
image(allpvals, col=c(mycols,'white'), xaxt="n", yaxt="n", bty="n", breaks=c(0,0.001,0.01,0.05,0.1,1), xlab = "PC")
axis(1, at=seq(0,1, length=22), label=F)
#mysigpcs = (1:36)[apply(myqvals, 1, min) < 0.1]
mysigpcs = c(1,5,10,15,20)
axis(1, at=seq(0,1, length=pcnum)[mysigpcs], label=mysigpcs, las=2)
axis(2, at=(0:21)/21, labels = mytraitsnice, las=2)
legend(-.55,-0.2, levels(mysig2), fill=c(mycols,'white'), bty="n", horiz=F, cex=1.5, ncol=3)

#for x in range
ycords = seq(0,1, length=22)
xcords = seq(0,1,length=pcnum)

for (x in 1:pcnum){ #read through PCs
  for (y in 1:22){ # read through traits
  if (myqvals[x,y] < 0.05){
    points(xcords[x],ycords[y], pch=16, col="white")
  }
  }
}

mtext('A', side=3, adj=0, cex=2, line=2)

### B  example of kernel pc1
par(mar=c(6,5,5,2), xpd=FALSE)
plot(fgmerge$X1, fgmerge$TotalKernelNo, bty="n", xlab = "PC 1", ylab = "", yaxt="n", col=as.factor(fgmerge$Subpopulation), lwd=2, ylim = c(130,420))
axis(2, las=2)
mtext("Kernel Number", side=2, line=4)
abline(lm(fgmerge$TotalKernelNo ~ fgmerge$X1), col='navy', lwd=2)
abline(a=mean(fgmerge$TotalKernelNo), b = 1.96*myCIsTKN[1], lty=2, col="#56B4E9", lwd=2)
abline(a=mean(fgmerge$TotalKernelNo), b = -1.96*myCIsTKN[1], lty=2, col="#56B4E9", lwd=2)
#par(xpd=T)
#legend(0,80, nicepops, bty="n", pch=1, pt.lwd=2, col = palette(), cex=1.5, ncol=2)
mtext('B', side=3, adj=0, cex=2, line=0)

### C example of flowering time pc 7
par(xpd=F)
plot(fgmerge$X1, fgmerge$DaysToSilk, bty="n", xlab = "PC 1", ylab = "Days to Silk", yaxt="n", col=as.factor(fgmerge$Subpopulation), lwd=2)
axis(2, las=2)
abline(lm(fgmerge$DaysToSilk ~ fgmerge$X1), col='navy', lwd=2)
abline(a=mean(fgmerge$DaysToSilk), b = 1.96*myCIsDTS[1], lty=2, col="#56B4E9", lwd=2)
abline(a=mean(fgmerge$DaysToSilk), b = -1.96*myCIsDTS[1], lty=2, col="#56B4E9", lwd=2)
#par(xpd=T)
#legend(-.13,45, nicepops, bty="n", pch=1, pt.lwd=2, col = palette(), cex=1.5, ncol=2)

mtext('C', side=3, adj=0, cex=2, line=0)

par(mar=c(0,0,0,0))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend('top', nicepops, bty="n", pch=1, pt.lwd=2, col = palette(), cex=1.5, ncol=3)

dev.off()



