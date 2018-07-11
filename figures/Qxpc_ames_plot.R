
##read in data from Qxpc-ames.

library(viridis)
library(qpctools)
library(qvalue)
load("../data/ames_qpc_data.rda")
ncpvals = sapply(ncamesOut, function(x) {x$pprime}) #matrix, rows are pvals, columns are traits
ncqvals = get_q_values(ncpvals)
pcpvals = sapply(qpcamesOut, function(x) {x$pprime}) #matrix, rows are pvals, columns are traits
qvals = get_q_values(pcpvals)


mycol = c(viridis(5)[2:5], 'white')
mysig2 =  cut((1:1000/1000), c(0,0.001,0.01,0.05,0.1,1)) #for legend

ycords = seq(0,1, length=dim(qvals)[2])
xcords = seq(0,1,length=dim(qvals)[1])


#plot stuff
postscript("Qxpc_ames_plot.eps",height=7,width=10,paper="special",horizontal=FALSE,colormodel="cymk")

layout(matrix(c(1,2,1,2,1,2,1,2,1,2,1,2,3,3),7,2,byrow=T))
par(mar=c(8,14,4,0), xpd=FALSE)

###non-conditional ames heatmap
image(ncpvals, col=mycol, xaxt="n", yaxt="n", breaks=c(0,0.001,0.01,0.05,0.1,1))
axis(1, at = c(0,0.2,0.4,0.6,0.8,1), labels=round(c(0,0.2,0.4,0.6,0.8,1)*nrow(allqvals)))
axis(2, at=(0:21)/21, labels = niceTraitnames, las=2, cex = 5)
mtext('A', side=3, adj=0, cex=2, line=1, at=0)

for (x in 1:dim(ncqvals)[1]){ #read through PCs
  for (y in 1:dim(ncqvals)[2]){ # read through traits
    if (ncqvals[x,y] < 0.05){
      points(xcords[x],ycords[y], pch=16, col="black")
    }
  }
}

####conditional ames heatmap
par(mar=c(8,4,4,10))
image(pcpvals, col=mycol, xaxt="n", yaxt="n", breaks=c(0,0.001,0.01,0.05,0.1,1))
axis(1, at = c(0,0.2,0.4,0.6,0.8,1), labels=round(c(0,0.2,0.4,0.6,0.8,1)*nrow(allqvals)))
axis(2, at=(0:21)/21, labels = FALSE)
#legend(-0.2,-0.15, levels(mysig2), fill=mycol, bty="n", horiz=T)

mtext('B', side=3, adj=0, cex=2, line=1, at = 0)

for (x in 1:dim(ncqvals)[1]){ #read through PCs
  for (y in 1:dim(ncqvals)[2]){ # read through traits
    if (qvals[x,y] < 0.05){
      points(xcords[x],ycords[y], pch=16, col="black")
    }
  }
}

par(mar=c(1,1,1,1))
plot.new()
legend("bottomleft", c('FDR < 0.05'), pch=16, col =c("black"), cex=2, bty="n", horiz=T)
legend("topleft", levels(mysig2), fill=mycol, bty="n", horiz=T, cex=2)

dev.off()


