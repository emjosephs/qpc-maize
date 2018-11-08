#make a figure showing how sensitive results are to choice of Va

library(viridis)
load('../data/changingVaPCresults.rda')


mycols = c(viridis(5)[2:5])
pcnum = 22

ycords = seq(0,1, length=22)
xcords = seq(0,1,length=pcnum)

mysig2 =  cut((1:1000/1000), c(0,0.001,0.01,0.05,0.1,1)) #for legend


postscript("Va-pc-sensitivity.eps",height=6,width=12,paper="special",horizontal=FALSE,colormodel="cymk")
layout(matrix(c(1,1,1,1,1,1,1,2,3,3,3,3,3,3,3,2,4,4,4,4,4,4,4,2), 8,3,byrow=F), heights = c(0.15,0.85))
par(mar=c(8,10,4,2), xpd=TRUE)

image(allpvalsAllpcs, col=c(mycols,'white'), xaxt="n", yaxt="n", bty="n", breaks=c(0,0.001,0.01,0.05,0.1,1), xlab = "PC")
axis(1, at=seq(0,1, length=pcmax), las=2, label=1:pcmax)
axis(2, at=(0:21)/21, labels = mytraitsnice, las=2)
#legend(-0.2,-0.15, levels(mysig2), fill=mycol, bty="n", horiz=T)
mtext('A', side=3, adj=-0.5, cex=2, line=1)

for (x in 1:pcnum){ #read through PCs
  for (y in 1:22){ # read through traits
    if (myqvalsAllpcs[x,y] < 0.05){
      points(xcords[x],ycords[y], pch=16, col="white")
    }
  }
}

par(mar=c(0,0,0,0), xpd=TRUE)
plot(-10,-10, xlim = c(1,0), ylim = c(1,0), bty="n", xaxt="n", yaxt = "n", xlab = "", ylab = "")
legend('center', levels(mysig2), fill=c(mycols,"white"), bty="n", horiz=F, ncol=3, cex=2)


par(mar=c(8,10,4,2), xpd=TRUE)
image(allpvals, col=c(mycols, 'white'), xaxt="n", yaxt="n", bty="n", breaks=c(0,0.001,0.01,0.05,0.1,1), xlab="PC")
axis(1, at=seq(0,1, length=pcmax), las=2, label=1:pcmax)
axis(2, at=(0:21)/21, labels = mytraitsnice, las=2)
#legend(-0.2,-0.15, levels(mysig2), fill=mycol, bty="n", horiz=T)
mtext('B', side=3, adj=-0.5, cex=2, line=1)

for (x in 1:pcnum){ #read through PCs
  for (y in 1:22){ # read through traits
    if (myqvals[x,y] < 0.05){
      points(xcords[x],ycords[y], pch=16, col="white")
    }
  }
}

image(allpvals50pcs, col=c(mycols, 'white'), xaxt="n", yaxt="n", bty="n", breaks=c(0,0.001,0.01,0.05,0.1,1), xlab = "PC")
axis(1, at=seq(0,1, length=pcmax), las=2, label=1:pcmax)
axis(2, at=(0:21)/21, labels = mytraitsnice, las=2)
#legend(-0.2,-0.15, levels(mysig2), fill=mycol, bty="n", horiz=T)
mtext('C', side=3, adj=-0.50, cex=2, line=1)


for (x in 1:pcnum){ #read through PCs
  for (y in 1:22){ # read through traits
    if (myqvals50pcs[x,y] < 0.05){
      points(xcords[x],ycords[y], pch=16, col="white")
    }
  }
}
dev.off()
