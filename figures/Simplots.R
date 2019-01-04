
##read in data from Simulations.Rmd
library(qpctools)
library(LaCroixColoR)

##load data
load('../data/simFiles/qxpc_nonconditional_ames_200') #ncamesOut
load('../data/simFiles/qxpc_nonconditional_euro_200') #nceuroOut
load('../data/simFiles/qpc_euro_200') #ceuroOut
load('../data/simFiles/qpc_ames_200') #camesOut

####compare nc and c for ames and euro (with inflation factor)
nap = sapply(ncamesOut, function(x) {x$pprime}) #matrix, rows are pvals, columns are traits
cap = sapply(camesOut, function(x) {x$pprime}) #matrix, rows are pvals, columns are traits
nep = sapply(nceuroOut, function(x) {x$pprime}) #matrix, rows are pvals, columns are traits
cep = sapply(ceuroOut, function(x) {x$pprime}) #matrix, rows are pvals, columns are traits

mycol = lacroix_palette('Lime')
prop05 <- function(pvals){apply(pvals, 1, function(x){sum(x< 0.05)/length(x)})}

#make the polygenic sore plot
load('../data/simFiles/simTraitsResults.rda')
load('../data/euro.282.E.rda')
myM=906
euro282 = myF
sigma11 = as.matrix(euro282[1:myM,1:myM])
myEig = eigen(sigma11)

#get all cors btw pc1 and sim trait
simcors = sapply(1:200, function(x){
  simTraits = myTraits[1,x]$simT[1:906]
  simTraits = simTraits - mean(simTraits)
  myCm = (simTraits %*% myEig$vectors)/sqrt(myEig$values)
  myVa = var0(myCm[(length(myEig$values)/2):length(myEig$values)])
  simcor = lm(simTraits/sqrt(myVa) ~ myEig$vectors[,1])
  return(simcor)
})

gwascors = sapply(1:200, function(x){
  gwasTraits = myTraits[2,x]$gwasT[1:906]
  gwasTraits = gwasTraits- mean(gwasTraits)
  myCm = (gwasTraits %*% myEig$vectors)/sqrt(myEig$values)
  myVa = var0(myCm[(length(myEig$values)/2):length(myEig$values)])
   gwascor = lm(gwasTraits/sqrt(myVa) ~ myEig$vectors[,1])
  return(gwascor)
})

myCm1Gwas = ((myTraits[2,1]$gwasT[1:906] - mean(myTraits[2,1]$gwasT[1:906])) %*% myEig$vectors)/sqrt(myEig$values)
myGwasTrait1 = (myTraits[2,1]$gwasT[1:906] - mean(myTraits[2,1]$gwasT[1:906]))/sqrt(var0(myCm1Gwas[(765:815)]))
myCm1Sim = ((myTraits[1,1]$simT[1:906] - mean(myTraits[1,1]$simT[1:906])) %*% myEig$vectors)/sqrt(myEig$values)
mySimTrait1 = (myTraits[1,1]$simT[1:906] - mean(myTraits[1,1]$simT[1:906]))/sqrt(var0(myCm1Sim[(765:815)]))

simVa = sapply(1:200, function(x){
  simTraits = myTraits[1,x]$simT[1:906]
  simTraits = simTraits- mean(simTraits)
  myCm = (simTraits %*% myEig$vectors)/sqrt(myEig$values)
  myVa = var0(myCm[(length(myEig$values)/2):length(myEig$values)])
  return(myVa)
  })

gwasVa = sapply(1:200, function(x){
  gwasTraits = myTraits[2,x]$gwasT[1:906]
  gwasTraits = gwasTraits- mean(gwasTraits)
  myCm = (gwasTraits %*% myEig$vectors)/sqrt(myEig$values)
  myVa = var0(myCm[(length(myEig$values)/2):length(myEig$values)])
  return(myVa)
})


postscript("Simplot.eps",height=8,width=8,paper="special",horizontal=FALSE,colormodel="cymk")
par(mfrow = c(2,2), mar=c(5,5,2,2), cex.lab=1.5)

plot(myEig$vectors[,1],mySimTrait1, bty="n", col = mycol[4], lwd=1, ylim = c(-4.5,4.5),
     xlab = "PC 1", ylab = "Simulated Trait", yaxt="n", xlim = c(-0.06,0.06))
#points(myEig$vectors[,1], mySimTrait1, col = mycol[4], lwd=1)
axis(2, las=2)
mtext('A', side=3, adj=0, cex=1.5, line=0)
#legend('bottomleft', c('Simulated traits','Polygenic scores'), bty="n", pch=1,pt.lwd=2, col = mycol[c(4,6)], cex=1.3)


plot(myEig$vectors[,1],myGwasTrait1, bty="n", col = mycol[6], lwd=1, ylim = c(-4.5,4.5),
     xlab = "PC 1", ylab = "Polygenic Score", yaxt="n", xlim = c(-0.06,0.06))
#points(myEig$vectors[,1], mySimTrait1, col = mycol[4], lwd=1)
axis(2, las=2)
mtext('B', side=3, adj=0, cex=1.5, line=0)
#legend('bottomleft', c('Simulated traits','Polygenic scores'), bty="n", pch=1,pt.lwd=2, col = mycol[c(4,6)], cex=1.3)



plot(-1,-1, ylim = c(0,1), xlim = c(1,30), bty="n", xlab = "PC", ylab = "Proportion significant tests", xaxt="n", yaxt = "n")
abline(h=0.05, col = mycol[6], lwd=2)
test = barplot(rbind(prop05(nap)[1:10], prop05(cap)[1:10]), beside=T, border=NA, col = mycol[c(2,6)], ylim=c(0,1), add=T, yaxt="n")
axis(1, at = test[1,]+ 0.5, lab = 1:10, cex=1.5)
axis(2, las=2)
legend('topleft', c('Ames Non-conditional test', 'Ames Conditional test'), fill = mycol[c(2,6)], border="white", bty="n", cex=1.3)
mtext('C', side=3, adj=0, cex=1.5, line=0)



plot(-1,-1, ylim = c(0,1), xlim = c(1,30), bty="n", xlab = "PC", ylab = "Proportion significant tests", xaxt="n", yaxt = "n")
abline(h=0.05, col = mycol[6], lwd=2)
test = barplot(rbind(prop05(nep)[1:10], prop05(cep)[1:10]), beside=T, border=NA, col = mycol[c(2,6)], ylim=c(0,1), add=T, yaxt="n")
axis(1, at = test[1,]+ 0.5, lab = 1:10, cex=1.5)
axis(2, las=2)
legend('topleft', c('Europe Non-conditional test', 'Europe Conditional test'), fill = mycol[c(2,6)], border="white", bty="n", cex=1.3)
mtext('D', side=3, adj=0, cex=1.5, line=0)

#plot(-1,-1, ylim = c(-5,5), xlim = c(-0.06,0.06), bty="n", xlab = "PC 1", yla = "Trait", col = mycol[1], yaxt="n")
#test = sapply(1:200, function(x) {abline(gwascors[1,x], col = mycol[3], lwd=1)})
#test = sapply(1:200, function(x) {abline(simcors[1,x], col = mycol[6], lwd=1)})
#axis(2, las=2)
#legend('topright',c('Simulated traits','Ascertained traits'), lwd = 2, col = mycol[c(6,3)], bty="n")
#mtext('D', side=3, adj=0, cex=2, line=0)

dev.off()

##how often is the gwas slope > the sim slope?
#gwasGreater = sapply(1:200, function(i){
#gwasslope = abs(gwascors[1,i][[1]][2])
#simslope = abs(simcors[1,i][[1]][2])
#x=0
#if(gwasslope > simslope){x = 1}
#return(x)
#})
#sum(gwasGreater)

postscript("Simplot50.eps",height=8,width=10,paper="special",horizontal=FALSE,colormodel="cymk")
par(xpd=F, mfrow = c(2,1), mar=c(5,5,2,2))

### Supp figure with the 50 snp sims
load('../data/simFiles/qxpc_nc_ames50_200.rda') #ncamesOut
load('../data/simFiles/qxpc_nc_euro50_200.rda') #nceuroOut
load('../data/simFiles/qxpc_euro50_200.rda') #ceuroOut
load('../data/simFiles/qxpc_ames50_200.rda') #camesOut

nap50 = sapply(ncamesOut50, function(x) {x$pprime}) #matrix, rows are pvals, columns are traits
cap50 = sapply(qxpcames_sims50, function(x) {x$pprime}) #matrix, rows are pvals, columns are traits
nep50 = sapply(nceuroOut50, function(x) {x$pprime}) #matrix, rows are pvals, columns are traits
cep50 = sapply(euroOut50, function(x) {x$pprime}) #matrix, rows are pvals, columns are traits

plot(-1,-1, ylim = c(0,1), xlim = c(1,30), bty="n", xlab = "PC", ylab = "Proportion significant tests", xaxt="n", yaxt = "n")
abline(h=0.05, col = mycol[6], lwd=2)
test = barplot(rbind(prop05(nap50)[1:10], prop05(cap50)[1:10]), beside=T, border=NA, col = mycol[c(2,5)], ylim=c(0,1), add=T)
axis(1, at = test[1,]+ 0.5, lab = 1:10, cex=1.5)
legend('topleft', c('Ames Non-conditional test', 'Ames Conditional test'), fill = mycol[c(2,5)], border="white", bty="n")
mtext('A', side=3, adj=0, cex=2, line=0)

plot(-1,-1, ylim = c(0,1), xlim = c(1,30), bty="n", xlab = "PC", ylab = "Proportion significant tests", xaxt="n", yaxt = "n")
abline(h=0.05, col = mycol[6], lwd=2)
test = barplot(rbind(prop05(nep50)[1:10], prop05(cep50)[1:10]), beside=T, border=NA, col = mycol[c(2,5)], ylim=c(0,1), add=T)
axis(1, at = test[1,]+ 0.5, lab = 1:10, cex=1.5)
legend('topleft', c('Europe Non-conditional test', 'Europe Conditional test'), fill = mycol[c(2,5)], border="white", bty="n")
mtext('B', side=3, adj=0, cex=2, line=0)

dev.off()



#for talks
# plot(-1,-1, ylim = c(0,1), xlim = c(1,30), bty="n", xlab = "PC", ylab = "Proportion significant tests", xaxt="n", yaxt = "n", cex.lab = 2, cex.axis=2)
# test = barplot(rbind(prop05(nep)[1:10], prop05(cep)[1:10]), beside=T, border=NA, col = c(mycol[2],'white'), ylim=c(0,1), add=T, yaxt="n")
# abline(h=0.05, col = mycol[6], lwd=2)
# axis(1, at = test[1,]+ 0.5, lab = 1:10, cex.axis=1.5)
# axis(2, las=2, cex.axis=1.5)
# 
# 
# plot(-1,-1, ylim = c(0,1), xlim = c(1,30), bty="n", xlab = "PC", ylab = "Proportion significant tests", xaxt="n", yaxt = "n", cex.lab = 2, cex.axis=2)
# test = barplot(rbind(prop05(nep)[1:10], prop05(cep)[1:10]), beside=T, border=NA, col = c(mycol[2],mycol[6]), ylim=c(0,1), add=T, yaxt="n")
# abline(h=0.05, col = mycol[6], lwd=2)
# axis(1, at = test[1,]+ 0.5, lab = 1:10, cex.axis=1.5)
# axis(2, las=2, cex.axis=1.5)
# 
# legend('topleft', c('Non-conditional test', 'Conditional test'), fill = mycol[c(2,6)], border="white", bty="n", cex=1.3)
# 

