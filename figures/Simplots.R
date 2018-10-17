
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


postscript("Simplot.eps",height=8,width=10,paper="special",horizontal=FALSE,colormodel="cymk")
par(xpd=F, mfrow = c(2,1), mar=c(5,5,2,2))


plot(-1,-1, ylim = c(0,1), xlim = c(1,30), bty="n", xlab = "PC", ylab = "Proportion significant tests", xaxt="n", yaxt = "n")
abline(h=0.05, col = mycol[6], lwd=2)
test = barplot(rbind(prop05(nap)[1:10], prop05(cap)[1:10]), beside=T, border=NA, col = mycol[c(2,5)], ylim=c(0,1), add=T, yaxt="n")
axis(1, at = test[1,]+ 0.5, lab = 1:10, cex=1.5)
axis(2, las=2)
legend('topleft', c('Ames Non-conditional test', 'Ames Conditional test'), fill = mycol[c(2,5)], border="white", bty="n")
mtext('A', side=3, adj=0, cex=2, line=0)



plot(-1,-1, ylim = c(0,1), xlim = c(1,30), bty="n", xlab = "PC", ylab = "Proportion significant tests", xaxt="n", yaxt = "n")
abline(h=0.05, col = mycol[6], lwd=2)
test = barplot(rbind(prop05(nep)[1:10], prop05(cep)[1:10]), beside=T, border=NA, col = mycol[c(2,5)], ylim=c(0,1), add=T, yaxt="n")
axis(1, at = test[1,]+ 0.5, lab = 1:10, cex=1.5)
axis(2, las=2)
legend('topleft', c('Europe Non-conditional test', 'Europe Conditional test'), fill = mycol[c(2,5)], border="white", bty="n")
mtext('B', side=3, adj=0, cex=2, line=0)


dev.off()

### Supp figure with the 50 snp sims