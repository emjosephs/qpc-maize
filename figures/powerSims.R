
##read in data from Simulations.Rmd
library(qpctools)
library(LaCroixColoR)

mycol = lacroix_palette('Lime')
prop05 <- function(pvals){apply(pvals, 1, function(x){sum(x< 0.05)/length(x)})}

load('../data/simFiles/ceuroOut_200_power')
cpvalsprime = sapply(ceuroOut_power, function(x) {x$pprime}) #matrix, rows are pvals, columns are traits

load('../data/simFiles/ceuroOut_200_power005')
cpvalsprime005 = sapply(ceuroOut_power005, function(x) {x$pprime}) #matrix, rows are pvals, columns are traits

load('../data/simFiles/ceuroOut_200_power05')
cpvalsprime05 = sapply(ceuroOut_power05, function(x) {x$pprime}) #matrix, rows are pvals, columns are traits



load('../data/simFiles/powerSimLatCors.rda')
latpcols = ifelse(latps<0.05, mycol[2], mycol[5])

postscript("powerSims.eps",height=8,width=10,paper="special",horizontal=FALSE,colormodel="cymk")


layout(matrix(c(1,1,1,1,1,1,2,2,2), 3, 3, byrow = TRUE))
par(cex.lab = 2, mar = c(6,6,3,3), cex.axis=1.5)

plot(-1,-1, ylim = c(0,1.2), xlim = c(1,20), bty="n", xlab = "PC", ylab = "Proportion significant tests", xaxt="n", yaxt = "n")
test = barplot(rbind(prop05(cpvalsprime005)[1:5],prop05(cpvalsprime)[1:5], prop05(cpvalsprime05)[1:5]), beside=T, col = mycol[c(2,4,5)], border = mycol[c(2,4,5)], add=T, yaxt="n")
axis(1, at = test[1,]+ 1, lab = 1:10)
axis(2, at = (0:5)/5,las=2)
mtext('A', side=3, adj=0, cex=2, line=0)
legend('topleft', c(expression(paste(alpha, '= 0.005')), expression(paste(alpha, '= 0.01')), expression(paste(alpha, '= 0.05'))),fill = mycol[c(2,4,5)], 
       border=NA,cex=2, bty="n", horiz=T, text.width=4)


plot(abs(latcors)[1:5], bty="n", xlab = "PC", ylab = "abs(Correlation)", xaxt="n", yaxt = "n", col = mycol[5], cex=2, lwd=4, ylim = c(0,0.8))
axis(1, lab = 1:21, at=1:21)
axis(2, las=2)
mtext('B', side=3, adj=0, cex=2, line=0)


dev.off()

