
library(LaCroixColoR)
library(qpctools)
mycol = lacroix_palette('Lime')

load('../data/simFiles/Qpc200.rda')


## read and sort data
allCms = sapply(neutQpcVe0, function(x) {x$cm}) #matrix where columns are the Cm values for a given simulation (from PC1 to PC 239)
allCmVars = apply(allCms, 1, var0) #get variance across all 200 sims for each pc
allPs = sapply(neutQpcVe0, function(x){x$pvals})
allCmsVe.1 = sapply(neutQpcVe.1, function(x) {x$cm})
allCmVarsVe.1 = apply(allCmsVe.1, 1, var0) #get variance across all 200 sims for each pc
allPs.1 = sapply(neutQpcVe.1, function(x){x$pvals})
allCmsVe.5 = sapply(neutQpcVe.5, function(x) {x$cm})
allCmVarsVe.5 = apply(allCmsVe.5, 1, var0) #get variance across all 200 sims for each pc
allPs.5 = sapply(neutQpcVe.5, function(x){x$pvals})

#get prop that are < 0.05
prop05 <- function(pvals){apply(pvals, 1, function(x){sum(x< 0.05)/length(x)})}
propPs = prop05(allPs)

#estimate Va
allVa = apply(allCms, 2, function(x){var0(x[(tailCutoff-50):tailCutoff])})
allVa.1 = apply(allCmsVe.1, 2, function(x){var0(x[(tailCutoff-50):tailCutoff])})
allVa.5 = apply(allCmsVe.5, 2, function(x){var0(x[(tailCutoff-50):tailCutoff])})

### MAKE PLOT
postscript("Ve_sims.eps",height=8,width=16,paper="special",horizontal=FALSE,colormodel="cymk")
par(mfrow = c(1,2), cex.lab = 1.5, cex.axis=1.5, mar = c(6,6,4,2))
#var of cms
plot(-10,-10, bty="n", xlab = 'PC', ylab = "", xlim = c(0, tailCutoff), ylim=c(0,3.5), yaxt = "n")
rect(xleft = tailCutoff-50, ybottom = 0, xright= tailCutoff, ytop = 4,col = mycol[3], border="NA")
axis(2, las=2)
mtext("Var(Cm)", side=2, line=4, cex=1.5)
points(allCmVarsVe.5[1:tailCutoff], col = mycol[4], lwd=2)
points(allCmVarsVe.1[1:tailCutoff], col = mycol[5], lwd=2)
points(allCmVars[1:tailCutoff], col = mycol[6], lwd=2)
legend('topleft', c('Ve=0','Ve=Va/10','Ve=Va/2'), col = mycol[c(6,5,4)], pch=1, bty="n", pt.lwd=3, cex=1.5)
mtext('A', side=3, adj=-.15, cex=2, line=0)



#sig tests
plot(-1,-1, ylim = c(0,0.1), xlim = c(1,20), bty="n", xlab = "PC", ylab = "", xaxt="n", yaxt = "n")
axis(2, las=2)
abline(h=0.05, col = mycol[6], lwd=2)
test = barplot(rbind(prop05(allPs)[1:5], prop05(allPs.1)[1:5], prop05(allPs.5)[1:5]), beside=T, border = mycol[c(6,5,4)], 
               col =  mycol[c(6,5,4)], yaxt = "n",ylim=c(0,1), add=T)
mtext("Proportion significant tests", side=2, line=4, cex=1.5)
axis(1, at = test[2,], lab = 1:5, cex=1.5)
legend('topleft', c('Ve=0','Ve=Va/10','Ve=Va/2'), fill = mycol[c(6,5,4)], border="white", bty="n", cex=1.5)
mtext('B', side=3, adj=-.15, cex=2, line=0)

dev.off()




