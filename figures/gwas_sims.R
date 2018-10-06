
library(LaCroixColoR)
library(qpctools)
mycol = lacroix_palette('Mango')

load('../data/gwasstats_euro.rda')
load('../data/gwasstats_ames.rda')

eurodf = data.frame(t(myGwasStatsEuro))
names(eurodf) = c('truepos','falsepos','trueneg','falseneg','cor','pval')

amesdf = data.frame(t(myGwasStats))
names(amesdf) = c('truepos','falsepos','trueneg','falseneg','cor','pval')



### MAKE PLOT
postscript("gwas_sims.eps",height=8,width=16,paper="special",horizontal=FALSE,colormodel="cymk")
par(cex.lab = 1.5, mar=c(7,9,3,3), cex.axis=2, mfrow=c(1,3))

plot(jitter(c(rep(0.3,200),rep(0.7,200)), factor=0.5),c(amesdf$truepos+amesdf$falsepos, eurodf$truepos+amesdf$falsepos), xlim = c(0,1), bty="n", 
     ylab = "", xlab = "", yaxt = "n", xaxt="n", ylim = c(0,1500), col = mycol[2])
points(c(0.3,0.7), c(mean(amesdf$truepos+amesdf$falsepos), mean(eurodf$truepos+amesdf$falsepos)), cex=2, lwd=4, pch=16)
axis(1, at = c(0.3,0.7), labels = c('Ames','Euro'), cex=2)
axis(2, las=2)
mtext('number of significant SNPs found in GWAS', side=2, line=5, cex=1.5)
mtext('A', side=3, adj=-0.1, cex=2, line=0)


plot(jitter(c(rep(0.3,200),rep(0.7,200)), factor=0.5),c(amesdf$truepos, eurodf$truepos), xlim = c(0,1), bty="n", 
     ylab = "", xlab = "", yaxt = "n", xaxt="n", ylim = c(0,25),col = mycol[2])
points(c(0.3,0.7), c(mean(amesdf$truepos), mean(eurodf$truepos)), cex=2, lwd=4, pch=16)
axis(1, at = c(0.3,0.7), labels = c('Ames','Euro'))
axis(2, las=2)
mtext('number of causal SNPs found in GWAS', side=2, line=5, cex=1.5)
mtext('B', side=3, adj=-0.1, cex=2, line=0)

plot(jitter(c(rep(0.3,200),rep(0.7,200)), factor=0.5),c(amesdf$cor, eurodf$cor), xlim = c(0,1), bty="n", 
     ylab = "", xlab = "", yaxt = "n", xaxt="n",,col = mycol[2])
points(c(0.3,0.7), c(mean(amesdf$cor), mean(eurodf$cor)), cex=2, lwd=4, pch=16)
axis(1, at = c(0.3,0.7), labels = c('Ames','Euro'))
axis(2, las=2)
mtext('Correlation between simulated and estimated effects', side=2, line=5, cex=1.5)
mtext('C', side=3, adj=-0.1, cex=2, line=0)



dev.off()




