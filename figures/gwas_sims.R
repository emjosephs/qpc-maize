
library(LaCroixColoR)
library(qpctools)
mycol = lacroix_palette('Mango')

load('../data/gwasstats_euro.rda')
load('../data/gwasstats_ames.rda')
load('../data/gwasstats50_ames.rda')
load('../data/gwasstats50bonf_ames.rda')
load('../data/gwasstats_euro50.rda')
#load('../data/gwasstats_euro50bonf_ames.rda')


##getting # of sig sites after ld filtering (this is hacky, fix later?)
#counts = read.table('../../polygenic-maize/data/simFiles/sigSnpsCounts', stringsAsFactors = F)
#counts$type = sapply(counts$V2, function(x){strsplit(x,'[.]')[[1]][2]})
#counts$num = sapply(counts$V2, function(x){strsplit(x,'[.]')[[1]][3]})


eurodf = data.frame(t(myGwasStatsEuro))
names(eurodf) = c('truepos','falsepos','falseneg', 'tot.hits')

eurodf50 = data.frame(t(myGwasStatsEuro50))
names(eurodf50) = c('truepos','falsepos','falseneg', 'tot.hits')

amesdf = data.frame(t(myGwasStats))
names(amesdf) = c('truepos','falsepos','falseneg', 'tot.hits')

amesdf50 = data.frame(t(myGwasStats50))
names(amesdf50) = c('truepos','falsepos','falseneg', 'tot.hits')
summary(amesdf50)

amesdf50b = data.frame(t(myGwasStats50bonf))
names(amesdf50b) = c('truepos','falsepos','falseneg', 'tot.hits')


### MAKE PLOT
postscript("gwas_sims.eps",height=8,width=16,paper="special",horizontal=FALSE,colormodel="cymk")
par(cex.lab = 1.5, mar=c(10,9,3,3), cex.axis=2, mfrow=c(1,3))

plot(jitter(c(rep(0.2,200),rep(0.4,200), rep(0.6,200),rep(0.8,200)), factor=0.5),c(amesdf$tot.hits,eurodf$tot.hits, 
                                                                                   amesdf50$tot.hits, eurodf50$tot.hits), 
     xlim = c(0,1), bty="n", ylim = c(0,1500),
     ylab = "", xlab = "", yaxt = "n", xaxt="n", col = mycol[2])
points(c(0.2, 0.4,0.6,0.8), c(mean(amesdf$tot.hits), mean(eurodf$tot.hits), mean(amesdf50$tot.hits),mean(eurodf50$tot.hits)), cex=2, lwd=4, pch=16)
axis(2, las=2)
axis(1, at = c(0.2,0.4,0.6,0.8), labels = c('Ames 500', 'Euro 500','Ames 50','Euro 50'), cex=2, las=2)
lines(x = c(0,.5), y = c(500,500), col = mycol[3], lty=2, lwd=3)
lines(x = c(.5,1), y = c(50,50), col = mycol[3], lty=2, lwd=3)
mtext('number of significant SNPs found in GWAS (with pruning)', side=2, line=6, cex=1.5)
mtext('A', side=3, adj=-0.3, cex=2, line=0)

plot(jitter(c(rep(0.2,200),rep(0.4,200), rep(0.6,200),rep(0.8,200)), factor=0.5),c(amesdf$truepos,eurodf$truepos, 
                                                                                   amesdf50$truepos, eurodf50$truepos), 
     xlim = c(0,1), bty="n", ylim = c(0,100),
     ylab = "", xlab = "", yaxt = "n", xaxt="n", col = mycol[2])
points(c(0.2, 0.4,0.6,0.8), c(mean(amesdf$truepos), mean(eurodf$truepos), mean(amesdf50$truepos),mean(eurodf50$truepos)), cex=2, lwd=4, pch=16)
axis(1, at = c(0.2,0.4,0.6,0.8), labels = c('Ames 500', 'Euro 500','Ames 50','Euro 50'), cex=2, las=2)
axis(2, las=2)
mtext('number of causal SNPs found in GWAS', side=2, line=5, cex=1.5)
lines(x = c(0,.5), y = c(500,500), col = mycol[3], lty=2, lwd=3)
lines(x = c(.5,1), y = c(50,50), col = mycol[3], lty=2, lwd=3)
mtext('B', side=3, adj=-0.3, cex=2, line=0)

plot(jitter(c(rep(0.2,200),rep(0.4,200), rep(0.6,200),rep(0.8,200)), factor=0.5),c(amesdf$truepos/amesdf$tot.hits,
                                                                                   eurodf$truepos/eurodf$tot.hits, 
                                                                                   amesdf50$truepos/amesdf50$tot.hits, 
                                                                                   eurodf50$truepos/eurodf50$tot.hits), 
     xlim = c(0,1), bty="n",
     ylab = "", xlab = "", yaxt = "n", xaxt="n", col = mycol[2])
points(c(0.2, 0.4,0.6,0.8), c(mean(amesdf$truepos/(amesdf$tot.hits)),
                                                             mean(eurodf$truepos/(eurodf$tot.hits)), 
                                                             mean(amesdf50$truepos/(amesdf50$tot.hits)), 
                                                             mean(eurodf50$truepos/(eurodf50$tot.hits))), cex=2, lwd=4, pch=16)
axis(1, at = c(0.2,0.4,0.6,0.8), labels = c('Ames 500', 'Euro 500','Ames 50','Euro 50'), cex=2, las=2)
axis(2, las=2)
mtext('Prop. of true positives that captured by a GWAS hit', side=2, line=6, cex=1.5)
mtext('C', side=3, adj=-0.3, cex=2, line=0)

dev.off()

#plot(jitter(c(rep(0.2,200),rep(0.4,200), rep(0.6,200),rep(0.8,200)), factor=0.5),c(amesdf$cor, eurodf$cor, amesdf50$cor, eurodf50$cor), xlim = c(0,1), bty="n", 
#     ylab = "", xlab = "", yaxt = "n", xaxt="n",,col = mycol[2])
#points(c(0.2, 0.4,0.6,0.8), c(mean(amesdf$cor), mean(eurodf$cor), mean(amesdf50$cor),mean(eurodf50$cor)), cex=2, lwd=4, pch=16)
#axis(1, at = c(0.2,0.4,0.6,0.8), labels = c('Ames 500', 'Euro 500','Ames 50','Euro 50'), cex=2, las=2)
#axis(2, las=2)
#mtext('Correlation between simulated and estimated effects', side=2, line=5, cex=1.5)
#mtext('C', side=3, adj=-0.1, cex=2, line=0)








