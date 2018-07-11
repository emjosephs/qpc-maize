
library(LaCroixColoR)
library(qpctools)
library(qvalue)

## AMES DATA
load("../data/ames_qpc_data.rda")
ncpvals = sapply(ncamesOut, function(x) {x$pprime}) #matrix, rows are pvals, columns are traits
ncqvals = get_q_values(ncpvals)
pcpvals = sapply(qpcamesOut, function(x) {x$pprime}) #matrix, rows are pvals, columns are traits
qvals = get_q_values(pcpvals)

## EURO DATA
load('../data/qpc_euro_output.rda')
pcpvalseuro = sapply(qpceuroOut, function(x) {x$pprime})[1:17,] #matrix, rows are pvals, columns are traits
qvalseuro = get_q_values(pcpvalseuro)

load('../data/qpc-euro-nc.rda')
ncpvalseuro = sapply(ncEuroOut, function(x) {x$pprime})[1:17,] #matrix, rows are pvals, columns are traits
ncqvalseuro = get_q_values(ncpvalseuro)


mycol = lacroix_palette('Lime')
postscript("Qxpc_pvals.eps",height=8,width=8,paper="special",horizontal=FALSE,colormodel="cymk")
par(mfrow=c(2,2), mar=c(6,6,2,2), cex.lab=1.2)
hist(c(ncpvals), main="", xlab = "Ames -- Non-conditional test p values", border="white", col = mycol[1], breaks=30)
mtext('A', side=3, adj=4, cex=2, line=0, at=0)
#abline(v = max(c(ncpvals)[c(ncqvals) < 0.05]), col = mycol[6], lwd=2)

hist(c(pcpvals), main="", xlab = "Ames -- Conditional test p values", border="white", col = mycol[1], breaks=30)
mtext('B', side=3, adj=4, cex=2, line=0, at=0)

hist(c(ncpvalseuro), main="", xlab = "Euro -- Non-conditional test p values", border="white", col = mycol[1], breaks=30)
mtext('D', side=3, adj=4, cex=2, line=0, at=0)
#abline(v = max(c(ncpvalseuro)[c(ncqvalseuro) < 0.05]), col = mycol[6], lwd=2)


hist(c(pcpvalseuro), main="", xlab = "Euro -- Conditional test p values", border="white", col = mycol[1], breaks=30)
mtext('D', side=3, adj=4, cex=2, line=0, at=0)


dev.off()
