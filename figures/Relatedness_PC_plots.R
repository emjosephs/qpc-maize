
#read in kinship matrix for 240 and eigen decompose
myK = read.table('../data/All_240E.nomaf.nomissing.K')
myKnames = read.table('../data/240.names', stringsAsFactors = F)$V1
row.names(myK) = myKnames[1:dim(myK)[1]]
eigF = eigen(myK)

#read in 240 pop info and combine with eigenvectors of K
fgtable = read.table('../data/FlintGarciaTableS1.csv', sep=',', header=T, stringsAsFactors = F)
mydf = data.frame(line = myKnames[1:dim(myK)[1]], eigF$vectors)
mydf$Inbred = sapply(as.character(mydf$line), function(x){strsplit(x,'_')[[1]][2]})
fgmerge = dplyr::inner_join(fgtable, mydf, by='Inbred')

#read in eigen decomp of kinship matrix for ames + 281
load('../data/ames.281E.regeig.rda')

#read in eigen decomp of euro matrix
load('../data/euro.282.eig.rda')

#read in euro pop info
load('../data/euro_qpc_data.rda')

#plot stuff
# mycol = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", '#D55E00', '#CC79A7')
# palette(mycol[c(1,4,5,3,7,8)])
# col282 = mycol[2]
# colAmes = mycol[6]
# colFlint = mycol[8] 
# colDent = mycol[4]
library(LaCroixColoR)
#palette(c('darkgray',lacroix_palette('Lime')[c(3,1,4,5,6)]))
palette(c('darkgray',lacroix_palette('Pamplemousse')[c(3,1,4,5,6)]))

 col282 = lacroix_palette('Lime')[2]
 colAmes = lacroix_palette('Lime')[6]
 colFlint = lacroix_palette('Lime')[4]
 colDent = lacroix_palette('Lime')[5]


nicepops = c('mixed','non-stiff stalk','popcorn','stiff-stalk','sweet','tropical') #readable names
postscript("Relatedness_PC_plots.eps",height=5,width=12,paper="special",horizontal=FALSE,colormodel="cymk")

par(mfrow=c(1,3), par(mar=c(5,7,3,2)))

##Panel A, PCs within the GWAS Panel
plot(fgmerge$X1, fgmerge$X2, col = as.factor(fgmerge$Subpopulation),
     lwd=1.5, bty="n", xlab = "", ylab = "", cex.lab=1.5, cex=1.5, yaxt="n", 
     xaxt = "n", xlim = c(-0.35,0.2), ylim = c(-0.17,0.17))
axis(1, lwd=2, cex.axis=1.5 )
axis(2, lwd=2, cex.axis=1.5, las=2)
legend("bottomleft", nicepops, bty="n", pch=1, pt.lwd=2, col = palette(), cex = 1.5)
mtext('A', side=3, adj=0, cex=2, line=-2)
mtext('GWAS panel PC 1', side=1, line = 3)
mtext('GWAS panel PC 2', side=2, line = 5)


###Panel B, PCs in 281 + Ames panel
plot(regEig$vectors[,1], regEig$vectors[,2], bty="n", xlab = "", ylab="", 
     col=c(rep(colAmes, 2704), rep(col282, 280)), lwd=1.5, cex = 1.5,cex.lab=1.5,
     yaxt="n", xaxt="n", ylim = c(-0.045, 0.09))
axis(1, lwd=2, cex.axis=1.5)
axis(2, lwd=2, cex.axis=1.5, las=2)
legend('bottomleft',c('Ames panel','GWAS panel'),pch=1, pt.lwd=2, col = c(colAmes, col282), cex=1.5, bty="n")
mtext('B', side=3, adj=0, cex=2, line=-2)
mtext('Combined PC 1', side=1, line = 3)
mtext('Combined PC 2', side=2, line = 5)


###Panel C, PCs in Euro + 282 panel
#palette(c(colDent,colFlint))
#plot(-euro282eigen$vectors[1:906,1], -euro282eigen$vectors[1:906,2], bty="n", xlab = "", ylab="", 
#     col=as.factor(mymerge$Type), lwd=1.5, cex.lab=1.5,
#     xaxt="n", yaxt="n", cex=1.5, ylim = c(-0.07, 0.077), xlim = c(-0.055, 0.05))
plot(euro282eigen$vectors[1:906,1], -euro282eigen$vectors[1:906,2], bty="n", xlab = "", ylab="", 
     col=colDent, lwd=1.5, cex.lab=1.5,
     xaxt="n", yaxt="n", cex=1.5, ylim = c(-0.077, 0.077), xlim = c(-0.055, 0.05))
points(-euro282eigen$vectors[906:1168,1], -euro282eigen$vectors[906:1168,2], 
       col = col282, lwd=1, cex=1.5)
axis(1, lwd=2, cex.axis=1.5)
axis(2, lwd=2, cex.axis=1.5, las=2, at = seq(-0.06, 0.06, by = 0.02))
#legend('bottomright',c('Dent','Flint','GWAS panel'),pch=1, pt.lwd=2, col = c(colDent,colFlint,col282), cex=1.5, bty="n")
legend('bottomleft',c('European Landraces','GWAS panel'),pch=1, pt.lwd=2, col = c(colDent,col282), cex=1.5, bty="n")

mtext('C', side=3, adj=0, cex=2, line=-2)
mtext('Combined PC 1', side=1, line = 3)
mtext('Combined PC 2', side=2, line = 5)

dev.off()


postscript("Supp_fig_pc10.eps",height=7,width=10,paper="special",horizontal=FALSE,colormodel="cymk")
par(mfrow=c(1,1), par(mar=c(5,7,3,2)))
plot(fgmerge$X10, fgmerge$X2, col = as.factor(fgmerge$Subpopulation),
     lwd=2, bty="n", xlab = "", ylab = "", cex.lab=1.5, cex=1.5, yaxt="n", 
     xaxt = "n", xlim = c(-0.25, 0.25), ylim = c(-0.2, 0.15))
axis(1, lwd=2, cex.axis=1.5 )
axis(2, lwd=2, cex.axis=1.5, las=2)
#legend("bottomright", levels(as.factor(fgmerge$Subpopulation)), bty="n", pch=1, pt.lwd=2, col = mycol)
legend("bottomleft", nicepops, bty="n", pch=1, pt.lwd=2, col = palette(), cex = 1.5)
mtext('GWAS panel PC 10', side=1, line = 3)
mtext('GWAS panel PC 2', side=2, line = 5)
dev.off()

