library(qpctools)
library(qqman)
args = commandArgs(trailingOnly=TRUE)


#get trait names
traitNames = read.table('../data/blup.names', stringsAsFactors=F)$V1
niceTraitnames = sapply(traitNames, function(x){
  strsplit(x, '0607')[[1]][1]
})

#plots for 600K only european lines
make_euro_plots <- function(myI){
  traitName = niceTraitnames[as.numeric(myI)-1]
  mydata = read.table(paste('/home/emjo/euro-maize/results/chip263.',myI,'.assoc.txt', sep=""), header=T, stringsAsFactors=F)
  mandata = data.frame(CHR = mydata$chr, BP = mydata$ps, P = mydata$p_lrt)
  
  
  #png(paste('gwas-plots/euro.',traitName,'.manhattan.png',sep=""))
  manhattan(mandata, main = paste(traitName, '263'), suggestiveline=F, genomewideline = -log10(0.005), cex.axis=1.5, cex.lab = 1.5)
  mtext('A', side=3, adj=0, cex=2, line=2)
  
  #dev.off()
  
  #png(paste('gwas-plots/euro.',traitName,'.qq.png',sep=""))
  qq(mandata$P, main = paste(traitName, '263'), cex.lab = 1.5, cex.axis = 1.5, bty="n", las=1)
  mtext('B', side=3, adj=0, cex=2, line=2)
  
  #dev.off()
}


#plots for ames gwas
make_ames_plots <- function(myI){
  traitName = traitNames[as.numeric(myI)-1]
  mydata = processGemmaOutput(paste('/home/emjo/282/gemma_output/',traitName,'.281.assoc.txt', sep=""))
  mandata = data.frame(CHR = mydata$chr, BP = mydata$ps, P = mydata$p_lrt)
  
  # png(paste('gwas-plots/ames.',traitName,'.manhattan.png',sep=""))
  manhattan(mandata, suggestiveline=F, main = paste(traitName, '281'), genomewideline = -log10(0.005), cex.axis=1.5, cex.lab = 1.5)
  mtext('C', side=3, adj=0, cex=2, line=2)
  
  # dev.off()
  
  # png(paste('gwas-plots/ames.',traitName,'.qq.png',sep=""))
  qq(mandata$P,  cex.lab = 1.5, cex.axis = 1.5, bty="n", las=1,main = paste(traitName, '281'))
  mtext('D', side=3, adj=0, cex=2, line=2)
  
  # dev.off()
  
}

postscript(paste('gwas-plots/',args[1],'.gwas.eps', sep=""), height=10, width=10, paper="special", horizontal=FALSE, colormodel="cymk")
par(mfrow=c(2,2), mar=c(6,6,4,2))
make_euro_plots(args[1])
make_ames_plots(args[1])
dev.off()
               

