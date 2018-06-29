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

png(paste('gwas-plots/euro.',traitName,'.manhattan.png',sep=""))
manhattan(mandata)

dev.off()

png(paste('gwas-plots/euro.',traitName,'.qq.png',sep=""))
qq(mandata$P)
dev.off()
}


#plots for ames gwas
make_ames_plots <- function(myI){
  traitName = traitNames[as.numeric(myI)-1]
  mydata = processGemmaOutput(paste('/home/emjo/282/gemma_output/',traitName,'.281.assoc.txt', sep=""))
  mandata = data.frame(CHR = mydata$chr, BP = mydata$ps, P = mydata$p_lrt)
  
  png(paste('gwas-plots/ames.',traitName,'.manhattan.png',sep=""))
  manhattan(mandata)
  
  dev.off()
  
  png(paste('gwas-plots/ames.',traitName,'.qq.png',sep=""))
  qq(mandata$P)
  dev.off()
  
}


make_euro_plots(args[1])
make_ames_plots(args[1])

