#rm(list=ls())

#wd <- "~/scratch5/mslimp/LakeElsinore/HumanAln"
#setwd(wd)

Files <- list.files(pattern=".nozeroes.cov.txt",full.names=TRUE)

Blast <- data.frame()
for(File in Files){
  BlastTMP <- read.table(File, header=F, comment.char="?", sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
  BlastTMP$filename <- gsub("./","",File)
  Blast <- rbind(Blast,BlastTMP)
}

write.table(Blast,"ConcatCoverage_SG_AlnENT_NM3.txt",quote=FALSE,sep="\t",row.names = FALSE)
