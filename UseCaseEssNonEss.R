###Use case 1 - recall of essential genes###
library(here)
source("EssNonEssfunctions.R")

source("Combat_HKfunctions.R")

library(ggplot2)
library(gridExtra)
library(magrittr)
library(tidyverse)
library(dplyr)

dir.Results<-"./ResultsFilter"
dir.Input<-"/path/to/downloaded/figshare/"
curveColours<-c("#3375A2","#E1822C","#3B9144","#C13E3F")
names(curveColours)<-c("ComBat",'ComBat+QN','ComBat+PC1','ComBat+PC2')
PipelineColours<-c("#F89952","#7ECBB2","#90C73E")
names(PipelineColours)<-c("CCR","CCRJ","CERES")

#see additional python scripts in github for generation of the NNMD statistic values:
nnmdvals<-read.csv(paste0(dir.Input,"/nnmd.csv"),header=T,stringsAsFactors = F)
cColours<-curveColours
names(cColours)<-unique(nnmdvals$Algorithm)[c(2,1,3,4)]
PipelineColours<-c(PipelineColours,"CCR-JACKS"="#4DFFC9FF" )
nnmdvals$Dataset<-factor(nnmdvals$Dataset,levels=c("CRISPRCleanR","CCR-JACKS","CERES"))

NNMDplot<-ggplot(nnmdvals,aes(x=Dataset,y=nnmd,fill=Algorithm))+geom_boxplot()+theme(axis.text.x = element_text(angle=45))+ylab("NNMD")+xlab("")+theme_bw()+scale_fill_manual(values=cColours)
pdf(paste0(dir.Results,"/NNMD_merge.pdf"))
print(NNMDplot)
dev.off()

wilcox.test(nnmdvals[nnmdvals$Dataset=="CRISPRCleanR"&nnmdvals$Algorithm=="ComBat","nnmd"],nnmdvals[nnmdvals$Dataset=="CCR-JACKS"&nnmdvals$Algorithm=="ComBat","nnmd"])
wilcox.test(nnmdvals[nnmdvals$Dataset=="CRISPRCleanR"&nnmdvals$Algorithm=="ComBat+QN","nnmd"],nnmdvals[nnmdvals$Dataset=="CCR-JACKS"&nnmdvals$Algorithm=="ComBat+QN","nnmd"])
wilcox.test(nnmdvals[nnmdvals$Dataset=="CRISPRCleanR"&nnmdvals$Algorithm=="ComBat+QN+PC1","nnmd"],nnmdvals[nnmdvals$Dataset=="CCR-JACKS"&nnmdvals$Algorithm=="ComBat+QN+PC1","nnmd"])
wilcox.test(nnmdvals[nnmdvals$Dataset=="CRISPRCleanR"&nnmdvals$Algorithm=="ComBat+QN+PC1-2","nnmd"],nnmdvals[nnmdvals$Dataset=="CCR-JACKS"&nnmdvals$Algorithm=="ComBat+QN+PC1-2","nnmd"])

wilcox.test(nnmdvals[nnmdvals$Dataset=="CRISPRCleanR"&nnmdvals$Algorithm=="ComBat","nnmd"],nnmdvals[nnmdvals$Dataset=="CERES"&nnmdvals$Algorithm=="ComBat","nnmd"])
wilcox.test(nnmdvals[nnmdvals$Dataset=="CRISPRCleanR"&nnmdvals$Algorithm=="ComBat+QN","nnmd"],nnmdvals[nnmdvals$Dataset=="CERES"&nnmdvals$Algorithm=="ComBat+QN","nnmd"])
wilcox.test(nnmdvals[nnmdvals$Dataset=="CRISPRCleanR"&nnmdvals$Algorithm=="ComBat+QN+PC1","nnmd"],nnmdvals[nnmdvals$Dataset=="CERES"&nnmdvals$Algorithm=="ComBat+QN+PC1","nnmd"])
wilcox.test(nnmdvals[nnmdvals$Dataset=="CRISPRCleanR"&nnmdvals$Algorithm=="ComBat+QN+PC1-2","nnmd"],nnmdvals[nnmdvals$Dataset=="CERES"&nnmdvals$Algorithm=="ComBat+QN+PC1-2","nnmd"])

wilcox.test(nnmdvals[nnmdvals$Dataset=="CERES"&nnmdvals$Algorithm=="ComBat","nnmd"],nnmdvals[nnmdvals$Dataset=="CCR-JACKS"&nnmdvals$Algorithm=="ComBat","nnmd"])
wilcox.test(nnmdvals[nnmdvals$Dataset=="CERES"&nnmdvals$Algorithm=="ComBat+QN","nnmd"],nnmdvals[nnmdvals$Dataset=="CCR-JACKS"&nnmdvals$Algorithm=="ComBat+QN","nnmd"])
wilcox.test(nnmdvals[nnmdvals$Dataset=="CERES"&nnmdvals$Algorithm=="ComBat+QN+PC1","nnmd"],nnmdvals[nnmdvals$Dataset=="CCR-JACKS"&nnmdvals$Algorithm=="ComBat+QN+PC1","nnmd"])
wilcox.test(nnmdvals[nnmdvals$Dataset=="CERES"&nnmdvals$Algorithm=="ComBat+QN+PC1-2","nnmd"],nnmdvals[nnmdvals$Dataset=="CCR-JACKS"&nnmdvals$Algorithm=="ComBat+QN+PC1-2","nnmd"])


mediannnmd<-nnmdvals%>%group_by(Algorithm,Dataset)%>%summarize(m=median(nnmd))

write.table(mediannnmd,file=paste0(dir.Results,"/mediannnmd.tsv"),sep="\t")





