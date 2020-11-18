
library(here)
library(pROC)
source("Combat_HKfunctions.R")


dir.FilterMerge<-"./Results"
dir.Results<-"./ResultsFilter"
dir.Input<-"/path/to/downloaded/figshare/"


CCR1<-readRDS(file=paste0(dir.FilterMerge,"/CCR_SQ_Combat_All_merge_F.Rds"))
CCR2<-readRDS(file=paste0(dir.FilterMerge,"/CCR_SQ_Combat_All_NoNorm_merge_F.Rds"))
CCR3<-readRDS(file=paste0(dir.FilterMerge,"/CCR_SQ_Combat_PC1_All_merge_F.Rds"))
CCR4<-readRDS(file=paste0(dir.FilterMerge,"/CCR_SQ_Combat_PC2_All_merge_F.Rds"))

CERES1<-readRDS(file=paste0(dir.FilterMerge,"/CERES_SQ_Combat_All_merge_F.Rds"))
CERES2<-readRDS(file=paste0(dir.FilterMerge,"/CERES_SQ_Combat_All_NoNorm_merge_F.Rds"))
CERES3<-readRDS(file=paste0(dir.FilterMerge,"/CERES_SQ_Combat_PC1_All_merge_F.Rds"))
CERES4<-readRDS(file=paste0(dir.FilterMerge,"/CERES_SQ_Combat_PC2_All_merge_F.Rds"))

CCRJ1<-readRDS(file=paste0(dir.FilterMerge,"/CCR_SQ_Combat_All_JACKS_merge_F.Rds"))
CCRJ2<-readRDS(file=paste0(dir.FilterMerge,"/CCR_SQ_Combat_All_NoNorm_JACKS_merge_F.Rds"))
CCRJ3<-readRDS(file=paste0(dir.FilterMerge,"/CCR_SQ_Combat_PC1_All_JACKS_merge_F.Rds"))
CCRJ4<-readRDS(file=paste0(dir.FilterMerge,"/CCR_SQ_Combat_PC2_All_JACKS_merge_F.Rds"))

load(paste0(dir.Input,"/BAGEL_essential.RData"))
load(paste0(dir.Input,"/BAGEL_nonEssential.RData"))
CCR_perfTH1<-list()
CCR_perfTH2<-list()
CCR_perfTH3<-list()
CCR_perfTH4<-list()


CCR_auc4<-list()
CCR_auc3<-list()
CCR_auc2<-list()
CCR_auc1<-list()
for(i in 1:ncol(CCR1)){
  output<-bangmag.roc(BFs=as.matrix(CCR1[,i]),ess_genes=BAGEL_essential,non_ess_genes = BAGEL_nonEssential,CL=colnames(CCR1)[i],th=0.05)
  CCR_perfTH1[[i]]<-output[[2]][1]
  CCR_auc1[[i]]<-output[[1]][1]
  output<-bangmag.roc(BFs=as.matrix(CCR2[,i]),ess_genes=BAGEL_essential,non_ess_genes = BAGEL_nonEssential,CL=colnames(CCR2)[i],th=0.05)
  CCR_perfTH2[[i]]<-output[[2]][1]
  CCR_auc2[[i]]<-output[[1]][1]
  output<-bangmag.roc(BFs=as.matrix(CCR3[,i]),ess_genes=BAGEL_essential,non_ess_genes = BAGEL_nonEssential,CL=colnames(CCR3)[i],th=0.05)
  CCR_perfTH3[[i]]<-output[[2]][1]
  CCR_auc3[[i]]<-output[[1]][1]
  output<-bangmag.roc(BFs=as.matrix(CCR4[,i]),ess_genes=BAGEL_essential,non_ess_genes = BAGEL_nonEssential,CL=colnames(CCR4)[i],th=0.05)
  CCR_perfTH4[[i]]<-output[[2]][1]
  CCR_auc4[[i]]<-output[[1]][1]
  
  
}

names(CCR_perfTH1)<-colnames(CCR1)
CCR1_FCs<-
  OT15_assembleDepletionMatrix(cellLines = colnames(CCR1),BFths = unlist(CCR_perfTH1),inputFolder = NULL,BFvalues = FALSE,BFmatrix=CCR1)
names(CCR_perfTH2)<-colnames(CCR2)
CCR2_FCs<-
  OT15_assembleDepletionMatrix(cellLines = colnames(CCR2),BFths = unlist(CCR_perfTH2),inputFolder = NULL,BFvalues = FALSE,BFmatrix=CCR2)
names(CCR_perfTH3)<-colnames(CCR3)
CCR3_FCs<-
  OT15_assembleDepletionMatrix(cellLines = colnames(CCR3),BFths = unlist(CCR_perfTH3),inputFolder = NULL,BFvalues = FALSE,BFmatrix=CCR3)
names(CCR_perfTH4)<-colnames(CCR4)
CCR4_FCs<-
  OT15_assembleDepletionMatrix(cellLines = colnames(CCR4),BFths = unlist(CCR_perfTH4),inputFolder = NULL,BFvalues = FALSE,BFmatrix=CCR4)


CERES_perfTH1<-list()
CERES_perfTH2<-list()
CERES_perfTH3<-list()
CERES_perfTH4<-list()

CERES_auc4<-list()
CERES_auc3<-list()
CERES_auc2<-list()
CERES_auc1<-list()
for(i in 1:ncol(CERES1)){
  output<-bangmag.roc(BFs=as.matrix(CERES1[,i]),ess_genes=BAGEL_essential,non_ess_genes = BAGEL_nonEssential,CL=colnames(CERES1)[i],th=0.05)
  CERES_perfTH1[[i]]<-output[[2]][1]
  CERES_auc1[[i]]<-output[[1]][1]
  output<-bangmag.roc(BFs=as.matrix(CERES2[,i]),ess_genes=BAGEL_essential,non_ess_genes = BAGEL_nonEssential,CL=colnames(CERES2)[i],th=0.05)
  CERES_perfTH2[[i]]<-output[[2]][1]
  CERES_auc2[[i]]<-output[[1]][1]
  output<-bangmag.roc(BFs=as.matrix(CERES3[,i]),ess_genes=BAGEL_essential,non_ess_genes = BAGEL_nonEssential,CL=colnames(CERES3)[i],th=0.05)
  CERES_perfTH3[[i]]<-output[[2]][1]
  CERES_auc3[[i]]<-output[[1]][1]
  output<-bangmag.roc(BFs=as.matrix(CERES4[,i]),ess_genes=BAGEL_essential,non_ess_genes = BAGEL_nonEssential,CL=colnames(CERES4)[i],th=0.05)
  CERES_perfTH4[[i]]<-output[[2]][1]
  CERES_auc4[[i]]<-output[[1]][1]
  
}

names(CERES_perfTH1)<-colnames(CERES1)
CERES1_FCs<-
  OT15_assembleDepletionMatrix(cellLines = colnames(CERES1),BFths = unlist(CERES_perfTH1),inputFolder = NULL,BFvalues = FALSE,BFmatrix=CERES1)
names(CERES_perfTH2)<-colnames(CERES2)
CERES2_FCs<-
  OT15_assembleDepletionMatrix(cellLines = colnames(CERES2),BFths = unlist(CERES_perfTH2),inputFolder = NULL,BFvalues = FALSE,BFmatrix=CERES2)
names(CERES_perfTH3)<-colnames(CERES3)
CERES3_FCs<-
  OT15_assembleDepletionMatrix(cellLines = colnames(CERES3),BFths = unlist(CERES_perfTH3),inputFolder = NULL,BFvalues = FALSE,BFmatrix=CERES3)
names(CERES_perfTH4)<-colnames(CERES4)
CERES4_FCs<-
  OT15_assembleDepletionMatrix(cellLines = colnames(CERES4),BFths = unlist(CERES_perfTH4),inputFolder = NULL,BFvalues = FALSE,BFmatrix=CERES4)

CCRJ_perfTH1<-list()
CCRJ_perfTH2<-list()
CCRJ_perfTH3<-list()
CCRJ_perfTH4<-list()


CCRJ_auc4<-list()
CCRJ_auc3<-list()
CCRJ_auc2<-list()
CCRJ_auc1<-list()
for(i in 1:ncol(CCRJ1)){
  output<-bangmag.roc(BFs=as.matrix(CCRJ1[,i]),ess_genes=BAGEL_essential,non_ess_genes = BAGEL_nonEssential,CL=colnames(CCRJ1)[i],th=0.05)
  CCRJ_perfTH1[[i]]<-output[[2]][1]
  CCRJ_auc1[[i]]<-output[[1]][1]
  output<-bangmag.roc(BFs=as.matrix(CCRJ2[,i]),ess_genes=BAGEL_essential,non_ess_genes = BAGEL_nonEssential,CL=colnames(CCRJ2)[i],th=0.05)
  CCRJ_perfTH2[[i]]<-output[[2]][1]
  CCRJ_auc2[[i]]<-output[[1]][1]
  output<-bangmag.roc(BFs=as.matrix(CCRJ3[,i]),ess_genes=BAGEL_essential,non_ess_genes = BAGEL_nonEssential,CL=colnames(CCRJ3)[i],th=0.05)
  CCRJ_perfTH3[[i]]<-output[[2]][1]
  CCRJ_auc3[[i]]<-output[[1]][1]
  output<-bangmag.roc(BFs=as.matrix(CCRJ4[,i]),ess_genes=BAGEL_essential,non_ess_genes = BAGEL_nonEssential,CL=colnames(CCRJ4)[i],th=0.05)
  CCRJ_perfTH4[[i]]<-output[[2]][1]
  CCRJ_auc4[[i]]<-output[[1]][1]
  
  
}

names(CCRJ_perfTH1)<-colnames(CCRJ1)
CCRJ1_FCs<-
  OT15_assembleDepletionMatrix(cellLines = colnames(CCRJ1),BFths = unlist(CCRJ_perfTH1),inputFolder = NULL,BFvalues = FALSE,BFmatrix=CCRJ1)
names(CCRJ_perfTH2)<-colnames(CCRJ2)
CCRJ2_FCs<-
  OT15_assembleDepletionMatrix(cellLines = colnames(CCRJ2),BFths = unlist(CCRJ_perfTH2),inputFolder = NULL,BFvalues = FALSE,BFmatrix=CCRJ2)
names(CCRJ_perfTH3)<-colnames(CCRJ3)
CCRJ3_FCs<-
  OT15_assembleDepletionMatrix(cellLines = colnames(CCRJ3),BFths = unlist(CCRJ_perfTH3),inputFolder = NULL,BFvalues = FALSE,BFmatrix=CCRJ3)
names(CCRJ_perfTH4)<-colnames(CCRJ4)
CCRJ4_FCs<-
  OT15_assembleDepletionMatrix(cellLines = colnames(CCRJ4),BFths = unlist(CCRJ_perfTH4),inputFolder = NULL,BFvalues = FALSE,BFmatrix=CCRJ4)


#can do assessment number of depleted genes identified versus number of unique depleted genes etc.
BinaryCCR1<-(CCR1_FCs<0)+0
BinaryCERES1<-(CERES1_FCs<0)+0
BinaryCCRJ1<-(CCRJ1_FCs<0)+0
BinaryCCR2<-(CCR2_FCs<0)+0
BinaryCERES2<-(CERES2_FCs<0)+0
BinaryCCRJ2<-(CCRJ2_FCs<0)+0
BinaryCCR3<-(CCR3_FCs<0)+0
BinaryCERES3<-(CERES3_FCs<0)+0
BinaryCCRJ3<-(CCRJ3_FCs<0)+0
BinaryCCR4<-(CCR4_FCs<0)+0
BinaryCERES4<-(CERES4_FCs<0)+0
BinaryCCRJ4<-(CCRJ4_FCs<0)+0

save(BinaryCCR1,file=paste0(dir.Results,"/BinaryCCR1m.RData"))
save(BinaryCCR2,file=paste0(dir.Results,"/BinaryCCR2m.RData"))
save(BinaryCCR3,file=paste0(dir.Results,"/BinaryCCR3m.RData"))
save(BinaryCCR4,file=paste0(dir.Results,"/BinaryCCR4m.RData"))

save(BinaryCERES1,file=paste0(dir.Results,"/BinaryCERES1m.RData"))
save(BinaryCERES2,file=paste0(dir.Results,"/BinaryCERES2m.RData"))
save(BinaryCERES3,file=paste0(dir.Results,"/BinaryCERES3m.RData"))
save(BinaryCERES4,file=paste0(dir.Results,"/BinaryCERES4m.RData"))

save(BinaryCCRJ1,file=paste0(dir.Results,"/BinaryCCRJ1m.RData"))
save(BinaryCCRJ2,file=paste0(dir.Results,"/BinaryCCRJ2m.RData"))
save(BinaryCCRJ3,file=paste0(dir.Results,"/BinaryCCRJ3m.RData"))
save(BinaryCCRJ4,file=paste0(dir.Results,"/BinaryCCRJ4m.RData"))

