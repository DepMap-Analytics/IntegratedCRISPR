
library(here)
library(MASS)
library(sn)
source("Combat_HKfunctions.R")


dir.MergeFile<-"./Results"
dir.Results<-"./ResultsFilter"
dir.Input<-"/path/to/downloaded/figshare/"
CCR_correctedPC2<-readRDS(file=paste0(dir.MergeFile,"/CCR_SQ_Combat_PC2_All_merge_F.Rds"))
CERES_correctedPC2<-readRDS(file=paste0(dir.MergeFile,"/CERES_SQ_Combat_PC2_All_merge_F.Rds"))
CCRJ_correctedPC2<-readRDS(file=paste0(dir.MergeFile,"/CCR_SQ_Combat_PC2_All_JACKS_merge_F.Rds"))

CCR_correctedPC1<-readRDS(file=paste0(dir.MergeFile,"/CCR_SQ_Combat_PC1_All_merge_F.Rds"))
CERES_correctedPC1<-readRDS(file=paste0(dir.MergeFile,"/CERES_SQ_Combat_PC1_All_merge_F.Rds"))
CCRJ_correctedPC1<-readRDS(file=paste0(dir.MergeFile,"/CCR_SQ_Combat_PC1_All_JACKS_merge_F.Rds"))

CCR_correctedQN<-readRDS(file=paste0(dir.MergeFile,"/CCR_SQ_Combat_All_merge_F.Rds"))
CERES_correctedQN<-readRDS(file=paste0(dir.MergeFile,"/CERES_SQ_Combat_All_merge_F.Rds"))
CCRJ_correctedQN<-readRDS(file=paste0(dir.MergeFile,"/CCR_SQ_Combat_All_JACKS_merge_F.Rds"))

CCR_corrected<-readRDS(file=paste0(dir.MergeFile,"/CCR_SQ_Combat_All_NoNorm_merge_F.Rds"))
CERES_corrected<-readRDS(file=paste0(dir.MergeFile,"/CERES_SQ_Combat_All_NoNorm_merge_F.Rds"))
CCRJ_corrected<-readRDS(file=paste0(dir.MergeFile,"/CCR_SQ_Combat_All_NoNorm_JACKS_merge_F.Rds"))

normLRTCCR2<-normLRT2(CCR_correctedPC2)
normLRTCERES2<-normLRT2(CERES_correctedPC2)
normLRTCCRJ2<-normLRT2(CCRJ_correctedPC2)
PCnumber<-2
save(normLRTCCR2,file=paste0(dir.Results,"/normLRTCCRmerge",PCnumber,".RData"))
save(normLRTCERES2,file=paste0(dir.Results,"/normLRTCERESmerge",PCnumber,".RData"))
save(normLRTCCRJ2,file=paste0(dir.Results,"/normLRTCCRJmerge",PCnumber,".RData"))
normLRTCCR1<-normLRT2(CCR_correctedPC1)
normLRTCERES1<-normLRT2(CERES_correctedPC1)
normLRTCCRJ1<-normLRT2(CCRJ_correctedPC1)
PCnumber<-1
save(normLRTCCR1,file=paste0(dir.Results,"/normLRTCCRmerge",PCnumber,".RData"))
save(normLRTCERES1,file=paste0(dir.Results,"/normLRTCERESmerge",PCnumber,".RData"))
save(normLRTCCRJ1,file=paste0(dir.Results,"/normLRTCCRJmerge",PCnumber,".RData"))
normLRTCCRQN<-normLRT2(CCR_correctedQN)
normLRTCERESQN<-normLRT2(CERES_correctedQN)
normLRTCCRJQN<-normLRT2(CCRJ_correctedQN)
PCnumber<-"QN"
save(normLRTCCRQN,file=paste0(dir.Results,"/normLRTCCRmerge",PCnumber,".RData"))
save(normLRTCERESQN,file=paste0(dir.Results,"/normLRTCERESmerge",PCnumber,".RData"))
save(normLRTCCRJQN,file=paste0(dir.Results,"/normLRTCCRJmerge",PCnumber,".RData"))
normLRTCCR<-normLRT2(CCR_corrected)
normLRTCERES<-normLRT2(CERES_corrected)
normLRTCCRJ<-normLRT2(CCRJ_corrected)
PCnumber<-"CT"
save(normLRTCCR,file=paste0(dir.Results,"/normLRTCCRmerge",PCnumber,".RData"))
save(normLRTCERES,file=paste0(dir.Results,"/normLRTCERESmerge",PCnumber,".RData"))
save(normLRTCCRJ,file=paste0(dir.Results,"/normLRTCCRJmerge",PCnumber,".RData"))