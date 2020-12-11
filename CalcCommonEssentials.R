
library(here)
source("EssNonEssfunctions.R")

#load a reference set of essential genes
library(CRISPRcleanR)
data(BAGEL_essential)
data(BAGEL_nonEssential)

source("Combat_HKfunctions.R")
library(ADaM2)
library(ggplot2)
library(gridExtra)
library(magrittr)
library(tidyverse)
library(dplyr)
library(magrittr)

dir.MergeFile<-"./Results"
dir.Results<-"./ResultsFilter"
dir.Input<-"/path/to/downloaded/figshare/"

cmp<-read.csv(paste0(dir.Input,"/model_list_latest.csv"),header=T,stringsAsFactors = F)
cmp2<-cmp
cmp2$model_id<-cmp2$BROAD_ID
cmp<-rbind(cmp,cmp2)
###use Broad 90th depletion method ###
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

B90CCR<-Broad90(CCR_corrected,BAGEL_essential,BAGEL_nonEssential)
B90CCRJ<-Broad90(CCRJ_corrected,BAGEL_essential,BAGEL_nonEssential)
B90CERES<-Broad90(CERES_corrected,BAGEL_essential,BAGEL_nonEssential)

B90CCRQN<-Broad90(CCR_correctedQN,BAGEL_essential,BAGEL_nonEssential)
B90CCRJQN<-Broad90(CCRJ_correctedQN,BAGEL_essential,BAGEL_nonEssential)
B90CERESQN<-Broad90(CERES_correctedQN,BAGEL_essential,BAGEL_nonEssential)

B90CCRPC1<-Broad90(CCR_correctedPC1,BAGEL_essential,BAGEL_nonEssential)
B90CCRJPC1<-Broad90(CCRJ_correctedPC1,BAGEL_essential,BAGEL_nonEssential)
B90CERESPC1<-Broad90(CERES_correctedPC1,BAGEL_essential,BAGEL_nonEssential)

B90CCRPC2<-Broad90(CCR_correctedPC2,BAGEL_essential,BAGEL_nonEssential)
B90CCRJPC2<-Broad90(CCRJ_correctedPC2,BAGEL_essential,BAGEL_nonEssential)
B90CERESPC2<-Broad90(CERES_correctedPC2,BAGEL_essential,BAGEL_nonEssential)

save(B90CCR,file=paste0(dir.Results,"/CE_IntCCR90.Rdata"))
save(B90CCRJ,file=paste0(dir.Results,"/CE_IntCCRJ90.Rdata"))
save(B90CERES,file=paste0(dir.Results,"/CE_IntCERES90.Rdata"))

save(B90CCRQN,file=paste0(dir.Results,"/CE_IntCCRQN90.Rdata"))
save(B90CCRJQN,file=paste0(dir.Results,"/CE_IntCCRJQN90.Rdata"))
save(B90CERESQN,file=paste0(dir.Results,"/CE_IntCERESQN90.Rdata"))

save(B90CCRPC1,file=paste0(dir.Results,"/CE_IntCCRPC190.Rdata"))
save(B90CCRJPC1,file=paste0(dir.Results,"/CE_IntCCRJPC190.Rdata"))
save(B90CERESPC1,file=paste0(dir.Results,"/CE_IntCERESPC190.Rdata"))

save(B90CCRPC2,file=paste0(dir.Results,"/CE_IntCCRPC290.Rdata"))
save(B90CCRJPC2,file=paste0(dir.Results,"/CE_IntCCRJPC290.Rdata"))
save(B90CERESPC2,file=paste0(dir.Results,"/CE_IntCERESPC290.Rdata"))

load(paste0(dir.Results,"/BinaryCCRJ4m.RData"))
load(paste0(dir.Results,"/BinaryCCR4m.RData"))
load(paste0(dir.Results,"/BinaryCERES4m.RData"))

PCcoreCERESPC2<-AdamTissue(BinaryCERES4,cmp)
PCcoreCCRPC2<-AdamTissue(BinaryCCR4,cmp)
PCcoreCCRJPC2<-AdamTissue(BinaryCCRJ4,cmp)

save(PCcoreCCRJPC2,file=paste0(dir.Results,"/AdamCoreCCRJPC2.RData"))
save(PCcoreCCRPC2,file=paste0(dir.Results,"/AdamCoreCCRPC2.RData"))
save(PCcoreCERESPC2,file=paste0(dir.Results,"/AdamCoreCERESPC2.RData"))


load(paste0(dir.Results,"/BinaryCCRJ3m.RData"))
load(paste0(dir.Results,"/BinaryCCR3m.RData"))
load(paste0(dir.Results,"/BinaryCERES3m.RData"))

PCcoreCERESPC1<-AdamTissue(BinaryCERES3,cmp)
PCcoreCCRPC1<-AdamTissue(BinaryCCR3,cmp)
PCcoreCCRJPC1<-AdamTissue(BinaryCCRJ3,cmp)

save(PCcoreCCRJPC1,file=paste0(dir.Results,"/AdamCoreCCRJPC1.RData"))
save(PCcoreCCRPC1,file=paste0(dir.Results,"/AdamCoreCCRPC1.RData"))
save(PCcoreCERESPC1,file=paste0(dir.Results,"/AdamCoreCERESPC1.RData"))

load(paste0(dir.Results,"/BinaryCCRJ2m.RData"))
load(paste0(dir.Results,"/BinaryCCR2m.RData"))
load(paste0(dir.Results,"/BinaryCERES2m.RData"))

PCcoreCERESQN<-AdamTissue(BinaryCERES2,cmp)
PCcoreCCRQN<-AdamTissue(BinaryCCR2,cmp)
PCcoreCCRJQN<-AdamTissue(BinaryCCRJ2,cmp)

save(PCcoreCCRJQN,file=paste0(dir.Results,"/AdamCoreCCRJQN.RData"))
save(PCcoreCCRQN,file=paste0(dir.Results,"/AdamCoreCCRQN.RData"))
save(PCcoreCERESQN,file=paste0(dir.Results,"/AdamCoreCERESQN.RData"))

load(paste0(dir.Results,"/BinaryCCRJ1m.RData"))
load(paste0(dir.Results,"/BinaryCCR1m.RData"))
load(paste0(dir.Results,"/BinaryCERES1m.RData"))

PCcoreCERES<-AdamTissue(BinaryCERES1,cmp)
PCcoreCCR<-AdamTissue(BinaryCCR1,cmp)
PCcoreCCRJ<-AdamTissue(BinaryCCRJ1,cmp)

save(PCcoreCCRJ,file=paste0(dir.Results,"/AdamCoreCCRJ.RData"))
save(PCcoreCCR,file=paste0(dir.Results,"/AdamCoreCCR.RData"))
save(PCcoreCERES,file=paste0(dir.Results,"/AdamCoreCERES.RData"))

#load in original Broad and Sanger datasets:
load(paste0(dir.MergeFile,"/BroadDataCERES.Rdata"))
load(paste0(dir.MergeFile,"/SangerDataCERES.Rdata"))
load(paste0(dir.MergeFile,"/BroadData.Rdata"))
load(paste0(dir.MergeFile,"/SangerData.Rdata"))
load(paste0(dir.MergeFile,"/JACKSsanger.Rdata"))
load(paste0(dir.MergeFile,"/JACKSbroad.Rdata"))

dn<-dimnames(BroadDataCERES)
BroadDataCERES<-normalize.quantiles(BroadDataCERES)
dimnames(BroadDataCERES)<-dn

dn<-dimnames(SangerDataCERES)
SangerDataCERES<-normalize.quantiles(SangerDataCERES)
dimnames(SangerDataCERES)<-dn


B90CCR_sanger<-Broad90(SangerData,BAGEL_essential,BAGEL_nonEssential)
B90CCR_broad<-Broad90(BroadData,BAGEL_essential,BAGEL_nonEssential)
B90CERES_sanger<-Broad90(SangerDataCERES,BAGEL_essential,BAGEL_nonEssential)
B90CERES_broad<-Broad90(BroadDataCERES,BAGEL_essential,BAGEL_nonEssential)

save(B90CCR_sanger,file=paste0(dir.Results,"/B90CCR_sanger.RData"))
save(B90CCR_broad,file=paste0(dir.Results,"/B90CCR_broad.RData"))
save(B90CERES_sanger,file=paste0(dir.Results,"/B90CERES_sanger.RData"))
save(B90CERES_broad,file=paste0(dir.Results,"/B90CERES_broad.RData"))

load(paste0(dir.Results,"/BinaryCCRSanger.RData"))
load(paste0(dir.Results,"/BinaryCCRBroad.RData"))
load(paste0(dir.Results,"/BinaryCERESSanger.RData"))
load(paste0(dir.Results,"/BinaryCERESBroad.RData"))


PCcoreCERES_Sanger<-AdamTissue(BinaryCERESSanger,cmp)
PCcoreCERES_Broad<-AdamTissue(BinaryCERESBroad,cmp)
PCcoreCCR_Sanger<-AdamTissue(BinaryCCRSanger,cmp)
PCcoreCCR_Broad<-AdamTissue(BinaryCCRBroad,cmp)

save(PCcoreCERES_Sanger,file=paste0(dir.Results,"/AdamCoreCERES_Sanger.RData"))
save(PCcoreCERES_Broad,file=paste0(dir.Results,"/AdamCoreCERES_Broad.RData"))
save(PCcoreCCR_Sanger,file=paste0(dir.Results,"/AdamCoreCCR_Sanger.RData"))
save(PCcoreCCR_Broad,file=paste0(dir.Results,"/AdamCoreCCR_Broad.RData"))
