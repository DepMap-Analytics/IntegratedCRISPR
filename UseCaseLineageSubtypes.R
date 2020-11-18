###Use case 2 - identify lineage subtypes###


#load a reference set of essential genes

library(here)
source("Lineagefunctions.R")
source("Combat_HKfunctions.R")

library(ggplot2)

library(RColorBrewer)

dir.MergeFile<-"./Results"
dir.Results<-"./ResultsFilter"
dir.Input<-"/path/to/downloaded/figshare/"

cmp<-read.csv(paste0(dir.Input,"/model_list_latest.csv"),header=T,stringsAsFactors = F)
cmp2<-cmp
cmp2$model_id<-cmp2$BROAD_ID
cmp<-rbind(cmp,cmp2)


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


curveColours<-c("#3375A2","#E1822C","#3B9144","#C13E3F")
names(curveColours)<-c("ComBat","ComBatQN","ComBatPC1","ComBatPC2")

#Get all Breast cell lines:

BreastCL<-cmp[cmp$tissue=="Breast",]
BreastCL<-BreastCL[BreastCL$model_id%in%colnames(CCR_corrected),]


BreastSubtypes<-read.delim(file=paste0(dir.Input,"/BreastSubtypes.txt"),header=T,stringsAsFactors = F,sep="\t")
BreastCL$pam50<-BreastSubtypes[match(BreastCL$COSMIC_ID,BreastSubtypes$COSMIC.ID),"pam50"]
BreastCLu<-BreastCL[BreastCL$pam50%in%c("Her2","Basal","LumA","Normal","LumB"),]
Blabels<-BreastCLu$pam50
names(Blabels)<-BreastCLu$model_id


#show subtypes according to crispr profiles
PCnumber<-2
load(file=paste0(dir.Results,"/normLRTCCRmerge",PCnumber,".RData"))
load(file=paste0(dir.Results,"/normLRTCERESmerge",PCnumber,".RData"))
load(file=paste0(dir.Results,"/normLRTCCRJmerge",PCnumber,".RData"))

PCnumber<-1
load(file=paste0(dir.Results,"/normLRTCCRmerge",PCnumber,".RData"))
load(file=paste0(dir.Results,"/normLRTCERESmerge",PCnumber,".RData"))
load(file=paste0(dir.Results,"/normLRTCCRJmerge",PCnumber,".RData"))

PCnumber<-"QN"
load(file=paste0(dir.Results,"/normLRTCCRmerge",PCnumber,".RData"))
load(file=paste0(dir.Results,"/normLRTCERESmerge",PCnumber,".RData"))
load(file=paste0(dir.Results,"/normLRTCCRJmerge",PCnumber,".RData"))

PCnumber<-"CT"
load(file=paste0(dir.Results,"/normLRTCCRmerge",PCnumber,".RData"))
load(file=paste0(dir.Results,"/normLRTCERESmerge",PCnumber,".RData"))
load(file=paste0(dir.Results,"/normLRTCCRJmerge",PCnumber,".RData"))

ccrlrt<-list(ComBat=normLRTCCR,ComBatQN=normLRTCCRQN,ComBat1=normLRTCCR1,ComBat2=normLRTCCR2)
ccrgenes<-getLRTgenes(ccrlrt)
ComBatPam50<-classPerfLineage(CCR_corrected[ccrgenes,],lineagelabel=Blabels)
ComBatQNPam50<-classPerfLineage(CCR_correctedQN[ccrgenes,],lineagelabel=Blabels)
ComBatPC1Pam50<-classPerfLineage(CCR_correctedPC1[ccrgenes,],lineagelabel=Blabels)
ComBatPC2Pam50<-classPerfLineage(CCR_correctedPC2[ccrgenes,],lineagelabel=Blabels)
plotKNN(list(ComBatPam50$CURV,ComBatQNPam50$CURV,ComBatPC1Pam50$CURV,ComBatPC2Pam50$CURV),labels=c("ComBat","ComBat+QN","ComBat+PC1","Combat+PC2"),filename=paste0(dir.Results,"/Figure_Pam50KNN_CCR.pdf"),COLS=curveColours)

ComBatPam50w<-classPerfLineage(CCR_corrected,lineagelabel=Blabels,weights=normLRTCCR[[3]])
ComBatQNPam50w<-classPerfLineage(CCR_correctedQN,lineagelabel=Blabels,weights=normLRTCCRQN[[3]])
ComBatPC1Pam50w<-classPerfLineage(CCR_correctedPC1,lineagelabel=Blabels,weights=normLRTCCR1[[3]])
ComBatPC2Pam50w<-classPerfLineage(CCR_correctedPC2,lineagelabel=Blabels,weights=normLRTCCR2[[3]])
plotKNN(list(ComBatPam50w$CURV,ComBatQNPam50w$CURV,ComBatPC1Pam50w$CURV,ComBatPC2Pam50w$CURV),labels=c("ComBat","ComBat+QN","ComBat+PC1","Combat+PC2"),filename=paste0(dir.Results,"/Figure_Pam50KNN_CCRweighted.pdf"),COLS=curveColours)



cereslrt<-list(ComBat=normLRTCERES,ComBatQN=normLRTCERESQN,ComBat1=normLRTCERES1,ComBat2=normLRTCERES2)
ceresgenes<-getLRTgenes(cereslrt)
ComBatPam50C<-classPerfLineage(CERES_corrected[ceresgenes,],lineagelabel=Blabels)
ComBatQNPam50C<-classPerfLineage(CERES_correctedQN[ceresgenes,],lineagelabel=Blabels)
ComBatPC1Pam50C<-classPerfLineage(CERES_correctedPC1[ceresgenes,],lineagelabel=Blabels)
ComBatPC2Pam50C<-classPerfLineage(CERES_correctedPC2[ceresgenes,],lineagelabel=Blabels)
plotKNN(list(ComBatPam50C$CURV,ComBatQNPam50C$CURV,ComBatPC1Pam50C$CURV,ComBatPC2Pam50C$CURV),labels=c("ComBat","ComBat+QN","ComBat+PC1","Combat+PC2"),filename=paste0(dir.Results,"/Figure_Pam50KNN_CERES.pdf"),COLS=curveColours)

ComBatPam50Cw<-classPerfLineage(CERES_corrected[ceresgenes,],lineagelabel=Blabels,weights=normLRTCERES[[3]])
ComBatQNPam50Cw<-classPerfLineage(CERES_correctedQN[ceresgenes,],lineagelabel=Blabels,weights=normLRTCERESQN[[3]])
ComBatPC1Pam50Cw<-classPerfLineage(CERES_correctedPC1[ceresgenes,],lineagelabel=Blabels,weights=normLRTCERES1[[3]])
ComBatPC2Pam50Cw<-classPerfLineage(CERES_correctedPC2[ceresgenes,],lineagelabel=Blabels,weights=normLRTCERES2[[3]])
plotKNN(list(ComBatPam50Cw$CURV,ComBatQNPam50Cw$CURV,ComBatPC1Pam50Cw$CURV,ComBatPC2Pam50Cw$CURV),labels=c("ComBat","ComBat+QN","ComBat+PC1","Combat+PC2"),filename=paste0(dir.Results,"/Figure_Pam50KNN_CERESweighted.pdf"),COLS=curveColours)


ccrjlrt<-list(ComBat=normLRTCCRJ,ComBatQN=normLRTCCRJQN,ComBat1=normLRTCCRJ1,ComBat2=normLRTCCRJ2)
ccrjgenes<-getLRTgenes(ccrjlrt)
ComBatPam50CJw<-classPerfLineage(CCRJ_corrected,lineagelabel=Blabels,weights=normLRTCCRJ[[3]])
ComBatQNPam50CJw<-classPerfLineage(CCRJ_correctedQN,lineagelabel=Blabels,weights=normLRTCCRJQN[[3]])
ComBatPC1Pam50CJw<-classPerfLineage(CCRJ_correctedPC1,lineagelabel=Blabels,weights=normLRTCCRJ1[[3]])
ComBatPC2Pam50CJw<-classPerfLineage(CCRJ_correctedPC2,lineagelabel=Blabels,weights=normLRTCCRJ2[[3]])
plotKNN(list(ComBatPam50CJw$CURV,ComBatQNPam50CJw$CURV,ComBatPC1Pam50CJw$CURV,ComBatPC2Pam50CJw$CURV),labels=c("ComBat","ComBat+QN","ComBat+PC1","Combat+PC2"),filename=paste0(dir.Results,"/Figure_Pam50KNN_CCRJweighted.pdf"),COLS=curveColours)

ComBatPam50CJ<-classPerfLineage(CCRJ_corrected[ccrjgenes,],lineagelabel=Blabels)
ComBatQNPam50CJ<-classPerfLineage(CCRJ_correctedQN[ccrjgenes,],lineagelabel=Blabels)
ComBatPC1Pam50CJ<-classPerfLineage(CCRJ_correctedPC1[ccrjgenes,],lineagelabel=Blabels)
ComBatPC2Pam50CJ<-classPerfLineage(CCRJ_correctedPC2[ccrjgenes,],lineagelabel=Blabels)
plotKNN(list(ComBatPam50CJ$CURV,ComBatQNPam50CJ$CURV,ComBatPC1Pam50CJ$CURV,ComBatPC2Pam50CJ$CURV),labels=c("ComBat","ComBat+QN","ComBat+PC1","Combat+PC2"),filename=paste0(dir.Results,"/Figure_Pam50KNN_CCRJ.pdf"),COLS=curveColours)


LungCL<-cmp[cmp$tissue=='Lung',]
LungCL<-LungCL[LungCL$model_id%in%colnames(CCR_corrected),]
Llabels<-LungCL$cancer_type
rownames(LungCL)<-LungCL$model_id
names(Llabels)<-LungCL$model_id

ComBatLungw<-classPerfLineage(CCR_corrected,lineagelabel=Llabels,weights=normLRTCCR[[3]])
ComBatQNLungw<-classPerfLineage(CCR_correctedQN,lineagelabel=Llabels,weights=normLRTCCRQN[[3]])
ComBatPC1Lungw<-classPerfLineage(CCR_correctedPC1,lineagelabel=Llabels,weights=normLRTCCR1[[3]])
ComBatPC2Lungw<-classPerfLineage(CCR_correctedPC2,lineagelabel=Llabels,weights=normLRTCCR2[[3]])
plotKNN(list(ComBatLungw$CURV,ComBatQNLungw$CURV,ComBatPC1Lungw$CURV,ComBatPC2Lungw$CURV),labels=c("ComBat","ComBat+QN","ComBat+PC1","Combat+PC2"),filename=paste0(dir.Results,"/Figure_LungKNN_CCRweighted.pdf"),COLS=curveColours)

ComBatLungCw<-classPerfLineage(CERES_corrected,lineagelabel=Llabels,weights=normLRTCERES[[3]])
ComBatQNLungCw<-classPerfLineage(CERES_correctedQN,lineagelabel=Llabels,weights=normLRTCERESQN[[3]])
ComBatPC1LungCw<-classPerfLineage(CERES_correctedPC1,lineagelabel=Llabels,weights=normLRTCERES1[[3]])
ComBatPC2LungCw<-classPerfLineage(CERES_correctedPC2,lineagelabel=Llabels,weights=normLRTCERES2[[3]])
plotKNN(list(ComBatLungCw$CURV,ComBatQNLungCw$CURV,ComBatPC1LungCw$CURV,ComBatPC2LungCw$CURV),labels=c("ComBat","ComBat+QN","ComBat+PC1","Combat+PC2"),filename=paste0(dir.Results,"/Figure_LungKNN_CERESweighted.pdf"),COLS=curveColours)

ComBatLungCJw<-classPerfLineage(CCRJ_corrected,lineagelabel=Llabels,weights=normLRTCCRJ[[3]])
ComBatQNLungCJw<-classPerfLineage(CCRJ_correctedQN,lineagelabel=Llabels,weights=normLRTCCRJQN[[3]])
ComBatPC1LungCJw<-classPerfLineage(CCRJ_correctedPC1,lineagelabel=Llabels,weights=normLRTCCRJ1[[3]])
ComBatPC2LungCJw<-classPerfLineage(CCRJ_correctedPC2,lineagelabel=Llabels,weights=normLRTCCRJ2[[3]])
plotKNN(list(ComBatLungCJw$CURV,ComBatQNLungCJw$CURV,ComBatPC1LungCJw$CURV,ComBatPC2LungCJw$CURV),labels=c("ComBat","ComBat+QN","ComBat+PC1","Combat+PC2"),filename=paste0(dir.Results,"/Figure_LungKNN_CCRJweighted.pdf"),COLS=curveColours)

