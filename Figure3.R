


library(here)

library(splines)
library(psych)
library(sva)
library(ggfortify)
library(e1071)
library(scales)
library(tsne)
library(RColorBrewer)
library(colorspace)
library(preprocessCore)
library(caTools)
seed <- 1809
set.seed(seed)

source("Combat_HKfunctions.R")
dir.Results<-"./Results"
dir.Input<-"/Volumes/GoogleDrive/My Drive/IntegratedCRISPR/ToFigshare/"

load(paste0(dir.Input,"/InstituteColours.Rdata"))
load(file=paste0(dir.Results,"/AllCrisprInfo.Rdata"))


AllCrispr$mn<-AllCrispr$stripped_cell_line_name
OverlapData<-AllCrispr[AllCrispr$source=="Both",c("Sanger_Model_ID","DepMap_ID","stripped_cell_line_name","mn")]
AllCrispr$InBatchCorrection<-"No"

AllCrispr[AllCrispr$Sanger_Model_ID=="SIDM00111","stripped_cell_line_name"]<-"U251"
###Load CRISPRcleanR processed data ###
load(file=paste0(dir.Results,"/SangerOverlap.Rdata"))
load(file=paste0(dir.Results,"/BroadOverlapBroadvC.Rdata"))

inOverlap<-unlist(sapply(colnames(SangerOverlap),function(x) strsplit(x,"---")[[1]][1]))
inOverlap<-toupper(gsub("\\.","",inOverlap))
AllCrispr[make.names(AllCrispr$stripped_cell_line_name)%in%inOverlap,"InBatchCorrection"]<-"Yes"

### Original overlap between CRISPRcleanR processed data ####
OrigCorrection<-combatadj(BroadOverlapBroadvC,SangerOverlap)
CCRcombat<-OrigCorrection$correctNq
OrigCorrectionNeb<-ComBatCP(OrigCorrection$mergedData,batch=unlist(sapply(colnames(OrigCorrection$mergedData),function(x)strsplit(x,"---")[[1]][2])),empBayes=FALSE)

saveRDS(OrigCorrection$corrected,file=paste0(dir.Results,"/OverlapCL_ComBat.Rds"))

PerfJoint<-classPerfCP(OrigCorrection$corrected,bycellline=TRUE)
mergedData<-OrigCorrection$mergedData


SkewWeightsAvg<-CalculateSkew(BroadOverlapBroadvC,SangerOverlap,method="mean")

### weighted correlation argument with ComBat correction only #####
OriginalWeightA<-classPerfCP(OrigCorrection$corrected,weights=abs(SkewWeightsAvg),bycellline=TRUE)

OrigWeightV<-classPerfCP(OrigCorrection$correctV,weights=abs(SkewWeightsAvg),bycellline=TRUE)
OrigWeightNQ<-classPerfCP(OrigCorrection$correctNq,weights=abs(SkewWeightsAvg),bycellline=TRUE)


###Load CCR + JACKS processed data###
load(file=paste0(dir.Results,"/JACKSoverlapBroad.Rdata"))
load(file=paste0(dir.Results,"/JACKSoverlapSanger.Rdata"))

OrigCorrectionJacksD<-combatadj(JBoverlap,JSoverlap)
OrigCorrectionJacksNeb<-ComBatCP(OrigCorrectionJacksD$mergedData,batch=unlist(sapply(colnames(OrigCorrectionJacksD$mergedData),function(x)strsplit(x,"---")[[1]][2])),empBayes=FALSE)
PerfJointJacks<-classPerfCP(OrigCorrectionJacksD$corrected,bycellline = TRUE)

SkewOverJacksAvg<-CalculateSkew(JBoverlap,JSoverlap,method="mean")
Jacksoweight<-classPerfCP(OrigCorrectionJacksD$corrected,weights=abs(SkewOverJacksAvg),bycellline = TRUE)


###Load CERES processed data ###
load(file=paste0(dir.Results,"/CeresBover.Rdata"))
load(file=paste0(dir.Results,"/CeresSover.Rdata"))

OrigCorrectionCERES<-combatadj(CeresBover,CeresSover)
saveRDS(OrigCorrectionCERES$corrected,file=paste0(dir.Results,"/CERES_Overlap_Combat.Rds"))
OrigCorrectionCERESNeb<-ComBatCP(OrigCorrectionCERES$mergedData,batch=unlist(sapply(colnames(OrigCorrectionCERES$mergedData),function(x)strsplit(x,"---")[[1]][2])),empBayes=FALSE)
CEREScombat<-OrigCorrectionCERES$correctNq


PerfJointCERES<-classPerfCP(OrigCorrectionCERES$corrected,bycellline = TRUE)


SkewOverCeres<-CalculateSkew(CeresBover,CeresSover)
SkewOverCeresAvg<-CalculateSkew(CeresBover,CeresSover,method="mean")
PerfCERESNeb<-classPerfCP(OrigCorrectionCERESNeb$correctedData,bycellline = TRUE,weights=abs(SkewOverCeresAvg))

CERESoweight<-classPerfCP(OrigCorrectionCERES$corrected,weights=abs(SkewOverCeresAvg),bycellline = TRUE)


###Screen quality + ComBat + Principal component
###CCR:

dm<-dimnames(SangerOverlap)
SangerOverlapq<-normalize.quantiles(SangerOverlap)
dimnames(SangerOverlapq)<-dm
dm<-dimnames(BroadOverlapBroadvC)
BroadOverlapBroadvCq<-normalize.quantiles(BroadOverlapBroadvC)
dimnames(BroadOverlapBroadvCq)<-dm
SangerMeans<-rowMeans(SangerOverlapq)
BroadMeans<-rowMeans(BroadOverlapBroadvC)

# Build the model
knots <- quantile(SangerMeans, p = c(0.25,0.5,0.75))
for(i in 1:ncol(SangerOverlap)){
  model <- lm (SangerOverlapq[,i] ~ bs(SangerMeans, knots = knots),na.action=na.exclude)

  if(i==1){
    SangerScreenC<-SangerOverlapq[,i]-fitted.values(model)+SangerMeans
  }else{
    SangerScreenC<-cbind(SangerScreenC,SangerOverlapq[,i]-fitted.values(model)+SangerMeans)
  }
  
}

dimnames(SangerScreenC)<-dimnames(SangerOverlapq)

knots <- quantile(BroadMeans, p = c(0.25,0.5,0.75))
for(i in 1:ncol(BroadOverlapBroadvCq)){
  model <- lm (BroadOverlapBroadvCq[,i] ~ bs(BroadMeans, knots = knots),na.action=na.exclude)

  if(i==1){
    BroadScreenC<-BroadOverlapBroadvCq[,i]-fitted.values(model)+BroadMeans
  }else{
    BroadScreenC<-cbind(BroadScreenC,BroadOverlapBroadvCq[,i]-fitted.values(model)+BroadMeans)
  }
  
}

dimnames(BroadScreenC)<-dimnames(BroadOverlapBroadvC)


save(SangerScreenC,file=paste0(dir.Results,"/SQsangerCCRoverlap.Rdata"))
save(BroadScreenC,file=paste0(dir.Results,"/SQbroadCCRoverlap.Rdata"))
SkewWeightsAvgSQ<-CalculateSkew(BroadScreenC,SangerScreenC,method="mean")
OrigCorrectionSC<-combatadj(BroadScreenC,SangerScreenC)
PerfJointSC<-classPerfCP(OrigCorrectionSC$corrected,bycellline=TRUE,weights=abs(SkewWeightsAvgSQ))
PerfJointSCnq<-classPerfCP(OrigCorrectionSC$correctNq,bycellline = TRUE,weights=abs(SkewWeightsAvgSQ))


CcrSQPC<-RemovePC(OrigCorrectionSC$corrected,1:2)

###CERES
dn<-dimnames(CeresSover)
CeresSoverq<-normalize.quantiles(CeresSover)
dimnames(CeresSoverq)<-dn
dn<-dimnames(CeresBover)
CeresBoverq<-normalize.quantiles(CeresBover)
dimnames(CeresBoverq)<-dn
SangerMeansCeres<-rowMeans(CeresSoverq)
BroadMeansCeres<-rowMeans(CeresBoverq)


# Build the model
knots <- quantile(SangerMeansCeres, p = c(0.25,0.5,0.75))
for(i in 1:ncol(CeresSover)){
  model <- lm (CeresSoverq[,i] ~ bs(SangerMeansCeres, knots = knots),na.action=na.exclude)

  if(i==1){
    
    SangerScreenCeres<-CeresSoverq[,i]-fitted.values(model)+SangerMeansCeres
  }else{
    SangerScreenCeres<-cbind(SangerScreenCeres,CeresSoverq[,i]-fitted.values(model)+SangerMeansCeres)
  }
  
}

dimnames(SangerScreenCeres)<-dimnames(CeresSover)

knots <- quantile(BroadMeansCeres, p = c(0.25,0.5,0.75))
for(i in 1:ncol(CeresBover)){
  model <- lm (CeresBoverq[,i] ~ bs(BroadMeansCeres, knots = knots),na.action=na.exclude)
  #model2<-smooth.spline(BroadMeansCeres,CeresBoverq[,i],cv=TRUE)
  if(i==1){
    BroadScreenCeres<-CeresBoverq[,i]-fitted.values(model)+BroadMeansCeres
  }else{
    BroadScreenCeres<-cbind(BroadScreenCeres,CeresBoverq[,i]-fitted.values(model)+BroadMeansCeres)
  }
  
}

dimnames(BroadScreenCeres)<-dimnames(CeresBover)
save(SangerScreenCeres,file=paste0(dir.Results,"/SQsangerCERESoverlap.Rdata"))
save(BroadScreenCeres,file=paste0(dir.Results,"/SQbroadCERESoverlap.Rdata"))
SkewWeightsAvgSQC<-CalculateSkew(BroadScreenCeres,SangerScreenCeres,method="mean")
OrigCorrectionSCeres<-combatadj(BroadScreenCeres,SangerScreenCeres)
PerfJointSCeres<-classPerfCP(OrigCorrectionSCeres$corrected,bycellline=TRUE,weights=abs(SkewWeightsAvgSQC))
OrigWeightCERESSQ<-classPerfCP(OrigCorrectionSCeres$correctV,bycellline=TRUE,weights=abs(SkewWeightsAvgSQC))
OrigWeightCERESSQnq<-classPerfCP(OrigCorrectionSCeres$correctNq,bycellline=TRUE,weights=abs(SkewWeightsAvgSQC))

###JACKS
dn<-dimnames(JSoverlap)
JacksSoverq<-normalize.quantiles(JSoverlap)
dimnames(JacksSoverq)<-dn
dn<-dimnames(JBoverlap)
JacksBoverq<-normalize.quantiles(JBoverlap)
dimnames(JacksBoverq)<-dn
SangerMeansJacks<-rowMeans(JacksSoverq)
BroadMeansJacks<-rowMeans(JacksBoverq)


# Build the model
knots <- quantile(SangerMeansJacks, p = c(0.25,0.5,0.75))
for(i in 1:ncol(JSoverlap)){
  model <- lm (JacksSoverq[,i] ~ bs(SangerMeansJacks, knots = knots),na.action=na.exclude)

  if(i==1){
    
    SangerScreenJacks<-JacksSoverq[,i]-fitted.values(model)+SangerMeansJacks
  }else{
    SangerScreenJacks<-cbind(SangerScreenJacks,JacksSoverq[,i]-fitted.values(model)+SangerMeansJacks)
  }
  
}

dimnames(SangerScreenJacks)<-dimnames(JSoverlap)

knots <- quantile(BroadMeansJacks, p = c(0.25,0.5,0.75))
for(i in 1:ncol(JBoverlap)){
  model <- lm (JacksBoverq[,i] ~ bs(BroadMeansJacks, knots = knots),na.action=na.exclude)

  if(i==1){
    BroadScreenJacks<-JacksBoverq[,i]-fitted.values(model)+BroadMeansJacks
  }else{
    BroadScreenJacks<-cbind(BroadScreenJacks,JacksBoverq[,i]-fitted.values(model)+BroadMeansJacks)
  }
  
}

dimnames(BroadScreenJacks)<-dimnames(JBoverlap)
save(SangerScreenJacks,file=paste0(dir.Results,"/SQsangerJACKSoverlap.Rdata"))
save(BroadScreenJacks,file=paste0(dir.Results,"/SQbroadJACKSoverlap.Rdata"))
SkewWeightsAvgSQJ<-CalculateSkew(BroadScreenJacks,SangerScreenJacks,method="mean")
OrigCorrectionSJacks<-combatadj(BroadScreenJacks,SangerScreenJacks)
PerfJointSJacks<-classPerfCP(OrigCorrectionSJacks$corrected,bycellline=TRUE,weights=abs(SkewWeightsAvgSQJ))
PerfJointSJacksNq<-classPerfCP(OrigCorrectionSJacks$correctNq,bycellline=TRUE,weights=abs(SkewWeightsAvgSQJ))
OrigWeightJACKSSQ<-classPerfCP(OrigCorrectionSJacks$correctV,bycellline=TRUE,weights=abs(SkewWeightsAvgSQJ))



###Remove first 1 or 2 principal components after screen quality and each batch correction method:
##CCR:
#ComBat:
CCRSQctPC4<-RemovePC(OrigCorrectionSC$corrected,1)
PerfCCRSQctPC4<-classPerfCP(CCRSQctPC4$correctedData,weights=abs(SkewWeightsAvgSQ),bycellline = TRUE)
CCRSQctPC2<-RemovePC(OrigCorrectionSC$corrected,1:2)
PerfCCRSQctPC2<-classPerfCP(CCRSQctPC2$correctedData,weights=abs(SkewWeightsAvgSQ),bycellline = TRUE)


##CERES:
#ComBat
CeresSQctPC4<-RemovePC(OrigCorrectionSCeres$corrected,1)
PerfCeresSQctPC4<-classPerfCP(CeresSQctPC4$correctedData,weights=abs(SkewWeightsAvgSQC),bycellline = TRUE)
CeresSQctPC2<-RemovePC(OrigCorrectionSCeres$corrected,1:2)
PerfCeresSQctPC2<-classPerfCP(CeresSQctPC2$correctedData,weights=abs(SkewWeightsAvgSQC),bycellline = TRUE)

##CCR + JACKS:
CCRSQctPC4J<-RemovePC(OrigCorrectionSJacks$corrected,1)
PerfCCRSQctPC4J<-classPerfCP(CCRSQctPC4J$correctedData,weights=abs(SkewWeightsAvgSQJ),bycellline = TRUE)
CCRSQctPC2J<-RemovePC(OrigCorrectionSJacks$corrected,1:2)
PerfCCRSQctPC2J<-classPerfCP(CCRSQctPC2J$correctedData,weights=abs(SkewWeightsAvgSQJ),bycellline = TRUE)


PerfCCRnoCorrect<-classPerfCP(OrigCorrection$mergedData,bycellline=TRUE,weights=abs(SkewWeightsAvg))
PerfCERESnoCorrect<-classPerfCP(OrigCorrectionCERES$mergedData,bycellline = TRUE,weights=abs(SkewOverCeresAvg))

###Figure 3a ###
pcaplot(OrigCorrection$mergedData,filename=paste0(dir.Results,"/Over_CCR_origData"),colours=InstituteColours,width=15,scalefig = FALSE)
pcaplot(OrigCorrectionJacksD$mergedData[rowSums(is.na(OrigCorrectionJacksD$mergedData))==0,],filename=paste0(dir.Results,"/Over_Jacks_origData"),colours=InstituteColours,width=15,scalefig = FALSE)
pcaplot(OrigCorrectionCERES$mergedData[rowSums(is.na(OrigCorrectionCERES$mergedData))==0,],filename=paste0(dir.Results,"/Over_CERES_origData"),colours=InstituteColours,width=15,scalefig=FALSE)

###FIGURE 3B ###

curveColours<-c("#3375A2","#E1822C","#3B9144","#C13E3F")
#CERES
plotKNN(list(OrigWeightCERESSQnq$CURVjoint,PerfJointSCeres$CURVjoint,PerfCeresSQctPC4$CURVjoint,PerfCeresSQctPC2$CURVjoint),xlim=c(1,length(OrigWeightCERESSQnq$CURVjoint)),labels=c("ComBat","ComBat+QN","ComBat+QN+PC1","ComBat+QN+PC1-2"),filename=paste0(dir.Results,"/Figure2b_CERES.pdf"),COLS=curveColours[1:4],ltylist=rep(1,4))
#CRISPRcleanR
plotKNN(list(PerfJointSCnq$CURVjoint,PerfJointSC$CURVjoint,PerfCCRSQctPC4$CURVjoint,PerfCCRSQctPC2$CURVjoint),xlim=c(1,length(PerfJointSCnq$CURVjoint)),labels=c("ComBat","ComBat+QN","ComBat+QN+PC1","ComBat+QN+PC1-2"),filename=paste0(dir.Results,"/Figure2b_CRISPRcleanR.pdf"),COLS=curveColours[1:4],ltylist=rep(1,4))
#CRISPRcleanR + JACKS
plotKNN(list(PerfJointSJacksNq$CURVjoint,PerfJointSJacks$CURVjoint,PerfCCRSQctPC4J$CURVjoint,PerfCCRSQctPC2J$CURVjoint),xlim=c(1,length(PerfJointSJacksNq$CURVjoint)),labels=c("ComBat","ComBat+QN","ComBat+QN+PC1","ComBat+QN+PC1-2"),filename=paste0(dir.Results,"/Figure2b_CRISPRcleanR_JACKS.pdf"),COLS=curveColours[1:4],ltylist=rep(1,4))


plotKNN(list(OrigWeightCERESSQnq$CURVjoint,PerfJointSCeres$CURVjoint,PerfCeresSQctPC4$CURVjoint,PerfCeresSQctPC2$CURVjoint),xlim=c(1,100),labels=c("ComBat","ComBat+QN","ComBat+QN+PC1","ComBat+QN+PC1-2"),filename=paste0(dir.Results,"/Figure2b_CERES_lim.pdf"),COLS=curveColours[1:4],ltylist=rep(1,4))
#CRISPRcleanR
plotKNN(list(PerfJointSCnq$CURVjoint,PerfJointSC$CURVjoint,PerfCCRSQctPC4$CURVjoint,PerfCCRSQctPC2$CURVjoint),xlim=c(1,100),labels=c("ComBat","ComBat+QN","ComBat+QN+PC1","ComBat+QN+PC1-2"),filename=paste0(dir.Results,"/Figure2b_CRISPRcleanR_lim.pdf"),COLS=curveColours[1:4],ltylist=rep(1,4))
#CRISPRcleanR + JACKS
plotKNN(list(PerfJointSJacksNq$CURVjoint,PerfJointSJacks$CURVjoint,PerfCCRSQctPC4J$CURVjoint,PerfCCRSQctPC2J$CURVjoint),xlim=c(1,100),labels=c("ComBat","ComBat+QN","ComBat+QN+PC1","ComBat+QN+PC1-2"),filename=paste0(dir.Results,"/Figure2b_CRISPRcleanR_JACKS_lim.pdf"),COLS=curveColours[1:4],ltylist=rep(1,4))



##Making up Supplementary Fig :
pcaplot(OrigCorrectionSC$corrected,filename=paste0(dir.Results,"/CombatOver_CCR_SQ"),colours=InstituteColours,width=15,scalefig = FALSE)
pcaplot(OrigCorrectionSCeres$corrected,filename=paste0(dir.Results,"/CombatOver_Ceres_SQ"),colours=InstituteColours,width=15,scalefig = FALSE)
pcaplot(OrigCorrectionSJacks$corrected,filename=paste0(dir.Results,"/ComabtOver_Jacks_SQ"),colours=InstituteColours,width=15,scalefig = FALSE)

saveRDS(OrigCorrectionSC$corrected,file=paste0(dir.Results,"/CombatOverlap_CCR_SQ.Rds"))
saveRDS(OrigCorrectionSCeres$corrected,file=paste0(dir.Results,"/CombatOverlap_CERES_SQ.Rds"))
saveRDS(OrigCorrectionSJacks$corrected,file=paste0(dir.Results,"/CombatOverlap_Jacks_SQ.Rds"))

OrigCorrectionSC_corrected<-readRDS(file=paste0(dir.Results,"/CombatOverlap_CCR_SQ.Rds"))
OrigCorrectionSCeres_corrected<-readRDS(file=paste0(dir.Results,"/CombatOverlap_CERES_SQ.Rds"))
OrigCorrectionSJacks_corrected<-readRDS(file=paste0(dir.Results,"/CombatOverlap_Jacks_SQ.Rds"))

pcaplot(OrigCorrectionSC_corrected,filename=paste0(dir.Results,"/CombatOver_CCR_SQ"),colours=InstituteColours,width=15,scalefig = FALSE)
pcaplot(OrigCorrectionSCeres_corrected,filename=paste0(dir.Results,"/CombatOver_Ceres_SQ"),colours=InstituteColours,width=15,scalefig = FALSE)
pcaplot(OrigCorrectionSJacks_corrected,filename=paste0(dir.Results,"/ComabtOver_Jacks_SQ"),colours=InstituteColours,width=15,scalefig = FALSE)


###Generate the fully batch corrected data sets ####

load(paste0(dir.Results,"/BroadDataCERES.Rdata"))
load(paste0(dir.Results,"/SangerDataCERES.Rdata"))
site<-sapply(colnames(OrigCorrectionCERES$mergedData),function(x) strsplit(x,"---",fixed=TRUE)[[1]][2])
CombatOrigCERES<-ComBatCP(as.matrix(OrigCorrectionCERES$mergedData),batch = site[colnames(OrigCorrectionCERES$mergedData)])
CeresAll<-BatchCorrection(data1=BroadDataCERES,data2=SangerDataCERES,CombatRes=CombatOrigCERES)
saveRDS(CeresAll$qNorm,file=paste0(dir.Results,"/CERES_CombatOnly_All_qN.Rds"))
saveRDS(CeresAll$NoNorm,file=paste0(dir.Results,"/CERES_CombatOnly_All_NoNorm.Rds"))


load(paste0(dir.Results,"/BroadData.Rdata"))
load(paste0(dir.Results,"/SangerData.Rdata"))
AllCrispr$source<-"Broad"

AllCrispr[AllCrispr$Sanger_Model_ID%in%colnames(SangerData),"source"]<-"Sanger"
AllCrispr[AllCrispr$DepMap_ID%in%colnames(BroadData)&AllCrispr$Sanger_Model_ID%in%colnames(SangerData),"source"]<-"Both"
###For CCR all data with ComBat###
site<-sapply(colnames(OrigCorrection$mergedData),function(x) strsplit(x,"---",fixed=TRUE)[[1]][2])
CombatOrig<-ComBatCP(as.matrix(OrigCorrection$mergedData),batch = site[colnames(OrigCorrection$mergedData)])
OrigAll<-BatchCorrection(data1=BroadData,data2=SangerData,CombatRes=CombatOrig)

saveRDS(OrigAll$qNorm,file=paste0(dir.Results,"/CCR_CombatOnly_All_qN.Rds"))
saveRDS(OrigAll$NoNorm,file=paste0(dir.Results,"/CCR_CombatOnly_All_NoNorm.Rds"))


##CERES
dn<-dimnames(BroadDataCERES)
BroadDataCERESq<-normalize.quantiles(BroadDataCERES)
dimnames(BroadDataCERESq)<-dn
dn<-dimnames(SangerDataCERES)
SangerDataCERESq<-normalize.quantiles(SangerDataCERES)
dimnames(SangerDataCERESq)<-dn
SangerMeansCeresA<-rowMeans(SangerDataCERESq,na.rm=TRUE)
BroadMeansCeresA<-rowMeans(BroadDataCERESq,na.rm=TRUE)
SangerMeansCeresnq<-rowMeans(SangerDataCERES,na.rm=TRUE)
BroadMeansCeresnq<-rowMeans(BroadDataCERES,na.rm=TRUE)

# Build the model
knots <- quantile(SangerMeansCeresA, p = c(0.25,0.5,0.75),na.rm=T)
knotsnq<-quantile(SangerMeansCeresnq,p=c(0.25,0.5,0.75),na.rm=T)
for(i in 1:ncol(SangerDataCERES)){
  model <- lm (SangerDataCERESq[,i] ~ bs(SangerMeansCeresA, knots = knots),na.action=na.exclude)
  modelnq<-lm(SangerDataCERES[,i]~bs(SangerMeansCeresnq,knots=knotsnq),na.action=na.exclude)

  if(i==1){
    
    SangerScreenCeresA<-SangerDataCERESq[,i]-fitted.values(model)+SangerMeansCeresA
    SangerScreenCeresAnq<-SangerDataCERES[,i]-fitted.values(modelnq)+SangerMeansCeresnq
  }else{
    SangerScreenCeresA<-cbind(SangerScreenCeresA,SangerDataCERESq[,i]-fitted.values(model)+SangerMeansCeresA)
    SangerScreenCeresAnq<-cbind(SangerScreenCeresAnq,SangerDataCERES[,i]-fitted.values(modelnq)+SangerMeansCeresnq)
  }
  
}

dimnames(SangerScreenCeresA)<-dimnames(SangerDataCERESq)
dimnames(SangerScreenCeresAnq)<-dimnames(SangerDataCERES)

knots <- quantile(BroadMeansCeresA, p = c(0.25,0.5,0.75),na.rm=T)
knotsnq<-quantile(BroadMeansCeresnq,p=c(0.25,0.5,0.75),na.rm=T)
for(i in 1:ncol(BroadDataCERES)){
  model <- lm (BroadDataCERESq[,i] ~ bs(BroadMeansCeresA, knots = knots),na.action=na.exclude)
  modelnq<-lm(BroadDataCERES[,i]~bs(BroadMeansCeresnq,knots=knotsnq),na.action=na.exclude)

  if(i==1){
    BroadScreenCeresA<-BroadDataCERESq[,i]-fitted.values(model)+BroadMeansCeresA
    BroadScreenCeresAnq<-BroadDataCERES[,i]-fitted.values(modelnq)+BroadMeansCeresnq
  }else{
    BroadScreenCeresA<-cbind(BroadScreenCeresA,BroadDataCERESq[,i]-fitted.values(model)+BroadMeansCeresA)
    BroadScreenCeresAnq<-cbind(BroadScreenCeresAnq,BroadDataCERES[,i]-fitted.values(modelnq)+BroadMeansCeresnq)
  }
  
}


dimnames(BroadScreenCeresA)<-dimnames(BroadDataCERESq)
dimnames(BroadScreenCeresAnq)<-dimnames(BroadDataCERES)


save(SangerScreenCeresA,file=paste0(dir.Results,"/SQsangerCERESall.Rdata"))
save(BroadScreenCeresA,file=paste0(dir.Results,"/SQbroadCERESall.Rdata"))
qnormScreenQ_ceres<-cbind(BroadScreenCeresA,SangerScreenCeresA[rownames(BroadScreenCeresA),])
site<-sapply(colnames(OrigCorrectionSCeres$mergedData),function(x) strsplit(x,"---",fixed=TRUE)[[1]][2])
CombatOrigCERESsq<-ComBatCP(as.matrix(OrigCorrectionSCeres$mergedData),batch = site[colnames(OrigCorrectionSCeres$mergedData)])
CeresAllsq<-BatchCorrection(data1=BroadScreenCeresA,data2=SangerScreenCeresA,CombatRes=CombatOrigCERESsq)
saveRDS(CeresAllsq$qNorm,file=paste0(dir.Results,"/CERES_SQ_Combat_All.Rds"))
saveRDS(CeresAllsq$NoNorm,file=paste0(dir.Results,"/CERES_SQ_Combat_All_NoNorm.Rds"))



ScreenQ_ceres<-cbind(BroadScreenCeresAnq,SangerScreenCeresAnq[rownames(BroadScreenCeresAnq),])

CeresAllsqnq<-BatchCorrection(data1=BroadScreenCeresAnq,data2=SangerScreenCeresAnq,CombatRes=CombatOrigCERESsq)
saveRDS(CeresAllsqnq$qNorm,file=paste0(dir.Results,"/CERES_SQ_Combat_All_NoPreq.Rds"))
saveRDS(CeresAllsqnq$NoNorm,file=paste0(dir.Results,"/CERES_SQ_Combat_All_NoNorm_NoPreq.Rds"))




CorrectedPCA_allCsq<-RemovePC(CeresAllsq$qNorm,1,perfCheck=FALSE)
saveRDS(CorrectedPCA_allCsq,file=paste0(dir.Results,"/CERES_SQ_Combat_PC1_All.Rds"))
CorrectedPCA_allCsq2<-RemovePC(CeresAllsq$qNorm,1:2,perfCheck=FALSE)
saveRDS(CorrectedPCA_allCsq2,file=paste0(dir.Results,"/CERES_SQ_Combat_PC2_All.Rds"))
CorrectedPCA_allCsqnq<-RemovePC(CeresAllsq$NoNorm,1,perfCheck=FALSE)
saveRDS(CorrectedPCA_allCsqnq,file=paste0(dir.Results,"/CERES_SQ_Combat_PC1_All_NoNorm.Rds"))

CorrectedPCA_allCsqnp<-RemovePC(CeresAllsqnq$qNorm,1,perfCheck=FALSE)
saveRDS(CorrectedPCA_allCsqnp,file=paste0(dir.Results,"/CERES_SQ_Combat_PC1_All_NoPre.Rds"))
CorrectedPCA_allCsqnqnp<-RemovePC(CeresAllsqnq$NoNorm,1,perfCheck=FALSE)
saveRDS(CorrectedPCA_allCsqnqnp,file=paste0(dir.Results,"/CERES_SQ_Combat_PC1_All_NoNorm_NoPreq.Rds"))

##CCR
dn<-dimnames(BroadData)
BroadDataq<-normalize.quantiles(BroadData)
dimnames(BroadDataq)<-dn
dn<-dimnames(SangerData)
SangerDataq<-normalize.quantiles(SangerData)
dimnames(SangerDataq)<-dn
SangerMeansA<-rowMeans(SangerDataq)
BroadMeansA<-rowMeans(BroadDataq)
SangerMeansAnq<-rowMeans(SangerData)
BroadMeansAnq<-rowMeans(BroadData)


# Build the model
knots <- quantile(SangerMeansA, p = c(0.25,0.5,0.75))
knotsnq<-quantile(SangerMeansAnq,p=c(0.25,0.5,0.75))
for(i in 1:ncol(SangerData)){
  model <- lm (SangerDataq[,i] ~ bs(SangerMeansA, knots = knots),na.action=na.exclude)
  modelnq<-lm(SangerData[,i]~bs(SangerMeansAnq,knots=knotsnq),na.action=na.exclude)

  if(i==1){
    
    SangerScreenA<-SangerDataq[,i]-fitted.values(model)+SangerMeansA
    SangerScreenAnq<-SangerData[,i]-fitted.values(modelnq)+SangerMeansAnq
  }else{
    SangerScreenA<-cbind(SangerScreenA,SangerDataq[,i]-fitted.values(model)+SangerMeansA)
    SangerScreenAnq<-cbind(SangerScreenAnq,SangerData[,i]-fitted.values(modelnq)+SangerMeansAnq)
  }
  
}

dimnames(SangerScreenA)<-dimnames(SangerDataq)
dimnames(SangerScreenAnq)<-dimnames(SangerData)

knots <- quantile(BroadMeansA, p = c(0.25,0.5,0.75))
knotsnq<-quantile(BroadMeansAnq,p=c(0.25,0.5,0.75))
for(i in 1:ncol(BroadData)){
  model <- lm (BroadDataq[,i] ~ bs(BroadMeansA, knots = knots),na.action=na.exclude)
  modelnq<-lm(BroadData[,i]~bs(BroadMeansAnq,knots=knotsnq),na.action=na.exclude)

  if(i==1){
    BroadScreenA<-BroadDataq[,i]-fitted.values(model)+BroadMeansA
    BroadScreenAnq<-BroadData[,i]-fitted.values(modelnq)+BroadMeansAnq
  }else{
    BroadScreenA<-cbind(BroadScreenA,BroadDataq[,i]-fitted.values(model)+BroadMeansA)
    BroadScreenAnq<-cbind(BroadScreenAnq,BroadData[,i]-fitted.values(modelnq)+BroadMeansAnq)
  }
  
}

dimnames(BroadScreenA)<-dimnames(BroadDataq)
dimnames(BroadScreenAnq)<-dimnames(BroadData)
save(SangerScreenA,file=paste0(dir.Results,"/SQsangerCCRall.Rdata"))
save(BroadScreenA,file=paste0(dir.Results,"/SQbroadCCRall.Rdata"))
site<-sapply(colnames(OrigCorrectionSC$mergedData),function(x) strsplit(x,"---",fixed=TRUE)[[1]][2])
CombatOrigsq<-ComBatCP(as.matrix(OrigCorrectionSC$mergedData),batch = site[colnames(OrigCorrectionSC$mergedData)])
CcrAllsq<-BatchCorrection(data1=BroadScreenA,data2=SangerScreenA,CombatRes=CombatOrigsq)


saveRDS(CcrAllsq$qNorm,file=paste0(dir.Results,"/CCR_SQ_Combat_All.Rds"))
saveRDS(CcrAllsq$NoNorm,file=paste0(dir.Results,"/CCR_SQ_Combat_All_NoNorm.Rds"))
CcrAllsqnq<-BatchCorrection(data1=BroadScreenAnq,data2=SangerScreenAnq,CombatRes=CombatOrigsq)


saveRDS(CcrAllsqnq$qNorm,file=paste0(dir.Results,"/CCR_SQ_Combat_All_NoPreq.Rds"))
saveRDS(CcrAllsqnq$NoNorm,file=paste0(dir.Results,"/CCR_SQ_Combat_All_NoNorm_NoPreq.Rds"))



CorrectedPCA_allsq<-RemovePC(CcrAllsq$qNorm,1,perfCheck=FALSE)

saveRDS(CorrectedPCA_allsq,file=paste0(dir.Results,"/CCR_SQ_Combat_PC1_All.Rds"))

CorrectedPCA_allsq2<-RemovePC(CcrAllsq$qNorm,1:2,perfCheck=FALSE)

saveRDS(CorrectedPCA_allsq2,file=paste0(dir.Results,"/CCR_SQ_Combat_PC2_All.Rds"))

CorrectedPCA_allsqnq<-RemovePC(CcrAllsq$NoNorm,1,perfCheck=FALSE)

saveRDS(CorrectedPCA_allsqnq,file=paste0(dir.Results,"/CCR_SQ_Combat_PC1_All_NoNorm.Rds"))


CorrectedPCA_allsqnp<-RemovePC(CcrAllsqnq$qNorm,1,perfCheck=FALSE)

saveRDS(CorrectedPCA_allsqnp,file=paste0(dir.Results,"/CCR_SQ_Combat_PC1_All_NoPre.Rds"))

CorrectedPCA_allsqnqnp<-RemovePC(CcrAllsqnq$NoNorm,1,perfCheck=FALSE)

saveRDS(CorrectedPCA_allsqnqnp,file=paste0(dir.Results,"/CCR_SQ_Combat_PC1_All_NoNorm_NoPre.Rds"))

###Load CCR + JACKS processed data###
#JACKSbroadD
load(file=paste0(dir.Results,"/JACKSbroad.Rdata"))
#JACKSsangerD
load(file=paste0(dir.Results,"/JACKSsanger.Rdata"))


site<-sapply(colnames(OrigCorrectionJacksD$mergedData),function(x) strsplit(x,"---",fixed=TRUE)[[1]][2])
CombatOrigJACKS<-ComBatCP(as.matrix(OrigCorrectionJacksD$mergedData),batch = site[colnames(OrigCorrectionJacksD$mergedData)])
JacksAll<-BatchCorrection(data1=JACKSbroadD,data2=JACKSsangerD,CombatRes=CombatOrigJACKS)

JBroad<-as.matrix(JACKSbroadD)
JSanger<-as.matrix(JACKSsangerD)

dn<-dimnames(JBroad)
JACKSBroadDq<-normalize.quantiles(JBroad)
dimnames(JACKSBroadDq)<-dn
dn<-dimnames(JSanger)
JACKSSangerDq<-normalize.quantiles(JSanger)
dimnames(JACKSSangerDq)<-dn
JSangerMeansA<-rowMeans(JACKSSangerDq)
JBroadMeansA<-rowMeans(JACKSBroadDq)



# Build the model
knots <- quantile(JSangerMeansA, p = c(0.25,0.5,0.75))

for(i in 1:ncol(JACKSSangerDq)){
  model <- lm (JACKSSangerDq[,i] ~ bs(JSangerMeansA, knots = knots),na.action=na.exclude)

  if(i==1){
    
    JSangerScreenA<-JACKSSangerDq[,i]-fitted.values(model)+JSangerMeansA

  }else{
    JSangerScreenA<-cbind(JSangerScreenA,JACKSSangerDq[,i]-fitted.values(model)+JSangerMeansA)
  }
  
}

dimnames(JSangerScreenA)<-dimnames(JACKSSangerDq)

knots <- quantile(JBroadMeansA, p = c(0.25,0.5,0.75))

for(i in 1:ncol(JACKSBroadDq)){
  model <- lm (JACKSBroadDq[,i] ~ bs(JBroadMeansA, knots = knots),na.action=na.exclude)

  if(i==1){
    JBroadScreenA<-JACKSBroadDq[,i]-fitted.values(model)+JBroadMeansA
  
  }else{
    JBroadScreenA<-cbind(JBroadScreenA,JACKSBroadDq[,i]-fitted.values(model)+JBroadMeansA)
  }
  
}

dimnames(JBroadScreenA)<-dimnames(JACKSBroadDq)

save(JSangerScreenA,file=paste0(dir.Results,"/SQsangerJACKSall.Rdata"))
save(JBroadScreenA,file=paste0(dir.Results,"/SQbroadJACKSall.Rdata"))
site<-sapply(colnames(OrigCorrectionJacksD$mergedData),function(x) strsplit(x,"---",fixed=TRUE)[[1]][2])
CombatOrigsqJacks<-ComBatCP(as.matrix(OrigCorrectionJacksD$mergedData),batch = site[colnames(OrigCorrectionJacksD$mergedData)])
CcrAllsqJ<-BatchCorrection(data1=JBroadScreenA,data2=JSangerScreenA,CombatRes=CombatOrigsqJacks)
saveRDS(CcrAllsqJ$qNorm,file=paste0(dir.Results,"/CCR_SQ_Combat_All_JACKS.Rds"))
saveRDS(CcrAllsqJ$NoNorm,file=paste0(dir.Results,"/CCR_SQ_Combat_All_NoNorm_JACKS.Rds"))
CorrectedPCA_allsqJ<-RemovePC(CcrAllsqJ$qNorm,1,perfCheck=FALSE)
saveRDS(CorrectedPCA_allsqJ,file=paste0(dir.Results,"/CCR_SQ_Combat_PC1_All_JACKS.Rds"))
CorrectedPCA_allsqJ2<-RemovePC(CcrAllsqJ$qNorm,1:2,perfCheck=FALSE)
saveRDS(CorrectedPCA_allsqJ2,file=paste0(dir.Results,"/CCR_SQ_Combat_PC2_All_JACKS.Rds"))


###All main datasets generated by this point####

###Checking for discordant CRISPR profiles:
CCR_Check<-classPerfCP(OrigCorrectionSC$corrected,bycellline=FALSE,weights=abs(SkewWeightsAvgSQ))
Jacks_Check<-classPerfCP(OrigCorrectionSJacks$corrected,bycellline=FALSE,weights=abs(SkewWeightsAvgSQJ))
Ceres_Check<-classPerfCP(OrigCorrectionSCeres$corrected,bycellline=FALSE,weights=abs(SkewWeightsAvgSQC))

#require screen quality metrics for each cell line:
load(paste0(dir.Input,"/BAGEL_essential.Rdata"))
load(paste0(dir.Input,"/BAGEL_nonEssential.Rdata"))
CCR_nnmdS<-nnmd(SangerScreenC,BAGEL_essential,BAGEL_nonEssential)
CCR_nnmdB<-nnmd(BroadScreenC,BAGEL_essential,BAGEL_nonEssential)
names(CCR_nnmdS)<-sapply(names(CCR_nnmdS),function(x) strsplit(x,'---',fixed=TRUE)[[1]][1])
names(CCR_nnmdB)<-sapply(names(CCR_nnmdB),function(x) strsplit(x,'---',fixed=TRUE)[[1]][1])
CCR_nnmd<-colMeans(rbind(CCR_nnmdS,CCR_nnmdB[names(CCR_nnmdS)]))
CCR_nnmdM<-apply(rbind(CCR_nnmdS,CCR_nnmdB[names(CCR_nnmdS)]),2,max)
CCR_nnmdBest<-apply(rbind(CCR_nnmdS,CCR_nnmdB[names(CCR_nnmdS)]),2,which.min)
CERES_nnmdS<-nnmd(SangerScreenCeres,BAGEL_essential,BAGEL_nonEssential)
CERES_nnmdB<-nnmd(BroadScreenCeres,BAGEL_essential,BAGEL_nonEssential)
names(CERES_nnmdS)<-sapply(names(CERES_nnmdS),function(x) strsplit(x,'---',fixed=TRUE)[[1]][1])
names(CERES_nnmdB)<-sapply(names(CERES_nnmdB),function(x) strsplit(x,'---',fixed=TRUE)[[1]][1])
CERES_nnmd<-colMeans(rbind(CERES_nnmdS,CERES_nnmdB[names(CERES_nnmdS)]))
CERES_nnmdM<-apply(rbind(CERES_nnmdS,CERES_nnmdB[names(CERES_nnmdS)]),2,max)
CERES_nnmdBest<-apply(rbind(CERES_nnmdS,CERES_nnmdB[names(CERES_nnmdS)]),2,which.min)
Jacks_nnmdS<-nnmd(SangerScreenJacks,BAGEL_essential,BAGEL_nonEssential)
Jacks_nnmdB<-nnmd(BroadScreenJacks,BAGEL_essential,BAGEL_nonEssential)
names(Jacks_nnmdS)<-sapply(names(Jacks_nnmdS),function(x) strsplit(x,'---',fixed=TRUE)[[1]][1])
names(Jacks_nnmdB)<-sapply(names(Jacks_nnmdB),function(x) strsplit(x,'---',fixed=TRUE)[[1]][1])
Jacks_nnmd<-colMeans(rbind(Jacks_nnmdS,Jacks_nnmdB[names(Jacks_nnmdS)]))
Jacks_nnmdM<-apply(rbind(Jacks_nnmdS,Jacks_nnmdB[names(Jacks_nnmdS)]),2,max)
Jacks_nnmdBest<-apply(rbind(Jacks_nnmdS,Jacks_nnmdB[names(Jacks_nnmdS)]),2,which.min)
CCR_dist<-CCR_Check$pairwiseCor[names(CCR_nnmd)]
CCR_data<-cbind(CCR_dist,CCR_nnmdM,CCR_nnmdBest)
CERES_dist<-Ceres_Check$pairwiseCor[names(CERES_nnmd)]
CERES_data<-cbind(CERES_dist,CERES_nnmdM,CERES_nnmdBest)
Jacks_dist<-Jacks_Check$pairwiseCor[names(Jacks_nnmd)]
Jacks_data<-cbind(Jacks_dist,Jacks_nnmdM,Jacks_nnmdBest)

matchStatsCCR<-CCR_Check$piestatss
matchCCR<-unique(unlist(sapply(names(matchStatsCCR)[matchStatsCCR=="match"],function(x) strsplit(x,"---",fixed=TRUE)[[1]][1])))

matchStatsCERES<-Ceres_Check$piestatss
matchCERES<-unique(unlist(sapply(names(matchStatsCERES)[matchStatsCERES=="match"],function(x) strsplit(x,"---",fixed=TRUE)[[1]][1])))

matchStatsJacks<-Jacks_Check$piestatss
matchJacks<-unique(unlist(sapply(names(matchStatsJacks)[matchStatsJacks=="match"],function(x) strsplit(x,"---",fixed=TRUE)[[1]][1])))

CCRcorquantile<-quantile(CCR_data[matchCCR,1],0.95)

CEREScorquantile<-quantile(CERES_data[matchCERES,1],0.95)

JACKScorquantile<-quantile(Jacks_data[matchJacks,1],0.95)


CCRnm<-CCR_data[CCR_data[,1]>CCRcorquantile,]
CERESnm<-CERES_data[CERES_data[,1]>CEREScorquantile,]
JACKSnm<-Jacks_data[Jacks_data[,1]>JACKScorquantile,]
write.table(CCR_data,quote=FALSE,file=paste0(dir.Results,"/CRISPRcleanR_DivergentData.tsv"),sep="\t",col.names = NA,row.names = T)
write.table(CERES_data,quote=FALSE,file=paste0(dir.Results,"/CERES_DivergentData.tsv"),sep="\t",col.names = NA,row.names = T)
write.table(Jacks_data,quote=FALSE,file=paste0(dir.Results,"/CRISPRcleanR_JACKS_DivergentData.tsv"),sep="\t",col.names = NA,row.names = T)

CCR_data<-read.table(file=paste0(dir.Results,"/CRISPRcleanR_DivergentData.tsv"),sep="\t",row.names = 1,header=T)
CERES_data<-read.table(file=paste0(dir.Results,"/CERES_DivergentData.tsv"),sep="\t",row.names = 1,header=T)
Jacks_data<-read.table(file=paste0(dir.Results,"/CRISPRcleanR_JACKS_DivergentData.tsv"),sep="\t",row.names = 1,header=T)


rownames(CERESnm)<-AllCrispr[match(rownames(CERESnm),AllCrispr$DepMap_ID),"mn"]
rownames(CERES_data)<-AllCrispr[match(rownames(CERES_data),AllCrispr$DepMap_ID),"mn"]
rownames(CCRnm)<-gsub("\\.","",rownames(CCRnm))
rownames(JACKSnm)<-gsub("\\.","",rownames(JACKSnm))
rownames(Jacks_data)<-gsub("\\.","",rownames(Jacks_data))
rownames(CCR_data)<-gsub("\\.","",rownames(CCR_data))
OverlapData$mn<-gsub("\\.","",OverlapData$mn)
divergentCL<-intersect(rownames(CCRnm),intersect(rownames(CERESnm),rownames(JACKSnm)))
write.table(divergentCL,quote=FALSE,file=paste0(dir.Results,"/DivergentCellLines.tsv"),sep="\t",row.names=F,col.names=F)
divergentCL<-read.table(file=paste0(dir.Results,"/DivergentCellLines.tsv"),sep="\t",stringsAsFactors = F)
divergentCL<-unlist(divergentCL)
JacksDiverge<-Jacks_data[divergentCL,]
CCRDiverge<-CCR_data[divergentCL,]
CERESDiverge<-CERES_data[divergentCL,]

AllQCsel<-cbind(JacksDiverge[,3],CCRDiverge[,3],CERESDiverge[,3])
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
UseCL_qc<-apply(AllQCsel,1,getmode)
JacksDiverge[,3]<-UseCL_qc
CCRDiverge[,3]<-UseCL_qc
CERESDiverge[,3]<-UseCL_qc

###Merged data sets ###
rm12<-read.csv(file=paste0(dir.Results,'Drop12NABroadCL.csv'))

CcrAllsqJm<-avgOverlapBroadAnnot(CcrAllsqJ$qNorm,OverlapData,JacksDiverge)
rm12<-intersect(colnames(CcrAllsqJm),rm12[,1])
submat<-CcrAllsqJm[,rm12]
submat[is.na(BroadDataCERES[rownames(submat),rm12])]<-NA
CcrAllsqJm[,rm12]<-submat
saveRDS(CcrAllsqJm,file=paste0(dir.Results,"/CCR_SQ_Combat_All_JACKS_merge_F.Rds"))

CcrAllsqnJm<-avgOverlapBroadAnnot(CcrAllsqJ$NoNorm,OverlapData,JacksDiverge)
submat<-CcrAllsqnJm[,rm12]
submat[is.na(BroadDataCERES[rownames(submat),rm12])]<-NA
CcrAllsqnJm[,rm12]<-submat
saveRDS(CcrAllsqnJm,file=paste0(dir.Results,"/CCR_SQ_Combat_All_NoNorm_JACKS_merge_F.Rds"))

CorrectedPCA_allsqJm<-avgOverlapBroadAnnot(CorrectedPCA_allsqJ,OverlapData,JacksDiverge)
submat<-CorrectedPCA_allsqJm[,rm12]
submat[is.na(BroadDataCERES[rownames(submat),rm12])]<-NA
CorrectedPCA_allsqJm[,rm12]<-submat
saveRDS(CorrectedPCA_allsqJm,file=paste0(dir.Results,"/CCR_SQ_Combat_PC1_All_JACKS_merge_F.Rds"))

CorrectedPCA_allsqJ2m<-avgOverlapBroadAnnot(CorrectedPCA_allsqJ2,OverlapData,JacksDiverge)
submat<-CorrectedPCA_allsqJ2m[,rm12]
submat[is.na(BroadDataCERES[rownames(submat),rm12])]<-NA
CorrectedPCA_allsqJ2m[,rm12]<-submat
saveRDS(CorrectedPCA_allsqJ2m,file=paste0(dir.Results,"/CCR_SQ_Combat_PC2_All_JACKS_merge_F.Rds"))

CcrAllsqm<-avgOverlapBroadAnnot(CcrAllsq$qNorm,OverlapData,CCRDiverge)
submat<-CcrAllsqm[,rm12]
submat[is.na(BroadDataCERES[rownames(submat),rm12])]<-NA
CcrAllsqm[,rm12]<-submat
saveRDS(CcrAllsqm,file=paste0(dir.Results,"/CCR_SQ_Combat_All_merge_F.Rds"))

CcrAllsqnm<-avgOverlapBroadAnnot(CcrAllsq$NoNorm,OverlapData,CCRDiverge)
submat<-CcrAllsqnm[,rm12]
submat[is.na(BroadDataCERES[rownames(submat),rm12])]<-NA
CcrAllsqnm[,rm12]<-submat
saveRDS(CcrAllsqnm,file=paste0(dir.Results,"/CCR_SQ_Combat_All_NoNorm_merge_F.Rds"))

CorrectedPCA_allsqm<-avgOverlapBroadAnnot(CorrectedPCA_allsq,OverlapData,CCRDiverge)
submat<-CorrectedPCA_allsqm[,rm12]
submat[is.na(BroadDataCERES[rownames(submat),rm12])]<-NA
CorrectedPCA_allsqm[,rm12]<-submat
saveRDS(CorrectedPCA_allsqm,file=paste0(dir.Results,"/CCR_SQ_Combat_PC1_All_merge_F.Rds"))

CorrectedPCA_allsq2m<-avgOverlapBroadAnnot(CorrectedPCA_allsq2,OverlapData,CCRDiverge)
submat<-CorrectedPCA_allsq2m[,rm12]
submat[is.na(BroadDataCERES[rownames(submat),rm12])]<-NA
CorrectedPCA_allsq2m[,rm12]<-submat
saveRDS(CorrectedPCA_allsq2m,file=paste0(dir.Results,"/CCR_SQ_Combat_PC2_All_merge_F.Rds"))


CeresAllsqm<-avgOverlapBroadAnnot(CeresAllsq$qNorm,OverlapData,CERESDiverge)
BroadC<-intersect(colnames(CeresAllsqm),colnames(BroadDataCERES))
submat<-CeresAllsqm[,BroadC]
submat[is.na(BroadDataCERES[rownames(submat),BroadC])]<-NA
CeresAllsqm[,BroadC]<-submat
saveRDS(CeresAllsqm,file=paste0(dir.Results,"/CERES_SQ_Combat_All_merge_F.Rds"))

CeresAllsqnm<-avgOverlapBroadAnnot(CeresAllsq$NoNorm,OverlapData,CERESDiverge)
submat<-CeresAllsqnm[,BroadC]
submat[is.na(BroadDataCERES[rownames(submat),BroadC])]<-NA
CeresAllsqnm[,BroadC]<-submat
saveRDS(CeresAllsqnm,file=paste0(dir.Results,"/CERES_SQ_Combat_All_NoNorm_merge_F.Rds"))

CorrectedPCA_allCsqm<-avgOverlapBroadAnnot(CorrectedPCA_allCsq,OverlapData,CERESDiverge)
submat<-CorrectedPCA_allCsqm[,BroadC]
submat[is.na(BroadDataCERES[rownames(submat),BroadC])]<-NA
CorrectedPCA_allCsqm[,BroadC]<-submat
saveRDS(CorrectedPCA_allCsqm,file=paste0(dir.Results,"/CERES_SQ_Combat_PC1_All_merge_F.Rds"))

CorrectedPCA_allCsq2m<-avgOverlapBroadAnnot(CorrectedPCA_allCsq2,OverlapData,CERESDiverge)
submat<-CorrectedPCA_allCsq2m[,BroadC]
submat[is.na(BroadDataCERES[rownames(submat),BroadC])]<-NA
CorrectedPCA_allCsq2m[,BroadC]<-submat
saveRDS(CorrectedPCA_allCsq2m,file=paste0(dir.Results,"/CERES_SQ_Combat_PC2_All_merge_F.Rds"))


##Supplementary Data 1:
##From CMP:
cmp<-read.csv(paste0(dir.Input,"/model_list_latest.csv"),header=T,stringsAsFactors = F)
rownames(cmp)<-cmp$model_id

temp<-sapply(CCRDiverge[,3],function(x) ifelse(x==1,"Sanger","Broad"))

write.table(CCRDiverge,file=paste0(dir.Results,"/DivergentCLandStats.tsv"),row.names=F,quote=FALSE,sep="\t")

AllCrispr[match(rownames(CCRDiverge),AllCrispr$mn),"source"]<-temp
write.table(AllCrispr,file="AllCrispr.tsv",sep="\t",quote=F,row.names = T,col.names=NA)


