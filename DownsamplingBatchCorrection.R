library(here)
source("Combat_HKfunctions.R")
source("DownSampleFunctions.R")
source("LineageFunctions.R")
library(psych)
library(sva)
library(ggfortify)
library(e1071)
library(scales)
library(stringr)
library(splines2)
library(RColorBrewer)
library(colorspace)
library(preprocessCore)
library(aricode)
library(cluster)
library(doParallel)
registerDoParallel(cores=24)
seed <- 1809
set.seed(seed)

dir.Results<-"./Results"
dir.Input<-"/path/to/Figshare"


saveObj<-function(objTS,filePrefix=dir.Results){
  save(objTS,file=paste0(filePrefix,'/',deparse(substitute(objTS)),".Rdata"))
  
}

cmp<-read.csv(paste0(dir.Input,"/model_list_latest.csv"),header=T,stringsAsFactors = F)
cmp$mn<-make.names(cmp$model_name)
cmpAnnot<-cmp
cmp2<-cmp
cmp2$model_id<-cmp2$BROAD_ID
cmpAnnot<-rbind(cmpAnnot,cmp2)
cmpAnnot<-cmpAnnot[!duplicated(cmpAnnot$model_id),]
rownames(cmpAnnot)<-cmpAnnot$model_id
B20Q2<-read.csv(paste0(dir.Input,"/Achilles_replicate_map.csv"),header=T,stringsAsFactors = F)

#this is fine with cmp incomplete as cmp does include Broad_ids for all crispr'd cell lines
OverlapData<-cmp[cmp$BROAD_ID%in%B20Q2$DepMap_ID,]
OverlapData<-OverlapData[OverlapData$BROAD_ID!="",]
OverlapData<-OverlapData[OverlapData$crispr_ko_data=="True",]
OverlapData$bmn<-make.names(OverlapData$BROAD_ID)
OverlapData$Bmn<-paste0(OverlapData$mn,"---Broad")
OverlapData$Smn<-paste0(OverlapData$mn,"---Sanger")
rmcl<-read.csv(file=paste0(dir.Results,'/Drop12NABroadCL.csv'),stringsAsFactors = FALSE)

OverlapData<-OverlapData[!OverlapData$BROAD_ID%in%rmcl[,1],]

###Load CRISPRcleanR and SQ processed data ###
load(file=paste0(dir.Results,"/SQsangerCCRoverlap.Rdata"))
load(file=paste0(dir.Results,"/SQbroadCCRoverlap.Rdata"))

OverlapData<-OverlapData[OverlapData$Smn%in%colnames(SangerScreenC),]
### Original overlap between CRISPRcleanR processed data ####
OrigCorrection<-combatadj(BroadScreenC,SangerScreenC)

mergedData<-OrigCorrection$mergedData

#load in the full dataset SQ processed with CCR:
load(file=paste0(dir.Results,"/SQsangerCCRall.Rdata"))
load(file=paste0(dir.Results,"/SQbroadCCRall.Rdata"))

#all data corrected:
CcrAllsq<-readRDS(file=paste0(dir.Results,"/CCR_SQ_Combat_All.Rds"))
DownSamples=seq(0.05,0.8,by=0.05)
corCCR1<-list()
eucCCR1<-list()
LineageCCR1<-list()
BatchCCR1<-list()

multiResultClass <- function(corR=NULL,euc=NULL,Lineage=NULL,Batch=NULL)
{
  me <- list(
    corR = corR,
    euc = euc,
    Lineage=Lineage,
    Batch=Batch
  )
  
  ## Set the name for the class
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}


###Load CCR + JACKS processed data###

load(file=paste0(dir.Results,"/SQsangerJACKSoverlap.Rdata"))
load(file=paste0(dir.Results,"/SQbroadJACKSoverlap.Rdata"))

OrigCorrectionJacksD<-combatadj(BroadScreenJacks,SangerScreenJacks)

#all data:
load(file=paste0(dir.Results,"/SQsangerJACKSall.Rdata"))
load(file=paste0(dir.Results,"/SQbroadJACKSall.Rdata"))

CcrAllsqJ<-readRDS(file=paste0(dir.Results,"/CCR_SQ_Combat_All_JACKS.Rds"))



###Load CERES processed data ###
load(file=paste0(dir.Results,"/SQsangerCERESoverlap.Rdata"))
load(file=paste0(dir.Results,"/SQbroadCERESoverlap.Rdata"))

#change so that each correction is applied to all data and then there will also be no difference in numbers of correlations across the down samples
OrigCorrectionCERES<-combatadj(BroadScreenCeres,SangerScreenCeres)

#for CERES need different annotation file:
OverlapDataCeres<-OverlapData
OverlapDataCeres$Smn<-paste0(OverlapData$BROAD_ID,"---Sanger")
OverlapDataCeres$Bmn<-paste0(OverlapData$BROAD_ID,"---Broad")
OverlapDataCeres$model_id<-OverlapDataCeres$BROAD_ID

#All data SQ corrected and CERES processed
load(file=paste0(dir.Results,"/SQsangerCERESall.Rdata"))
load(file=paste0(dir.Results,"/SQbroadCERESall.Rdata"))

#result with all data:
CeresAllsq<-readRDS(file=paste0(dir.Results,"/CERES_SQ_Combat_All.Rds"))


DownSamplesL=c(0.3,0.6,0.9)
multiResultClass2 <- function(corR=NULL,euc=NULL,Batch=NULL)
{
  me <- list(
    corR = corR,
    euc = euc,
    Batch=Batch
  )
  
  ## Set the name for the class
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}
corCCR2<-list()
eucCCR2<-list()


OverlapDataSubtype<-OverlapData[OverlapData$tissue%in%c("Lung"),]
BatchCCR2<-foreach(i=1:length(DownSamplesL))%dopar%{
  print(paste("CCR2",i))
  result<-multiResultClass2()
  CombatDS_CCR<-batch_downsample(mergedData,empBayes=TRUE,annotation=OverlapDataSubtype,numberSamples=50,DownRange=DownSamplesL[i],allBroad=BroadScreenA,allSanger=SangerScreenA,OrigCorrection=CcrAllsq,cmpAnnot)
  
  result$corR<-CombatDS_CCR$corRes
  result$euc<-CombatDS_CCR$eucRes

  result$Batch<-CombatDS_CCR$clusterB
  return(result)
}
corCCR2<-lapply(BatchCCR2,function(x) x$corR)
eucCCR2<-lapply(BatchCCR2,function(x) x$euc)
BatchCCR2<-lapply(BatchCCR2,function(x) x$Batch)
saveObj(corCCR2)
saveObj(eucCCR2)

saveObj(BatchCCR2)


corCCRJ2<-list()
eucCCRJ2<-list()

BatchCCRJ2<-list()
OverlapDataSubtype<-OverlapData[OverlapData$tissue%in%c("Lung"),]
BatchCCRJ2<-foreach(i=1:length(DownSamplesL))%dopar%{
print(paste("CCRJ2",i))
  result<-multiResultClass2()
 CombatDS_CCRJ<-batch_downsample(OrigCorrectionJacksD$mergedData,empBayes=TRUE,annotation=OverlapDataSubtype,numberSamples=50,DownRange=DownSamplesL[i],allBroad=JBroadScreenA,allSanger=JSangerScreenA,OrigCorrection=CcrAllsqJ,cmpAnnot)
  
  result$corR<-CombatDS_CCRJ$corRes
  result$euc<-CombatDS_CCRJ$eucRes
  
  result$Batch<-CombatDS_CCRJ$clusterB
  return(result)
}
corCCRJ2<-lapply(BatchCCRJ2,function(x) x$corR)
eucCCRJ2<-lapply(BatchCCRJ2,function(x) x$euc)
BatchCCRJ2<-lapply(BatchCCRJ2,function(x) x$Batch)
saveObj(corCCRJ2)
saveObj(eucCCRJ2)

saveObj(BatchCCRJ2)

corCERES2<-list()
eucCERES2<-list()

OverlapDataCeresSub<-OverlapDataCeres[OverlapDataCeres$tissue=="Lung",]
BatchCERES2<-foreach(i=1:length(DownSamplesL))%dopar%{
  print(paste("CERES2",i))
  result<-multiResultClass2()
  CombatDS_CERES<-batch_downsample(OrigCorrectionCERES$mergedData,empBayes=TRUE,annotation=OverlapDataCeresSub,numberSamples=50,DownRange=DownSamplesL[i],allBroad=BroadScreenCeresA,allSanger=SangerScreenCeresA,OrigCorrection=CeresAllsq,cmpAnnot)
  
  result$corR<-CombatDS_CERES$corRes
  result$euc<-CombatDS_CERES$eucRes

  result$Batch<-CombatDS_CERES$clusterB
  return(result)
}
corCERES2<-lapply(BatchCERES2,function(x) x$corR)
eucCERES2<-lapply(BatchCERES2,function(x) x$euc)
BatchCERES2<-lapply(BatchCERES2,function(x) x$Batch)
saveObj(corCERES2)
saveObj(eucCERES2)

saveObj(BatchCERES2)


