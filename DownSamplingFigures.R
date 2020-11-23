library(here)
source("Combat_HKfunctions.R")
source("DownSampleFunctions.R")
source("LineageFunctions.R")

library(ggfortify)
library(RColorBrewer)
library(colorspace)

seed <- 1809
set.seed(seed)

dir.Results<-"./Results"
dir.Input<-"/path/to/Figshare"


loadObj<-function(objTS,filePrefix=dir.Results){
  load(file=paste0(filePrefix,'/',objTS,".Rdata"))
  return(objTS)
}
load(paste0(dir.Input,"/PipelineColours.Rdata"))

corCCR1<-loadObj("corCCR1")

DownSamples=seq(0.05,0.8,by=0.05)
plotDS(corCCR1,nsamples=168,plotname=paste0(dir.Results,"/CompareDS_CCR_ComBat.pdf"),ylim=c(0.95,1),DownRange = DownSamples)

corCCRJ1<-loadObj("corCCRJ1")

plotDS(corCCRJ1,nsamples=168,plotname=paste0(dir.Results,"/CompareDS_CCRJ_ComBat.pdf"),ylim=c(0.95,1),DownRange = DownSamples)


corCERES1<-loadObj("corCERES1")

plotDS(corCERES1,nsamples=168,plotname=paste0(dir.Results,"/CompareDS_CERES_ComBat.pdf"),ylim=c(0.95,1),DownRange = DownSamples)

allCor<-list(corCCR1,corCCRJ1,corCERES1)
setnames<-c("CRISPRcleanR","CCR-JACKS","CERES")
names(PipelineColours)<-setnames
noCombat<-list()
load(file=paste0(dir.Results,"/SQsangerCCRall.Rdata"))
load(file=paste0(dir.Results,"/SQbroadCCRall.Rdata"))
AllCCR<-cbind(BroadScreenA,SangerScreenA[rownames(BroadScreenA),])
#all data corrected:
CcrAllsq<-readRDS(file=paste0(dir.Results,"/CCR_SQ_Combat_All.Rds"))

noCombat[["CCR"]]<-sapply(colnames(AllCCR),function(x) cor(AllCCR[,x],CcrAllsq[,x],use="pairwise"))
load(file=paste0(dir.Results,"/SQsangerJACKSall.Rdata"))
load(file=paste0(dir.Results,"/SQbroadJACKSall.Rdata"))
AllCCRJ<-cbind(JBroadScreenA,JSangerScreenA[rownames(JBroadScreenA),])
CcrAllsqJ<-readRDS(file=paste0(dir.Results,"/CCR_SQ_Combat_All_JACKS.Rds"))

noCombat[["CCRJ"]]<-sapply(colnames(AllCCRJ),function(x) cor(AllCCRJ[,x],CcrAllsqJ[,x],use="pairwise"))

load(file=paste0(dir.Results,"/SQsangerCERESall.Rdata"))
load(file=paste0(dir.Results,"/SQbroadCERESall.Rdata"))
AllCERES<-cbind(BroadScreenCeresA,SangerScreenCeresA[rownames(BroadScreenCeresA),])
#result with all data:
CeresAllsq<-readRDS(file=paste0(dir.Results,"/CERES_SQ_Combat_All.Rds"))
noCombat[["CERES"]]<-sapply(colnames(AllCERES),function(x) cor(AllCERES[,x],CeresAllsq[,x],use="pairwise"))

plotDS(allCor,nsamples=168,plotname=paste0(dir.Results,"/CompareDS_All_ComBat_withN0.pdf"),ylim=c(0.94,1),DownRange = DownSamples,setNames = setnames,listRes=TRUE,noAdj=noCombat,ylab="Pearson Correlation",plotcolours=PipelineColours,pwidth=15,pheight=12)
plotDS(allCor,nsamples=168,plotname=paste0(dir.Results,"/CompareDS_All_ComBat.pdf"),ylim=c(0.98,1),DownRange = DownSamples,setNames = setnames,listRes=TRUE,ylab="Pearson Correlation",plotcolours=PipelineColours,pwidth=15,pheight=12)



noCombatBatch<-list()
ll<-c(rep(1,ncol(BroadScreenA)),rep(2,ncol(SangerScreenA)))
names(ll)<-colnames(AllCCR)
noCombatBatch[["CCR"]]<-ASWpc(AllCCR,numberPCs = 20,lineagelabels=ll)
noCombatBatch[["CCRJ"]]<-ASWpc(AllCCRJ,numberPCs = 20,lineagelabels=ll)
noCombatBatch[["CERES"]]<-ASWpc(AllCERES,numberPCs = 20,lineagelabels=ll)



BatchCCR1<-loadObj("BatchCCR1")
BatchCCRJ1<-loadObj("BatchCCRJ1")
BatchCERES1<-loadObj("BatchCERES1")
allBatch<-list(BatchCCR1,BatchCCRJ1,BatchCERES1)
setnames<-c("CRISPRcleanR","CCR-JACKS","CERES")
plotDS(allBatch,nsamples=168,plotname=paste0(dir.Results,"/CompareDS_All_ComBat_Batch.pdf"),ylim=c(0,0.3),DownRange = DownSamples,setNames = setnames,listRes=TRUE,ylab="Average Silhouette width",plotcolours=PipelineColours,pwidth=15,pheight=12)
plotDS(allBatch,nsamples=168,plotname=paste0(dir.Results,"/CompareDS_All_ComBat_Batch_noCombat.pdf"),ylim=c(0,0.75),DownRange = DownSamples,setNames = setnames,listRes=TRUE,noAdj=noCombatBatch,ylab="Average Silhouette width",plotcolours=PipelineColours,pwidth=15,pheight=12)


#add in results using mixture of cell lines as well:
DownSamplesL=c(0.3,0.6,0.9)
BatchCCR2<-loadObj("BatchCCR2")
BatchCCRJ2<-loadObj("BatchCCRJ2")
BatchCERES2<-loadObj("BatchCERES2")
allBatch2<-list(BatchCCR2,BatchCCRJ2,BatchCERES2)
setnames<-c("CRISPRcleanR","CCR-JACKS","CERES")
plotDS(allBatch2,nsamples=28,plotname=paste0(dir.Results,"/CompareDS_Lung_ComBat_Batch.pdf"),ylim=c(0,0.3),DownRange = DownSamplesL,setNames = setnames,listRes=TRUE,plotcolours = PipelineColours,pwidth=10,pheight=10)

corCCR2<-loadObj("corCCR2")
corCCRJ2<-loadObj("corCCRJ2")
corCERES2<-loadObj("corCERES2")
allCor2<-list(corCCR2,corCCRJ2,corCERES2)
setnames<-c("CRISPRcleanR","CCR-JACKS","CERES")
plotDS(allCor2,nsamples=28,plotname=paste0(dir.Results,"/CompareDS_Lung_ComBat_Cor.pdf"),ylim=c(0.98,1),DownRange = DownSamplesL,setNames = setnames,listRes=TRUE,plotcolours = PipelineColours,pwidth=10,pheight=10)

#Do a subset plot for all cell lines and then put them side by side
allCorSub<-list(corCCR1[1:3],corCCRJ1[1:3],corCERES1[1:3])
plotDS(allCorSub,nsamples=168,plotname=paste0(dir.Results,"/CompareDS_All_ComBat_sub.pdf"),ylim=c(0.98,1),DownRange = DownSamples[1:3],setNames = setnames,listRes=TRUE,ylab="Pearson Correlation",plotcolours=PipelineColours,pwidth=10,pheight=10)

allBatchSub<-list(BatchCCR1[1:3],BatchCCRJ1[1:3],BatchCERES1[1:3])
plotDS(allBatchSub,nsamples=168,plotname=paste0(dir.Results,"/CompareDS_All_ComBat_Batch_sub.pdf"),ylim=c(0,0.3),DownRange = DownSamples[1:3],setNames = setnames,listRes=TRUE,ylab="Average Silhouette width",plotcolours=PipelineColours,pwidth=10,pheight=10)


