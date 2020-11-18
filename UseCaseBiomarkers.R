##Run CalculateNormLRT script before running this script.
library(here)
library(ggplot2)

library(MASS)
library(sn)
source("Combat_HKfunctions.R")
source("BiomarkerFunctions.R")

dir.MergeFile<-"./Results"
dir.Results<-"./ResultsFilter"
#dir.Input<-"/path/to/downloaded/figshare/"

load(paste0(dir.Input,"/InstituteColours.Rdata"))
load(paste0(dir.Input,"/PipelineColours.Rdata"))
cmp<-read.csv(paste0(dir.Input,"/model_list_latest.csv"),header=T,stringsAsFactors = F)
cmp2<-cmp
cmp2$model_id<-cmp2$BROAD_ID
cmp<-rbind(cmp,cmp2)
PCnumber<-2
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

ccrnt<-normLRTCCR[[3]]
ceresnt<-normLRTCERES[[3]]
ccrjnt<-normLRTCCRJ[[3]]
CTlist<-list(ccrnt,ceresnt,ccrjnt)
MakeNormLRTplot(CTlist,"CT")

ccrnt<-normLRTCCRQN[[3]]
ceresnt<-normLRTCERESQN[[3]]
ccrjnt<-normLRTCCRJQN[[3]]
CTlist<-list(ccrnt,ceresnt,ccrjnt)
MakeNormLRTplot(CTlist,"QN")

ccrnt<-normLRTCCR1[[3]]
ceresnt<-normLRTCERES1[[3]]
ccrjnt<-normLRTCCRJ1[[3]]
CTlist<-list(ccrnt,ceresnt,ccrjnt)
MakeNormLRTplot(CTlist,"PC1")

ccrnt<-normLRTCCR2[[3]]
ceresnt<-normLRTCERES2[[3]]
ccrjnt<-normLRTCCRJ2[[3]]
CTlist<-list(ccrnt,ceresnt,ccrjnt)
MakeNormLRTplot(CTlist,"PC2")


## loading Multiomic Cancer Functional Event binary matrix (MoBEM)
## and selecting columns corresponding to the cell lines in the inventory
load(paste0(dir.Input,'/MoBEM.RData'))
MSS<-cmp$msi_status[match(colnames(MoBEM),cmp$COSMIC_ID)]
MSS[MSS=='MSI']<-1
MSS[MSS!='1']<-0
MSS<-as.numeric(MSS)
MSS[is.na(MSS)]<-0
MoBEM<-rbind(MoBEM,MSS)
rownames(MoBEM)[nrow(MoBEM)]<-'MSI'

## decode CNA identifiers in the MoBEM

load(paste0(dir.Input,"/CNAdecode.RData"))

MoBEM<-decodeCNA_cp(MoBEM)
MoBEM<-unique(MoBEM)
rownames(MoBEM)<-make.names(rownames(MoBEM),unique=TRUE)
PCnumber<-"CT"
load(file=paste0(dir.Results,"/normLRTCCRmerge",PCnumber,".RData"))
load(file=paste0(dir.Results,"/normLRTCERESmerge",PCnumber,".RData"))
load(file=paste0(dir.Results,"/normLRTCCRJmerge",PCnumber,".RData"))
PCnumber<-"QN"
load(file=paste0(dir.Results,"/normLRTCCRmerge",PCnumber,".RData"))
load(file=paste0(dir.Results,"/normLRTCERESmerge",PCnumber,".RData"))
load(file=paste0(dir.Results,"/normLRTCCRJmerge",PCnumber,".RData"))
PCnumber<-"1"
load(file=paste0(dir.Results,"/normLRTCCRmerge",PCnumber,".RData"))
load(file=paste0(dir.Results,"/normLRTCERESmerge",PCnumber,".RData"))
load(file=paste0(dir.Results,"/normLRTCCRJmerge",PCnumber,".RData"))
PCnumber<-"2"
load(file=paste0(dir.Results,"/normLRTCCRmerge",PCnumber,".RData"))
load(file=paste0(dir.Results,"/normLRTCERESmerge",PCnumber,".RData"))
load(file=paste0(dir.Results,"/normLRTCCRJmerge",PCnumber,".RData"))

normLRTlist<-list(normLRTCCR,normLRTCERES,normLRTCCRJ)
ssdCT<-getSDDgenes(normLRTlist)
normLRTlist<-list(normLRTCCRQN,normLRTCERESQN,normLRTCCRJQN)
ssdQN<-getSDDgenes(normLRTlist = normLRTlist)
normLRTlist<-list(normLRTCCR1,normLRTCERES1,normLRTCCRJ1)
ssdPC1<-getSDDgenes(normLRTlist = normLRTlist)
normLRTlist<-list(normLRTCCR2,normLRTCERES2,normLRTCCRJ2)
ssdPC2<-getSDDgenes(normLRTlist = normLRTlist)

ssdgenes<-unique(c(ssdCT,ssdQN,ssdPC1,ssdPC2))

load(paste0(dir.Input,"/PipelineColours.Rdata"))
curveColours<-c("#3375A2","#E1822C","#3B9144","#C13E3F")
names(curveColours)<-c("ComBat",'ComBatQN','ComBatPC1','ComBatPC2')

preProc<-c("CCR","CERES","CCRJ")

#tissue specific biomarkers:

FClist<-list(CCR_corrected,CERES_corrected,CCRJ_corrected)
AnovaResCTt<-BiomarkerTissue(MoBEM,FClist,cmp,"CT",ssdgenes)


FClist<-list(CCR_correctedQN,CERES_correctedQN,CCRJ_correctedQN)
AnovaResQNt<-BiomarkerTissue(MoBEM,FClist,cmp,"QN",ssdgenes)


FClist<-list(CCR_correctedPC1,CERES_correctedPC1,CCRJ_correctedPC1)
AnovaResPC1t<-BiomarkerTissue(MoBEM,FClist,cmp,"PC1",ssdgenes)


FClist<-list(CCR_correctedPC2,CERES_correctedPC2,CCRJ_correctedPC2)
AnovaResPC2t<-BiomarkerTissue(MoBEM,FClist,cmp,"PC2",ssdgenes)

save(AnovaResCTt,file=paste0(dir.Results,"/AnovaResCTt.Rdata"))
save(AnovaResQNt,file=paste0(dir.Results,"/AnovaResQNt.Rdata"))
save(AnovaResPC1t,file=paste0(dir.Results,"/AnovaResPC1t.Rdata"))
save(AnovaResPC2t,file=paste0(dir.Results,"/AnovaResPC2t.Rdata"))
load(file=paste0(dir.Results,"/AnovaResCTt.Rdata"))
load(file=paste0(dir.Results,"/AnovaResQNt.Rdata"))
load(file=paste0(dir.Results,"/AnovaResPC1t.Rdata"))
load(file=paste0(dir.Results,"/AnovaResPC2t.Rdata"))



AnovaResCTtr<-lapply(AnovaResCTt,function(x) Reduce(rbind,x))
AnovaResQNtr<-lapply(AnovaResQNt,function(x) Reduce(rbind,x))
AnovaResPC1tr<-lapply(AnovaResPC1t,function(x) Reduce(rbind,x))
AnovaResPC2tr<-lapply(AnovaResPC2t,function(x) Reduce(rbind,x))

AnovaProp(AnovaResCTtr,"CT",dir.Results,fdr=0.05)
AnovaProp(AnovaResQNtr,"QN",dir.Results,fdr=0.05)
AnovaProp(AnovaResPC1tr,"PC1",dir.Results,fdr=0.05)
AnovaProp(AnovaResPC2tr,"PC2",dir.Results,fdr=0.05)

write.table(AnovaResCTtr[[1]][AnovaResCTtr[[1]]$fdr<0.05,],file=paste0(dir.Results,"/AnovaResultsCT_tissue_CCR.tsv"),sep="\t",quote=F)
write.table(AnovaResCTtr[[2]][AnovaResCTtr[[2]]$fdr<0.05,],file=paste0(dir.Results,"/AnovaResultsCT_tissue_CERES.tsv"),sep="\t",quote=F)
write.table(AnovaResCTtr[[3]][AnovaResCTtr[[3]]$fdr<0.05,],file=paste0(dir.Results,"/AnovaResultsCT_tissue_CCRJ.tsv"),sep="\t",quote=F)
write.table(AnovaResQNtr[[1]][AnovaResQNtr[[1]]$fdr<0.05,],file=paste0(dir.Results,"/AnovaResultsQN_tissue_CCR.tsv"),sep="\t",quote=F)
write.table(AnovaResQNtr[[2]][AnovaResQNtr[[2]]$fdr<0.05,],file=paste0(dir.Results,"/AnovaResultsQN_tissue_CERES.tsv"),sep="\t",quote=F)
write.table(AnovaResQNtr[[3]][AnovaResQNtr[[3]]$fdr<0.05,],file=paste0(dir.Results,"/AnovaResultsQN_tissue_CCRJ.tsv"),sep="\t",quote=F)
write.table(AnovaResPC1tr[[1]][AnovaResPC1tr[[1]]$fdr<0.05,],file=paste0(dir.Results,"/AnovaResultsPC1_tissue_CCR.tsv"),sep="\t",quote=F)
write.table(AnovaResPC1tr[[2]][AnovaResPC1tr[[2]]$fdr<0.05,],file=paste0(dir.Results,"/AnovaResultsPC1_tissue_CERES.tsv"),sep="\t",quote=F)
write.table(AnovaResPC1tr[[3]][AnovaResPC1tr[[3]]$fdr<0.05,],file=paste0(dir.Results,"/AnovaResultsPC1_tissue_CCRJ.tsv"),sep="\t",quote=F)
write.table(AnovaResPC2tr[[1]][AnovaResPC2tr[[1]]$fdr<0.05,],file=paste0(dir.Results,"/AnovaResultsPC2_tissue_CCR.tsv"),sep="\t",quote=F)
write.table(AnovaResPC2tr[[2]][AnovaResPC2tr[[2]]$fdr<0.05,],file=paste0(dir.Results,"/AnovaResultsPC2_tissue_CERES.tsv"),sep="\t",quote=F)
write.table(AnovaResPC2tr[[3]][AnovaResPC2tr[[3]]$fdr<0.05,],file=paste0(dir.Results,"/AnovaResultsPC2_tissue_CCRJ.tsv"),sep="\t",quote=F)



ResListT<-list(ComBat=AnovaResCTtr,ComBatQN=AnovaResQNtr,ComBatPC1=AnovaResPC1tr,ComBatPC2=AnovaResPC2tr)
plotBMres(ResListT,curveColours,fdr=0.05,plotName="tissue")
cfelistT<-plotBMresByCFE(ResListT,curveColours,fdr=0.05)
AnovaFDRcurve(ResListT,plotcolours = PipelineColours,plotName="tissue")
AnovaFDRcurve(ResListT,plotcolours = PipelineColours,noCFE=TRUE,plotName="tissue")




