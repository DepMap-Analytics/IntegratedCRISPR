
library(here)
library(preprocessCore)

library(stringr)
dir.Results<-"./Results"
dir.Input<-"/Director/ToFigshare/"
source("./Combat_HKfunctions.R")
#create folders for output if they don't exist already:
if(!dir.exists("./Results")){dir.create("./Results")}
if(!dir.exists("./ResultsFilter")){dir.create("./ResultsFilter")}
###Load CRISPRCleanR processed only data ####
### load Broad data ####

#Broad Data CRISPRcleanR 20Q2:
geneFCmatrixv<-readRDS(file=paste0(dir.Input,"/Broad_CCR_20Q2.Rds"))
### load Sanger data. ####
S_ccr<-read.table(paste0(dir.Input,"/Sanger_corrected_logFCs.tsv"),sep="\t",stringsAsFactors = F,header=T,row.names=1)
S_ccr<-as.matrix(S_ccr)
### load meta info of overlapping cell lines ####
cmp<-read.csv(paste0(dir.Input,"/model_list_latest.csv"),header=T,stringsAsFactors = F)
cmp$mn<-make.names(cmp$model_name)
###load in the Broad meta data
BroadSI<-read.csv(paste0(dir.Input,"/sample_info.csv"),header=T,stringsAsFactors = F)

B20Q2<-read.csv(paste0(dir.Input,"/Achilles_replicate_map.csv"),header=T,stringsAsFactors = F)

#this is fine with cmp incomplete as cmp does include Broad_ids for all crispr'd cell lines
OverlapData<-cmp[cmp$BROAD_ID%in%B20Q2$DepMap_ID,]
OverlapData<-OverlapData[OverlapData$BROAD_ID!="",]
OverlapData<-OverlapData[OverlapData$crispr_ko_data=="True",]
OverlapData$bmn<-make.names(OverlapData$BROAD_ID)

OverlapData$Bmn<-paste0(OverlapData$mn,"---Broad")
OverlapData$Smn<-paste0(OverlapData$mn,"---Sanger")

###Generate All data objects for CRISPRcleanR only data ###

SangerMatrix<-S_ccr

###Read in CERES data ###
#this is the Sanger data processed with CERES
geneeffect<-read.csv(paste0(dir.Input,"Sanger_gene_effect_unscaled.csv"),header=T,stringsAsFactors = F,row.names=1)
geneeffect<-t(geneeffect)
geneeffect<-geneeffect[rowSums(is.na(geneeffect))==0,]
gnamesB<-sapply(rownames(geneeffect),function(x) strsplit(x,"..",fixed=TRUE)[[1]][1])
gidB<-sapply(rownames(geneeffect),function(x) strsplit(x,"..",fixed=TRUE)[[1]][1])
gnames<-sapply(rownames(geneeffect),function(x) strsplit(x,"..",fixed=TRUE)[[1]][1])
rownames(geneeffect)<-unlist(gnames)

dn<-dimnames(geneeffect)
qnorm_corrected_logFCs<-normalize.quantiles(geneeffect)
dimnames(qnorm_corrected_logFCs)<-dn

dn<-dimnames(S_ccr)
qnorm_corrected_logFCs_t<-normalize.quantiles(S_ccr)
dimnames(qnorm_corrected_logFCs_t)<-dn

### Broad data processed with CERES ####

geneeffectB<-read.csv(file=paste0(dir.Input,"/Achilles_gene_effect_unscaled.csv"),header=T,stringsAsFactors = F,row.names=1)
geneeffectB<-t(geneeffectB)
nacelllines<-which(colSums(is.na(geneeffectB))!=0)
minNAgenes<-min(colSums(is.na(geneeffectB[,nacelllines])))
rmcl<-which(colSums(is.na(geneeffectB))>minNAgenes)
geneeffectBall<-geneeffectB
totalRemove<-intersect(colnames(geneeffectB)[rmcl],OverlapData$BROAD_ID)
geneeffectBall<-geneeffectBall[,!colnames(geneeffectBall)%in%totalRemove]
geneeffectB<-geneeffectB[,-rmcl]

write.csv(rmcl,file=paste0(dir.Results,'Drop12NABroadCL.csv'))
OverlapDataAll<-OverlapData
OverlapData<-OverlapData[!OverlapData$BROAD_ID%in%names(rmcl),]

geneeffectB2<-geneeffectB[rowSums(is.na(geneeffectB))==0,]
droppedgenes<-setdiff(rownames(geneeffectB),rownames(geneeffectB2))
numberNAgenes<-rowSums(is.na(geneeffectB[droppedgenes,]))
CLnas<-which(colSums(is.na(geneeffectB[droppedgenes,]))!=0)
write.csv(droppedgenes,file=paste0(dir.Results,"DroppedGenes.csv"))
CLnas<-which(colSums(is.na(geneeffectB[droppedgenes,]))!=0)
write.csv(CLnas,file=paste0(dir.Results,"/CLwithNA.csv"))
geneeffectB<-geneeffectB[,-CLnas]
totalRemove<-intersect(colnames(geneeffectB)[CLnas],OverlapData$BROAD_ID)
geneeffectBall<-geneeffectBall[,!colnames(geneeffectBall)%in%totalRemove]
OverlapData<-OverlapData[!OverlapData$BROAD_ID%in%names(CLnas),]
gnames<-sapply(rownames(geneeffectB),function(x) strsplit(x,"..",fixed=TRUE)[[1]][1])
gids<-sapply(rownames(geneeffectB),function(x) strsplit(x,"..",fixed=TRUE)[[1]][2])
gids<-sapply(gids,function(x) gsub(".","",x,fixed=TRUE))
rownames(geneeffectB)<-unlist(gnames)
rownames(geneeffectBall)<-unlist(gnames)
#need to remove any that are also screened at the sanger and in the CLna list:

dn<-dimnames(geneeffectB)
qnorm_corrected_logFCs<-normalize.quantiles(geneeffectB)
dimnames(qnorm_corrected_logFCs)<-dn

OverlapCeres<-OverlapData[which(OverlapData$BROAD_ID%in%colnames(geneeffectB)),]

CeresSover<-geneeffect[,OverlapCeres$BROAD_ID]
CeresBover<-geneeffectB[,OverlapCeres$BROAD_ID]

### Get Broad data overlapping between two institutes ####
BroadOverlapBroadvC<-geneFCmatrixv[,OverlapData$BROAD_ID]
colnames(BroadOverlapBroadvC)<-OverlapData$mn

### Get Sanger data overlapping between two institutes ####

SangerOverlap<-SangerMatrix[,OverlapData$model_id]

colnames(SangerOverlap)<-paste(OverlapData$mn,'---Sanger',sep='')

colnames(BroadOverlapBroadvC)<-paste(colnames(BroadOverlapBroadvC),'---Broad',sep='')

###Update everything to new Hugo gene symbols consistently ###

#hierarchically rematch:
load(paste0(dir.Input,"/AllSymbolsHNGCmap.Rdata"))
symbolS1<-allSymbol[allSymbol$Match.type=="Approved symbol",]
symbolS2<-allSymbol[allSymbol$Match.type=="Previous symbol",]
symbolS3<-allSymbol[allSymbol$Match.type=="Synonyms",]


allgenes<-c(rownames(SangerMatrix),rownames(geneeffect),rownames(geneeffectB),rownames(geneFCmatrixv))
allgenes<-unique(allgenes)
allMap<-symbolS1[match(allgenes,symbolS1$Input),"Approved.symbol"]
names(allMap)<-allgenes
m1<-names(allMap)[is.na(allMap)]
m1M<-symbolS2[match(m1,symbolS2$Input),"Approved.symbol"]
names(m1M)<-m1
m2<-names(m1M)[is.na(m1M)]
m2M<-symbolS3[match(m2,symbolS3$Input),"Approved.symbol"]
names(m2M)<-m2

allMap[names(m1M)]<-m1M
allMap[names(m2M)]<-m2M


allMap[is.na(allMap)]<-names(allMap)[is.na(allMap)]

CeresBover<-updateRownames(CeresBover,allMap)
CeresSover<-updateRownames(CeresSover,allMap)
SangerMatrix<-updateRownames(SangerMatrix,allMap)
SangerOverlap<-updateRownames(SangerOverlap,allMap)
BroadOverlapBroadvC<-updateRownames(BroadOverlapBroadvC,allMap)
geneeffect<-updateRownames(geneeffect,allMap)
geneeffectB<-updateRownames(geneeffectB,allMap)
geneeffectBall<-updateRownames(geneeffectBall,allMap)
geneFCmatrixv<-updateRownames(geneFCmatrixv,allMap)


#read in JACKS default normalisation#
JACKSsangerD<-read.table(paste0(dir.Input,"/SangerJACKSd_gene_JACKS_results.txt"),header=T,sep="\t",stringsAsFactors = F,row.names=1)

JACKSbroadD<-read.table(paste0(dir.Input,"/BroadJACKS_gene_JACKS_results.txt"),header=T,sep="\t",row.names=1,stringsAsFactors = FALSE)

JACKSbroadD<-updateRownames(JACKSbroadD,allMap)
JACKSsangerD<-updateRownames(JACKSsangerD,allMap)




###Get all data for same cell lines and same genes ###
genesinc<-intersect(intersect(rownames(geneeffect),rownames(geneeffectB)),intersect(rownames(SangerMatrix),rownames(geneFCmatrixv)))
save(genesinc,file=paste0(dir.Results,"/GeneIncluded.Rdata"))
clinc<-intersect(which(OverlapData$Smn%in%colnames(SangerOverlap)),which(OverlapData$BROAD_ID%in%colnames(CeresBover)))
OverlapData<-OverlapData[clinc,]
save(OverlapData,file=paste0(dir.Results,"/OverlapData.Rdata"))

SangerOverlap<-SangerOverlap[genesinc,OverlapData[,"Smn"]]
BroadOverlapBroadvC<-BroadOverlapBroadvC[genesinc,OverlapData[,"Bmn"]]



CeresBover<-CeresBover[genesinc,OverlapData[,"BROAD_ID"]]
CeresSover<-CeresSover[genesinc,OverlapData[,"BROAD_ID"]]
colnames(CeresSover)<-paste0(colnames(CeresSover),"---Sanger")
colnames(CeresBover)<-paste0(colnames(CeresBover),"---Broad")

save(CeresBover,file=paste0(dir.Results,"/CeresBover.Rdata"))
save(CeresSover,file=paste0(dir.Results,"/CeresSover.Rdata"))

save(SangerOverlap,file=paste0(dir.Results,"/SangerOverlap.Rdata"))
save(BroadOverlapBroadvC,file=paste0(dir.Results,"/BroadOverlapBroadvC.Rdata"))


#convert to cell model passport ID so can distinguish Sanger screens from Broad screens and column names in integrated data set are unique

#colnames(SangerMatrix)<-cmp[match(colnames(SangerMatrix),cmp$mn),"model_id"]
colnames(geneeffect)<-cmp[match(colnames(geneeffect),cmp$BROAD_ID),"model_id"]
BroadData<-geneFCmatrixv[genesinc,]
IncSangerLines<-cmp[cmp$crispr_ko_data=="True","model_id"]
SangerData<-SangerMatrix[genesinc,colnames(SangerMatrix)%in%IncSangerLines]

SangerDataCERES<-geneeffect[genesinc,colnames(geneeffect)%in%IncSangerLines]
BroadDataCERES<-geneeffectBall[genesinc,]
clincB<-intersect(colnames(BroadDataCERES),colnames(BroadData))
clincS<-intersect(colnames(SangerDataCERES),colnames(SangerData))


BroadDataCERES<-BroadDataCERES[,clincB]
SangerDataCERES<-SangerDataCERES[,clincS]

save(BroadDataCERES,file=paste0(dir.Results,"/BroadDataCERES.Rdata"))
save(SangerDataCERES,file=paste0(dir.Results,"/SangerDataCERES.Rdata"))

BroadData<-BroadData[,clincB]
SangerData<-SangerData[,clincS]
save(BroadData,file=paste0(dir.Results,"/BroadData.Rdata"))
save(SangerData,file=paste0(dir.Results,"/SangerData.Rdata"))

JSoverlap<-as.matrix(JACKSsangerD[genesinc,OverlapData$model_id])
JBoverlap<-as.matrix(JACKSbroadD[genesinc,make.names(OverlapData$BROAD_ID)])

colnames(JSoverlap)<-paste0(OverlapData[,"mn"],"---Sanger")
colnames(JBoverlap)<-paste0(OverlapData[,"mn"],"---Broad")

save(JSoverlap,file=paste0(dir.Results,"/JACKSoverlapBroad.Rdata"))
save(JBoverlap,file=paste0(dir.Results,"/JACKSoverlapSanger.Rdata"))

JACKSsangerD<-JACKSsangerD[genesinc,clincS]
colnames(JACKSbroadD)<-gsub(".","-",colnames(JACKSbroadD),fixed=TRUE)
JACKSbroadD<-JACKSbroadD[genesinc,clincB]

save(JACKSsangerD,file=paste0(dir.Results,"/JACKSsanger.Rdata"))
save(JACKSbroadD,file=paste0(dir.Results,"/JACKSbroad.Rdata"))


