library(here)
library(CRISPRcleanR)
library(ggplot2)
library(beeswarm)
library(qusage)

dir.MergeFile<-"./Results"
dir.Results<-"./ResultsFilter"
dir.Input<-"/Location/Of/Figshare/Directory"

source("Combat_HKfunctions.R")
load(paste0(dir.Input,"/InstituteColours.Rdata"))
load(paste0(dir.Input,"/PipelineColours.Rdata"))
load(paste0(dir.Input,"/BAGEL_essential"))
load(paste0(dir.Input,"/BAGEL_nonEssential"))

cmp<-read.csv(paste0(dir.Input,"/model_list_latest.csv"),header=T,stringsAsFactors = F)
cmp2<-cmp
cmp2$model_id<-cmp2$BROAD_ID
cmp<-rbind(cmp,cmp2)

load(file=paste0(dir.MergeFile,"/AllSymbolsHNGCmap.Rdata"))

#start with the SQ+ComBat+PC1 merged all data

CCR_corrected<-readRDS(file=paste0(dir.MergeFile,"/CCR_SQ_Combat_PC1_All_merge_F.Rds"))
CERES_corrected<-readRDS(file=paste0(dir.MergeFile,"/CERES_SQ_Combat_PC1_All_merge_F.Rds"))

###Load the binary depletion matrices
###Note make sure the items in the Figshare folder Binary Depletion matrices are available from dir.Input

load(paste0(dir.Input,"/BinaryCCRSanger.RData"))
load(paste0(dir.Input,"/BinaryCCRBroad.RData"))
load(paste0(dir.Input,"/BinaryCERESSanger.RData"))
load(paste0(dir.Input,"/BinaryCERESBroad.RData"))

###binary depletion for merged data PC1, see CalculateBinaryDepletions.R code for source:

load(paste0(dir.Results,"/BinaryCCR3m.RData"))
load(paste0(dir.Results,"/BinaryCERES3m.RData"))


###number depletions:
depCERESsanger<-rowSums(BinaryCERESSanger,na.rm=T)
depCERESbroad<-rowSums(BinaryCERESBroad,na.rm=T)

depCERESint<-rowSums(BinaryCERES3,na.rm=T)

onedep<-c(sum(depCERESsanger==1,na.rm=T),sum(depCERESbroad==1,na.rm=T),sum(depCERESint==1,na.rm=T))
threedep<-c(sum(depCERESsanger==3,na.rm=T),sum(depCERESbroad==3,na.rm=T),sum(depCERESint==3,na.rm=T))
fivedep<-c(sum(depCERESsanger==5,na.rm=T),sum(depCERESbroad==5,na.rm=T),sum(depCERESint==5,na.rm=T))
gtendep<-c(sum(depCERESsanger>=10,na.rm=T),sum(depCERESbroad>=10,na.rm=T),sum(depCERESint>=10,na.rm=T))

depDataCeres<-rbind(data.frame(deps=onedep,numDep="1",data=c("Sanger","Broad","Integrated")),
               data.frame(deps=threedep,numDep="3",data=c("Sanger","Broad","Integrated")),
               data.frame(deps=fivedep,numDep="5",data=c("Sanger","Broad","Integrated")),
               data.frame(deps=gtendep,numDep=">=10",data=c("Sanger","Broad","Integrated")))
names(InstituteColours)[3]<-"Integrated"
depDataCeres$data<-factor(depDataCeres$data,levels=c("Integrated","Broad","Sanger"))
depplot<-ggplot(depDataCeres,aes(x=numDep,y=deps,fill=data))+geom_bar(stat="identity",position="dodge")+scale_fill_manual(values=InstituteColours)+theme(axis.text.x=element_text(angle=45))+ylab("Number Genes")+xlab("Number depletions")+theme_bw()+theme(legend.position = c(0.2,0.8),legend.background = element_blank())

pdf(paste0(dir.Results,"/BinaryComparison_CERES_IntegratedFig5_PC1.pdf"))
print(depplot)
dev.off()



###use Broad 90th depletion method ###
load(paste0(dir.MergeFile,"/BroadDataCERES.Rdata"))
load(paste0(dir.MergeFile,"/SangerDataCERES.Rdata"))
load(paste0(dir.MergeFile,"/BroadData.Rdata"))
load(paste0(dir.MergeFile,"/SangerData.Rdata"))


library(ADaM2)


##CERES Broad method###
load(paste0(dir.Results,"/CE_IntCERESPC190.Rdata"))
CeresCE<-B90CERESPC1
CE_IntCeres90g<-CeresCE



###CERES ADaM method ###
load(paste0(dir.Results,"/AdamCoreCERESPC1.RData"))
#PCcoreCERESPC2
PCcoreCERES<-PCcoreCERESPC1


allMethods<-intersect(CeresCE,PCcoreCERES)

###Figure 5b###
##venn diagram of genes found across four different methods
library(eulerr)
A<-setdiff(CeresCE,allMethods)
B<-setdiff(PCcoreCERES,allMethods)

venn<-euler(c("A&B"=length(allMethods),A=length(A),B=length(B)))


###use the rank of the normLRT values to compare across datasets of different numbers of cell lines.
PCnumber=1
load(file=paste0(dir.Results,"/normLRTCCRmerge",PCnumber,".RData"))
load(file=paste0(dir.Results,"/normLRTCERESmerge",PCnumber,".RData"))
load(file=paste0(dir.Results,"/normLRTCCRJmerge",PCnumber,".Rdata"))


load(file=paste0(dir.Input,"/Kegg.DNArep.Rdata"))
load(file=paste0(dir.Input,"/Kegg.Ribosome.Rdata"))
load(file=paste0(dir.Input,"/Kegg.Proteasome.Rdata"))
load(file=paste0(dir.Input,"/Kegg.Spliceosome.Rdata"))
load(file=paste0(dir.Input,"/Kegg.RNApoly.Rdata"))
load(file=paste0(dir.Input,"/Histones.Rdata"))

allRefEss<-unique(c(Kegg.DNArep,Kegg.Proteasome,Kegg.Ribosome,Kegg.RNApoly,Kegg.Spliceosome,Histones))
allRefEss<-unique(allSymbol[na.omit(match(allRefEss,allSymbol$Input)),"Approved.symbol"])
getAllMap<-function(allSymbol,allgenes){
  symbolS1<-allSymbol[allSymbol$Match.type=="Approved symbol",]
  symbolS2<-allSymbol[allSymbol$Match.type=="Previous symbol",]
  symbolS3<-allSymbol[allSymbol$Match.type=="Synonyms",]
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
  return(allMap)
}

                       
#PanCancerCoreFitnessGenes, Behan 2019:
load(paste0(dir.Input,'/10_PANCANCER_coreFitness_genes.RData'))
Behan<-PanCancerCoreFitnessGenes
#load in the latest Traver Hart pan cancer genes:
load(paste0(dir.Input,'/BAGEL_v2_ESSENTIAL_GENES.Rdata'))


#CERES common essentials:
allCEsC<-unique(c(CeresCE,PCcoreCERES))

CEmatrixC<-matrix(0,nrow=length(allCEsC),ncol=2)
rownames(CEmatrixC)<-allCEsC
colnames(CEmatrixC)<-c("Ceres_CE","Ceres_ADaM")

CEmatrixC[CeresCE,1]<-1

CEmatrixC[PCcoreCERES,2]<-1

save(CEmatrixC,file=paste0(dir.Results,"/CEmatrixCERES.Rdata"))
NumberMethods<-rowSums(CEmatrixC)
Tier1<-sum(NumberMethods==2)
Tier2<-sum(NumberMethods==1)

Tier1g<-allCEsC[NumberMethods==2]
Tier2g<-allCEsC[NumberMethods==1]

TierCEC<-data.frame(gene=c(Tier1g,Tier2g),tier=c(rep("Tier1",Tier1),rep("Tier2",Tier2)),stringsAsFactors = F)
save(TierCEC,file=paste0(dir.Results,"/TierCE_CERES_PC1.Rdata"))
write.csv(TierCEC,file=paste0(dir.Results,"/TierCE_CERES_PC1.csv"),quote=FALSE)

load(file=paste0(dir.Results,"/TierCE_CERES_PC1.Rdata"))
Tier1<-as.character(TierCEC[which(TierCEC$tier=="Tier1"),1])
#Tier1<-as.character(TierCE[,1])
Tier2<-as.character(TierCEC[which(TierCEC$tier=="Tier1"|TierCEC$tier=="Tier2"),1])

allMap<-getAllMap(allSymbol = allSymbol,allgenes)
BehanM<-as.matrix(Behan,ncol=1)
rownames(BehanM)<-Behan
Behan<-names(updateRownames(BehanM,allMap))
Bagel<-as.matrix(BAGEL_essential,ncol=1)
rownames(Bagel)<-BAGEL_essential
BAGELessential<-names(updateRownames(Bagel,allMap))


RecallSpliceosome<-lapply(list(Behan=Behan,Hart=BAGELessential,Broad=BroadCE,Sanger=SangerCE,Int=Tier2),function(x) length(intersect(x,Kegg.Spliceosome))/length(Kegg.Spliceosome))
RecallRibosome<-lapply(list(Behan=Behan,Hart=BAGELessential,Broad=BroadCE,Sanger=SangerCE,Int=Tier2),function(x) length(intersect(x,Kegg.Ribosome))/length(Kegg.Ribosome))
RecallProteasome<-lapply(list(Behan=Behan,Hart=BAGELessential,Broad=BroadCE,Sanger=SangerCE,Int=Tier2),function(x) length(intersect(x,Kegg.Proteasome))/length(Kegg.Proteasome))
RecallRNApoly<-lapply(list(Behan=Behan,Hart=BAGELessential,Broad=BroadCE,Sanger=SangerCE,Int=Tier2),function(x) length(intersect(x,Kegg.RNApoly))/length(Kegg.RNApoly))
RecallDNARep<-lapply(list(Behan=Behan,Hart=BAGELessential,Broad=BroadCE,Sanger=SangerCE,Int=Tier2),function(x) length(intersect(x,Kegg.DNArep))/length(Kegg.DNArep))
RecallHistones<-lapply(list(Behan=Behan,Hart=BAGELessential,Broad=BroadCE,Sanger=SangerCE,Int=Tier2),function(x) length(intersect(x,Histones))/length(Histones))

RecallDataC<-rbind(data.frame(deps=unlist(RecallSpliceosome),set="Spliceosome",data=unlist(names(RecallSpliceosome))),
                   data.frame(deps=unlist(RecallRibosome),set="Ribosome",data=unlist(names(RecallRibosome))),
                   data.frame(deps=unlist(RecallProteasome),set="Proteasome",data=unlist(names(RecallProteasome))),
                   data.frame(deps=unlist(RecallRNApoly),set="RNA polymerase",data=unlist(names(RecallRNApoly))),
                   data.frame(deps=unlist(RecallDNARep),set="DNA Replication",data=unlist(names(RecallDNARep))),
                   data.frame(deps=unlist(RecallHistones),set="Histones",data=unlist(names(RecallHistones)))
                   )

RecallDataC$data<-factor(RecallDataC$data,levels=c("Int","Broad","Sanger","Behan","Hart"))
Pcolours<-c("blue","orange","lightblue","pink","purple")
names(Pcolours)<-c("Int","Broad","Sanger","Behan","Hart")
RecallDataC$set<-factor(RecallDataC$set,levels=c("Proteasome","Ribosome","RNA polymerase","DNA Replication","Spliceosome","Histones"))
Recallplot<-ggplot(RecallDataC,aes(x=set,y=deps,fill=data))+geom_bar(stat="identity",position="dodge")+scale_fill_manual(values=Pcolours)+theme(axis.text.x=element_text(angle=45))+ylab("Recall")+xlab("")+theme_bw()+theme(legend.position = c(0.25,0.925),legend.background = element_blank(),legend.title=element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                                                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

pdf(paste0(dir.Results,"/RecallKnownSets_IntegratedFig5_ConsensusVScurrent_CERES_PC1_All.pdf"))
print(Recallplot)
dev.off()

                       #load gene expression data to get negative controls
load(paste0(dir.Input,"/CCLEexpression.RData"))
allMap<-getAllMap(allSymbol,rownames(CCLEexpression))
CCLEexpression<-updateRownames(CCLEexpression,allMap)  
NotExpr<-ADAM2.PercentileCF(1/CCLEexpression,display=FALSE)$cfgenes
NotExpr<-as.matrix(NotExpr,ncol=1)
rownames(NotExpr)<-NotExpr
allMap<-getAllMap(allSymbol,NotExpr[,1])
NotExpr<-names(updateRownames(NotExpr,allMap))
NotExpr<-setdiff(NotExpr,allRefEss)

FDR_NotExpr<-c(
  length(intersect(NotExpr,Tier2))/length(Tier2),
  length(intersect(NotExpr,Behan))/length(Behan),
  length(intersect(NotExpr,BAGELessential))/length(BAGELessential),
  length(intersect(NotExpr,BroadCE))/length(BroadCE),
  length(intersect(NotExpr,SangerCE))/length(SangerCE))
FDRdata<-data.frame(FDR=FDR_NotExpr,Set=c("Int","Behan","Hart","Broad","Sanger"))
FDRdata$Set<-factor(FDRdata$Set,levels=c("Hart","Int","Broad","Behan","Sanger"))
Pcolours<-c("red","blue","orange","lightblue","pink","purple")
names(Pcolours)<-c("Tier1","Int","Broad","Sanger","Behan","Hart")
FDRplot<-ggplot(FDRdata,aes(x=Set,y=FDR,fill=Set))+geom_bar(stat="identity",position="dodge")+scale_fill_manual(values=Pcolours)+theme(axis.text.x=element_text(angle=45))+ylab("Estimated FDR")+xlab("")+theme_bw()+theme(legend.position = c(0.25,0.925),legend.background = element_blank(),legend.title=element_blank())


pdf(paste0(dir.Results,"/FDR_IntegratedFig5_CERES_PC1_notExpr.pdf"))
print(FDRplot)
dev.off()

                       
                       
###Comparison of biomarker analysis###                       

load(paste0(dir.Input,'/MoBEM.RData'))



CCR_sanger<-SangerData
CCR_broad<-BroadData
dn<-dimnames(CCR_sanger)
CCR_sanger<-normalize.quantiles(CCR_sanger)
dimnames(CCR_sanger)<-dn

dn<-dimnames(CCR_broad)
CCR_broad<-normalize.quantiles(CCR_broad)
dimnames(CCR_broad)<-dn
#load the relevant Binary matrices and show how many more genes are depleted at least once in integrated versus
#individual data sets


## creating MSS_status dummy variables and attaching them to the MoBEM
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
##individual not batch corrected data
colnames(CCR_sanger)<-cmp[match(colnames(CCR_sanger),cmp$model_id),"COSMIC_ID"]
colnames(CCR_broad)<-cmp[match(colnames(CCR_broad),cmp$BROAD_ID),"COSMIC_ID"]

CCR_allCosmic<-CCR_corrected

colnames(CCR_allCosmic)<-cmp[match(colnames(CCR_allCosmic),cmp$model_id),"COSMIC_ID"]


#filter by highest normLRT genes

n1<-names(normLRTCCR1[[3]])[normLRTCCR1[[3]]>200]

usegenessd<-n1
utissues<-unique(cmp$tissue)
##for the integrated data##:
ccr_anova<-list()


j=1
for(i in 1:length(utissues)){
  cosmiccl<-unique(cmp[cmp$tissue==utissues[i],"COSMIC_ID"])
  cosmicavail<-unique(colnames(CCR_allCosmic)[colnames(CCR_allCosmic)%in%cosmiccl])

  if(length(cosmicavail)>9){
    ccr_anova[[j]]<-Allanova(MoBEM,CCR_allCosmic[usegenessd,colnames(CCR_allCosmic)%in%cosmiccl],utissues[i])
  j=j+1
  }
}

cclnames<-c()
for(i in 1:length(utissues)){
  cosmiccl<-unique(cmp[cmp$tissue==utissues[i],"COSMIC_ID"])
  cosmicavail<-unique(colnames(CCR_allCosmic)[colnames(CCR_allCosmic)%in%cosmiccl])

  if(length(cosmicavail)>9){
    cclnames<-c(cclnames,utissues[i])
  }
}
names(ccr_anova)<-cclnames
sig_listCCR<-list()

match_listCCR<-list()

for(i in 1:length(ccr_anova)){
  
  temp<-ccr_anova[[i]]$p
  if(length(temp)>0){
    new<-ccr_anova[[i]]
    new$fdr<-p.adjust(temp,method="fdr")
    ccr_anova[[i]]<-new
    sigmat<-new[new$fdr<0.01,]
    sig_listCCR[[i]]<-sigmat
    mut_list<-sigmat[grep("_mut",sigmat$CFE),]
    
    mut_list$CFE<-unlist(sapply(mut_list$CFE,function(x) strsplit(x,"_mut",fixed=TRUE)[[1]][1]))
    
    match_listCCR[[i]]<-mut_list[which(mut_list$CFE==mut_list$GENE),]
  }else{
    sig_listCCR[[i]]<-NULL
  }
}


##for the single data sets:
ccr_sanger<-list()

ccr_broad<-list()

j=1
for(i in 1:length(utissues)){
  cosmiccl<-unique(cmp[cmp$tissue==utissues[i],"COSMIC_ID"])
  cosmicavail<-unique(colnames(CCR_sanger)[colnames(CCR_sanger)%in%cosmiccl])

  if(length(cosmicavail)>9){
    ccr_sanger[[j]]<-Allanova(MoBEM,CCR_sanger[usegenessd,colnames(CCR_sanger)%in%cosmiccl],utissues[i])

    j=j+1
  }
}
cclnamesS<-c()
for(i in 1:length(utissues)){
  cosmiccl<-unique(cmp[cmp$tissue==utissues[i],"COSMIC_ID"])
  cosmicavail<-unique(colnames(CCR_sanger)[colnames(CCR_sanger)%in%cosmiccl])

  if(length(cosmicavail)>9){
    cclnamesS<-c(cclnamesS,utissues[i])
  }
}
names(ccr_sanger)<-cclnamesS
j=1
for(i in 1:length(utissues)){
  cosmiccl<-unique(cmp[cmp$tissue==utissues[i],"COSMIC_ID"])
  cosmicavail<-unique(colnames(CCR_broad)[colnames(CCR_broad)%in%cosmiccl])

  if(length(cosmicavail)>9){

    ccr_broad[[j]]<-Allanova(MoBEM,CCR_broad[usegenessd,colnames(CCR_broad)%in%cosmiccl],utissues[i])
   j=j+1
  }
}

cclnamesB<-c()
for(i in 1:length(utissues)){
  cosmiccl<-unique(cmp[cmp$tissue==utissues[i],"COSMIC_ID"])
  cosmicavail<-unique(colnames(CCR_broad)[colnames(CCR_broad)%in%cosmiccl])

  if(length(cosmicavail)>9){
    cclnamesB<-c(cclnamesB,utissues[i])
  }
}
names(ccr_broad)<-cclnamesB
for(i in 1:length(ccr_sanger)){
  
  temp<-ccr_sanger[[i]]$p
  if(length(temp)>0){
    new<-ccr_sanger[[i]]
    new$fdr<-p.adjust(temp,method="fdr")
    ccr_sanger[[i]]<-new

  }
}
for(i in 1:length(ccr_broad)){
  
  temp<-ccr_broad[[i]]$p
  if(length(temp)>0){
    new<-ccr_broad[[i]]
    new$fdr<-p.adjust(temp,method="fdr")
    ccr_broad[[i]]<-new
    
  }
}
###set test names:
ccrN<-lapply(ccr_anova,function(x) {rownames(x)<-paste(x$CFE,paste(x$GENE,x$tissue,sep="-"),sep="-")
return(x)})

ccrS<-lapply(ccr_sanger,function(x) {rownames(x)<-paste(x$CFE,paste(x$GENE,x$tissue,sep="-"),sep="-")
return(x)})

ccrB<-lapply(ccr_broad,function(x) {rownames(x)<-paste(x$CFE,paste(x$GENE,x$tissue,sep="-"),sep="-")
return(x)})

names(ccrN)<-names(ccr_anova)
names(ccrS)<-names(ccr_sanger)
names(ccrB)<-names(ccr_broad)

save(ccrS,file=paste0(dir.Results,"/ccrSangerBiomarkersPC1.Rdata"))
save(ccrB,file=paste0(dir.Results,"/ccrBroadBiomarkersPC1.Rdata"))
save(ccrN,file=paste0(dir.Results,"/ccrIntegratedBiomarkersPC1.Rdata"))
load(file=paste0(dir.Results,"/ccrSangerBiomarkersPC1.Rdata"))
load(file=paste0(dir.Results,"/ccrBroadBiomarkersPC1.Rdata"))
load(file=paste0(dir.Results,"/ccrIntegratedBiomarkersPC1.Rdata"))
NewSanger<-FindSignif(ccrS,ccrN,fdr=0.01)
NewBroad<-FindSignif(ccrB,ccrN,fdr=0.01)

SangerNew<-Reduce(rbind,NewSanger$criteria)
BroadNew<-Reduce(rbind,NewBroad$criteria)
Newboth<-intersect(rownames(SangerNew),rownames(BroadNew))

cat(paste("Number new biomarker associations CCR Integrated versus Sanger data PC1:",nrow(SangerNew)),file=paste0(dir.Results,"/Figure5Log.txt"),sep="\n",append=TRUE)
cat(paste("Number new biomarker associations CCR Integrated versus Broad data PC1:",nrow(BroadNew)),file=paste0(dir.Results,"/Figure5Log.txt"),sep="\n",append=TRUE)


save(SangerNew,file=paste0(dir.Results,"/SangerNewPC1.Rdata"))
save(BroadNew,file=paste0(dir.Results,"/BroadNewPC1.Rdata"))

write.csv(SangerNew,file=paste0(dir.Results,"/SangerNewBiomarkersPC1.csv"))
write.csv(BroadNew,file=paste0(dir.Results,"/BroadNewBiomarkersPC1.csv"))
dim(SangerNew)
dim(BroadNew)
length(Newboth)

#Sanger TP53mut MDM2 Lung
#Broad STAG2_mut-STAG1-Central Nervous System
pdf(paste0(dir.Results,"/ExampleNewBiomarkersPC1.pdf"),width=6.5,height=2.25)
par(mfrow=c(1,4))
par(mar=c(0.3,2,0.5,0.1)+0.1,mgp=c(1.25,0.5,0))
plotAssociationExamplesCP(CCR_sanger,CCR_allCosmic,20,MoBEM,cmp,SangerNew)

plotAssociationExamplesCP(CCR_broad,CCR_allCosmic,60,MoBEM,cmp,BroadNew)
dev.off()

ccrS[["Lung"]]["TP53_mut-MDM2-Lung",]
ccrB[["Central Nervous System"]]["STAG2_mut-STAG1-Central Nervous System",]
pdf(paste0(dir.Results,"/Figure6BioMarkersPC1.pdf"),width=6.5,height=2.25)
par(mfrow=c(1,2))
par(mar=c(0.3,2,0.5,0.1)+0.1,mgp=c(1.25,0.5,0))
#plot 1 CDKN2C SMAD7 Haem
#plot 2 PPM1D FOXA1 Breast
plotSingleAssociation(CCR_allCosmic,10,MoBEM,cmp,SangerNew)
plotSingleAssociation(CCR_allCosmic,38,MoBEM,cmp,SangerNew)
dev.off()

