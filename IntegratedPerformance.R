library(here)
library(CRISPRcleanR)
library(ggplot2)
library(beeswarm)


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



###Figure 5b - recall of essential genes###
CEbroad<-read.csv(paste0(dir.Input,"/priorcommonessentials.csv"),header=T,stringsAsFactors = F)
CEbroad<-apply(CEbroad,1,function(x) strsplit(x," ",fixed=TRUE)[[1]][1])
CEbroad<-intersect(rownames(BinaryCCR3),CEbroad)
NCEbroad<-read.csv(paste0(dir.Input,"/priornonessentials.csv"),header=T,stringsAsFactors = F)
NCEbroad<-apply(NCEbroad,1,function(x) strsplit(x," ",fixed=TRUE)[[1]][1])
NCEbroad<-intersect(rownames(BinaryCCR3),NCEbroad)

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


load(paste0(dir.Input,'/EssGenes.DNA_REPLICATION_cons.RData'))
load(paste0(dir.Input,'EssGenes.KEGG_rna_polymerase.RData'))
load(paste0(dir.Input,'EssGenes.PROTEASOME_cons.RData'))
load(paste0(dir.Input,'EssGenes.ribosomalProteins.RData'))
load(paste0(dir.Input,'EssGenes.SPLICEOSOME_cons.RData'))

#PanCancerCoreFitnessGenes, Behan 2019:
load(paste0(dir.Input,'/10_PANCANCER_coreFitness_genes.RData'))
Behan<-PanCancerCoreFitnessGenes
#load in the latest Traver Hart pan cancer genes:
load(paste0(dir.Input,'/BAGEL_v2_ESSENTIAL_GENES.Rdata'))


Recall_BehanD9<-length(intersect(EssGenes.DNA_REPLICATION_cons,Behan))/length(EssGenes.DNA_REPLICATION_cons)
Recall_HartD9<-length(intersect(EssGenes.DNA_REPLICATION_cons,BAGEL_essential))/length(EssGenes.DNA_REPLICATION_cons)


Recall_BehanK9<-length(intersect(EssGenes.KEGG_rna_polymerase,Behan))/length(EssGenes.KEGG_rna_polymerase)
Recall_HartK9<-length(intersect(EssGenes.KEGG_rna_polymerase,BAGEL_essential))/length(EssGenes.KEGG_rna_polymerase)


Recall_BehanP9<-length(intersect(EssGenes.PROTEASOME_cons,Behan))/length(EssGenes.PROTEASOME_cons)
Recall_HartP9<-length(intersect(EssGenes.PROTEASOME_cons,BAGEL_essential))/length(EssGenes.PROTEASOME_cons)

Recall_BehanR9<-length(intersect(EssGenes.ribosomalProteins,Behan))/length(EssGenes.ribosomalProteins)
Recall_HartR9<-length(intersect(EssGenes.ribosomalProteins,BAGEL_essential))/length(EssGenes.ribosomalProteins)

Recall_BehanS9<-length(intersect(EssGenes.SPLICEOSOME_cons,Behan))/length(EssGenes.SPLICEOSOME_cons)
Recall_HartS9<-length(intersect(EssGenes.SPLICEOSOME_cons,BAGEL_essential))/length(EssGenes.SPLICEOSOME_cons)

#option 2 for common essentials, just use CERES:
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

Recall_IntCeresD9c<-length(intersect(EssGenes.DNA_REPLICATION_cons,Tier1))/length(EssGenes.DNA_REPLICATION_cons)
Recall_IntCeresD9t2c<-length(intersect(EssGenes.DNA_REPLICATION_cons,Tier2))/length(EssGenes.DNA_REPLICATION_cons)


Recall_IntCeresDc<-length(intersect(EssGenes.DNA_REPLICATION_cons,Tier1))/length(EssGenes.DNA_REPLICATION_cons)
Recall_IntCeresDt2c<-length(intersect(EssGenes.DNA_REPLICATION_cons,Tier2))/length(EssGenes.DNA_REPLICATION_cons)


Recall_IntCeresK9c<-length(intersect(EssGenes.KEGG_rna_polymerase,Tier1))/length(EssGenes.KEGG_rna_polymerase)
Recall_IntCeresK9t2c<-length(intersect(EssGenes.KEGG_rna_polymerase,Tier2))/length(EssGenes.KEGG_rna_polymerase)

Recall_IntCeresP9c<-length(intersect(EssGenes.PROTEASOME_cons,Tier1))/length(EssGenes.PROTEASOME_cons)
Recall_IntCeresP9t2c<-length(intersect(EssGenes.PROTEASOME_cons,Tier2))/length(EssGenes.PROTEASOME_cons)

Recall_IntCeresR9c<-length(intersect(EssGenes.ribosomalProteins,Tier1))/length(EssGenes.ribosomalProteins)
Recall_IntCeresR9t2c<-length(intersect(EssGenes.ribosomalProteins,Tier2))/length(EssGenes.ribosomalProteins)


Recall_IntCeresS9c<-length(intersect(EssGenes.SPLICEOSOME_cons,Tier1))/length(EssGenes.SPLICEOSOME_cons)
Recall_IntCeresS9t2c<-length(intersect(EssGenes.SPLICEOSOME_cons,Tier2))/length(EssGenes.SPLICEOSOME_cons)

RecallDataC<-rbind(data.frame(deps=Recall_BehanS9,set="Spliceosome",data="Behan"),
                   data.frame(deps=Recall_HartS9,set="Spliceosome",data="Hart"),
                   data.frame(deps=Recall_IntCeresS9c,set="Spliceosome",data="Integrated Tier 1"),
                   data.frame(deps=Recall_IntCeresS9t2c,set="Spliceosome",data="Integrated Tier 2"),
                   
                   data.frame(deps=Recall_BehanR9,set="Ribosomal",data="Behan"),
                   data.frame(deps=Recall_HartR9,set="Ribosomal",data="Hart"),
                   data.frame(deps=Recall_IntCeresR9c,set="Ribosomal",data="Integrated Tier 1"),
                   data.frame(deps=Recall_IntCeresR9t2c,set="Ribosomal",data="Integrated Tier 2"),
                   
                   data.frame(deps=Recall_BehanP9,set="Proteasome",data="Behan"),
                   data.frame(deps=Recall_HartP9,set="Proteasome",data="Hart"),
                   data.frame(deps=Recall_IntCeresP9c,set="Proteasome",data="Integrated Tier 1"),
                   data.frame(deps=Recall_IntCeresP9t2c,set="Proteasome",data="Integrated Tier 2"),
                   
                   data.frame(deps=Recall_BehanK9,set="RNA_polymerase",data="Behan"),
                   data.frame(deps=Recall_HartK9,set="RNA_polymerase",data="Hart"),
                   data.frame(deps=Recall_IntCeresK9c,set="RNA_polymerase",data="Integrated Tier 1"),
                   data.frame(deps=Recall_IntCeresK9t2c,set="RNA_polymerase",data="Integrated Tier 2"),
                   
                   data.frame(deps=Recall_BehanD9,set="DNA_Replication",data="Behan"),
                   data.frame(deps=Recall_HartD9,set="DNA_Replication",data="Hart"),
                   data.frame(deps=Recall_IntCeresD9c,set="DNA_Replication",data="Integrated Tier 1"),
                   data.frame(deps=Recall_IntCeresD9t2c,set="DNA_Replication",data="Integrated Tier 2"))


RecallDataC$data<-factor(RecallDataC$data,levels=c("Integrated Tier 1","Integrated Tier 2","Behan","Hart"))
Pcolours<-c("red","blue","pink","purple")
names(Pcolours)<-c("Integrated Tier 1","Integrated Tier 2","Behan","Hart")
RecallDataC$set<-factor(RecallDataC$set,levels=c("Proteasome","Ribosomal","RNA_polymerase","DNA_Replication","Spliceosome"))
Recallplot<-ggplot(RecallDataC,aes(x=set,y=deps,fill=data))+geom_bar(stat="identity",position="dodge")+scale_fill_manual(values=Pcolours)+theme(axis.text.x=element_text(angle=45))+ylab("Number Genes")+xlab("")+theme_bw()+theme(legend.position = c(0.25,0.925),legend.background = element_blank(),legend.title=element_blank())

pdf(paste0(dir.Results,"/RecallKnownSets_IntegratedFig5_ConsensusVScurrent_CERES_PC1.pdf"))
print(Recallplot)
dev.off()


###load and run the normLRT recall of DEMETER genes to assess false positives.
ssc<-read.csv(file=paste0(dir.Input,"/skewed_tdist.csv"),header=T,stringsAsFactors = F)
SSC20<-ssc[ssc[,2]>=20,1]
SSC50<-ssc[ssc[,2]>=50,1]
SSC100<-ssc[ssc[,2]>=100,1]
SSC200<-ssc[ssc[,2]>=200,1]

FDR_SC20<-c(length(intersect(SSC20,Tier1))/length(Tier1),
            length(intersect(SSC20,Tier2))/length(Tier2),
            length(intersect(SSC20,Behan))/length(Behan),
            length(intersect(SSC20,BAGEL_essential))/length(BAGEL_essential))
FDR_SC50<-c(length(intersect(SSC50,Tier1))/length(Tier1),
            length(intersect(SSC50,Tier2))/length(Tier2),
            length(intersect(SSC50,Behan))/length(Behan),
            length(intersect(SSC50,BAGEL_essential))/length(BAGEL_essential))
FDR_SC100<-c(length(intersect(SSC100,Tier1))/length(Tier1),
             length(intersect(SSC100,Tier2))/length(Tier2),
             length(intersect(SSC100,Behan))/length(Behan),
             length(intersect(SSC100,BAGEL_essential))/length(BAGEL_essential))
FDR_SC200<-c(length(intersect(SSC200,Tier1))/length(Tier1),
             length(intersect(SSC200,Tier2))/length(Tier2),
             length(intersect(SSC200,Behan))/length(Behan),
             length(intersect(SSC200,BAGEL_essential))/length(BAGEL_essential))
FDRs<-rbind(FDR_SC20,FDR_SC50,FDR_SC100,FDR_SC200)
pdf(paste0(dir.Results,"/skewLRT_demeter_Tier_CERES+PC1.pdf"))
barplot(t(FDRs),beside = TRUE,col=Pcolours,border=FALSE,ylab='Estimated FDR',xlab = 'Log-likelihood of skew-t distribution',names.arg = c(20,50,100,200))
dev.off()

allPriorSets<-c(EssGenes.DNA_REPLICATION_cons,EssGenes.KEGG_rna_polymerase,EssGenes.PROTEASOME_cons,EssGenes.ribosomalProteins,EssGenes.SPLICEOSOME_cons)
allPriorSets<-unique(allPriorSets)

PrecisionT1<-length(intersect(Tier1,allPriorSets))/length(Tier1)
PrecisionT2<-length(intersect(Tier2,allPriorSets))/length(Tier2)
PrecisionBehan<-length(intersect(Behan,allPriorSets))/length(Behan)
PrecisionHart<-length(intersect(BAGEL_essential,allPriorSets))/length(BAGEL_essential)

PrecisionData<-data.frame(c("Tier1"=PrecisionT1,"Tier2"=PrecisionT2,"Behan"=PrecisionBehan,"Hart"=PrecisionHart))
PrecisionData$set<-rownames(PrecisionData)
colnames(PrecisionData)<-c("Precision","Set")
Pcolours<-c("red","blue","pink","purple")
names(Pcolours)<-c("Tier1","Tier2","Behan","Hart")
Precisionplot<-ggplot(PrecisionData,aes(x=Set,y=Precision,fill=Set))+geom_bar(stat="identity",position="dodge")+scale_fill_manual(values=Pcolours)+theme(axis.text.x=element_text(angle=45))+ylab("Precision")+xlab("")+theme_bw()+theme(legend.position = c(0.25,0.925),legend.background = element_blank(),legend.title=element_blank())

pdf(paste0(dir.Results,"/PrecisionKnownSets_IntegratedFig5_CERES_PC1.pdf"))
print(Precisionplot)
dev.off()

load(paste0(dir.Input,'/MoBEM.RData'))



CCR_sanger<-SangerData
CCR_broad<-BroadData

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

CCRint_Broad<-CCR_corrected[,grep("ACH",colnames(CCR_corrected))]
CCRint_Sanger<-CCR_corrected[,grep("SIDM",colnames(CCR_corrected))]
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
plotAssociationExamplesCP(CCR_sanger,CCR_allCosmic,16,MoBEM,cmp,SangerNew)

plotAssociationExamplesCP(CCR_broad,CCR_allCosmic,52,MoBEM,cmp,BroadNew)
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

