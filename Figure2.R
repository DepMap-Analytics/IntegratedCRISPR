
library(here)
source("Combat_HKfunctions.R")

dir.Results<-"./Results"
dir.Input<-"/path/to/downloaded/figshare/"

library(RColorBrewer)
library(colorspace)
library(scales)
library(psych)
library(sva)
library(ggfortify)
library(sunburstR)
library(VennDiagram)

#Note this annotation file adds filler CMPIDs to file to allow for lineage data to be read for all Broad samples
cmp<-read.csv(paste0(dir.Input,"/model_list_latest.csv"),header=T,stringsAsFactors = FALSE)


#latest 20Q2 sample information from Broad:
rep_map<-read.csv(paste0(dir.Input,"/Achilles_replicate_map.csv"),header=T,stringsAsFactors = F)
sample_info<-read.csv(paste0(dir.Input,"/sample_info.csv"),header=T,stringsAsFactors = F)


##20Q2 broad samples processed by CERES:
load(file=paste0(dir.Results,"/BroadDataCERES.Rdata"))
load(file=paste0(dir.Results,"/SangerData.Rdata"))


load(file=paste0(dir.Results,"/OverlapData.Rdata"))
load(file=paste0(dir.Results,"/SangerOverlap.Rdata"))

load(file=paste0(dir.Results,"/BroadOverlapBroadvC.Rdata"))


BroadSamples<-colnames(BroadDataCERES)
BCrisprcmp<-cmp[cmp$BROAD_ID%in%colnames(BroadDataCERES),]
BCrispr<-sample_info[sample_info$DepMap_ID%in%colnames(BroadDataCERES),]
MB<-setdiff(BroadSamples,cmp$BROAD_ID)

SCrisprcmp<-cmp[cmp$model_id%in%colnames(SangerData),]
SCrispr<-sample_info[sample_info$Sanger_Model_ID%in%colnames(SangerData),]


BCrispr$source<-"Broad"
SCrispr$source<-"Sanger"
BCrisprcmp$source<-"Broad"
SCrisprcmp$source<-"Sanger"
BCrispr$source[BCrispr$DepMap_ID%in%SCrispr$DepMap_ID]<-"Both"
SCrispr<-SCrispr[!SCrispr$DepMap_ID%in%BCrispr$DepMap_ID,]
AllCrispr<-rbind(BCrispr,SCrispr)
BCrisprcmp$source[BCrisprcmp$BROAD_ID%in%SCrisprcmp$BROAD_ID]<-"Both"
SCrisprcmp<-SCrisprcmp[!SCrisprcmp$BROAD_ID%in%BCrisprcmp$BROAD_ID,]
AllCrisprcmp<-rbind(BCrisprcmp,SCrisprcmp)

AllCrispr$SBtype<-paste(AllCrispr$lineage,AllCrispr$source,sep="-")
AllCrispr$mn<-AllCrispr$stripped_cell_line_name
OverlapData$Smn<-paste(OverlapData$model_name,"---Sanger")
OverlapData$Bmn<-paste(OverlapData$model_name,"---Broad")
AllCrispr$InBatchCorrection<-"No"
AllCrispr[AllCrispr$DepMap_ID%in%OverlapData$BROAD_ID,"InBatchCorrection"]<-"Yes"
AllCrispr[AllCrispr$Sanger_Model_ID%in%OverlapData$model_ID,"InBatchCorrection"]<-"Yes"

save(AllCrispr,file=paste0(dir.Results,"/AllCrisprInfo.Rdata"))
AllCrisprcmp$SBtype<-paste(AllCrisprcmp$tissue,AllCrisprcmp$source,sep="-")
AllCrisprcmp$mn<-make.names(AllCrisprcmp$model_name)
save(AllCrisprcmp,file=paste0(dir.Results,"/AllCrisprInfoCMP.Rdata"))
SBdata<-data.frame(cbind(names(table(AllCrispr$SBtype)),table(AllCrispr$SBtype)))
AllCrispr$Sanger_Model_ID<-cmp[match(AllCrispr$DepMap_ID,cmp$BROAD_ID),"model_id"]
write.table(AllCrispr[,c(1,21)],paste0(dir.Results,"/AllCrispr.txt"),quote=FALSE,sep="\t")


load(paste0(dir.Input,"/PipelineColours.Rdata"))
load(paste0(dir.Input,"/InstituteColours.Rdata"))
load(paste0(dir.Input,"/TissueColours.Rdata"))
Piedata<-data.frame(lineage=names(table(AllCrispr$lineage)),value=as.numeric(table(AllCrispr$lineage)))

Piedata<-Piedata[Piedata$value!=0,]

Piedata$value<-Piedata$value/sum(Piedata$value)
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "none",
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
bp<- ggplot(Piedata, aes(x="", y=value, fill=lineage))+
  geom_bar(width = 1, stat = "identity")
#+scale_fill_manual(values=TissueColours)
pie<-bp + coord_polar("y", start=0)+  blank_theme +
  theme(axis.text.x=element_blank()) +
  geom_text(aes(label = paste(lineage,percent(value)), x = 1.6), position = position_fill(),size=2)
#geom_text(aes( y = value/3 + c(0, cumsum(value)[-length(value)]), 
#              label = percent(value)), size=3)
pdf(paste0(dir.Results,"/piechart_lineage.pdf"))
print(pie)
dev.off()

library(ggplot2)

g<-ggplot(AllCrisprcmp,aes(tissue))
g<-g+geom_bar(aes(fill=source))+
  coord_flip() +
  theme(legend.position = "top",
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
        legend.key = element_rect(fill = "transparent", colour = NA), # get rid of key legend fill, and of the surrounding
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold")) + 
  ylab("Number of Cell lines")+
  xlab("tissue")+
  scale_fill_manual("legend", values = InstituteColours)
pdf(paste0(dir.Results,"/Figure1_Barcharttissue.pdf"))
print(g)
dev.off()


#lung split by subtype - should be new groups to test
lungall<-rbind(BCrisprcmp[BCrisprcmp$tissue=="Lung",],SCrisprcmp[SCrisprcmp$tissue=="Lung",])


g<-ggplot(lungall,aes(cancer_type))
g<-g+geom_bar(aes(fill=source))+
  coord_flip() +
  theme(legend.position = "top",
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
        legend.key = element_rect(fill = "transparent", colour = NA), # get rid of key legend fill, and of the surrounding
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold")) + 
  ylab("Number of Cell lines")+
  xlab("Cancer Type")+
  scale_fill_manual("legend", values = InstituteColours)
pdf(paste0(dir.Results,"/Figure1_BarchartLung.pdf"))
print(g)
dev.off()


BreastSubtypes<-read.delim(file=paste0(dir.Input,"/BreastSubtypes.txt"),header=T,stringsAsFactors = F,sep="\t")
BreastSubtypes$source<-AllCrispr[match(BreastSubtypes[,1],AllCrispr$COSMICID),"source"]
BreastSubtypes<-BreastSubtypes[!is.na(BreastSubtypes$source),]
g<-ggplot(BreastSubtypes,aes(pam50))
g<-g+geom_bar(aes(fill=source))+
  coord_flip() +
  theme(legend.position = "top",
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
        legend.key = element_rect(fill = "transparent", colour = NA), # get rid of key legend fill, and of the surrounding
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold")) + 
  ylab("Number of Cell lines")+
  xlab("Cancer Type")+
  scale_fill_manual("legend", values = InstituteColours)
pdf(paste0(dir.Results,"/Figure1_BarchartPAM50.pdf"))
print(g)
dev.off()


