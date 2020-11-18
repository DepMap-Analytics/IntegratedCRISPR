


Broad90<-function(FCs,Essgenes,NonEssgenes){
  perfTH<-list()
  for(i in 1:ncol(FCs)){
    perfTH[[i]]<-bangmag.roc(BFs=as.matrix(FCs[,i]),ess_genes=Essgenes,non_ess_genes =NonEssgenes,CL=colnames(FCs)[i],th=0.05)[[2]][1]
  }
  names(perfTH)<-colnames(FCs)
  sFCs<-
    assembleDepletionMatrix(cellLines = colnames(FCs),
                            BFths = unlist(perfTH),
                            inputFolder = NULL,
                            BFvalues = FALSE,BFmatrix=FCs)
  
  CE_Int90<-ADAM2.PercentileCF(sFCs*(-1),display=FALSE)
  return(CE_Int90$cfgenes)
}

Broad90_LD<-function(FCs,Essgenes,NonEssgenes){
  perfTH<-list()
  for(i in 1:ncol(FCs)){
    perfTH[[i]]<-bangmag.roc(BFs=as.matrix(FCs[,i]),ess_genes=Essgenes,non_ess_genes =NonEssgenes,CL=colnames(FCs)[i],th=0.05)[[2]][1]
  }
  names(perfTH)<-colnames(FCs)
  sFCs<-
    assembleDepletionMatrix(cellLines = colnames(FCs),
                            BFths = unlist(perfTH),
                            inputFolder = NULL,
                            BFvalues = FALSE,BFmatrix=FCs)
  
  CE_Int90<-ADAM2.PercentileCF(sFCs*(-1),display=FALSE)
  return(CE_Int90$LeastDependent)
}

AdamTissue<-function(Binary,annot){
  CCRlabels<-data.frame(name=colnames(Binary),sidm=annot[match(colnames(Binary),annot$model_id),"tissue"],bid=annot[match(colnames(Binary),annot$BROAD_ID),"tissue"],stringsAsFactors = FALSE)
  
  CCRlabels[is.na(CCRlabels[,2]),2]<-0
  CCRlabels[is.na(CCRlabels[,3]),3]<-0
  utissues<-unique(annot$tissue)
  ntissues<-length(utissues)
  print(ntissues)
  coreCCR<-list()
  j=1
  tissuenames<-c()
  for(i in 1:ntissues){
    print(utissues[i])
    selCL<-CCRlabels[CCRlabels[,2]==utissues[i]|CCRlabels[,3]==utissues[i],1]
    selCL<-intersect(selCL,colnames(Binary))
    if(length(selCL)>9){
      coreCCR[[j]]<-ADAMwrapper(Binary[,selCL])
      j=j+1
      tissuenames<-c(tissuenames,utissues[i])
    }
  }
  names(coreCCR)<-tissuenames
  
  BinaryTissue<-matrix(0,nrow=nrow(Binary),ncol=length(coreCCR))
  rownames(BinaryTissue)<-rownames(Binary)
  colnames(BinaryTissue)<-names(coreCCR)
  for(i in 1:length(coreCCR)){
    BinaryTissue[coreCCR[[i]],i]<-1
  }
  PCcoreCCR<-ADAMwrapper(BinaryTissue)
  return(PCcoreCCR)
}

CoreStats<-function(Coregenes90,CoregenesAdam,allPriorSets,FDRlist,ssc){
  CoreGenes<-union(Coregenes90,CoregenesAdam)
  NumberCG<-length(CoreGenes)
  Recall_IntCeresD9<-length(intersect(allPriorSets,CoreGenes))/length(allPriorSets)
  FDRs<-lapply(FDRlist,function(x) length(intersect(x,CoreGenes))/length(CoreGenes))
  ssc<-ssc[order(ssc[,2],decreasing=F),]
  ssc<-ssc[ssc[,2]>=20&ssc[,2]<=200,]
  FDRcurve<-sapply(1:max(ssc[,2]),function(x) length(intersect(ssc[ssc[,2]>=x,1],CoreGenes))/length(CoreGenes))
  skewvals<-1:max(ssc[,2])
  PrecisionT1<-length(intersect(CoreGenes,allPriorSets))/length(CoreGenes)
  return(list(Recall=Recall_IntCeresD9,FDRs=FDRs,Precision=PrecisionT1,FDRcurve=FDRcurve,NumberCG=NumberCG,skewvals=skewvals))
}

RecallPlot<-function(inputlist,plotcolours,plotAnnotation){
  datanames<-names(inputlist)
  Recalls<-lapply(inputlist,function(x) x$Recall)
  RecallData<-data.frame(Recall=unlist(Recalls),data=datanames,preProc=plotAnnotation[datanames,"preProc"],batch=plotAnnotation[datanames,"batch"])
   
  RecallData$preProc<-factor(RecallData$preProc,levels=c("CCR","CCRJ","CERES"))

  RecallData$batch<-factor(RecallData$batch,levels=c("ComBat","ComBat+QN","ComBat+PC1","ComBat+PC2"))
  Recallplot<-ggplot(RecallData,aes(x=preProc,y=Recall,fill=batch))+geom_bar(stat="identity",position="dodge")+scale_fill_manual(values=plotcolours)+theme(axis.text.x=element_text(angle=45))+ylab("Recall")+xlab("")+theme_bw()+theme(legend.position = c(0.25,0.925),legend.background = element_blank(),legend.title=element_blank())
  
  pdf(paste0(dir.Results,"/RecallKnownSets.pdf"))
  print(Recallplot)
  dev.off()
  return(RecallData)
}

FDRPlot<-function(inputlist,plotcolours,plotAnnotation){
  datanames<-names(inputlist)
  FDRs1<-lapply(inputlist[1:3],function(x) x$FDRs)
  FDRs2<-lapply(inputlist[4:6],function(x) x$FDRs)
  FDRs3<-lapply(inputlist[7:9],function(x) x$FDRs)
  FDRs4<-lapply(inputlist[10:12],function(x) x$FDRs)
  FDRData1<-data.frame(FDR=unlist(FDRs1),preProc=c(rep("CCR",4),rep('CCRJ',4),rep('CERES',4)),fdrlevel=factor(rep(c("20",'50','100','200'),3),levels=c("20",'50','100','200')))
  FDRData2<-data.frame(FDR=unlist(FDRs2),preProc=c(rep("CCR",4),rep('CCRJ',4),rep('CERES',4)),fdrlevel=factor(rep(c("20",'50','100','200'),3),levels=c("20",'50','100','200')))
  FDRData3<-data.frame(FDR=unlist(FDRs3),preProc=c(rep("CCR",4),rep('CCRJ',4),rep('CERES',4)),fdrlevel=factor(rep(c("20",'50','100','200'),3),levels=c("20",'50','100','200')))
  FDRData4<-data.frame(FDR=unlist(FDRs4),preProc=c(rep("CCR",4),rep('CCRJ',4),rep('CERES',4)),fdrlevel=factor(rep(c("20",'50','100','200'),3),levels=c("20",'50','100','200')))
  
  FDRData1$preProc<-factor(FDRData1$preProc,levels=c("CCR","CCRJ","CERES"))
  FDRData2$preProc<-factor(FDRData2$preProc,levels=c("CCR","CCRJ","CERES"))
  FDRData3$preProc<-factor(FDRData3$preProc,levels=c("CCR","CCRJ","CERES"))
  FDRData4$preProc<-factor(FDRData4$preProc,levels=c("CCR","CCRJ","CERES"))
  

  
  FDRplot1<-ggplot(FDRData1,aes(x=fdrlevel,y=FDR,fill=preProc))+geom_bar(stat="identity",position="dodge")+scale_fill_manual(values=plotcolours)+theme(axis.text.x=element_text(angle=45))+ylab("Estimated FDR")+xlab("")+ggtitle("ComBat")+theme_bw()+theme(legend.position = c(0.25,0.925),legend.background = element_blank(),legend.title=element_blank())
  FDRplot2<-ggplot(FDRData2,aes(x=fdrlevel,y=FDR,fill=preProc))+geom_bar(stat="identity",position="dodge")+scale_fill_manual(values=plotcolours)+theme(axis.text.x=element_text(angle=45))+ylab("Estimated FDR")+xlab("")+ggtitle("ComBat+QN")+theme_bw()+theme(legend.position = c(0.25,0.925),legend.background = element_blank(),legend.title=element_blank())
  FDRplot3<-ggplot(FDRData3,aes(x=fdrlevel,y=FDR,fill=preProc))+geom_bar(stat="identity",position="dodge")+scale_fill_manual(values=plotcolours)+theme(axis.text.x=element_text(angle=45))+ylab("Estimated FDR")+xlab("")+ggtitle("ComBat+PC1")+theme_bw()+theme(legend.position = c(0.25,0.925),legend.background = element_blank(),legend.title=element_blank())
  FDRplot4<-ggplot(FDRData4,aes(x=fdrlevel,y=FDR,fill=preProc))+geom_bar(stat="identity",position="dodge")+scale_fill_manual(values=plotcolours)+theme(axis.text.x=element_text(angle=45))+ylab("Estimated FDR")+xlab("")+ggtitle("ComBat+PC2")+theme_bw()+theme(legend.position = c(0.25,0.925),legend.background = element_blank(),legend.title=element_blank())
  
  pdf(paste0(dir.Results,"/FDRrnai_CoreSets.pdf"))
    print(grid.arrange(FDRplot1,FDRplot2,FDRplot3,FDRplot4,nrow=2))
  dev.off()
  
}

PrecisionPlot<-function(inputlist,plotcolours,plotAnnotation){
  datanames<-names(inputlist)
  Precision<-lapply(inputlist,function(x) x$Precision)
  PData<-data.frame(Precision=unlist(Precision),data=datanames,preProc=plotAnnotation[datanames,"preProc"],batch=plotAnnotation[datanames,"batch"])
  
  PData$preProc<-factor(PData$preProc,levels=c("CCR","CCRJ","CERES"))
  
  PData$batch<-factor(PData$batch,levels=c("ComBat","ComBat+QN","ComBat+PC1","ComBat+PC2"))
  Pplot<-ggplot(PData,aes(x=preProc,y=Precision,fill=batch))+geom_bar(stat="identity",position="dodge")+scale_fill_manual(values=plotcolours)+theme(axis.text.x=element_text(angle=45))+ylab("Recall")+xlab("")+theme_bw()+theme(legend.position = c(0.25,0.925),legend.background = element_blank(),legend.title=element_blank())
  
  pdf(paste0(dir.Results,"/PrecisionKnownSets.pdf"))
  print(Pplot)
  dev.off()
  return(PData)
}

NCGPlot<-function(inputlist,plotcolours,plotAnnotation){
  datanames<-names(inputlist)
  Recalls<-lapply(inputlist,function(x) x$NumberCG)
  RecallData<-data.frame(Recall=unlist(Recalls),data=datanames,preProc=plotAnnotation[datanames,"preProc"],batch=plotAnnotation[datanames,"batch"])
  
  RecallData$preProc<-factor(RecallData$preProc,levels=c("CCR","CCRJ","CERES"))
  
  RecallData$batch<-factor(RecallData$batch,levels=c("ComBat","ComBat+QN","ComBat+PC1","ComBat+PC2"))
  Recallplot<-ggplot(RecallData,aes(x=preProc,y=Recall,fill=batch))+geom_bar(stat="identity",position="dodge")+scale_fill_manual(values=plotcolours)+theme(axis.text.x=element_text(angle=45))+ylab("Number Common Essentials")+xlab("")+theme_bw()+theme(legend.position = c(0.25,0.925),legend.background = element_blank(),legend.title=element_blank())
  
  pdf(paste0(dir.Results,"/NumberCoreEss.pdf"))
  print(Recallplot)
  dev.off()
  return(RecallData)
}

FDRcurvePlot<-function(inputlist,plotcolours,plotAnnotation){
  datanames<-names(inputlist)
  dfCT<-NULL
  dfCTQN<-NULL
  dfCTPC1<-NULL
  dfCTPC2<-NULL
  pP<-rep(c("CCR","CCRJ","CERES"),4)
  for(i in 1:3){
    temp<-unlist(inputlist[[i]]$FDRcurve)
    tempdf<-data.frame(FDR=temp,preProc=rep(pP[i],length(temp)),skewvals=unlist(inputlist[[i]]$skewvals))
    dfCT<-rbind(dfCT,tempdf)
  }
  for(i in 4:6){
    temp<-unlist(inputlist[[i]]$FDRcurve)
    tempdf<-data.frame(FDR=temp,preProc=rep(pP[i],length(temp)),skewvals=unlist(inputlist[[i]]$skewvals))
    dfCTQN<-rbind(dfCTQN,tempdf)
  }
  for(i in 7:9){
    temp<-unlist(inputlist[[i]]$FDRcurve)
    tempdf<-data.frame(FDR=temp,preProc=rep(pP[i],length(temp)),skewvals=unlist(inputlist[[i]]$skewvals))
    dfCTPC1<-rbind(dfCTPC1,tempdf)
  }
  for(i in 10:12){
    temp<-unlist(inputlist[[i]]$FDRcurve)
    tempdf<-data.frame(FDR=temp,preProc=rep(pP[i],length(temp)),skewvals=unlist(inputlist[[i]]$skewvals))
    dfCTPC2<-rbind(dfCTPC2,tempdf)
  }

  
  FDRplot1<-ggplot(dfCT,aes(x=skewvals,y=FDR,col=preProc))+geom_line()+scale_colour_manual(values=plotcolours)+theme(axis.text.x=element_text(angle=45))+ylab("Estimated FDR")+xlab("")+ggtitle("ComBat")+theme_bw()+theme(legend.position = "none",legend.background = element_blank(),legend.title=element_blank())
  FDRplot2<-ggplot(dfCTQN,aes(x=skewvals,y=FDR,col=preProc))+geom_line()+scale_colour_manual(values=plotcolours)+theme(axis.text.x=element_text(angle=45))+ylab("Estimated FDR")+xlab("")+ggtitle("ComBat+QN")+theme_bw()+theme(legend.position = "none",legend.background = element_blank(),legend.title=element_blank())
  FDRplot3<-ggplot(dfCTPC1,aes(x=skewvals,y=FDR,col=preProc))+geom_line()+scale_colour_manual(values=plotcolours)+theme(axis.text.x=element_text(angle=45))+ylab("Estimated FDR")+xlab("")+ggtitle("ComBat+PC1")+theme_bw()+theme(legend.position = "none",legend.background = element_blank(),legend.title=element_blank())
  FDRplot4<-ggplot(dfCTPC2,aes(x=skewvals,y=FDR,col=preProc))+geom_line()+scale_colour_manual(values=plotcolours)+theme(axis.text.x=element_text(angle=45))+ylab("Estimated FDR")+xlab("")+ggtitle("ComBat+PC2")+theme_bw()+theme(legend.position = "none",legend.background = element_blank(),legend.title=element_blank())
  
  pdf(paste0(dir.Results,"/FDRrnai_CoreSetsCurves.pdf"))
  print(grid.arrange(FDRplot1,FDRplot2,FDRplot3,FDRplot4,nrow=2))
  dev.off()
  return(list(ComBat=dfCT,ComBatQC=dfCTQN,ComBatPC1=dfCTPC1,ComBatPC2=dfCTPC2))
}

