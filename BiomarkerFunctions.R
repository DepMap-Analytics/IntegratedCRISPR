MakeNormLRTplot<-function(lrtobjs,filename){
  ccrnt<-lrtobjs[[1]]
  ceresnt<-lrtobjs[[2]]
  ccrjnt<-lrtobjs[[3]]
  ccr1<-sum(ccrnt<=100)
  ceres1<-sum(ceresnt<=100)
  ccrj1<-sum(ccrjnt<=100)
  ccr2<-sum(ccrnt>100&ccrnt<=200)
  ceres2<-sum(ceresnt>100&ceresnt<=200)
  ccrj2<-sum(ccrjnt>100&ccrjnt<=200)
  ccr3<-sum(ccrnt>200&ccrnt<=500)
  ceres3<-sum(ceresnt>200&ceresnt<=500)
  ccrj3<-sum(ccrjnt>200&ccrjnt<=500)
  ccr4<-sum(ccrnt>500)
  ceres4<-sum(ceresnt>500)
  ccrj4<-sum(ccrjnt>500)
  plotnt<-rbind(data.frame(x="<=100",v=ccr1,method="CCR"),data.frame(x="<=100",v=ceres1,method="CERES"),data.frame(x="<=100",v=ccrj1,method="CCRJ"),
                data.frame(x=">100<=200",v=ccr2,method="CCR"),data.frame(x=">100<=200",v=ceres2,method="CERES"),data.frame(x=">100<=200",v=ccrj2,method="CCRJ"),
                data.frame(x=">200<=500",v=ccr3,method="CCR"),data.frame(x=">200<=500",v=ceres3,method="CERES"),data.frame(x=">200<=500",v=ccrj3,method="CCRJ"),
                data.frame(x=">500",v=ccr4,method="CCR"),data.frame(x=">500",v=ceres4,method="CERES"),data.frame(x=">500",v=ccrj4,method="CCRJ"))
  
  normlrtplot<-ggplot(plotnt,aes(x=x,y=v,fill=method))+geom_bar(stat="identity",position="dodge")+scale_fill_manual(values=PipelineColours)+theme(axis.text.x=element_text(angle=45))+ylab("Number Genes")+xlab("Norm LRT")+theme_bw()+theme(legend.position = c(0.8,0.8),legend.background = element_blank())
  pdf(paste0(dir.Results,"/",filename,"NormLRTGene_plot.pdf"))
  print(normlrtplot)
  dev.off()
  plotnt$logcount<-log(plotnt$v)
  normloglrtplot<-ggplot(plotnt,aes(x=x,y=logcount,fill=method))+geom_bar(stat="identity",position="dodge")+scale_fill_manual(values=PipelineColours)+theme(axis.text.x=element_text(angle=45))+ylab("Number Genes")+xlab("Norm LRT")+theme_bw()+theme(legend.position = c(0.8,0.8),legend.background = element_blank())
  pdf(paste0(dir.Results,"/",filename,"NormLRTGene_plotLog.pdf"))
  print(normloglrtplot)
  dev.off()
  ngenes<-length(ccrnt)
  p3<-prop.test(x=c(ccr2,ceres2),n=c(ngenes,ngenes))
  p4<-prop.test(x=c(ccrj2,ceres2),n=c(ngenes,ngenes))
  p5<-prop.test(x=c(ccrj2,ccr2),n=c(ngenes,ngenes))
  
  cat(paste("CCR versus CERES proportion test, number normLRT:",filename,p3$estimate,"p-value",p3$p.value),file=paste0(dir.Results,"/Figure4Log.txt"),sep="\n",append=TRUE)
  cat(paste("CCR_JACKS versus CERES proportion test, number normLRT:",filename,p4$estimate,"p-value",p4$p.value),file=paste0(dir.Results,"/Figure4Log.txt"),sep="\n",append=TRUE)
  cat(paste("CCR_JACKS versus CCR proportion test, number normLRT:",filename,p5$estimate,"p-value",p5$p.value),file=paste0(dir.Results,"/Figure4Log.txt"),sep="\n",append=TRUE)
  
  
}
getSDDgenes<-function(normLRTlist){
  normLRTCCR<-normLRTlist[[1]]
  normLRTCERES<-normLRTlist[[2]]
  normLRTCCRJ<-normLRTlist[[3]]
  n1<-names(normLRTCCR[[3]])[normLRTCCR[[3]]>200]
  n2<-names(normLRTCERES[[3]])[normLRTCERES[[3]]>200]
  n3<-names(normLRTCCRJ[[3]])[normLRTCCRJ[[3]]>200]
  geneCCRceres<-union(n1,n2)
  geneCCRccrj<-union(n1,n3)
  geneCERESccrj<-union(n3,n2)
  usegenessd<-union(geneCCRceres,geneCCRccrj)
  return(usegenessd)
}
BiomarkerTests<-function(MoBEM,FCdata,cmp,name,usegenessd,fdr=0.05){

  
  CCR_allCosmic<-FCdata[[1]]
  CERES_allCosmic<-FCdata[[2]]
  CCRJ_allCosmic<-FCdata[[3]]
  colnames(CCR_allCosmic)<-cmp[match(colnames(CCR_allCosmic),cmp$model_id),"COSMIC_ID"]
  colnames(CERES_allCosmic)<-cmp[match(colnames(CERES_allCosmic),cmp$model_id),"COSMIC_ID"]
  colnames(CCRJ_allCosmic)<-cmp[match(colnames(CCRJ_allCosmic),cmp$model_id),"COSMIC_ID"]
  All_ccr<-Allanova(MoBEM,CCR_allCosmic[usegenessd,])
  All_ceres<-Allanova(MoBEM,CERES_allCosmic[usegenessd,])
  All_ccrJ<-Allanova(MoBEM,CCRJ_allCosmic[usegenessd,])
  cat(paste("Number of ssd genes tested in pancancer ANOVA",name,length(usegenessd)),file=paste0(dir.Results,"/Figure4Log.txt"),sep="\n",append=TRUE)
  rownames(All_ccr)<-paste(All_ccr$CFE,All_ccr$GENE,sep="-")
  rownames(All_ceres)<-paste(All_ceres$CFE,All_ceres$GENE,sep="-")
  rownames(All_ccrJ)<-paste(All_ccrJ$CFE,All_ccrJ$GENE,sep="-")
  alltested<-intersect(rownames(All_ccr),intersect(rownames(All_ceres),rownames(All_ccrJ)))
  All_ccr<-All_ccr[alltested,]
  All_ceres<-All_ceres[alltested,]
  All_ccrJ<-All_ccrJ[alltested,]
  cor(All_ccr$delta,All_ceres[rownames(All_ccr),"delta"],method="spearman")
  cor(All_ccrJ$delta,All_ceres[rownames(All_ccrJ),"delta"],method="spearman")
  
  cor(All_ccr$delta,All_ccrJ[rownames(All_ccr),"delta"],method="spearman")
  
  c1<-cor.test(All_ccr$delta,All_ceres[rownames(All_ccr),"delta"],method="spearman")
  c2<-cor.test(All_ccrJ$delta,All_ceres[rownames(All_ccrJ),"delta"],method="spearman")
  
  c3<-cor.test(All_ccr$delta,All_ccrJ[rownames(All_ccr),"delta"],method="spearman")
  
  cat(paste("CCR versus CERES biomarker deltas, spearman correlation:",name,c1$estimate,"p-value",c1$p.value),file=paste0(dir.Results,"/Figure4Log.txt"),sep="\n",append=TRUE)
  cat(paste("CCR_JACKS versus CERES biomarker deltas, spearman correlation:",name,c2$estimate,"p-value",c2$p.value),file=paste0(dir.Results,"/Figure4Log.txt"),sep="\n",append=TRUE)
  cat(paste("CCR versus CCR_JACKS biomarker deltas, spearman correlation:",name,c3$estimate,"p-value",c3$p.value),file=paste0(dir.Results,"/Figure4Log.txt"),sep="\n",append=TRUE)
  
  All_ccr$fdr<-p.adjust(All_ccr$p,method="fdr")
  All_ceres$fdr<-p.adjust(All_ceres$p,method="fdr")
  All_ccrJ$fdr<-p.adjust(All_ccrJ$p,method="fdr")
  
  nBM_CCR<-sum(All_ccr$fdr<fdr)
  nBM_CERES<-sum(All_ceres$fdr<fdr)
  nBM_CCRJ<-sum(All_ccrJ$fdr<fdr)
  cat(paste("Number of associations tested: ",name,nrow(All_ccr)),file=paste0(dir.Results,"/Figure4Log.txt"),sep="\n",append=TRUE)
  cat(paste("Number signif biomarkers assoc CCR",name,nBM_CCR),file=paste0(dir.Results,"/Figure4Log.txt"),sep="\n",append=TRUE)
  cat(paste("Number signif biomarkers assoc CERES",name,nBM_CERES),file=paste0(dir.Results,"/Figure4Log.txt"),sep="\n",append=TRUE)
  cat(paste("Number signif biomarkers assoc CCR_JACKS",name,nBM_CCRJ),file=paste0(dir.Results,"/Figure4Log.txt"),sep="\n",append=TRUE)
  
  prop.test(c(nBM_CCR,nBM_CERES),c(nrow(All_ccr),nrow(All_ccr)))
  prop.test(c(nBM_CCR,nBM_CCRJ),c(nrow(All_ccr),nrow(All_ccr)))
  prop.test(c(nBM_CERES,nBM_CCRJ),c(nrow(All_ccr),nrow(All_ccr)))
  write.table(All_ccr[All_ccr$fdr<fdr,],file=paste0(dir.Results,"/SuppTable_SigCCRbiomarkers",name,".csv"),quote=F)
  write.table(All_ceres[All_ceres$fdr<fdr,],file=paste0(dir.Results,"/SuppTable_SigCERESbiomarkers",name,".csv"),quote=F)
  write.table(All_ccrJ[All_ccrJ$fdr<fdr,],file=paste0(dir.Results,"/SuppTable_SigCCR_JACKSbiomarkers",name,".csv"),quote=F)
  return(list(All_ccr,All_ceres,All_ccrJ,alltested))
  
  
}
AnovaProp<-function(AnovaRes,SetName,dir.Results,fdr=0.05){
  All_ccr<-AnovaRes[[1]]
  All_ceres<-AnovaRes[[2]]
  All_ccrJ<-AnovaRes[[3]]
  nBM_CCR<-sum(All_ccr$fdr<fdr)
  nBM_CERES<-sum(All_ceres$fdr<fdr)
  nBM_CCRJ<-sum(All_ccrJ$fdr<fdr)
  pt1<-prop.test(c(nBM_CCR,nBM_CERES),c(nrow(All_ccr),nrow(All_ccr)))
  pt2<-prop.test(c(nBM_CCR,nBM_CCRJ),c(nrow(All_ccr),nrow(All_ccr)))
  pt3<-prop.test(c(nBM_CERES,nBM_CCRJ),c(nrow(All_ccr),nrow(All_ccr)))
  cat(paste("CCR versus CERES proportion test, number bm:",SetName,pt1$estimate,"p-value",pt1$p.value),file=paste0(dir.Results,"/Figure4Log.txt"),sep="\n",append=TRUE)
  cat(paste("CCR_JACKS versus CERES proportion test, number bm:",SetName,pt3$estimate,"p-value",pt3$p.value),file=paste0(dir.Results,"/Figure4Log.txt"),sep="\n",append=TRUE)
  cat(paste("CCR_JACKS versus CCR proportion test, number bm:",SetName,pt2$estimate,"p-value",pt2$p.value),file=paste0(dir.Results,"/Figure4Log.txt"),sep="\n",append=TRUE)
  
}
AnovaFDRcurve<-function(AnovaResList,ResultNames=c("CCR","CERES","CCRJ"),batchlev=c("ComBat","ComBatQN","ComBatPC1","ComBatPC2"),plotcolours=NULL,noCFE=FALSE,plotName=""){
  fdrlevels<-seq(0,0.1,by=0.0001)
  batchnames<-names(AnovaResList)
  preProc<-c("CCR","CERES","CCRJ")
  BMData<-NULL
  for(i in 1:length(AnovaResList)){
    for(j in 1:3){
      if(noCFE){
        temp<-data.frame(nmarkers=unlist(sapply(fdrlevels,function(x) length(unique(AnovaResList[[i]][[j]][AnovaResList[[i]][[j]]$fdr<x,"CFE"])))),preProc=preProc[j],batch=batchnames[i],fdr=fdrlevels)
        
      }else{
        temp<-data.frame(nmarkers=unlist(sapply(fdrlevels,function(x) sum(AnovaResList[[i]][[j]]$fdr<x))),preProc=preProc[j],batch=batchnames[i],fdr=fdrlevels)
      }
      BMData<-rbind(BMData,temp)
    }
    
  }
  BMData$preProc<-factor(BMData$preProc,levels=c("CCR","CCRJ","CERES"))
  
  BMData$batch<-factor(BMData$batch,levels=batchlev)
  for(i in 1:length(batchlev)){
    BMDatab<-BMData[BMData$batch==batchlev[i],]
    BMplot<-ggplot(BMDatab,aes(x=fdr,y=nmarkers))+geom_line(aes(color=preProc))+scale_color_manual(values=plotcolours)+ylab("Number Significant Associations")+xlab("FDR level")+theme_bw()+theme(legend.position = "none",legend.background = element_blank(),legend.title=element_blank())
    if(noCFE){
      pdf(paste0(dir.Results,"/NumberCFEs_",batchlev[i],plotName,".pdf"))
      print(BMplot)
      dev.off()
    }else{
    pdf(paste0(dir.Results,"/NumberBiomarkers_",batchlev[i],plotName,".pdf"))
    print(BMplot)
    dev.off()
    }
  }

  

}

plotBMres<-function(AnovaResList,plotcolours,fdr=0.05,preProc=c("CCR","CERES","CCRJ"),plotName="",CFEtype="All"){
  batchnames<-names(AnovaResList)

  BMData<-NULL
  for(i in 1:length(AnovaResList)){
    for(j in 1:3){
      if(CFEtype=="All"){
        temp<-data.frame(nmarkers=sum(AnovaResList[[i]][[j]]$fdr<fdr),preProc=preProc[j],batch=batchnames[i])
        BMData<-rbind(BMData,temp)
      }
      if(CFEtype%in%c("_mut","_HypMET")){
        Ares<-AnovaResList[[i]][[j]]
        Ares<-Ares[grep(CFEtype,Ares$CFE),]
        temp<-data.frame(nmarkers=sum(Ares$fdr<fdr),preProc=preProc[j],batch=batchnames[i])
        BMData<-rbind(BMData,temp)
      }
      if(CFEtype=="CNA"){
        Ares<-AnovaResList[[i]][[j]]
        Ares<-Ares[grep("_mut",Ares$CFE,invert=TRUE),]
        Ares<-Ares[grep("_HypMET",Ares$CFE,invert=TRUE),]
        Ares<-Ares[grep("MSI",Ares$CFE,invert=TRUE),]
        temp<-data.frame(nmarkers=sum(Ares$fdr<fdr),preProc=preProc[j],batch=batchnames[i])
        BMData<-rbind(BMData,temp)
      }
    }
    
  }
  
  BMData$preProc<-factor(BMData$preProc,levels=c("CCR","CCRJ","CERES"))
  
  BMData$batch<-factor(BMData$batch,levels=c("ComBat","ComBatQN","ComBatPC1","ComBatPC2"))

  BMplot<-ggplot(BMData,aes(x=preProc,y=nmarkers,fill=batch))+geom_bar(stat="identity",position="dodge")+scale_fill_manual(values=plotcolours)+theme(axis.text.x=element_text(angle=45))+ylab("Number Significant Associations")+xlab("")+theme_bw()+theme(legend.position = "none",legend.background = element_blank(),legend.title=element_blank())
  

  pdf(paste0(dir.Results,"/NumberBiomarkers",plotName,".pdf"))
  print(BMplot)
  dev.off()
}
plotBMresByCFE<-function(AnovaResList,plotcolours,fdr=0.05,preProc=c("CCR","CERES","CCRJ")){
  batchnames<-names(AnovaResList)

  BMData<-NULL
  cfe<-list()
  cfev<-list()
  for(i in 1:length(AnovaResList)){
    cfe[[i]]<-list()
    cfev[[i]]<-list()
    for(j in 1:3){
      siglist<-AnovaResList[[i]][[j]][AnovaResList[[i]][[j]]$fdr<fdr,]
      cfev[[i]][[j]]<-unique(siglist$CFE)
      cfe[[i]][[j]]<-length(unique(siglist$CFE))
      temp<-data.frame(nmarkers=length(unique(siglist$CFE)),preProc=preProc[j],batch=batchnames[i])
      BMData<-rbind(BMData,temp)
      
    }
    
  }
  
  BMData$preProc<-factor(BMData$preProc,levels=c("CCR","CCRJ","CERES"))
  
  BMData$batch<-factor(BMData$batch,levels=c("ComBat","ComBatQN","ComBatPC1","ComBatPC2"))
  
  BMplot<-ggplot(BMData,aes(x=preProc,y=nmarkers,fill=batch))+geom_bar(stat="identity",position="dodge")+scale_fill_manual(values=plotcolours)+theme(axis.text.x=element_text(angle=45))+ylab("Number Significant Associations")+xlab("")+theme_bw()+theme(legend.position = "none",legend.background = element_blank(),legend.title=element_blank())
  
  
  pdf(paste0(dir.Results,"/NumberBiomarkersCFEs.pdf"))
  print(BMplot)
  dev.off()
  cfeoverlap<-list()
  cfediff<-list()
  namesl<-c()
  for(i in 1:length(AnovaResList)){

    combinations<-cbind(combn(1:3,2),combn(3:1,2))
    cfeoverlap[[i]]<-apply(combinations,2,function(x) c(length(intersect(cfev[[i]][[x[1]]],cfev[[i]][[x[2]]])),paste(preProc[x[1]],preProc[x[2]])))
    cfediff[[i]]<-apply(combinations,2,function(x) setdiff(cfev[[i]][[x[1]]],cfev[[i]][[x[2]]]))
    namesl<-c(namesl,apply(combinations,2,function(x) paste(preProc[x[1]],preProc[x[2]])))
  }
 
  return(list(cfe=cfe,cfeoverlap=cfeoverlap,cfediff=cfediff,namesl=namesl))
}

CFEbm<-function(AnovaResList,cfes,fdr=0.1,preProc=c("CCR","CERES","CCRJ")){
  batchnames<-names(AnovaResList)

  BMData<-NULL
  cferes<-list()
  cfev<-list()
  for(i in 1:length(AnovaResList)){
    cferes[[i]]<-list()
    cfev[[i]]<-list()
    for(j in 1:3){
      siglist<-AnovaResList[[i]][[j]][AnovaResList[[i]][[j]]$CFE==cfes&AnovaResList[[i]][[j]]$fdr<fdr,]
      if(nrow(siglist)>0){
        cferes[[i]][[j]]<-siglist
      }else{
        cferes[[i]][[j]]<-1
      }

    }

  }
  return(cferes)
  
}

scatterBMres<-function(Res1,Res2,column,pname,preprocnames,fdr=0.05){
  sig1<-rownames(Res1)[Res1$fdr<fdr]
  sig2<-rownames(Res2)[Res2$fdr<fdr]
  allassoc<-union(sig1,sig2)
  both<-intersect(sig1,sig2)
  r1only<-setdiff(sig1,sig2)
  r2only<-setdiff(sig2,sig1)
  bm1only<-setdiff(unique(Res1[sig1,"CFE"]),unique(Res2[sig2,"CFE"]))
  bm2only<-setdiff(unique(Res2[sig2,"CFE"]),unique(Res1[sig1,"CFE"]))
  pointcols<-rep("grey",length(allassoc))
  names(pointcols)<-allassoc
  pointcols[r1only]<-"blue"
  pointcols[r2only]<-"red"
  pdf(paste0(dir.Results,"/",pname,".pdf"),useDingbats = FALSE)
    plot(Res1[allassoc,column],Res2[allassoc,column],col=pointcols,xlab=paste(preprocnames[1],column),ylab=paste(preprocnames[2],column))
    if(column=="effect_size"){
      abline(0,1)
    }
  dev.off()
  return(list(sig1=Res1[r1only,],sig2=Res2[r2only,],bm1only=bm1only,bm2only=bm2only))
}

BiomarkerTissue<-function(MoBEM,FCdata,cmp,name,usegenessd,fdr=0.05){
  CCR_allCosmic<-FCdata[[1]]
  CERES_allCosmic<-FCdata[[2]]
  CCRJ_allCosmic<-FCdata[[3]]
  colnames(CCR_allCosmic)<-cmp[match(colnames(CCR_allCosmic),cmp$model_id),"COSMIC_ID"]
  colnames(CERES_allCosmic)<-cmp[match(colnames(CERES_allCosmic),cmp$model_id),"COSMIC_ID"]
  colnames(CCRJ_allCosmic)<-cmp[match(colnames(CCRJ_allCosmic),cmp$model_id),"COSMIC_ID"]
  utissues<-unique(cmp$tissue)
  j=1
  ccr<-list()
  ceres<-list()
  ccrj<-list()
  for(i in 1:length(utissues)){
    cosmiccl<-unique(cmp[cmp$tissue==utissues[i],"COSMIC_ID"])
    cosmicavail1<-unique(colnames(CCR_allCosmic)[colnames(CCR_allCosmic)%in%cosmiccl])
    cosmicavail2<-unique(colnames(CERES_allCosmic)[colnames(CERES_allCosmic)%in%cosmiccl])
    cosmicavail3<-unique(colnames(CCRJ_allCosmic)[colnames(CCRJ_allCosmic)%in%cosmiccl])
    if(length(cosmicavail1)>9&length(cosmicavail2)>9&length(cosmicavail3)>9){
      
      ccr[[j]]<-Allanova(MoBEM,CCR_allCosmic[usegenessd,colnames(CCR_allCosmic)%in%cosmiccl],utissues[i])
      ceres[[j]]<-Allanova(MoBEM,CERES_allCosmic[usegenessd,colnames(CERES_allCosmic)%in%cosmiccl],utissues[i])
      ccrj[[j]]<-Allanova(MoBEM,CCRJ_allCosmic[usegenessd,colnames(CCRJ_allCosmic)%in%cosmiccl],utissues[i])
      j=j+1
    }
  }
  cclnames<-c()
  for(i in 1:length(utissues)){
    cosmiccl<-unique(cmp[cmp$tissue==utissues[i],"COSMIC_ID"])
    cosmicavail1<-unique(colnames(CCR_allCosmic)[colnames(CCR_allCosmic)%in%cosmiccl])
    cosmicavail2<-unique(colnames(CERES_allCosmic)[colnames(CERES_allCosmic)%in%cosmiccl])
    cosmicavail3<-unique(colnames(CCRJ_allCosmic)[colnames(CCRJ_allCosmic)%in%cosmiccl])
    if(length(cosmicavail1)>9&length(cosmicavail2)>9&length(cosmicavail3)>9){
      cclnames<-c(cclnames,utissues[i])
    }
  }
  names(ccr)<-cclnames
  names(ceres)<-cclnames
  names(ccrj)<-cclnames
  for(i in 1:length(ccr)){
    
    temp<-ccr[[i]]$p
    if(length(temp)>0){
      new<-ccr[[i]]
      new$fdr<-p.adjust(temp,method="fdr")
      ccr[[i]]<-new
      
    }
  }
  for(i in 1:length(ceres)){
    
    temp<-ceres[[i]]$p
    if(length(temp)>0){
      new<-ceres[[i]]
      new$fdr<-p.adjust(temp,method="fdr")
      ceres[[i]]<-new
      
    }
  }
  for(i in 1:length(ccrj)){
    
    temp<-ccrj[[i]]$p
    if(length(temp)>0){
      new<-ccrj[[i]]
      new$fdr<-p.adjust(temp,method="fdr")
      ccrj[[i]]<-new
      
    }
  }
  ###set test names:
  ccrN<-lapply(ccr,function(x) {rownames(x)<-paste(x$CFE,paste(x$GENE,x$tissue,sep="-"),sep="-")
  return(x)})
  
  ccrS<-lapply(ceres,function(x) {rownames(x)<-paste(x$CFE,paste(x$GENE,x$tissue,sep="-"),sep="-")
  return(x)})
  
  ccrB<-lapply(ccrj,function(x) {rownames(x)<-paste(x$CFE,paste(x$GENE,x$tissue,sep="-"),sep="-")
  return(x)})
  
  names(ccrN)<-names(ccr)
  names(ccrS)<-names(ceres)
  names(ccrB)<-names(ccrj)
  
  return(list(CCR=ccrN,CERES=ccrS,CCRJ=ccrB))
}
