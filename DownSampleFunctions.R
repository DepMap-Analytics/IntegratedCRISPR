batch_downsample<-function(mergeddata,empBayes=FALSE,numberSamples=100,annotation,DownRange=seq(0.1,0.9,by=0.05),allBroad,allSanger,OrigCorrection,cmpAnnot){
  nsamples<-nrow(annotation)
  annotAll<-rbind(annotation,annotation)
  annotAll$cn<-c(annotation$Smn,annotation$Bmn)
  annotAll$cn<-gsub(" ","",annotAll$cn)

  rownames(annotAll)<-annotAll$cn
  annottab<-table(annotation$tissue)
  lineageFreq<-as.vector(annottab)/nsamples
  names(lineageFreq)<-names(annottab)
  sampleprobs<-lineageFreq[annotation$tissue]
  names(sampleprobs)<-annotation$model_id
  lineagelabels<-cmpAnnot[c(colnames(allBroad),colnames(allSanger)),"tissue"]
  names(lineagelabels)<-c(colnames(allBroad),colnames(allSanger))
  batchlabel<-c(rep(1,ncol(allBroad)),rep(2,ncol(allSanger)))
  names(batchlabel)<-names(lineagelabels)
  cor_res<-list()
  euc_res<-list()
  cluster_resLineage<-list()
  cluster_resBatch<-list()
  for(i in 1:length(DownRange)){
      samples_gen<-sapply(1:numberSamples,function(y) annotAll[annotAll$model_id%in%sample(annotation$model_id,size=round(nsamples*DownRange[i]),prob=sampleprobs),"cn"])
      #samples_gen matrix nrow sample size, ncol numberSamples
      cor_res[[i]]<-vector("numeric",length = numberSamples)
      euc_res[[i]]<-vector("numeric",length = numberSamples)
      cluster_resLineage[[i]]<-vector("numeric",length = numberSamples)
      cluster_resBatch[[i]]<-vector("numeric",length = numberSamples)
      for(j in 1:numberSamples){
        result<-BatchCorrection(data1=allBroad,data2=allSanger,CombatRes=ComBatCP(mergeddata[,samples_gen[,j]],batch=unlist(sapply(samples_gen[,j],function(y)strsplit(y,"---")[[1]][2])),empBayes=empBayes))
        perf_res<-compareCD(OrigCorrection,result$qNorm,DownRange = DownRange)
        #cluster_resL<-ASWpc(result$qNorm,numberPCs = 20,lineagelabels)
       
        cluster_resB<-ASWpc(result$qNorm,numberPCs = 20,lineagelabels=batchlabel)
        cor_res[[i]][j]<-perf_res$corRes
        euc_res[[i]][j]<-perf_res$eucRes
        #cluster_resLineage[[i]][j]<-cluster_resL
        cluster_resLineage[[i]][j]<-i+j
        cluster_resBatch[[i]][j]<-cluster_resB
      }

      
  }
  #ideal output is a list of matrices. Each matrix is one parameter with Sample value x DownRange
  return(list(corRes=cor_res,eucRes=euc_res,clusterL=cluster_resLineage,clusterB=cluster_resBatch))
}

compareCD<-function(origCorrection,sampleCorrection,DownRange=seq(0.1,0.9,by=0.05),numberPCs=20,lineagelabels){
  corRes<-list()
  eucRes<-list()

  for(i in 1:length(DownRange)){
    sampleRes<-sampleCorrection
    corRes[[i]]<-unlist(sapply(colnames(sampleRes),function(y) cor(sampleRes[,y],origCorrection[,y],use="pairwise")))
    eucRes[[i]]<-unlist(sapply(colnames(sampleRes),function(y) dist(t(na.omit(cbind(sampleRes[,y],origCorrection[,y]))))))
  }
  return(list(corRes=corRes,eucRes=eucRes))
}


batch_downsampleAMI<-function(mergeddata,empBayes=FALSE,numberSamples=100,annotation,DownRange=seq(0.1,0.9,by=0.05),allBroad,allSanger,OrigCorrection,cmpAnnot){
  nsamples<-nrow(annotation)
  annotAll<-rbind(annotation,annotation)
  annotAll$cn<-c(annotation$Smn,annotation$Bmn)
  annotAll$cn<-gsub(" ","",annotAll$cn)
  
  rownames(annotAll)<-annotAll$cn
  annottab<-table(annotation$tissue)
  lineageFreq<-as.vector(annottab)/nsamples
  names(lineageFreq)<-names(annottab)
  sampleprobs<-lineageFreq[annotation$tissue]
  names(sampleprobs)<-annotation$model_id
  lineagelabels<-cmpAnnot[c(colnames(allBroad),colnames(allSanger)),"tissue"]
  names(lineagelabels)<-c(colnames(allBroad),colnames(allSanger))

  cluster_AMILineage<-list()

  for(i in 1:length(DownRange)){
    samples_gen<-sapply(1:numberSamples,function(y) annotAll[annotAll$model_id%in%sample(annotation$model_id,size=round(nsamples*DownRange[i]),prob=sampleprobs),"cn"])
    #samples_gen matrix nrow sample size, ncol numberSamples

    cluster_AMILineage[[i]]<-vector("numeric",length = numberSamples)
   
    for(j in 1:numberSamples){
      result<-BatchCorrection(data1=allBroad,data2=allSanger,CombatRes=ComBatCP(mergeddata[,samples_gen[,j]],batch=unlist(sapply(samples_gen[,j],function(y)strsplit(y,"---")[[1]][2])),empBayes=empBayes))
  
      cluster_resL<-RepeatKmeansAMI(list(result$qNorm),genes="var",nrepeat=20,length(lineageFreq),lineagelabels)

      cluster_AMILineage[[i]][j]<-median(unlist(cluster_resL[[1]]$ami))

    }
    
    
  }
  #ideal output is a list of matrices. Each matrix is one parameter with Sample value x DownRange
  return(cluster_AMILineage)
}

ASWpc<-function(inputdata,numberPCs=20,lineagelabels){

  inputdata<-inputdata[,names(lineagelabels)[!is.na(lineagelabels)]]
  lineagelabels<-lineagelabels[!is.na(lineagelabels)]
  if(sum(is.na(inputdata))!=0){
    #Have NAs and need to impute missing values
    #data is genes x cell lines
    meanVals<-rowMeans(inputdata,na.rm=TRUE)
    genesToimpute<-which(rowSums(is.na(inputdata))!=0)
    for(i in 1:length(genesToimpute)){
      selcl<-which(is.na(inputdata[genesToimpute[i],]))
      inputdata[genesToimpute[i],selcl]<-meanVals[genesToimpute[i]]
    }
  }
  estpca<-prcomp(t(inputdata),scale.=TRUE)
  
  subsetPCs <- estpca$x[,1:numberPCs]
  
  #do the clustering and silhouette values:
  distM<-dist(subsetPCs)
  res<-silhouetteScores(lineagelabels,distM)

  return(res)
}

plotDS<-function(compareCD,DownRange=seq(0.1,0.9,by=0.05),nsamples=NULL,plotname,ylim=c(0.98,1),setNames=NULL,listRes=FALSE,noAdj=NULL,ylab="",plotcolours=NULL,axisS=10,labelS=10,pwidth=5,pheight=5){
  xfactors<-paste0("Sample",DownRange)
  if(!is.null(nsamples)){
    xfactors<-paste0("N=",round(nsamples*DownRange))
  }
  if(!listRes){
    names(compareCD)<-xfactors
    corvals<-c()
    xfvals<-c()
    for(i in 1:length(xfactors)){
      corvals<-c(corvals,unlist(compareCD[[i]]))
      xfvals<-c(xfvals,rep(xfactors[i],length(unlist(compareCD[[i]]))))
    
    }
    if(!is.null(noAdj)){
      #also add in the values with correlations with no ComBat adjustment
      xfno<-"N=0"
      corvals<-c(corvals,unlist(noAdj))
      xfvals<-c(xfvals,rep(xfno,length(unlist(noAdj))))
    }
    plotdata<-data.frame(x=corvals,value=xfvals)
    if(!is.null(noAdj)){
      plotdata$value<-factor(plotdata$value,levels=c(xfno,xfactors))
    }else{
      plotdata$value<-factor(plotdata$value,levels=xfactors)
    }
    plotDS<-ggplot(aes(x=value,y=x),data=plotdata)+geom_boxplot()+theme_bw()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.title = element_text(size = axisS),axis.text = element_text(size = labelS))+ylim(ylim)+labs(x="",y=ylab)
    if(!is.null(plotcolours)){
      plotDS<-plotDS+scale_fill_manual(plotcolours)
    }
    savepdf(plotname,width=pwidth,height=pheight)
      print(plotDS)
    dev.off()
  }else{
    #have a set of correlation downsamples e.g. different pre-processing methods
    nSets<-length(compareCD)
    corvals<-c()
    xfvals<-c()
    setvals<-c()
    for(i in 1:nSets){
      CD<-compareCD[[i]]
      if(!is.null(noAdj)){
        NC<-noAdj[[i]]
      }
      cvals<-c()
      xvals<-c()
      for(j in 1:length(xfactors)){
        cvals<-c(cvals,unlist(CD[[j]]))
        xvals<-c(xvals,rep(xfactors[j],length(unlist(CD[[j]]))))
        
      }
      if(!is.null(noAdj)){
        cvals<-c(cvals,NC)
        xvals<-c(xvals,rep("N=0",length(NC)))
      }
      setvals<-c(setvals,rep(setNames[i],length(cvals)))
      corvals<-c(corvals,cvals)
      xfvals<-c(xfvals,xvals)
      
    }
    plotdata<-data.frame(x=corvals,value=xfvals,set=setvals)
    if(!is.null(noAdj)){
      plotdata$value<-factor(plotdata$value,levels=c("N=0",xfactors))
    }else{
    plotdata$value<-factor(plotdata$value,levels=xfactors)}
    plotDS<-ggplot(aes(x=value,y=x,fill=setvals),data=plotdata)+geom_boxplot()+theme_bw()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none",axis.title = element_text(size = axisS),axis.text = element_text(size = labelS))+ylim(ylim)+labs(x="",y=ylab)
    if(!is.null(plotcolours)){
      plotDS<-plotDS+scale_fill_manual(values=plotcolours)
    }
    savepdf(plotname,width=pwidth,height=pheight)
    print(plotDS)
    dev.off()
  }

}


plotASW<-function(ASWscores,DownRange=seq(0.1,0.9,by=0.05),nsamples=NULL,plotname,ylim=c(0.98,1)){
  xfactors<-paste0("Sample",DownRange)
  if(!is.null(nsamples)){
    xfactors<-paste("N=",round(nsamples*DownRange))
  }
  names(compareCD)<-xfactors
  corvals<-c()
  xfvals<-c()
  for(i in 1:length(xfactors)){
    corvals<-c(corvals,unlist(compareCD[[i]]))
    xfvals<-c(xfvals,rep(xfactors[i],length(unlist(compareCD[[i]]))))
    
  }
  plotdata<-data.frame(x=corvals,value=xfvals)
  plotdata$value<-factor(plotdata$value,levels=xfactors)
  plotDS<-ggplot(aes(x=value,y=x),data=plotdata)+geom_boxplot()+theme_bw()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ylim(ylim)
  pdf(plotname)
  print(plotDS)
  dev.off()
}