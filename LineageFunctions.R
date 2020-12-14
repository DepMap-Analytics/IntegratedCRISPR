lineageCluster<-function(InputData,genes){
  res<-list()
  resj<-list()
  for(i in 1:length(InputData)){
    InData<-InputData[[i]][genes,]
    
    InData<-scale(t(na.omit(InData)))
    for(j in 1:100){
      resj[[j]]<-NbClust(InData, distance = "euclidean", min.nc=2, max.nc=15, 
                 method = "kmeans", index = "silhouette")$All.index
    }
    res[[i]]<-matrix(unlist(resj),ncol=14,byrow=T)
  }
  
  names(res)<-names(InputData)
  return(res)
}

RepeatKmeansAMI<-function(InputData,genes,nrepeat=100,k,annot){
  res<-list()
  kres<-list()
  ami<-list()
  totwith<-list()
  betweens<-list()
  for(i in 1:length(InputData)){
    if(genes=="var"){
      rvs<-rowVars(InputData[[i]],na.rm=T)
      names(rvs)<-rownames(InputData[[i]])
      rvs<-rvs[order(rvs,decreasing=TRUE)]
      genes<-names(rvs)[1:500]
      InData<-InputData[[i]][genes,]
    }else{
    InData<-InputData[[i]][genes,]}
    
    InData<-scale(t(na.omit(InData)))
    rownames(InData)<-colnames(InputData[[i]])
    colnames(InData)<-genes
    incCL<-intersect(rownames(InData),names(annot))
    InData<-InData[incCL,]
    annot<-annot[incCL]
    InData<-InData[!is.na(annot[rownames(InData)]),]
    annot<-annot[!is.na(annot)]
    for(j in 1:nrepeat){
      temp<-kmeans(InData,k)
      kres[[j]]<-temp$cluster
      ami[[j]]<-AMI(temp$cluster,annot[names(temp$cluster)])
      totwith[[j]]<-temp$tot.withinss
      betweens[[j]]<-temp$betweenss
    }
    res[[i]]<-list(ami=ami,totwithss=totwith,betweenss=betweens)
  }
  names(res)<-names(InputData)
  return(res)
}

plotAMI<-function(AMIlist,plotcolours,pname){
  batchnames<-names(AMIlist)
  preProc<-names(AMIlist[[1]])
  AMIData<-NULL
  for(i in 1:length(AMIlist)){
    for(j in 1:3){
      temp<-data.frame(ami=unlist(AMIlist[[i]][[j]]$ami),preProc=rep(preProc[j],100),batch=rep(batchnames[i],100))
      AMIData<-rbind(AMIData,temp)
    }
    
  }

  AMIData$preProc<-factor(AMIData$preProc,levels=c("CCR","CCRJ","CERES"))
  
  AMIData$batch<-factor(AMIData$batch,levels=c("ComBat","ComBatQN","ComBatPC1","ComBatPC2"))
  dodge <- position_dodge(width = 1)
  AMIplot<-ggplot(AMIData,aes(x=preProc,y=ami,fill=batch))+scale_fill_manual(values=plotcolours)+theme(axis.text.x=element_text(angle=45))+ylab("AMI")+xlab("")+geom_boxplot(position = dodge,width=0.5)+theme_bw()+theme(legend.position = "none",legend.background = element_blank(),legend.title=element_blank())
  

  pdf(paste0(dir.Results,"/LineageAMI",pname,".pdf"))
  print(AMIplot)
  dev.off()
}




classPerfLineage<-function(dataset,qualityTH=Inf,QC=NULL,weights=NULL,geneset=NULL,distmetric="Cor",lineagelabel,withRange=1){
  if(!is.null(geneset)){
    dataset<-dataset[intersect(rownames(dataset),geneset),]
  }

  dataset<-na.omit(dataset)

  clnames<-colnames(dataset)
  labnames<-intersect(clnames,names(lineagelabel))
  dataset<-dataset[,labnames]
  lineagelabel<-lineagelabel[labnames]
  if(distmetric=="Cor"){
    cdist<-as.matrix(1-cor(dataset))
  }
  if(distmetric=="Euc"){
    cdist<-as.matrix(dist(t(dataset),method = 'euclidean'))
  }
  if(distmetric=="Jacc"){
    cdist<-as.matrix(dist(t(dataset),method="binary"))
  }
  if(!is.null(weights)){
    weights<-weights[rownames(dataset)]
    weights[weights<0]<-0
    weightmat<-matrix(weights,nrow=nrow(dataset),ncol=ncol(dataset))
    cdist<-as.matrix(1-cor.wt(dataset,w=weights)$r)
  }
  ncls<-ncol(cdist)
  cdist<-cdist[names(lineagelabel),names(lineagelabel)]
  MATCH<-NULL
  
  for(i in 1:ncls){
    
    lineage<-lineagelabel[i]
    neighbours<-colnames(cdist)[setdiff(order(cdist[i,]),i)]
    closestLineage<-which(lineagelabel[neighbours]==lineage)[1]
    mRow<-rep(0,length(neighbours))
    mRow[closestLineage]<-1
    MATCH<-rbind(MATCH,mRow)
    
  }
  
  
  
  MATCH<-MATCH+0
  rownames(MATCH)<-rownames(cdist)
  curv<-cumsum(colSums(MATCH))/ncol(cdist)
  
  return(list(MATCHmat=MATCH,CURV=curv,nb=neighbours))
}


getLRTgenes<-function(normLRTlist,thresh=200){
  glist<-lapply(normLRTlist,function(x) names(x[[3]])[x[[3]]>thresh])
  return(unique(unlist(glist)))
}
