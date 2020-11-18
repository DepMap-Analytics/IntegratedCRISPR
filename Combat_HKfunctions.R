
EssNonEss_plot<-function(DataMatrix,essgenes,nonessgenes,filename,fdr=0.05){
  if(is.list(DataMatrix)){
    DataVectors<-lapply(DataMatrix,function(x) as.vector(unlist(x)))
    nameMats<-lapply(DataMatrix,function(x) {
      temp<-matrix(rownames(x),nrow=nrow(x),ncol=ncol(x))
      colnames(temp)<-colnames(x)
      Nvec_un<-matrix(unlist(temp),ncol=1)
    })
    essNonEss<-list()
    for(i in 1:length(DataMatrix)){
      essNonEss[[i]]<-getEssNonEss(DataVectors[[i]],essgenes,nonessgenes,nameMats[[i]])
    }
    xlim_values<-c(min(unlist(lapply(essNonEss,function(x) min(x$ess)))),max(unlist(lapply(essNonEss,function(x) max(x$noness)))))
    color = grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), 
                                     invert = T)]
   
    COLS <- sample(color, length(DataMatrix))
    pdf(filename)
      plot(density(essNonEss[[1]]$noness),xlim=xlim_values,col=COLS[1],type="l",lty=1,main="")
      lines(density(essNonEss[[1]]$ess),col=COLS[1],type="l",lty=2)
      for(i in 2:length(DataMatrix)){
        lines(density(essNonEss[[i]]$noness),col=COLS[i],type="l",lty=1)
        lines(density(essNonEss[[i]]$ess),col=COLS[i],type="l",lty=2)
      }
      legend("topleft",legend=names(DataMatrix),col = COLS,lty=1)
      legend("left",legend=c("Essentials","Non Essentials"),lwd=1,lty=c(2,1))
    dev.off()
    
  
  }else{
    DataVector<-as.vector(unlist(DataMatrix))
    nameMat<-matrix(rownames(DataMatrix),nrow=nrow(DataMatrix),ncol=ncol(DataMatrix))
    colnames(nameMat)<-colnames(DataMatrix)
    Nvec_un<-matrix(unlist(nameMat),ncol=1)
    EssNess<-getEssNonEss(DataVector,essgenes,nonessgenes,Nvec_un)
    RESBun<-PrRc_Curve(DataVector,EssNess$ess,EssNess$noness,display = FALSE,FDRth = fdr)
    pdf(filename)
      plot(density(EssNess$noness,xlim=c(min(EssNess$ess),max(EssNess$noness))),col="black")
      lines(density(EssNess$ess,col="red"))
      
    dev.off()
  }
}

getEssNonEss<-function(DataVector,essgenes,nonessgenes,Datanames){
  names(DataVector)<-Datanames
  essuse<-intersect(essgenes,names(DataVector))
  nonessuse<-intersect(nonessgenes,names(DataVector))
  return(list(ess=DataVector[essuse],noness=DataVector[nonessuse]))
}

combatadj<-function(Data1,Data2,fdr=0.05,distmetric="Cor",subsetByQC=NULL,geneset=NULL,customBatch=NULL,extraBatch=NULL,std=FALSE,ess=NULL,noness=NULL,stdglobal=FALSE,postStd=FALSE,weights=NULL,qualityTH=Inf,QC=NULL,bycellline=FALSE){
  ingenes<-intersect(rownames(Data1),rownames(Data2))
  combinedData<-cbind(Data1[ingenes,],Data2[ingenes,])
  if(!is.null(subsetByQC)){
    clnames<-strsplit(colnames(combinedData),'---')
    clnames<-unlist(lapply(clnames,function(x){x[1]}))
    
    uc<-unique(clnames)
    for (i in 1:length(uc)){
      QC[which(clnames==uc[i])]<-rep(min(QC[which(clnames==uc[i])]),2)
    }
    
    combinedData<-combinedData[,QC>=qualityTH]
    QC<-NULL
  }
  mergeOnly<-combinedData
  dn<-dimnames(combinedData)
  combinedData<-normalize.quantiles(as.matrix(combinedData))
  dimnames(combinedData)<-dn
  
  if(std){
    if(stdglobal){
      if(!is.null(ess)){
        ess_in<-intersect(ingenes,ess)
        essCD<-combinedData[ess_in,]
        ess_meds<-apply(essCD,2,median)

        stdData<-(combinedData-median(ess_meds)-1)
        
        #stdData<-combinedData-ess_meds-1
      }
    }else{
      if(!is.null(ess)){
        ess_in<-intersect(ingenes,ess)
        essCD<-combinedData[ess_in,]
        ess_meds<-apply(essCD,2,median)
        scaleMat<-matrix(ess_meds,ncol=ncol(combinedData),nrow=nrow(combinedData),byrow=T)
        noness_in<-intersect(noness,ingenes)
        noness_meds<-apply(combinedData[noness_in,],2,median)
        shiftMat<-matrix(noness_meds,ncol=ncol(combinedData),nrow=nrow(combinedData),byrow = T)
        stdData<-(combinedData-shiftMat)/abs(scaleMat)

      #stdData<-combinedData-ess_meds-1
      }
    }
  }else{

    stdData<-combinedData
  }
  site<-strsplit(colnames(combinedData),'---')
  site<-unlist(site)
  site<-site[seq(2,length(site),2)]
  names(site)<-colnames(combinedData)
  if(!is.null(customBatch)){
    corrected<-ComBat(as.matrix(stdData),batch = customBatch)
  }else{
  if(!is.null(extraBatch)){
    #allbatch<-cbind(site[colnames(combinedData)],extraBatch)
    corrected<-ComBat(as.matrix(stdData),mod = extraBatch[colnames(combinedData)],batch=site[colnames(combinedData)])
  }else{
    corrected<-ComBat(as.matrix(stdData),batch = site[colnames(combinedData)])
  }
  if(postStd){
    if(!is.null(ess)){
      ess_in<-intersect(ingenes,ess)
      essCD<-corrected[ess_in,]
      ess_meds<-apply(essCD,2,median)
      shiftMat<-matrix(ess_meds+1,ncol=ncol(combinedData),nrow=nrow(combinedData),byrow = T)
      corrected<-(corrected-shiftMat)
      
      #stdData<-combinedData-ess_meds-1
    }
  }
  }
  correctNq<-corrected
  dn<-dimnames(corrected)
  corrected<-normalize.quantiles(as.matrix(corrected))
  dimnames(corrected)<-dn
  correctV<-ComBat(as.matrix(mergeOnly),batch=site[colnames(mergeOnly)])
  RES_allV<-classPerfCP(correctV,qualityTH=qualityTH,QC=QC,geneset=geneset,distmetric=distmetric,bycellline=bycellline)
  if(!is.null(weights)){
    
    RES_all<-classPerfCP(corrected,weights=weights,qualityTH=qualityTH,QC=QC,geneset=geneset,bycellline=bycellline)
    RES_allNq<-classPerfCP(correctNq,weights=weights,qualityTH=qualityTH,QC=QC,geneset=geneset,bycellline=bycellline)
  }else{
    weights=NULL
    if(distmetric=="Jacc"){
      #need to generate the per cell line binary significance matrix.
      #binarymatrix<-binarysignif(corrected,ess=BAGEL_essential,noness=BAGEL_nonEssential)
      binarymatrix<-corrected
      for(i in 1:ncol(corrected)){
        FCs<-corrected[,i]
        names(FCs)<-rownames(corrected)
        RES<-PrRc_Curve(FCs,ess,noness,display = FALSE,FDRth = fdr)
        
        binarymatrix[,i]<-(FCs<=RES$sigthreshold)*1
      }
      RES_all<-classPerfCP(binarymatrix,distmetric = distmetric,geneset=geneset,bycellline=bycellline)
    }else{
      RES_all<-classPerfCP(corrected,qualityTH=qualityTH,QC=QC,geneset=geneset,distmetric=distmetric,bycellline=bycellline)
      RES_allNq<-classPerfCP(correctNq,qualityTH=qualityTH,QC=QC,geneset=geneset,distmetric=distmetric,bycellline=bycellline)}
  }

  
  return(list(mergedData=stdData,corrected=corrected,RES_all=RES_all,mergeOnly=mergeOnly,RES_allNq=RES_allNq,correctNq=correctNq,correctV=correctV,RES_allV=RES_allV  ))
}

binarysignif<-function(corrected,ess=BAGEL_essential,noness=BAGEL_nonEssential,fdr=0.05){
  noness<-intersect(noness,rownames(corrected))
  nullvalues<-as.vector(corrected[noness,])
  meanvec<-mean(nullvalues)
  sdvec<-sd(nullvalues)

  pvalmat<-apply(corrected,2,function(x) pnorm(x,mean=meanvec,sd=sdvec,lower.tail=FALSE))
  padjmat<-apply(pvalmat,2,function(x) p.adjust(x,method="fdr"))
  binarymat<-(padjmat<fdr)*1
}

classPerfCP<-function(dataset,qualityTH=Inf,QC=NULL,weights=NULL,geneset=NULL,distmetric="Cor",pairmat=NULL,withRange=1,bycellline=FALSE){
  if(!is.null(geneset)){
    dataset<-dataset[intersect(rownames(dataset),geneset),]
  }
  
  if(length(QC)>0){
    clnames<-strsplit(colnames(dataset),'---')
    clnames<-unlist(lapply(clnames,function(x){x[1]}))
    
    uc<-unique(clnames)
    for (i in 1:length(uc)){
      QC[which(clnames==uc[i])]<-rep(min(QC[which(clnames==uc[i])]),2)
    }
    
    dataset<-dataset[,QC>=qualityTH]
  }
  
  clnames<-strsplit(colnames(dataset),'---')
  clnames<-unlist(lapply(clnames,function(x){x[1]}))
  
  site<-strsplit(colnames(dataset),'---')
  site<-unlist(lapply(site,function(x){x[[2]]}))
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
    weightmat<-matrix(weights,nrow=nrow(dataset),ncol=ncol(dataset))
    cdist<-as.matrix(1-cor.wt(dataset,w=weights)$r)
  }
  ncls<-ncol(cdist)
  pairwiseCor<-sapply(unique(clnames),function(x) cdist[clnames==x,clnames==x])
  pairwiseCor<-apply(pairwiseCor,2,max)
  MATCH<-NULL
  MATCHjoint<-NULL
  curvjoint<-NULL
  MATCHrank<-NULL
  CURVrank<-NULL
  if(!is.null(pairmat)){
    for(i in 1:ncls){
      pairmatch<-colnames(pairmat)[setdiff(order(pairmat[i,]),i)][withRange]
      neighbours<-colnames(cdist)[setdiff(order(cdist[i,]),i)]
      matchingfeatures<-which(neighbours%in%pairmatch)[1]
      mRow<-rep(0,length(neighbours))
      mRow[matchingfeatures]<-1
      MATCH<-rbind(MATCH,mRow)
      
    }
  }else{
    piestatss<-list()
    for (i in 1:ncls){
      neighbours<-clnames[setdiff(order(cdist[i,]),i)]
      
      MATCH<-rbind(MATCH,neighbours==clnames[i])
      if(neighbours[1]==clnames[i]){
      piestatss[[i]]<-"match"}else{
        closest<-colnames(cdist)[setdiff(order(cdist[i,]),i)]
        siteclosest<-strsplit(closest[1], "---")[[1]][2]
        site_i<-strsplit(colnames(cdist)[i],"---")[[1]][2]
        if(siteclosest==site_i){
          piestatss[[i]]<-"Batch"        }else{
            piestatss[[i]]<-"lineage"
          }
      }
    }
    names(piestatss)<-rownames(cdist)
    rankmatrix<-apply(cdist,1,function(x) rank(x))
    dimnames(rankmatrix)<-dimnames(cdist)
    pmean <- function(x,y) (x+y)/2
    symRank<- pmean(rankmatrix, matrix(rankmatrix, nrow(rankmatrix), byrow=TRUE))
    dimnames(symRank)<-dimnames(cdist)
    RankDist<-as.matrix(dist(symRank))
    dimnames(RankDist)<-dimnames(cdist)
    MATCHrd<-NULL
    piestats<-list()
    for (i in 1:ncls){
      neighbours<-clnames[setdiff(order(symRank[i,]),i)]
      neighbours2<-clnames[setdiff(order(RankDist[i,]),i)]
      MATCHrank<-rbind(MATCHrank,neighbours==clnames[i])
     
      if(neighbours[1]==clnames[i]){
        piestats[[i]]<-"match"}else{
          closest<-colnames(symRank)[setdiff(order(symRank[i,]),i)]
          siteclosest<-strsplit(closest[1], "---")[[1]][2]
          site_i<-strsplit(colnames(symRank)[i],"---")[[1]][2]
          if(siteclosest==site_i){
            piestats[[i]]<-"Batch"        }else{
              piestats[[i]]<-"lineage"
            }
        }
    }
    names(piestats)<-rownames(symRank)
  }
  piestatsJ<-list()
  if(bycellline){
    ucl<-unique(clnames)

    for(j in 1:length(ucl)){
 
      idx<-which(clnames==ucl[j])
      neighbours1<-clnames[setdiff(order(cdist[idx[1],]),idx[1])]
      neighbours2<-clnames[setdiff(order(cdist[idx[2],]),idx[2])]
      n1<-neighbours1==ucl[j]
      n2<-neighbours2==ucl[j]
      newvec1<-rep(0,ncol(cdist)-1)
      newvec1[min(which(n1==1),which(n2==1))]<-1
      
      MATCHjoint<-rbind(MATCHjoint,newvec1)
     
      if(newvec1[1]==1){
        piestatsJ[[j]]<-"Match"}else{
          closest1<-colnames(cdist)[setdiff(order(cdist[idx[1],]),idx[1])]
          closest2<-colnames(cdist)[setdiff(order(cdist[idx[2],]),idx[2])]
          siteclosest1<-strsplit(closest1[1], "---")[[1]][2]
          siteclosest2<-strsplit(closest2[1], "---")[[1]][2]
          site_i1<-strsplit(colnames(cdist)[idx[1]],"---")[[1]][2]
          site_i2<-strsplit(colnames(cdist)[idx[2]],"---")[[1]][2]
          #need to know which is min idx1 or idx2:
          minofpair<-which.min(c(which(n1==1),which(n2==1)))
          if(minofpair==1){
            if(siteclosest1==site_i1){
              piestatsJ[[j]]<-"Batch"
            }else{
              piestatsJ[[j]]<-"Lineage"
            }
          }else{
            if(siteclosest2==site_i2){
              piestatsJ[[j]]<-"Batch"
            }else{
              piestatsJ[[j]]<-"Lineage"
            }
            
          }
        }
    }
    rownames(MATCHjoint)<-ucl
    names(piestatsJ)<-ucl
    curvjoint<-cumsum(colSums(MATCHjoint))*2/ncol(cdist)
  }
  
  MATCH<-MATCH+0
  rownames(MATCH)<-rownames(cdist)
  curv<-cumsum(colSums(MATCH))/ncol(cdist)
  MATCHrank<-MATCHrank+0
  rownames(MATCHrank)<-rownames(cdist)
  CURVrank<-cumsum(colSums(MATCHrank))/ncol(cdist)
  
  
  return(list(pairwiseCor=pairwiseCor,MATCHmat=MATCH,CURV=curv,MATCHjoint=MATCHjoint,CURVjoint=curvjoint,piestatss=piestatss,piestats=piestats,piestatsJ=piestatsJ,CURVrank=CURVrank,MATCHrank=MATCHrank,cdist=cdist))
}
classPerfAll<-function(dataset,qualityTH=Inf,QC=NULL,weights=NULL,geneset=NULL,distmetric="Cor",pairmat,withRange=1){
  if(!is.null(geneset)){
    dataset<-dataset[intersect(rownames(dataset),geneset),]
  }
  
clnames<-colnames(dataset)

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
    weightmat<-matrix(weights,nrow=nrow(dataset),ncol=ncol(dataset))
    cdist<-as.matrix(1-cor.wt(dataset,w=weights)$r)
  }
  ncls<-ncol(cdist)
  cdist<-cdist[colnames(pairmat),colnames(pairmat)]
  MATCH<-NULL

    for(i in 1:ncls){
      pairmatch<-colnames(pairmat)[setdiff(order(pairmat[i,]),i)][withRange]
      neighbours<-colnames(cdist)[setdiff(order(cdist[i,]),i)]
      matchingfeatures<-which(neighbours%in%pairmatch)[1]
      mRow<-rep(0,length(neighbours))
      mRow[matchingfeatures]<-1
      MATCH<-rbind(MATCH,mRow)
      
    }

  
  
  MATCH<-MATCH+0
  rownames(MATCH)<-rownames(cdist)
  curv<-cumsum(colSums(MATCH))/ncol(cdist)
  
  return(list(MATCHmat=MATCH,CURV=curv,pm=pairmatch,nb=neighbours))
}
pcaplotAll<-function(combatoutput,filename,colours=NULL,width=10,sitedata,tissue=tissue){
  
  if(!is.null(colours)){metai<-data.frame(name=colnames(combatoutput),batch=sitedata,tissue=tissue)}else{
  metai<-data.frame(name=colnames(combatoutput),batch=sitedata)}
  pcs<-prcomp(t(combatoutput),scale. = TRUE)

  eigs<-(pcs$sdev)^2
  totaleigs<-sum(eigs)
  prop1<-eigs[1]/totaleigs
  prop2<-eigs[2]/totaleigs
  hplot<-(prop2/prop1)*width
  savepdf(paste0(filename,"_combatPCA"),width=width,height=hplot)
  if(is.null(colours)){
  print(autoplot(prcomp(t(combatoutput)),data=metai,colour='batch',alpha=0.5))}else{
    print(autoplot(prcomp(t(combatoutput),scale. = TRUE),data=metai,alpha=0.5,colour='tissue',shape='batch')+scale_colour_manual(values=colours)+ theme(legend.position="none"))
  }
  dev.off()
}

pcaplot<-function(combatoutput,filename,colours=NULL,width=10,scalefig=TRUE){
  sitedata<-sapply(colnames(combatoutput),function(x) strsplit(x,"---",fixed=T)[[1]][2])
  ds<-dimnames(combatoutput)
  combatoutput<-normalize.quantiles(as.matrix(combatoutput))
  
  dimnames(combatoutput)<-ds
  if(!is.null(colours)){metai<-data.frame(name=colnames(combatoutput),batch=sitedata,cfill=colours[sitedata])}else{
    metai<-data.frame(name=colnames(combatoutput),batch=sitedata)}
  pcs<-prcomp(t(combatoutput),scale. = TRUE)
  
  eigs<-(pcs$sdev)^2
  totaleigs<-sum(eigs)
  prop1<-eigs[1]/totaleigs
  prop2<-eigs[2]/totaleigs
  hplot<-(prop2/prop1)*width
  if(scalefig){
  savepdf(paste0(filename,"_combatPCA"),width=width,height=hplot)
  if(is.null(colours)){
    print(autoplot(prcomp(t(combatoutput)),data=metai,colour='batch',alpha=0.5))}else{
      print(autoplot(prcomp(t(combatoutput),scale. = TRUE),data=metai,alpha=0.5,colour='batch')+scale_colour_manual(values=colours))
    }
  dev.off()
  }else{
    savepdf(paste0(filename,"_combatPCA_square"),width=width)
    if(is.null(colours)){
      print(autoplot(prcomp(t(combatoutput)),data=metai,colour='batch',alpha=0.5))}else{
        print(autoplot(prcomp(t(combatoutput),scale. = TRUE),data=metai,alpha=0.5,colour='batch')+scale_colour_manual(values=colours))
      }
    dev.off()
  }
}
tsneplot_cp<-function(combatoutput,sitedata,filename,colourgroups,weights=NULL,COLS=NULL,perplexity=30,init_dims=30){

  ds<-dimnames(combatoutput)
  combatoutput<-normalize.quantiles(as.matrix(combatoutput))
  
  dimnames(combatoutput)<-ds
  #metai<-data.frame(name=colnames(combatoutput),batch=sitedata)
  #savepdf(paste0(filename,"_combatPCA"),width=10)
  #print(autoplot(prcomp(t(combatoutput)),data=metai,colour='batch',alpha=0.5))
  #dev.off()
  
  nit<-1000
  if(!is.null(weights)){
    weights<-weights[rownames(combatoutput)]
    weights[weights<0]<-0
    cdists<-as.matrix(1-cor.wt(combatoutput,w=weights)$r)
  }else{
  cdists <- as.dist(1 - cor(combatoutput))}
  set.seed(679661)
  
  
  tfits <- tsne(cdists, max_iter = nit,perplexity = perplexity,initial_dims = init_dims)
  
  SYMBOLS <- rep(21, length(sitedata))
  SYMBOLS[sitedata == "Broad"] <- 23
  names(SYMBOLS) <- sitedata
  CEX <- rep(1.3, length(sitedata))
  CEX[sitedata == "Broad"] <- 1.8
  uc <- unique(colourgroups)
  n <- length(uc)
  if(is.null(COLS)){
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    COLS = sample(col_vector,n)
  }
  
  names(COLS) <- uc
  COLS["BroadSpecific"]<-NA
  tCOLS <- makeTransparent(COLS, alpha = 140)
 
  
  colnames(tfits)<-c("posx","posy")
  pdf(paste0(filename,"_combatTSNE.pdf"))
  plot(tfits[, 1], tfits[, 2], col = "black", 
       cex = CEX, pch = SYMBOLS[sitedata], bg = tCOLS[colourgroups],frame.plot = FALSE, xaxt = "n", 
       yaxt = "n", xlab = "tSNE AU", ylab = "tSNE AU", main = filename)
  dev.off()
}
splineadjustment<-function(data,quantiles=c(0.25,0.5,0.75)){
  dm<-dimnames(data)
  dataq<-normalize.quantiles(data)
  dimnames(dataq)<-dm
  GeneMeans<-rowMeans(dataq,na.rm = TRUE)
  knots <- quantile(Means, p = quantiles)
  for(i in 1:ncol(dataq)){
    model <- lm (dataq[,i] ~ bs(GeneMeans, knots = knots),na.action=na.exclude)
    #model2<-smooth.spline(SangerMeans,SangerOverlap[,i],cv=TRUE)
    if(i==1){
      DataScreen<-dataq[,i]-fitted.values(model)+GeneMeans
    }else{
      DataScreen<-cbind(DataScreen,dataq[,i]-fitted.values(model)+GeneMeans)
    }
    
  }
  dimnames(DataScreen)<-dm
  return(DataScreen)
}
nnmd<-function(FCs,essential,nonessential){
  ess<-intersect(rownames(FCs),essential)
  noness<-intersect(rownames(FCs),nonessential)
  numerators<-apply(FCs,2,function(x) mean(x[ess])-mean(x[noness]))
  denom<-apply(FCs,2,function(x) sd(x[noness]))
  unlist(numerators)/unlist(denom)
}
ssmd <- function (a, b, verbose=TRUE,...) 
{
  if (length(a) < 2 | length(b) < 2) {
    stop(call. = FALSE, "Inputs need to be greater at least 2 elements long")
  }
  if (is.numeric(a) == FALSE | is.numeric(b) == FALSE) {
    stop(call. = FALSE, "Input needs to be numeric.")
  }
  
  mu_a <- mean(a, ...)
  mu_b <- mean(b, ...)
  var_a <- var(a, ...)
  var_b <- var(b, ...)
  
  # if lengths are equal assume correlation and calculate covariance
  if (length(a) == length(b)) {
    cov_ab <- cov(a, b)
    beta <- (mu_a - mu_b)/sqrt(var_a + var_b - 2 * cov_ab)
  } else{ # unequal lengths => not paired , cannot calc covariance
    beta <- (mu_a - mu_b) / sqrt(var_a + var_b)
    if(verbose == TRUE)
    {
      warning("a and b have different lengths. Calculations assumed no correlation.",
              call. = FALSE)
    }
  }
  

  
  return(beta)
}

CL_qualitySSMD<-function(FCmatrix,ess,noness){
  out<-vector("numeric",length=ncol(FCmatrix))
  for(i in 1:ncol(FCmatrix)){
    input<-FCmatrix[,i]
    names(input)<-rownames(FCmatrix)
    out[i]<-quality_SSMD(input,ess,noness)
  }
  return(out)
}
quality_SSMD<-function(fcs,ess,noness){
  essfc<-fcs[intersect(ess,names(fcs))]
  nonessfc<-fcs[intersect(noness,names(fcs))]
  ssmd(essfc,nonessfc)
  #mu_a<-mean(essfc)
  #mu_b<-mean(non)
  #beta <- (mu_a - mu_b) / sqrt(var_a + var_b)
}

quality_D2<-function(fcMat,ess,noness){
  med_pos<-apply(fcMat[intersect(ess,rownames(fcMat)),],1,mean)
  med_neg<-apply(fcMat[intersect(noness,rownames(fcMat)),],1,mean)
  diff_pos<-fcMat[intersect(ess,rownames(fcMat)),]-med_pos
  diff_neg<-fcMat[intersect(noness,rownames(fcMat)),]-med_neg
  qjs<-apply(diff_neg,2,median)-apply(diff_pos,2,median)
  qjs<-qjs-mean(qjs)+1
  return(qjs)
}

ComBatCP <- function (dat, batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE,
                      mean.only = FALSE, ref.batch = NULL, BPPARAM = bpparam("SerialParam"),empBayes=TRUE) {
  ## make batch a factor and make a set of indicators for batch
  if(mean.only==TRUE){
    message("Using the 'mean only' version of ComBat")
  }
  if(length(dim(batch))>1){
    stop("This version of ComBat only allows one batch variable")
  }  ## to be updated soon!
  batch <- as.factor(batch)
  batchmod <- model.matrix(~-1+batch)  
  if (!is.null(ref.batch)){
    ## check for reference batch, check value, and make appropriate changes
    if (!(ref.batch%in%levels(batch))) {
      stop("reference level ref.batch is not one of the levels of the batch variable")
    }
    cat("Using batch =",ref.batch, "as a reference batch (this batch won't change)\n")
    ref <- which(levels(as.factor(batch))==ref.batch) # find the reference
    batchmod[,ref] <- 1
  } else {
    ref <- NULL
  }
  message("Found", nlevels(batch), "batches")
  
  ## A few other characteristics on the batches
  n.batch <- nlevels(batch)
  batches <- list()
  for (i in 1:n.batch) {
    batches[[i]] <- which(batch == levels(batch)[i])
  } # list of samples in each batch  
  n.batches <- sapply(batches, length)
  if(any(n.batches==1)){
    mean.only=TRUE
    message("Note: one batch has only one sample, setting mean.only=TRUE")
  }
  n.array <- sum(n.batches)
  ## combine batch variable and covariates
  design <- cbind(batchmod,mod)
  
  ## check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  if(!is.null(ref)){
    check[ref] <- FALSE
  } ## except don't throw away the reference batch indicator
  design <- as.matrix(design[,!check])
  
  ## Number of covariates or covariate levels
  message("Adjusting for", ncol(design)-ncol(batchmod), 'covariate(s) or covariate level(s)')
  
  ## Check if the design is confounded
  if(qr(design)$rank < ncol(design)) {
    ## if(ncol(design)<=(n.batch)){stop("Batch variables are redundant! Remove one or more of the batch variables so they are no longer confounded")}
    if(ncol(design)==(n.batch+1)) {
      stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")
    }
    if(ncol(design)>(n.batch+1)) {
      if((qr(design[,-c(1:n.batch)])$rank<ncol(design[,-c(1:n.batch)]))){
        stop('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded')
      } else {
        stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")
      }
    }
  }
  
  ## Check for missing values
  NAs <- any(is.na(dat))
  if(NAs){
    message(c('Found',sum(is.na(dat)),'Missing Data Values'), sep=' ')}
  ## print(dat[1:2,])
  
  ##Standardize Data across genes
  cat('Standardizing Data across genes\n')
  if (!NAs){
    B.hat <- solve(crossprod(design), tcrossprod(t(design), as.matrix(dat)))
  } else { 
    B.hat <- apply(dat, 1, Beta.NA, design) # FIXME
  }
  
  ## change grand.mean for ref batch
  if(!is.null(ref.batch)){
    grand.mean <- t(B.hat[ref, ])
  } else {
    grand.mean <- crossprod(n.batches/n.array, B.hat[1:n.batch,])
  }
  
  ## change var.pooled for ref batch
  if (!NAs){
    if(!is.null(ref.batch)) {
      ref.dat <- dat[, batches[[ref]]]
      var.pooled <- ((ref.dat-t(design[batches[[ref]], ] %*% B.hat))^2) %*% rep(1/n.batches[ref],n.batches[ref]) # FIXME
    } else {
      var.pooled <- ((dat-t(design %*% B.hat))^2) %*% rep(1/n.array,n.array) # FIXME
    }
  } else {
    if(!is.null(ref.batch)) {
      ref.dat <- dat[, batches[[ref]]]
      var.pooled <- rowVars(ref.dat-t(design[batches[[ref]], ]%*%B.hat), na.rm=TRUE)
    } else {
      var.pooled <- rowVars(dat-t(design %*% B.hat), na.rm=TRUE)
    }
  }
  
  stand.mean <- t(grand.mean) %*% t(rep(1,n.array)) # FIXME
  if(!is.null(design)){
    tmp <- design
    tmp[,c(1:n.batch)] <- 0
    stand.mean <- stand.mean+t(tmp %*% B.hat) #FIXME
  }  
  s.data <- (dat-stand.mean)/(sqrt(var.pooled) %*% t(rep(1,n.array))) # FIXME
  
  ##Get regression batch effect parameters
  message("Fitting L/S model and finding priors")
  batch.design <- design[, 1:n.batch]
  if (!NAs){
    gamma.hat <- solve(crossprod(batch.design), tcrossprod(t(batch.design),
                                                           as.matrix(s.data)))
  } else{
    gamma.hat <- apply(s.data, 1, Beta.NA, batch.design) # FIXME
  }
  delta.hat <- NULL
  for (i in batches){
    if(mean.only==TRUE) {
      delta.hat <- rbind(delta.hat,rep(1,nrow(s.data))) 
    } else {
      delta.hat <- rbind(delta.hat, rowVars(s.data[,i], na.rm=TRUE))
    }
  }
  
  if(empBayes){
  ##Find Priors
  gamma.bar <- rowMeans(gamma.hat)
  t2 <- rowVars(gamma.hat)
  a.prior <- apply(delta.hat, 1, sva:::aprior) # FIXME 
  b.prior <- apply(delta.hat, 1, sva:::bprior) # FIXME
  
  ## Plot empirical and parametric priors
  
  if (prior.plots && par.prior) {
    par(mfrow=c(2,2))
    
    ## Top left
    tmp <- density(gamma.hat[1,])
    plot(tmp,  type='l', main=expression(paste("Density Plot of First Batch ",  hat(gamma))))
    xx <- seq(min(tmp$x), max(tmp$x), length=100)
    lines(xx,dnorm(xx,gamma.bar[1],sqrt(t2[1])), col=2)
    
    ## Top Right
    qqnorm(gamma.hat[1,], main=expression(paste("Normal Q-Q Plot of First Batch ", hat(gamma))))
    qqline(gamma.hat[1,], col=2)
    
    ## Bottom Left
    tmp <- density(delta.hat[1,])
    xx <- seq(min(tmp$x), max(tmp$x), length=100)
    tmp1 <- list(x=xx, y=dinvgamma(xx, a.prior[1], b.prior[1]))
    plot(tmp, typ="l", ylim=c(0, max(tmp$y, tmp1$y)),
         main=expression(paste("Density Plot of First Batch ", hat(delta))))
    lines(tmp1, col=2)
    
    ## Bottom Right
    invgam <- 1/qgamma(1-ppoints(ncol(delta.hat)), a.prior[1], b.prior[1])
    qqplot(invgam, delta.hat[1,],
           main=expression(paste("Inverse Gamma Q-Q Plot of First Batch ", hat(delta))),
           ylab="Sample Quantiles", xlab="Theoretical Quantiles")
    lines(c(0, max(invgam)), c(0, max(invgam)), col=2)
  }
  
  ## Find EB batch adjustments
  
  gamma.star <- delta.star <- matrix(NA, nrow=n.batch, ncol=nrow(s.data))
  if (par.prior) {
    message("Finding parametric adjustments")
    results <- bplapply(1:n.batch, function(i) {
      if (mean.only) {
        gamma.star <- postmean(gamma.hat[i,], gamma.bar[i], 1, 1, t2[i])
        delta.star <- rep(1, nrow(s.data))
      }
      else {
        temp <- sva:::it.sol(s.data[, batches[[i]]], gamma.hat[i, ],
                             delta.hat[i, ], gamma.bar[i], t2[i], a.prior[i],
                             b.prior[i])
        gamma.star <- temp[1, ]
        delta.star <- temp[2, ]
      }
      list(gamma.star=gamma.star, delta.star=delta.star)
    }, BPPARAM = BPPARAM)
    for (i in 1:n.batch) {
      gamma.star[i,] <- results[[i]]$gamma.star
      delta.star[i,] <- results[[i]]$delta.star
    }
  }
  else {
    message("Finding nonparametric adjustments")
    results <- bplapply(1:n.batch, function(i) {
      if (mean.only) {
        delta.hat[i, ] = 1
      }
      temp <- int.eprior(as.matrix(s.data[, batches[[i]]]),
                         gamma.hat[i, ], delta.hat[i, ])
      list(gamma.star=temp[1,], delta.star=temp[2,])
    }, BPPARAM = BPPARAM)
    for (i in 1:n.batch) {
      gamma.star[i,] <- results[[i]]$gamma.star
      delta.star[i,] <- results[[i]]$delta.star
    }
  }
  }else{
    #no empirical bayes adjustment:
    gamma.star<-gamma.hat
    delta.star<-delta.hat
  }
  if(!is.null(ref.batch)){
    gamma.star[ref,] <- 0  ## set reference batch mean equal to 0
    delta.star[ref,] <- 1  ## set reference batch variance equal to 1
  }
  
  ## Normalize the Data ###
  message("Adjusting the Data\n")
  
  bayesdata <- s.data
  j <- 1
  for (i in batches){
    bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j]))) # FIXME
    j <- j+1
  }
  
  bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean # FIXME
  
  ## tiny change still exist when tested on bladder data
  ## total sum of change within each batch around 1e-15 
  ## (could be computational system error).  
  ## Do not change ref batch at all in reference version
  if(!is.null(ref.batch)){
    bayesdata[, batches[[ref]]] <- dat[, batches[[ref]]]
  }
  
  return(list(correctedData=bayesdata,batchDesign=batch.design,gamma.star=gamma.star,delta.star=delta.star,varpool=var.pooled,stdmean=stand.mean))
}

ComBatCP_mixed <- function (dat, batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE,
                      mean.only = FALSE, ref.batch = NULL, BPPARAM = bpparam("SerialParam"),Z,empBayes=TRUE) {
  ## make batch a factor and make a set of indicators for batch
  if(mean.only==TRUE){
    message("Using the 'mean only' version of ComBat")
  }
  if(length(dim(batch))>1){
    stop("This version of ComBat only allows one batch variable")
  }  ## to be updated soon!
  batch <- as.factor(batch)
  batchmod <- model.matrix(~-1+batch)  
  if (!is.null(ref.batch)){
    ## check for reference batch, check value, and make appropriate changes
    if (!(ref.batch%in%levels(batch))) {
      stop("reference level ref.batch is not one of the levels of the batch variable")
    }
    cat("Using batch =",ref.batch, "as a reference batch (this batch won't change)\n")
    ref <- which(levels(as.factor(batch))==ref.batch) # find the reference
    batchmod[,ref] <- 1
  } else {
    ref <- NULL
  }
  message("Found", nlevels(batch), "batches")
  
  ## A few other characteristics on the batches
  n.batch <- nlevels(batch)
  batches <- list()
  for (i in 1:n.batch) {
    batches[[i]] <- which(batch == levels(batch)[i])
  } # list of samples in each batch  
  n.batches <- sapply(batches, length)
  if(any(n.batches==1)){
    mean.only=TRUE
    message("Note: one batch has only one sample, setting mean.only=TRUE")
  }
  n.array <- sum(n.batches)
  ## combine batch variable and covariates
  design <- cbind(batchmod,mod)
  
  ## check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  if(!is.null(ref)){
    check[ref] <- FALSE
  } ## except don't throw away the reference batch indicator
  design <- as.matrix(design[,!check])
  
  ## Number of covariates or covariate levels
  message("Adjusting for", ncol(design)-ncol(batchmod), 'covariate(s) or covariate level(s)')
  
  ## Check if the design is confounded
  if(qr(design)$rank < ncol(design)) {
    ## if(ncol(design)<=(n.batch)){stop("Batch variables are redundant! Remove one or more of the batch variables so they are no longer confounded")}
    if(ncol(design)==(n.batch+1)) {
      stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")
    }
    if(ncol(design)>(n.batch+1)) {
      if((qr(design[,-c(1:n.batch)])$rank<ncol(design[,-c(1:n.batch)]))){
        stop('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded')
      } else {
        stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")
      }
    }
  }
  
  ## Check for missing values
  NAs <- any(is.na(dat))
  if(NAs){
    message(c('Found',sum(is.na(dat)),'Missing Data Values'), sep=' ')}
  ## print(dat[1:2,])
  
  ##Standardize Data across genes
  cat('Standardizing Data across genes\n')
  if (!NAs){
    B.hat <- solve(crossprod(design), tcrossprod(t(design), as.matrix(dat)))
  } else { 
    B.hat <- apply(dat, 1, Beta.NA, design) # FIXME
  }
  
  ## change grand.mean for ref batch
  if(!is.null(ref.batch)){
    grand.mean <- t(B.hat[ref, ])
  } else {
    grand.mean <- crossprod(n.batches/n.array, B.hat[1:n.batch,])
  }
  
  ## change var.pooled for ref batch
  if (!NAs){
    if(!is.null(ref.batch)) {
      ref.dat <- dat[, batches[[ref]]]
      var.pooled <- ((ref.dat-t(design[batches[[ref]], ] %*% B.hat))^2) %*% rep(1/n.batches[ref],n.batches[ref]) # FIXME
    } else {
      var.pooled <- ((dat-t(design %*% B.hat))^2) %*% rep(1/n.array,n.array) # FIXME
    }
  } else {
    if(!is.null(ref.batch)) {
      ref.dat <- dat[, batches[[ref]]]
      var.pooled <- rowVars(ref.dat-t(design[batches[[ref]], ]%*%B.hat), na.rm=TRUE)
    } else {
      var.pooled <- rowVars(dat-t(design %*% B.hat), na.rm=TRUE)
    }
  }
  
  stand.mean <- t(grand.mean) %*% t(rep(1,n.array)) # FIXME
  if(!is.null(design)){
    tmp <- design
    tmp[,c(1:n.batch)] <- 0
    stand.mean <- stand.mean+t(tmp %*% B.hat) #FIXME
  }  
  s.data <- (dat-stand.mean)/(sqrt(var.pooled) %*% t(rep(1,n.array))) # FIXME
  
  ##Get regression batch effect parameters
  message("Fitting L/S model and finding priors")
  batch.design <- design[, 1:n.batch]
  if (!NAs){
    gamma.hat <- solve(crossprod(batch.design), tcrossprod(t(batch.design),
                                                           as.matrix(s.data)))
  } else{
    gamma.hat <- apply(s.data, 1, Beta.NA, batch.design) # FIXME
  }
  delta.hat <- NULL
  for (i in batches){
    if(mean.only==TRUE) {
      delta.hat <- rbind(delta.hat,rep(1,nrow(s.data))) 
    } else {
      delta.hat <- rbind(delta.hat, rowVars(s.data[,i], na.rm=TRUE))
    }
  }
  if(empBayes){
  ##Find Priors
  gamma.bar <- rowMeans(gamma.hat)
  t2 <- rowVars(gamma.hat)
  a.prior <- apply(delta.hat, 1, sva:::aprior) # FIXME 
  b.prior <- apply(delta.hat, 1, sva:::bprior) # FIXME
  
  ## Plot empirical and parametric priors
  
  if (prior.plots && par.prior) {
    par(mfrow=c(2,2))
    
    ## Top left
    tmp <- density(gamma.hat[1,])
    plot(tmp,  type='l', main=expression(paste("Density Plot of First Batch ",  hat(gamma))))
    xx <- seq(min(tmp$x), max(tmp$x), length=100)
    lines(xx,dnorm(xx,gamma.bar[1],sqrt(t2[1])), col=2)
    
    ## Top Right
    qqnorm(gamma.hat[1,], main=expression(paste("Normal Q-Q Plot of First Batch ", hat(gamma))))
    qqline(gamma.hat[1,], col=2)
    
    ## Bottom Left
    tmp <- density(delta.hat[1,])
    xx <- seq(min(tmp$x), max(tmp$x), length=100)
    tmp1 <- list(x=xx, y=dinvgamma(xx, a.prior[1], b.prior[1]))
    plot(tmp, typ="l", ylim=c(0, max(tmp$y, tmp1$y)),
         main=expression(paste("Density Plot of First Batch ", hat(delta))))
    lines(tmp1, col=2)
    
    ## Bottom Right
    invgam <- 1/qgamma(1-ppoints(ncol(delta.hat)), a.prior[1], b.prior[1])
    qqplot(invgam, delta.hat[1,],
           main=expression(paste("Inverse Gamma Q-Q Plot of First Batch ", hat(delta))),
           ylab="Sample Quantiles", xlab="Theoretical Quantiles")
    lines(c(0, max(invgam)), c(0, max(invgam)), col=2)
  }
  
  ## Find EB batch adjustments
  
  gamma.star <- delta.star <- matrix(NA, nrow=n.batch, ncol=nrow(s.data))
  if (par.prior) {
    message("Finding parametric adjustments")
    results <- bplapply(1:n.batch, function(i) {
      if (mean.only) {
        gamma.star <- postmean(gamma.hat[i,], gamma.bar[i], 1, 1, t2[i])
        delta.star <- rep(1, nrow(s.data))
      }
      else {
        temp <- sva:::it.sol(s.data[, batches[[i]]], gamma.hat[i, ],
                             delta.hat[i, ], gamma.bar[i], t2[i], a.prior[i],
                             b.prior[i])
        gamma.star <- temp[1, ]
        delta.star <- temp[2, ]
      }
      list(gamma.star=gamma.star, delta.star=delta.star)
    }, BPPARAM = BPPARAM)
    for (i in 1:n.batch) {
      gamma.star[i,] <- results[[i]]$gamma.star
      delta.star[i,] <- results[[i]]$delta.star
    }
  }
  else {
    message("Finding nonparametric adjustments")
    results <- bplapply(1:n.batch, function(i) {
      if (mean.only) {
        delta.hat[i, ] = 1
      }
      temp <- int.eprior(as.matrix(s.data[, batches[[i]]]),
                         gamma.hat[i, ], delta.hat[i, ])
      list(gamma.star=temp[1,], delta.star=temp[2,])
    }, BPPARAM = BPPARAM)
    for (i in 1:n.batch) {
      gamma.star[i,] <- results[[i]]$gamma.star
      delta.star[i,] <- results[[i]]$delta.star
    }
  }
  }else{
    #dont do the empirical bayes adjustment:
    gamma.star<-gamma.hat
    delta.star<-delta.hat
  }
  if(!is.null(ref.batch)){
    gamma.star[ref,] <- 0  ## set reference batch mean equal to 0
    delta.star[ref,] <- 1  ## set reference batch variance equal to 1
  }
  
  ## Normalize the Data ###
  message("Adjusting the Data\n")
  
  bayesdata <- s.data
  j <- 1
  for (i in batches){
    bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j]))) # FIXME
    j <- j+1
  }
  
  bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean # FIXME
  
  ## tiny change still exist when tested on bladder data
  ## total sum of change within each batch around 1e-15 
  ## (could be computational system error).  
  ## Do not change ref batch at all in reference version
  if(!is.null(ref.batch)){
    bayesdata[, batches[[ref]]] <- dat[, batches[[ref]]]
  }
  
  return(list(correctedData=bayesdata,batchDesign=batch.design,gamma.star=gamma.star,delta.star=delta.star,varpool=var.pooled,stdmean=stand.mean))
}


BatchCorrection<-function(data1,data2,site1="Broad",site2="Sanger",CombatRes,stdPrior=TRUE,qcThresh=NULL,qcvalues1=NULL,qcvalues2=NULL){
  #need to make sure it's broad first then sanger for passing to adjustnewdata
  if(!is.null(qcThresh)){
    data1<-data1[,qcvalues1>=qcThresh]
    data2<-data2[,qcvalues2>=qcThresh]
  }
  site=c(rep(site1,ncol(data1)),rep(site2,ncol(data2)))
  adjusted<-AdjustNewData(data1,data2,CombatRes,site,stdPrior)
  return(adjusted)
}

AdjustNewData<-function(data1,data2,CombatRes,site,stdPrior=TRUE){
  blevels<-levels(as.factor(site))
  batch.design2<-CombatRes$batchDesign
  mean.star<-CombatRes$gamma.star
  var.star<-CombatRes$delta.star
  ngenes<-nrow(data1)
  mergedata<-cbind(data1,data2)
  sd1<-ncol(data1)
  
  dn<-dimnames(mergedata)
  mergedata<-normalize.quantiles(as.matrix(mergedata))
  dimnames(mergedata)<-dn
  data1<-mergedata[,colnames(data1)]
  data2<-mergedata[,colnames(data2)]
  usegenes<-intersect(rownames(data1),rownames(CombatRes$stdmean))
  data1<-data1[usegenes,]
  data2<-data2[usegenes,]
  ngenes<-length(usegenes)
  if(stdPrior){
    B_meanStd2<-matrix(CombatRes$stdmean[usegenes,1],nrow=ngenes,ncol=ncol(data1))
    S_meanStd2<-matrix(CombatRes$stdmean[usegenes,1],nrow=ngenes,ncol=ncol(data2))
    B_varStd2<-matrix(sqrt(CombatRes$varpool[usegenes,]),nrow=ngenes,ncol=ncol(data1))
    S_varStd2<-matrix(sqrt(CombatRes$varpool[usegenes,]),nrow=ngenes,ncol=ncol(data2))
  
    D1_all_std2<-(data1-B_meanStd2)/B_varStd2
    D2_all_std2<-(data2-S_meanStd2)/S_varStd2
  }else{
    grand.mean1<-rowMeans(data1)
    grand.mean2<-rowMeans(data2)
    stand.mean1 <- t(grand.mean1) %*% t(rep(1,ncol(data1)))
    stand.mean2 <- t(grand.mean2) %*% t(rep(1,ncol(data2)))
    sd1<-rowSds(data1)
    sd2<-rowSds(data2)
    B_varStd2<-t(sd1)%*%t(rep(1,ncol(data1)))
    S_varStd2<-t(sd2)%*%t(rep(1,ncol(data2)))
    
    D1_all_std2<-(data1-stand.mean1)/B_varStd2
    D2_all_std2<-(data2-stand.mean2)/S_varStd2
    
  }
  colnames(mean.star)<-rownames(CombatRes$stdmean)
  B_meanAll2<-matrix(mean.star[1,usegenes],nrow=ngenes,ncol=ncol(data1))
  S_meanAll2<-matrix(mean.star[2,usegenes],nrow=ngenes,ncol=ncol(data2))
  colnames(var.star)<-rownames(CombatRes$stdmean)
  B_var2<-sqrt(var.star[which(blevels=="Broad"),usegenes])%*%t(rep(1,ncol(data1)))
  S_var2<-sqrt(var.star[which(blevels=="Sanger"),usegenes])%*%t(rep(1,ncol(data2)))
  Broad_all_adjust2<-(D1_all_std2-B_meanAll2)/B_var2
  Sanger_all_adjust2<-(D2_all_std2-S_meanAll2)/S_var2
  
  Broad_all_adjust2<-(Broad_all_adjust2)*B_varStd2+B_meanStd2
  Sanger_all_adjust2<-(Sanger_all_adjust2)*S_varStd2+S_meanStd2
  #need to plot this data:
  alldata2<-cbind(Broad_all_adjust2,Sanger_all_adjust2)
  AdjData<-cbind(Broad_all_adjust2,Sanger_all_adjust2)
  dn<-dimnames(alldata2)
  alldata2<-normalize.quantiles(as.matrix(alldata2))
  dimnames(alldata2)<-dn
  return(list(qNorm=alldata2,NoNorm=AdjData))
  
}




qualityScale<-function(inputdata,ess_genes){
  ess_genes<-intersect(ess_genes,rownames(inputdata))
  scaleddata<-inputdata
  for(i in 1:nrow(inputdata)){
    mv<-mean(inputdata[i,])
    sdv<-sd(inputdata[i,])
    scaleddata[i,]<-(inputdata[i,]-mv)/sdv
  }
  ess_mod<-apply(inputdata[ess_genes,],2,function(x) median(x)+1)
  weights<-apply(scaleddata,2,function(x) x/sum(abs(x)))
  weights2<-apply(weights,2,function(x) (x-min(x))/(max(x)-min(x)))
  adjusteddata<-inputdata
  for(i in 1:ncol(inputdata)){
    adjusteddata[,i]<-inputdata[,i]-ess_mod[i]-weights2[,i]*ess_mod[i]
  }
  dn<-dimnames(adjusteddata)
  newdata<-normalize.quantiles(as.matrix(adjusteddata))
  dimnames(newdata)<-dn
  return(list(adjdata=adjusteddata,normdata=newdata))
}
plotPCA_2batch<-function(data,col_vec,scaleF=TRUE){
  pcaOall<-prcomp(t(data),scale.=scaleF)
eigs<-(pcaOall$sdev)^2
totaleigs<-sum(eigs)

  plot(t(pcaOall$rotation[,1])%*%data,t(pcaOall$rotation[,2])%*%data,col=col_vec,xlab=paste("PC1",round((eigs[1]/totaleigs)*100,2),"% variance"),ylab=paste("PC2",round((eigs[2]/totaleigs)*100,2),"% variance"))
}
RemovePC<-function(data,droppcanumber=1,perfCheck=TRUE){
  if(sum(is.na(data))!=0){
    #Have NAs and need to impute missing values
    #data is genes x cell lines
    meanVals<-rowMeans(data,na.rm=TRUE)
    genesToimpute<-which(rowSums(is.na(data))!=0)
    for(i in 1:length(genesToimpute)){
      selcl<-which(is.na(data[genesToimpute[i],]))
      data[genesToimpute[i],selcl]<-meanVals[genesToimpute[i]]
    }
  }
  estpca<-prcomp(t(data),scale.=TRUE)
  npcas<-1:ncol(data)
  pcause<-npcas[!npcas%in%droppcanumber]
  df.denoised <- estpca$x[,pcause] %*% t(estpca$rotation[,pcause])
  df.denoised<-t(df.denoised)
  correctedData<-df.denoised*estpca$scale+estpca$center
  #test data
  if(perfCheck){
  res<-classPerfCP(correctedData)
  return(list(correctedData=correctedData,Res=res))
  }else{
    return(correctedData)
  }
}

silhouetteScores<-function(clusterlabels,distmat){
  if(!is.numeric(clusterlabels)){
    clustercodes<-as.integer(1:length(unique(clusterlabels)))
    names(clustercodes)<-unique(clusterlabels)
    cc<-clustercodes[clusterlabels]
    names(cc)<-names(clusterlabels)
    resS<-silhouette(x=round(cc),dist=distmat)
  }else{
    resS<-silhouette(x=clusterlabels,dist=distmat)
  }

  return(summary(resS)$avg.width)
}
CompareRNAseqProfiles<-function(BroadData,SangerData,topn=500,Sgenelist=NULL,Bgenelist=NULL){
  usegenes<-intersect(rownames(BroadData),rownames(SangerData))
  BroadData<-BroadData[usegenes,]
  SangerData<-SangerData[usegenes,]
  Broadvars<-rowVars(BroadData)
  Sangervars<-rowVars(SangerData)
  TopxB<-Broadvars[order(Broadvars,decreasing = T)][1:topn]
  TopxS<-Sangervars[order(Sangervars,decreasing=T)][1:topn]
  N<-nrow(SangerData)
  #BroadData = Broad RNA data
  #SangerData = Sanger RNA data

  if(!is.null(Sgenelist)){
    sublist<-Sgenelist[Sgenelist[,1]%in%names(TopxS),]
    rownames(sublist)<-sublist[,1]
    SangerData<-SangerData[sublist[names(TopxS),1],]
    rownames(SangerData)<-sublist[names(TopxS),2]
  }
  if(!is.null(Bgenelist)){
    sublist<-Bgenelist[Bgenelist[,1]%in%names(TopxS),]
    rownames(sublist)<-sublist[,1]
    BroadData<-BroadData[sublist[names(TopxB),1],]
    rownames(BroadData)<-sublist[names(TopxB),2]
  }
Broad_Var<-BroadData[TopxB,]
Sanger_Var<-SangerData[TopxS,]
  

  
  usegenes<-union(rownames(Broad_Var),rownames(Sanger_Var))
  Bdata<-BroadData[usegenes,]
  Sdata<-SangerData[usegenes,]
  
  distMat<-as.matrix(1-cor(cbind(Bdata,Sdata)))
  return(distMat)
}

ROCCurvesWvB<-function(dataset,features,site,weights=NULL){
  #first do the test based on all samples for the features to use
  #will have to remove some features if they are specific to one institute:
  dataset<-dataset[,names(features)]
  #first test features when they can match from either institute:
  resBoth<-classPerfFeature(dataset,weights=weights,feature=features)
  #then need to subset the features by site (prob best to do when doing the classPerf as otherwise need a subsetting loop here)
  groups<-unique(site)
  ResSplit<-list()
  for(i in 1:length(groups)){
    ResSplit[[i]]<-classPerfFeature(dataset,weights=weights,feature=features,site=site,groups[i])
  }
}

classPerfFeature<-function(dataset,qualityTH=Inf,QC=NULL,weights=NULL,geneset=NULL,distmetric="Cor",feature,site=NULL,subGroup=NULL){
  if(!is.null(geneset)){
    dataset<-dataset[intersect(rownames(dataset),geneset),]
  }
  
  clnames<-colnames(dataset)
  
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
    weightmat<-matrix(weights,nrow=nrow(dataset),ncol=ncol(dataset))
    cdist<-as.matrix(1-cor.wt(dataset,w=weights)$r)
  }
  ncls<-ncol(cdist)
  cdist<-cdist[names(feature),names(feature)]
  MATCH<-NULL
  outputDist<-list()
  if(!is.null(site)){
    Totest<-rownames(cdist)[site==subGroup]
    Toexc<-rownames(cdist)[site!=subGroup]
    featureTouse<-feature[names(site)[site==subGroup]]
    featuresTesting<-feature[names(site)[site!=subGroup]]
    checkNotSingleton<-sapply(featureTouse,function(x) sum(featuresTesting==x,na.rm=T))
    Totest<-Totest[checkNotSingleton>0]
    cdist<-cdist[c(Totest,Toexc),c(Totest,Toexc)]
 
    for(i in 1:length(Totest)){
      neighbours<-Toexc[order(cdist[i,Toexc])]
      feature_byneighbour<-feature[neighbours]
      featurematch<-feature[Totest[i]]
      matchingfeatures<-which(feature_byneighbour==featurematch)[1]
      mRow<-rep(0,length(neighbours))
      mRow[matchingfeatures]<-1
      MATCH<-rbind(MATCH,mRow)

      outputDist[[i]]<-feature_byneighbour
      
    }
  }else{
    for(i in 1:ncls){
    #needs to be changed to get the closest neighbour that matches the cell lines input feature
   
      neighbours<-colnames(cdist)[setdiff(order(cdist[i,]),i)]
      feature_byneighbour<-feature[neighbours]
      featurematch<-feature[colnames(cdist)[i]]
      matchingfeatures<-which(feature_byneighbour==featurematch)[1]
      mRow<-rep(0,length(neighbours))
      mRow[matchingfeatures]<-1
      MATCH<-rbind(MATCH,mRow)
      outputDist[[i]]<-feature_byneighbour
    }
    
  }
  outinfo<-feature[colnames(cdist)]
  MATCH<-MATCH+0
  if(is.null(site)){
    rownames(MATCH)<-rownames(cdist)
    curv<-cumsum(colSums(MATCH))/ncol(cdist)
  }else{
    rownames(MATCH)<-Totest
    curv<-cumsum(colSums(MATCH))/length(Totest)
  }
  return(list(MATCHmat=MATCH,CURV=curv,distinfo=outputDist,Tested=outinfo))
}

distPlotCP<-function(dataset,title,XLIM,YLIMS,extraDist=NULL,weights=TRUE){
  clnames<-strsplit(colnames(dataset),'---')
  clnames<-unlist(lapply(clnames,function(x){x[1]}))
  
  site<-strsplit(colnames(dataset),'---')
  site<-unlist(lapply(site,function(x){x[[2]]}))
  
  
  site1dataset<-dataset[,which(site=='Broad')]
  site2dataset<-dataset[,which(site=='Sanger')]
  if(!weights){
    cdist1<-c(as.dist(cor(site1dataset)))
    cdist2<-c(as.dist(cor(site2dataset)))
  
    cdist<-cor(dataset)}else{
      w1<-apply(site1dataset,1,e1071::skewness)
      w2<-apply(site2dataset,1,e1071::skewness)
      cdist1<-cor.wt(site1dataset,w=abs(w1))$r
      cdist2<-cor.wt(site2dataset,w=abs(w2))$r
      wa<-CalculateSkew(site1dataset,site2dataset,"mean")
      cdist<-cor.wt(dataset,w=abs(wa))$r
  }
  if(!is.null(extraDist)){
    if(weights){
      clnamese<-strsplit(colnames(extraDist),'---')
      clnamese<-unlist(lapply(clnamese,function(x){x[1]}))
      
      sitee<-strsplit(colnames(extraDist),'---')
      sitee<-unlist(lapply(sitee,function(x){x[2]}))
      
      
      ex1data<-extraDist[,which(site=='Broad')]
      ex2data<-extraDist[,which(site=='Sanger')]
      #wex<-CalculateSkew(ex1data,ex2data,"mean")
      w1e<-apply(ex1data,1,e1071::skewness)
      w2e<-apply(ex2data,1,e1071::skewness)
      cdist1e<-cor.wt(ex1data,w=abs(w1e))$r
      cdist2e<-cor.wt(ex2data,w=abs(w2e))$r
      #cdistex<-cor.wt(extraDist,w=abs(wex))$r
    }else{
      cdistex<-cor(extraDist)}
  }
  ucl<-unique(clnames)
  
  cdistSAME<-NULL
  
  for (i in 1:length(ucl)){
    id<-which(clnames==ucl[i])
    
    cdistSAME<-c(cdistSAME,cdist[id[1],id[2]])
    cdist[id[1],id[2]]<-NA
    cdist[id[2],id[1]]<-NA
    
    
  }
  
  cdistALL<-as.dist(cdist)
  if(!is.null(extraDist)){
    
    multDensPlot(list(
      density(cdist2e,na.rm=TRUE),
      density(cdist1e,na.rm = TRUE),
      density(cdist2,na.rm =TRUE),density(cdist1,na.rm=TRUE)),XLIMS=XLIM,TITLE = title,YLIMS=YLIMS,
      COLS = c('darkblue','blue','darkgreen','green'),LEGentries = c('Unprocessed Sanger','Unprocessed Broad',
        'Processed Sanger',
        'Processed Broad'))
  }else{
  multDensPlot(list(
                    density(cdistALL,na.rm = TRUE),
                    density(cdistSAME,na.rm =TRUE)),XLIMS=XLIM,TITLE = title,YLIMS=YLIMS,
               COLS = c('gray','darkgreen'),LEGentries = c(
                                                                           'Overall',
                                                                           'Same Cell Line'))
  }
  #print(paste('Expected R within screen (Broad)',median(cdist1)))
  #print(paste('Expected R within screen (Sanger)',median(cdist2)))
  #print(paste('Expected overall R',median(cdistALL,na.rm = TRUE)))
  #print(paste('R withn cell line across screens',median(cdistSAME)))
  
  #print(t.test(cdistSAME,cdist1)$p.value)
  #print(t.test(cdistSAME,cdist2)$p.value)
  
}



plotKNN<-function(Reslist,labels=NULL,filename,transparentScale=FALSE,COLS=NULL,ltylist=NULL,xlim=NULL){
  color = grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), 
                                   invert = T)]
  
  if(is.null(COLS)){
    if(is.list(Reslist)){
      COLS <- sample(color, length(Reslist))
    }else{
      COLS=sample(color,1)
    }
  }
  if(is.null(ltylist)){
    ltylist<-rep(1,length(Reslist))
  }
  if(transparentScale){
    alphaseq<-seq(200,50,length.out=length(Reslist))
    #COLS<-sample(color,1)
    collist<-sapply(alphaseq,function(x) makeTransparent(COLS,alpha=x))
  }else{
    collist <- makeTransparent(COLS, alpha = 150)
  }

if(is.null(labels)){labels=names(Reslist)}

  par(mfrow=c(1,1))
  par(mar=c(0.5,0.5,0.1,0.2)+0.1,mgp=c(1,0.25,0))
  if(is.list(Reslist)){
    
    Res1<-Reslist[[1]]
    npoints<-length(Res1)
  pdf(filename,width=4,height=4)
  par(mar=c(2,2,0.1,0.2)+0.1,mgp=c(1,0.25,0))
  if(is.null(xlim)){
  plot(100*Res1,frame.plot = FALSE,col=collist[1],lwd=5,type='l',
       xlab='neighbourhood',ylab='% cell lines matching counterpart',cex=0.8,cex.axis=0.8,cex.lab=0.8,lty=ltylist[1])}else{
         plot(100*Res1,frame.plot = FALSE,col=collist[1],lwd=5,type='l',
              xlab='neighbourhood',ylab='% cell lines matching counterpart',xlim=xlim,cex=0.8,cex.axis=0.8,cex.lab=0.8,lty=ltylist[1])
       }
  #lines(1:npoints,100*1:npoints*1/npoints,col=makeTransparent('black',150))
  labels[1]<-paste(labels[1],percent(Res1[1]))
  for(i in 2:length(Reslist)){
    lines(1:npoints,100*Reslist[[i]],col=collist[i],lwd=5,type='l',lty=ltylist[i])
    labels[i]<-paste(labels[i],percent(Reslist[[i]][1]))
  }
  

nAUClist<-lapply(Reslist,function(x) trapz(1:npoints,100*x)/(100*npoints))

  
  legend('bottomright',legend=paste(labels,
                                    ' (nAUC =',
                                    format(unlist(nAUClist)
                                    ,digits=2),
                                    ')',sep=''),
         col=collist,lwd=4,cex=0.7,lty=ltylist)
  
  dev.off()
  }else{
    npoints<-length(Reslist)
    pdf(filename,width=3,height=3)
    par(mar=c(2.0,2.0,0.1,0.2)+0.1,mgp=c(1,0.25,0))
    if(is.null(xlim)){
      plot(100*Reslist,frame.plot = FALSE,col=collist[1],lwd=5,type='l',
         xlab='neighbourhood',ylab='% cell lines matching counterpart',cex=0.8,cex.axis=0.8,cex.lab=0.8)
    }else{
           plot(100*Reslist,frame.plot = FALSE,col=collist[1],lwd=5,type='l',xlim=xlim,
                xlab='neighbourhood',ylab='% cell lines matching counterpart',cex=0.8,cex.axis=0.8,cex.lab=0.8)
    }
    lines(1:npoints,100*1:npoints*1/npoints,col=makeTransparent('black',150))
    nauc<-trapz(1:npoints,100*Reslist)/(100*npoints)
    legend('bottomright',legend=paste(labels,
                                      ' (nAUC =',
                                      format(nauc,digits=2),
                                      ')',sep=''),col=collist[1],lwd=4,cex=0.7)
    
    dev.off()
  }
  
}
plotKNNjoint<-function(Reslist,labels=NULL,filename,transparentScale=FALSE){
  color = grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), 
                                   invert = T)]
  if(is.list(Reslist)){
    COLS <- sample(color, length(Reslist))
  }else{
    COLS=sample(color,1)
  }
  if(transparentScale){
    alphaseq<-seq(200,50,length.out=length(Reslist))
    COLS<-sample(color,1)
    collist<-sapply(alphaseq,function(x) makeTransparent(COLS,alpha=x))
  }else{
    collist <- makeTransparent(COLS, alpha = 150)
  }
  
  if(is.null(labels)){labels=names(Reslist)}
  
  par(mfrow=c(1,1))
  par(mar=c(0.5,0.5,0.1,0.2)+0.1,mgp=c(1,0.25,0))
  if(is.list(Reslist)){
    
    Res1<-Reslist[[1]]
    npoints<-length(Res1$CURVjoint)
    pdf(filename,width=4,height=4)
    par(mar=c(2,2,0.1,0.2)+0.1,mgp=c(1,0.25,0))
    plot(100*Res1$CURVjoint,frame.plot = FALSE,col=collist[1],lwd=5,type='l',
         xlab='minimum neighbourhood',ylab='% matching cell lines',cex=0.8,cex.axis=0.8,cex.lab=0.8)
    lines(1:npoints,100*1:npoints*1/npoints,col=makeTransparent('black',150))
    for(i in 2:length(Reslist)){
      lines(1:npoints,100*Reslist[[i]]$CURVjoint,col=collist[i],lwd=5,type='l')
    }
    
    
    nAUClist<-lapply(Reslist,function(x) trapz(1:npoints,100*x$CURVjoint)/(100*npoints))
    
    
    legend('bottomright',legend=paste(labels,
                                      ' (nAUC =',
                                      format(unlist(nAUClist)
                                             ,digits=2),
                                      ')',sep=''),
           col=collist,lwd=4,cex=0.7)
    
    dev.off()
  }else{
    npoints<-length(Reslist)
    pdf(filename,width=3,height=3)
    par(mar=c(2.0,2.0,0.1,0.2)+0.1,mgp=c(1,0.25,0))
    plot(100*Reslist,frame.plot = FALSE,col=collist[1],lwd=5,type='l',
         xlab='minimum neighbourhood',ylab='% matching cell lines',cex=0.8,cex.axis=0.8,cex.lab=0.8)
    lines(1:npoints,100*1:npoints*1/npoints,col=makeTransparent('black',150))
    nauc<-trapz(1:npoints,100*Reslist)/(100*npoints)
    legend('bottomright',legend=paste(labels,
                                      ' (nAUC =',
                                      format(nauc,digits=2),
                                      ')',sep=''),col=collist[1],lwd=4,cex=0.7)
    
    dev.off()
  }
  
}

Sampledata<-function(dataset,feature,setsize){
  #proportions based on sanger data set
  proportions<-table(feature)
  prob<-proportions/sum(proportions)
  names(prob)<-names(proportions)
  
  #sample from broad, with proportions same as sanger
  #want to know how many sanger cell lines match to sanger by tissue. Then compare sanger to broad
  
  testsample<-sample(colnames(dataset),setsize,prob=prob[feature])
  return(dataset[,testsample])
  #want to know how many broad versus broad match.
}

plotPCApanel<-function(librarydfB,col=NULL,filename){

  BCas9activity<-librarydfB$cas9
  
  ##PC1##
  savepdf(paste0(filename,"PC1"),height=6)
  par(mfrow=c(1,3))
  plot(librarydfB$PC1,librarydfB$SSMD,xlab="PC 1",ylab="SSMD")
  abline(lm(librarydfB$SSMD~librarydfB$PC1))
  xp<-quantile(librarydfB$PC1,0.65)
  yp<-quantile(librarydfB$SSMD,0.9)
  text(xp,yp,labels=paste("p-value =",signif(summary(lm(librarydfB$SSMD~librarydfB$PC1))$coefficients[2,4],2)))
  
  plot(librarydfB$PC1[!is.na(librarydfB$cas9)],librarydfB$cas9[!is.na(librarydfB$cas9)],xlab="PC 1",ylab="Cas9")
  abline(lm(BCas9activity[!is.na(librarydfB$cas9)]~librarydfB$PC1[!is.na(librarydfB$cas9)]))
  xp<-quantile(librarydfB$PC1,0.65)
  yp<-quantile(BCas9activity[!is.na(librarydfB$cas9)],0.9)
  text(xp,yp,labels=paste("p-value =",signif(summary(lm(BCas9activity[!is.na(librarydfB$cas9)]~librarydfB$PC1[!is.na(librarydfB$cas9)]))$coefficients[2,4],2)))
  
  boxplot(librarydfB$PC1~librarydfB$libV,xlab="pDNA batch")
  dev.off()
  ##PC2##
  savepdf(paste0(filename,"PC2"),height=6)
  par(mfrow=c(1,3))
  plot(librarydfB$PC2,librarydfB$SSMD,xlab="PC 2",ylab="SSMD")
  abline(lm(librarydfB$SSMD~librarydfB$PC2))
  xp<-quantile(librarydfB$PC2,0.65)
  yp<-quantile(librarydfB$SSMD,0.9)
  text(xp,yp,labels=paste("p-value =",signif(summary(lm(librarydfB$SSMD~librarydfB$PC2))$coefficients[2,4],2)))
  
  plot(librarydfB$PC2[!is.na(librarydfB$cas9)],librarydfB$cas9[!is.na(librarydfB$cas9)],xlab="PC 2",ylab="Cas9")
  abline(lm(BCas9activity[!is.na(librarydfB$cas9)]~librarydfB$PC2[!is.na(librarydfB$cas9)]))
  xp<-quantile(librarydfB$PC2,0.65)
  yp<-quantile(BCas9activity[!is.na(librarydfB$cas9)],0.9)
  text(xp,yp,labels=paste("p-value =",signif(summary(lm(BCas9activity[!is.na(librarydfB$cas9)]~librarydfB$PC2[!is.na(librarydfB$cas9)]))$coefficients[2,4],2)))
  
  boxplot(librarydfB$PC2~librarydfB$libV,xlab="pDNA batch")
  dev.off()
  ###PC3###
  savepdf(paste0(filename,"PC3"),height=6)
  par(mfrow=c(1,3))
  plot(librarydfB$PC3,librarydfB$SSMD,xlab="PC 3",ylab="SSMD")
  abline(lm(librarydfB$SSMD~librarydfB$PC3))
  xp<-quantile(librarydfB$PC3,0.65)
  yp<-quantile(librarydfB$SSMD,0.9)
  text(xp,yp,labels=paste("p-value =",signif(summary(lm(librarydfB$SSMD~librarydfB$PC3))$coefficients[2,4],2)))
  
  plot(librarydfB$PC3[!is.na(librarydfB$cas9)],librarydfB$cas9[!is.na(librarydfB$cas9)],xlab="PC 3",ylab="Cas9")
  abline(lm(BCas9activity[!is.na(librarydfB$cas9)]~librarydfB$PC3[!is.na(librarydfB$cas9)]))
  xp<-quantile(librarydfB$PC3,0.65)
  yp<-quantile(BCas9activity[!is.na(librarydfB$cas9)],0.9)
  text(xp,yp,labels=paste("p-value =",signif(summary(lm(BCas9activity[!is.na(librarydfB$cas9)]~librarydfB$PC3[!is.na(librarydfB$cas9)]))$coefficients[2,4],2)))
  
  boxplot(librarydfB$PC3~librarydfB$libV,xlab="pDNA batch")
  dev.off()
  ###PC4 ###
  savepdf(paste0(filename,"PC4"),height=6)
  par(mfrow=c(1,3))
  plot(librarydfB$PC4,librarydfB$SSMD,xlab="PC 4",ylab="SSMD")
  abline(lm(librarydfB$SSMD~librarydfB$PC4))
  xp<-quantile(librarydfB$PC4,0.65)
  yp<-quantile(librarydfB$SSMD,0.9)
  text(xp,yp,labels=paste("p-value =",signif(summary(lm(librarydfB$SSMD~librarydfB$PC4))$coefficients[2,4],2)))
  
  plot(librarydfB$PC4[!is.na(librarydfB$cas9)],librarydfB$cas9[!is.na(librarydfB$cas9)],xlab="PC 4",ylab="Cas9")
  abline(lm(BCas9activity[!is.na(librarydfB$cas9)]~librarydfB$PC4[!is.na(librarydfB$cas9)]))
  xp<-quantile(librarydfB$PC4,0.65)
  yp<-quantile(BCas9activity[!is.na(librarydfB$cas9)],0.9)
  text(xp,yp,labels=paste("p-value =",signif(summary(lm(BCas9activity[!is.na(librarydfB$cas9)]~librarydfB$PC4[!is.na(librarydfB$cas9)]))$coefficients[2,4],2)))
  
  boxplot(librarydfB$PC4~librarydfB$libV,xlab="pDNA batch")
  dev.off()
  
  ###PC5###
  savepdf(paste0(filename,"PC5"),height=6)
  par(mfrow=c(1,3))
  plot(librarydfB$PC5,librarydfB$SSMD,xlab="PC 5",ylab="SSMD")
  abline(lm(librarydfB$SSMD~librarydfB$PC5))
  xp<-quantile(librarydfB$PC5,0.65)
  yp<-quantile(librarydfB$SSMD,0.9)
  text(xp,yp,labels=paste("p-value =",signif(summary(lm(librarydfB$SSMD~librarydfB$PC5))$coefficients[2,4],2)))
  
  plot(librarydfB$PC5[!is.na(librarydfB$cas9)],librarydfB$cas9[!is.na(librarydfB$cas9)],xlab="PC 5",ylab="Cas9")
  abline(lm(BCas9activity[!is.na(librarydfB$cas9)]~librarydfB$PC5[!is.na(librarydfB$cas9)]))
  xp<-quantile(librarydfB$PC5,0.65)
  yp<-quantile(BCas9activity[!is.na(librarydfB$cas9)],0.9)
  text(xp,yp,labels=paste("p-value =",signif(summary(lm(BCas9activity[!is.na(librarydfB$cas9)]~librarydfB$PC5[!is.na(librarydfB$cas9)]))$coefficients[2,4],2)))
  
  boxplot(librarydfB$PC5~librarydfB$libV,xlab="pDNA batch")
  dev.off()
  ###PC6###
  savepdf(paste0(filename,"PC6"),height=6)
  par(mfrow=c(1,3))
  plot(librarydfB$PC6,librarydfB$SSMD,xlab="PC 6",ylab="SSMD")
  abline(lm(librarydfB$SSMD~librarydfB$PC6))
  xp<-quantile(librarydfB$PC6,0.65)
  yp<-quantile(librarydfB$SSMD,0.9)
  text(xp,yp,labels=paste("p-value =",signif(summary(lm(librarydfB$SSMD~librarydfB$PC6))$coefficients[2,4],2)))
  
  plot(librarydfB$PC6[!is.na(librarydfB$cas9)],librarydfB$cas9[!is.na(librarydfB$cas9)],xlab="PC 6",ylab="Cas9")
  abline(lm(BCas9activity[!is.na(librarydfB$cas9)]~librarydfB$PC6[!is.na(librarydfB$cas9)]))
  xp<-quantile(librarydfB$PC6,0.65)
  yp<-quantile(BCas9activity[!is.na(librarydfB$cas9)],0.9)
  text(xp,yp,labels=paste("p-value =",signif(summary(lm(BCas9activity[!is.na(librarydfB$cas9)]~librarydfB$PC6[!is.na(librarydfB$cas9)]))$coefficients[2,4],2)))
  
  boxplot(librarydfB$PC6~librarydfB$libV,xlab="pDNA batch")
  dev.off()
}

distFeat<-function(distmatrix,featuresets){
 #assume featuresets is a BEM
    
    output<-list()
    for(i in 1:nrow(featuresets)){
      sel<-colnames(featuresets)[which(featuresets[i,]==1)]
      submat<-distmatrix[sel,sel]
      output[[i]]<-mean(submat[lower.tri(submat)])
    }

  return(output)
}

savepdf <- function(file, width=16, height=10)
{
  fname <- paste(file,".pdf",sep="")
  pdf(fname, width=width/2.54, height=height/2.54,
      pointsize=10)
  par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
}

plotPIE<-function(classCPres){
  piecounts<-as.vector(table(unlist(classCPres)))
  names(piecounts)<-names(table(unlist(classCPres)))
  piepercent<- round(100*piecounts/sum(piecounts), 1)
  pie(piecounts,labels=paste0(piepercent,"%"),col=c("white","blue","pink"))
  legend("bottomright",legend=names(piecounts),fill=c("white","blue","pink"),cex=0.8)
}
cordists<-function(group,fulldata,backgroundcor){
  if(length(group)>1){
  ingroup<-fulldata[,group]
  ingroupCor<-cor(ingroup)

  inCor<-ingroupCor[lower.tri(ingroupCor)]
  if(length(inCor)>5){
    DensDist<-ssmd(inCor,backgroundcor)
  }else{DensDist=NULL}
  return(list(ingroup=inCor,ssmd=DensDist))}else{
    NULL
  }
}
backdist<-function(data,labels){
  corres<-list()
  for(i in 1:length(labels)){
    if(length(labels[[i]])>0){
    corres[[i]]<-cor(data[,labels[[i]]],data[,!colnames(data)%in%labels[[i]]])}
  }
  return(unlist(corres))
}
CalculateSkew<-function(Bdata,Sdata,method="max"){
  bgenes<-intersect(rownames(Bdata),rownames(Sdata))
  BroadSkew<-apply(Bdata[bgenes,],1,e1071:::skewness)
  SangerSkew<-apply(Sdata[bgenes,],1,e1071:::skewness)
  BothSkew<-cbind(BroadSkew,SangerSkew)
  if(method=="max"){
    SkewWeights<-apply(BothSkew,1,max)}
  if(method=="mean"){
    SkewWeights<-apply(BothSkew,1,mean)
  }
  if(method=="min"){
    SkewWeights<-apply(BothSkew,1,min)
  }
  names(SkewWeights)<-bgenes
  return(SkewWeights)
  
}

ccr.geneMeanFCsGE<-function (sgRNA_FCprofile, libraryAnnotation,w) 
{
  GENES=libraryAnnotation[names(sgRNA_FCprofile),"GENES"]
  inputdf<-cbind(sgRNA_FCprofile,w[names(sgRNA_FCprofile)])
  splitdf<-split.data.frame(inputdf,GENES)
  FCsprofile<-lapply(splitdf,function(x) weighted.mean(x$FC,w=x$weights,na.rm=T))
  nn<-names(FCsprofile)
  FCsprofile<-unlist(FCsprofile)


  names(FCsprofile) <- as.character(nn)
  return(FCsprofile)
}


compareADaMoutput<-function(dir1,dir2,filesuffix){
  Scompare<-list.dirs(dir1)
  So<-list.dirs(dir2)
  temp1<-grep(paste0(".*OUTPUT"),Scompare,value=T)
  temp1<-unique(temp1)
  temp2<-grep(paste0(".*OUTPUT"),So,value=T)
  temp2<-unique(temp2)
  cTypes<-sapply(temp1,function(x) strsplit(x,"_")[[1]][1])
  cTypes<-sapply(cTypes,function(x) strsplit(x,"/")[[1]][4])
  
  cType2<-sapply(temp2,function(x) strsplit(x,"_")[[1]][1])
  cType2<-sapply(cType2,function(x) strsplit(x,"/")[[1]][4])
  cType2<-make.names(cType2)
  compareType<-intersect(cType2,cTypes)
  #compareType<-setdiff(compareType,"Gastric.Carcinoma")
  MissingMarkers<-c()
  comparisonTypes<-data.frame(cbind(compareType,compareType))
  k=1
  for(i in 1:length(compareType)){
    res1<-read.table(paste0(temp1[which(cTypes==compareType[i])],"/ANOVA_results.txt"),header=T,sep="\t",stringsAsFactors = F)
    res2<-read.table(paste0(temp2[which(cType2==compareType[i])],"/ANOVA_results.txt"),header=T,sep="\t",stringsAsFactors = F)
    
    rownames(res1)<-paste(res1[,2],res1[,3],sep="_")
    rownames(res2)<-paste(res2[,2],res2[,3],sep="_")
    
    fdrBH1<-p.adjust(res1[,"FEATURE_ANOVA_pval"],method="fdr")
    fdrBH2<-p.adjust(res2[,"FEATURE_ANOVA_pval"],method="fdr")
    names(fdrBH2)<-rownames(res2)
    names(fdrBH1)<-rownames(res1)
    print(compareType[i])
    print(length(fdrBH1))
    print(length(fdrBH2))
    #res 1 is corrected
    #res 2 in uncorrected
    sigUncorrected<-rownames(res2)[res2$ANOVA.FEATURE.FDR..<30]
    print(paste('number uncorrected <30',length(sigUncorrected)))
    resC<-res1[intersect(sigUncorrected,rownames(res1)),]
    sigUnOnly<-rownames(resC)[resC$ANOVA.FEATURE.FDR..>30]
    uncorrectedEffect<-res2[sigUnOnly,"FEATUREpos_Glass_delta"]
    correctedEffect<-res1[sigUnOnly,"FEATUREpos_Glass_delta"]
    uncorrectedpval<-res2[sigUnOnly,"FEATURE_ANOVA_pval"]
    correctedpval<-res1[sigUnOnly,"FEATURE_ANOVA_pval"]
    uncorrectedFDR<-res2[sigUnOnly,"ANOVA.FEATURE.FDR.."]
    correctedFDR<-res1[sigUnOnly,"ANOVA.FEATURE.FDR.."]

    feature<-res1[sigUnOnly,2]
    dependency<-res1[sigUnOnly,3]
    if(k==1){
      if(length(sigUnOnly)>0){
        MissingMarkers<-cbind(feature,dependency,compareType[i],uncorrectedEffect,correctedEffect,uncorrectedpval,correctedpval,uncorrectedFDR,correctedFDR)
        k=2}
    }else{
      if(length(sigUnOnly)>0){
        MissingMarkers<-rbind(MissingMarkers,cbind(feature,dependency,compareType[i],uncorrectedEffect,correctedEffect,uncorrectedpval,correctedpval,uncorrectedFDR,correctedFDR))}
    }

  }
  write.csv(MissingMarkers,file=paste0(dir.Results,"/MissingMarkers_",filesuffix,".csv"))
}


compareADaMplots<-function(dir1,dir2,plotprefix){
  Scompare<-list.dirs(dir1)
  So<-list.dirs(dir2)
  temp1<-grep(paste0(".*OUTPUT"),Scompare,value=T)
  temp2<-grep(paste0(".*OUTPUT"),So,value=T)
  cTypes<-sapply(temp1,function(x) strsplit(x,"_")[[1]][1])
  cTypes<-sapply(cTypes,function(x) strsplit(x,"/")[[1]][4])
  
  cType2<-sapply(temp2,function(x) strsplit(x,"_")[[1]][1])
  cType2<-sapply(cType2,function(x) strsplit(x,"/")[[1]][4])
  cType2<-make.names(cType2)
  compareType<-intersect(cType2,cTypes)
  #compareType<-setdiff(compareType,"Gastric.Carcinoma")

  k=1
  for(i in 1:length(compareType)){
    res1<-read.table(paste0(temp1[which(cTypes==compareType[i])],"/ANOVA_results.txt"),header=T,sep="\t",stringsAsFactors = F)
    res2<-read.table(paste0(temp2[which(cType2==compareType[i])],"/ANOVA_results.txt"),header=T,sep="\t",stringsAsFactors = F)
    
    rownames(res1)<-paste(res1[,2],res1[,3],sep="_")
    rownames(res2)<-paste(res2[,2],res2[,3],sep="_")
    
    sig25_1<-rownames(res1)[res1$ANOVA.FEATURE.FDR..<10]
    sig25_2<-rownames(res2)[res2$ANOVA.FEATURE.FDR..<10]
    #res 1 is corrected
    #res 2 in uncorrected

    testedboth<-intersect(rownames(res1),rownames(res2))
    plotboth<-intersect(testedboth,union(sig25_1,sig25_2))

  #plot(res1[plotboth,"ANOVA.FEATURE.FDR.."],res2[plotboth,"ANOVA.FEATURE.FDR.."],xlab="Corrected",ylab="Un corrected")
  
  cT<-cor(res1[testedboth,"ANOVA.FEATURE.FDR.."],res2[testedboth,"ANOVA.FEATURE.FDR.."])
  plotdata<-cbind(res1[testedboth,"ANOVA.FEATURE.FDR.."],res2[testedboth,"ANOVA.FEATURE.FDR.."])
  
  colnames(plotdata)<-c("Batch_Corrected","Site_specific")
  rownames(plotdata)<-testedboth
  label_data<-data.frame(Batch_Corrected=as.numeric(plotdata[plotboth,1]),Site_specific=as.numeric(plotdata[plotboth,2]),Label=plotboth,stringsAsFactors = FALSE)
  
  pdf(paste0(dir.Results,"/Corr_SangerOrigVsSangerCorrected_",compareType[i],".pdf"))
  
  #add labels according to plotboth i.e. signif at 25% FDR either set.
  plotout<-ggplot(plotdata, aes(x=Batch_Corrected,y=Site_specific))+geom_point(alpha=0.5,col='grey')+
    #xlim(-2.5,1)+
    #ylim(-2.5,1)+
    geom_hline(yintercept= 25,linetype="dashed",col="gray")+
    geom_vline(xintercept = 25,linetype="dashed",col="gray")+
    theme_bw()+
    geom_text_repel(data=label_data,aes(x=Batch_Corrected,y=Site_specific,label=Label))+
    geom_point(data=plotdata,aes(x=Batch_Corrected,y=Site_specific),col='purple')+
    labs(x="Sanger Batch Corrected FDR %",y="Sanger Un corrected FDR %",cex=0.8)+
    annotate("text", x = c(50), y=95, label = paste("Correlation:",round(cT,2)))
  
  print(plotout)
  #add labels according to plotboth i.e. signif at 25% FDR either set.
  dev.off()
  
  
  sig25_1<-rownames(res1)[res1$ANOVA.FEATURE.FDR..<10]
  sig25_2<-rownames(res2)[res2$ANOVA.FEATURE.FDR..<10]
  
  
  plotSingle<-intersect(setdiff(sig25_1,sig25_2),testedboth)
  if(length(plotSingle)>0){
    plotboth<-intersect(testedboth,union(sig25_1,sig25_2))
    
    label_data<-data.frame(Batch_Corrected=as.numeric(plotdata[plotSingle,1]),Site_specific=as.numeric(plotdata[plotSingle,2]),Label=plotSingle,stringsAsFactors = FALSE)
    
    pdf(paste0(dir.Results,"/CorrCombine_SangerOrigVsSangerCorrected_",compareType[i],".pdf"))
    
    #add labels according to plotboth i.e. signif at 25% FDR either set.
    plotout<-ggplot(plotdata, aes(x=Batch_Corrected,y=Site_specific))+geom_point(alpha=0.5,col='grey')+
      #xlim(-2.5,1)+
      #ylim(-2.5,1)+
      geom_hline(yintercept= 25,linetype="dashed",col="gray")+
      geom_vline(xintercept = 25,linetype="dashed",col="gray")+
      theme_bw()+
      geom_text_repel(data=label_data,aes(x=Batch_Corrected,y=Site_specific,label=Label))+
      geom_point(data=plotdata,aes(x=Batch_Corrected,y=Site_specific),col='purple')+
      labs(x="Sanger Batch Corrected FDR %",y="Sanger Un corrected FDR %",cex=0.8)+
      annotate("text", x = c(50), y=95, label = paste("Correlation:",round(cT,2)))
    
    print(plotout)
    #add labels according to plotboth i.e. signif at 25% FDR either set.
    dev.off()
  }
  plotCombine<-intersect(setdiff(sig25_2,sig25_1),testedboth)
  if(length(plotCombine)>0){
    label_data<-data.frame(Batch_Corrected=as.numeric(plotdata[plotCombine,1]),Site_specific=as.numeric(plotdata[plotCombine,2]),Label=plotCombine,stringsAsFactors = FALSE)
    
    pdf(paste0(dir.Results,"/CorrSingle_SangerOrigVsSangerCorrected_",compareType[i],".pdf"))
    
    #add labels according to plotboth i.e. signif at 25% FDR either set.
    plotout<-ggplot(plotdata, aes(x=Batch_Corrected,y=Site_specific))+geom_point(alpha=0.5,col='grey')+
      #xlim(-2.5,1)+
      #ylim(-2.5,1)+
      geom_hline(yintercept= 25,linetype="dashed",col="gray")+
      geom_vline(xintercept = 25,linetype="dashed",col="gray")+
      theme_bw()+
      geom_text_repel(data=label_data,aes(x=Batch_Corrected,y=Site_specific,label=Label))+
      geom_point(data=plotdata,aes(x=Batch_Corrected,y=Site_specific),col='purple')+
      labs(x="Sanger Batch Corrected FDR %",y="Sanger Un corrected FDR %",cex=0.8)+
      annotate("text", x = c(50), y=95, label = paste("Correlation:",round(cT,2)))
    
    print(plotout)
    #add labels according to plotboth i.e. signif at 25% FDR either set.
    dev.off()
  }
}
  
}

makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

OT15_assembleDepletionMatrix<-function(cellLines,BFths,inputFolder='../../RESULTS/BAGEL-R_output/',BFvalues=FALSE,BFmatrix=NULL){
  
  for (i in 1:length(cellLines)){
    print(i)
    if(!is.null(inputFolder)){
      if(file.exists(paste(inputFolder,cellLines[i],'_GeneLevelBF.rdata',sep=''))){
        load(paste(inputFolder,cellLines[i],'_GeneLevelBF.rdata',sep=''))
        if(BFvalues){
          tmp<-matrix(GeneLevelBFs[,1],nrow(GeneLevelBFs),1,
                      dimnames = list(rownames(GeneLevelBFs),cellLines[i]))
        }else{
          tmp<-matrix(GeneLevelBFs[,1]-BFths[cellLines[i]],nrow(GeneLevelBFs),1,
                      dimnames = list(rownames(GeneLevelBFs),cellLines[i])) 
        }
      }
    }else{
      if(!is.null(BFmatrix)){
        tmp<-matrix(BFmatrix[,cellLines[i]]-BFths[cellLines[i]],nrow(BFmatrix),1,
                    dimnames = list(rownames(BFmatrix),cellLines[i]))   
      }
    }
    #if(i==1){
    if(!exists("resMat")){
      resMat<-tmp
    }else{
      if(nrow(tmp)<nrow(resMat)){
        missingG<-setdiff(rownames(resMat),rownames(tmp))
        toAdd<-matrix(NA,length(missingG),dimnames = list(missingG,NULL))
        tmp<-rbind(tmp,toAdd)
        tmp<-matrix(tmp[rownames(resMat),],nrow(resMat),dimnames = list(rownames(resMat),cellLines[i]))
      }
      resMat<-cbind(resMat,tmp)
    }
    
    
  }
  
  resMat<-resMat+0
  return(resMat)
}

bangmag.roc<-function(BFs,mgkdata_pval,mgkdata_fdr,ess_genes,non_ess_genes,CL,indSet=FALSE,th,doplot=FALSE){
  Genes<-rownames(BFs)
  
  ess_genes<-intersect(ess_genes,Genes)
  non_ess_genes<-intersect(non_ess_genes,Genes)
  
  allg<-c(ess_genes,non_ess_genes)
  
  essentiality<-rep(0,length(allg))
  names(essentiality)<-allg
  
  essentiality[ess_genes]<-1
  predBFs<- BFs[allg,1]
  
  BF_roc<-roc(essentiality,predBFs-min(predBFs,na.rm=TRUE))
  
  if(indSet){
    TITLE<-'Independent gene sets'
  }else{
    TITLE<-'BAGEL reference gene sets'
  }
  bf_coords<-pROC::coords(roc=BF_roc,x='all',ret=c('threshold','sensitivity','specificity','ppv','precision','recall'),transpose=TRUE)
  bf_coords['threshold',]<-bf_coords['threshold',]+min(predBFs,na.rm=TRUE)
  if(doplot){
  plot(BF_roc,col='blue',lwd=3,main=TITLE)
  
  legendEntry<-c(paste('AUC = ',format(BF_roc$auc,digits=3),sep=''))
  legend('bottomright',legend=legendEntry,inset = c(0.1,0.1),bty = 'n',y.intersp = 2)
  
  }
  
  
  id<-min(which(bf_coords['threshold',] > 0))
  bf_0<-c(0,bf_coords['specificity',id],bf_coords['sensitivity',id])
  if(doplot){points(bf_coords['specificity',id],bf_coords['sensitivity',id],col='blue',lwd=3,cex=1.5)}
  
  id<-min(which(bf_coords['threshold',] > 1))
  bf_1<-c(0,bf_coords['specificity',id],bf_coords['sensitivity',id])
  if(doplot){points(bf_coords['specificity',id],bf_coords['sensitivity',id],col='blue',lwd=3,cex=1.5,pch=16)}
  
  id<-min(which(bf_coords['ppv',]>(1-th)))
  if(id=="Inf"){
    id<-min(which(bf_coords['ppv',]>=(1-th)))}
  if(id=="Inf"){
    id<-max(which(round(bf_coords['ppv',])>=(1-th)))
  }
  if(doplot){
  abline(h=bf_coords['sensitivity',id],col='blue',lwd=1)
  text(0.8,pos = 4,bf_coords['sensitivity',id]+0.015,
       paste('PPV > 95% (BF > ',format(bf_coords['threshold',id],digits=3),')',sep=''),col='blue',cex = 0.9)
  }
  bf_best_prec<-c(bf_coords['threshold',id],bf_coords['specificity',id],bf_coords['sensitivity',id],bf_coords['ppv',id])
  bestPrecisionTh<-rbind(bf_best_prec)
  colnames(bestPrecisionTh)<-c('thresholds','specificity','sensitivity','ppv')
  rownames(bestPrecisionTh)<-c('Ban_BF')
  
  return(list(bfAUC=BF_roc$auc,bestTh=bestPrecisionTh,RocCurve=bf_coords))
  
}

decodeCNA_cp<-function(MoBEM){
  
  rn <- rownames(MoBEM)
  ii <- grep("cna", rownames(MoBEM))
  cnaId <- rownames(MoBEM)[ii]
  containedGenes <- unlist(lapply(strsplit(cnaId, " "), function(x) {
    x[2]
  }))
  containedGenes[is.na(containedGenes)] <- ""
  segments <- unlist(lapply(strsplit(unlist(lapply(strsplit(cnaId, 
                                                              ":"), function(x) {
                                                                x[2]
                                                              })), " "), function(x) {
                                                                x[1]
                                                              }))
  loci <- as.character(CNAdecode$locus[match(segments, CNAdecode$Identifier)])
  altType <- as.character(CNAdecode$Recurrent.Amplification..Deletion[match(segments, 
                                                                            CNAdecode$Identifier)])
  altType[altType == "Amplification"] <- "G:"
  altType[altType == "Deletion"] <- "L:"
  rownames(MoBEM)[ii] <- paste(altType, loci, " ", containedGenes, 
                               sep = "")
  return(MoBEM)
}

Allanova<-function (MoBEM, dataset,tissue=NULL) {
usecl<-intersect(colnames(MoBEM),colnames(dataset))
if(length(usecl)>5){
MoBEM<-MoBEM[,usecl]
## Filtering out cancer functional events present in less than 2 and more than 144 cell lines
MoBEM<-MoBEM[rowSums(MoBEM,na.rm=T)>2,]
MoBEM<-MoBEM[rowSums(MoBEM,na.rm=T)<ncol(MoBEM)-2,]
nfet<-nrow(MoBEM)
if(nfet>0){
res<-NULL
#subdataset <- na.omit(dataset[, usecl])
subdataset<-dataset[,usecl]
for(j in 1:nfet){
  print(paste('feature ',j,'of',nfet))


  cfePattern <- MoBEM[j,]

  
  subSetWT <- subdataset[, names(cfePattern)[cfePattern == 0]]
  subSetALT <- subdataset[, names(cfePattern)[cfePattern == 1]]

  effect_size <- unlist(lapply(seq_len(nrow(subSetALT)), function(i) {
    cohens_d(na.exclude(subSetALT[i, ]), na.exclude(subSetWT[i, ]))
  }))
  delta <- unlist(lapply(seq_len(nrow(subSetALT)), function(i) {
    mean(subSetALT[i, ],na.rm=T) - mean(subSetWT[i, ],na.rm=T)
  }))
  pval <- unlist(lapply(seq_len(nrow(subSetALT)), function(i) {
    RES <- t.test(na.exclude(subSetALT[i, ]), na.exclude(subSetWT[i, ]), var.equal = TRUE)
    RES$p.value
  }))
  npos=sum(cfePattern)
  nneg=length(cfePattern)-npos
 
  #res <- data.frame(CFE = CFE, GENE = rownames(dataset), delta = delta, 
             #       effect_size = effect_size, p = pval, fdr = FDR, stringsAsFactors = FALSE)
  if(!is.null(tissue)){
  temp<- data.frame(CFE = rownames(MoBEM)[j], GENE = rownames(subdataset), delta = delta, 
                    effect_size = effect_size, p = pval, npos=npos,nneg=nneg, tissue=tissue,stringsAsFactors = FALSE)
  }else{
    temp<- data.frame(CFE = rownames(MoBEM)[j], GENE = rownames(subdataset), delta = delta, 
                      effect_size = effect_size, p = pval, npos=npos,nneg=nneg, stringsAsFactors = FALSE)
    
  }
  res<-rbind(res,temp)
}
  return(res)}else{
    #no features to run analysis with
    return(NULL)
  }
}else{
  #not enough cell lines with molecular data to run analysis
  return(NULL)
}
}



marker_groups<-function(subpt,Biomarker=c("A","B"),T1=4,T2=c(2,3,5,6,7),T3=1){
  subpt$GROUP<-"None"
  subpt[subpt$MARKERCLASS%in%Biomarker&subpt$TRACTABILITY%in%T1,"GROUP"]<-"High"
  subpt[(!subpt$MARKERCLASS%in%Biomarker)&subpt$TRACTABILITY%in%T1,"GROUP"]<-"Medium"
  subpt[subpt$MARKERCLASS%in%Biomarker&subpt$TRACTABILITY%in%T2,"GROUP"]<-"Medium"
  subpt[subpt$MARKERCLASS%in%Biomarker&subpt$TRACTABILITY%in%T3,"GROUP"]<-"None"
  subpt[subpt$MARKERCLASS%in%Biomarker&(!subpt$TRACTABILITY%in%c(T1,T2,T3)),"GROUP"]<-"Low"
  subpt[(!subpt$MARKERCLASS%in%Biomarker)&subpt$TRACTABILITY%in%T2,"GROUP"]<-"Medium_Low"
  subpt[(!subpt$MARKERCLASS%in%Biomarker)&subpt$TRACTABILITY%in%T3,"GROUP"]<-"None"
  subpt[(!subpt$MARKERCLASS%in%Biomarker)&(!subpt$TRACTABILITY%in%c(T1,T2,T3)),"GROUP"]<-"Low"
  return(subpt)
}

GenerateManifest<-function(origmanifest,feature,mobem,cancertype=NULL,tissue=NULL){
  if(!is.null(cancertype)){
    origmanifest<-origmanifest[origmanifest$CANCER_TYPE==cancertype,]
  }
  if(!is.null(tissue)){
    origmanifest<-origmanifest[origmanifest$TISSUE==tissue,]
  }
  celllineswithfeature<-colnames(mobem)[which(mobem[feature,]==1)]
  origmanifest$CL_GROUP<-NA
  origmanifest[origmanifest$COSMIC_ID%in%celllineswithfeature,"CL_GROUP"]<-feature
  newmanifest<-origmanifest
  newmanifest<-newmanifest[!is.na(newmanifest$CL_GROUP),]
  return(newmanifest)
}

compareDependencies<-function(testData,referenceData){
  overlapCL<-intersect(colnames(testData),colnames(referenceData))
  corCL<-sapply(overlapCL,function(x) mcc(testData[,x],referenceData[,x]))
  names(corCL)<-overlapCL
  return(corCL)
}

perCTdependencies<-function(binarymat,annotation,minCL=5){
  annotinc<-annotation[colnames(binarymat),]
  CTs<-table(annotinc$cancer_type)
  CTinc<-names(CTs)[which(CTs>4)]
  CTdependencies<-list()
  for(i in 1:length(CTinc)){
    selCL<-rownames(annotinc)[which(annotinc$cancer_type==CTinc[i])]
    CTmat<-binarymat[,selCL]
    CTexc<-binarymat[,colnames(binarymat)[!colnames(binarymat)%in%(selCL)]]
    CTdep<-rownames(CTmat)[rowSums(CTmat)>2]
    CTexc<-rownames(CTexc)[rowSums(CTexc)>0]
    CTdependencies[[i]]<-setdiff(CTdep,CTexc)
  }
  names(CTdependencies)<-CTinc
  return(CTdependencies)
}

normLRT<-function(data){

  normLike<-apply(data,1,function(x) fitdistr(na.omit(x),"normal")$loglik)
  
  tLike<-apply(data,1,function(x) fittdist(na.omit(x)))
  lrt<-2*(tLike-normLike)
  return(list(norm=normLike,t=tLike,LRT=lrt))
}
normLRT2<-function(data){
  
  normLike<-apply(data,1,function(x) fitdistr(na.omit(x),"normal")$loglik)
  
  tLike<-apply(data,1,function(x) fitskewtdist(na.omit(x)))
  lrt<-2*(tLike-normLike)
  return(list(norm=normLike,t=tLike,LRT=lrt))
}
fitskewtdist<-function(data){
  mm<-matrix(1,ncol=1,nrow=length(data))
  result=tryCatch({
    st.mple(x=mm,y=data)$logL
  
  },error=function(e){
    f1<-st.mple(x=mm,y=data,fixed.nu=2)
  
    res2=tryCatch({
      st.mple(x=mm,y=data,dp =f1$dp,fixed.nu=2)$logL
    }, error=function(e){
      f1<-st.mple(x=mm,y=data,fixed.nu=5)
      res3=tryCatch({
        st.mple(x=mm,y=data,dp=f1$dp,fixed.nu=5)$logL
      }, error=function(e){
        f1<-st.mple(x=mm,y=data,fixed.nu=10)
        res4=tryCatch({
          st.mple(x=mm,y=data,dp=f1$dp,fixed.nu=10)$logL
        }, error=function(e){
          f1<-st.mple(x=mm,y=data,fixed.nu=25)
          res5=tryCatch({
            st.mple(x=mm,y=data,dp=f1$dp,fixed.nu=25)$logL
          }, error=function(e){
            f1<-st.mple(x=mm,y=data,fixed.nu=50)
            res6=tryCatch({
              st.mple(x=mm,y=data,dp=f1$dp,fixed.nu=50)$logL
            }, error=function(e){
              f1<-st.mple(x=mm,y=data,fixed.nu=100)
              st.mple(x=mm,y=data,dp=f1$dp,fixed.nu=100)$logL
            })
          })
          
        })
        
      })
      
    })
  })
  
  return(result) 
}

fittdist<-function(data){
  result=tryCatch({
    fitdistr(data,"t")$loglik
  },error=function(e){
    f1<-fitdistr(data,"t",df=2)
    res2=tryCatch({
      fitdistr(data,"t",start =list(m=f1$estimate["m"],s=f1$estimate["s"],df=2))$loglik
    }, error=function(e){
      f1<-fitdistr(data,"t",df=5)
      res3=tryCatch({
        fitdistr(data,"t",start =list(m=f1$estimate["m"],s=f1$estimate["s"],df=5))$loglik
      }, error=function(e){
        f1<-fitdistr(data,"t",df=10)
        res4=tryCatch({
          fitdistr(data,"t",start =list(m=f1$estimate["m"],s=f1$estimate["s"],df=10))$loglik
        }, error=function(e){
          f1<-fitdistr(data,"t",df=25)
          res5=tryCatch({
            fitdistr(data,"t",start =list(m=f1$estimate["m"],s=f1$estimate["s"],df=25))$loglik
          }, error=function(e){
            f1<-fitdistr(data,"t",df=50)
            res6=tryCatch({
              fitdistr(data,"t",start =list(m=f1$estimate["m"],s=f1$estimate["s"],df=50))$loglik
            }, error=function(e){
              f1<-fitdistr(data,"t",df=100)
              fitdistr(data,"t",start =list(m=f1$estimate["m"],s=f1$estimate["s"],df=100))$loglik
            })
          })

        })
        
      })
      
    })
  })
  
 return(result) 
}

avgOverlap<-function(data,overlapAnnot,divergentCL=NULL){
  if(!is.null(divergentCL)){
    singledropSanger<-overlapAnnot[overlapAnnot$mn%in%rownames(divergentCL)[which(divergentCL[,3]==2)],"model_id"]
    singledropBroad<-overlapAnnot[overlapAnnot$mn%in%rownames(divergentCL)[which(divergentCL[,3]==1)],"BROAD_ID"]
    dropcl<-c(singledropBroad,singledropSanger)
    data<-data[,!colnames(data)%in%dropcl]
    overlapAnnot<-overlapAnnot[!overlapAnnot$mn%in%rownames(divergentCL),]
      
  }
  
  for(i in 1:nrow(overlapAnnot)){
    subdata<-data[,colnames(data)%in%unlist(overlapAnnot[i,c("model_id","BROAD_ID")])]
    print(i)
    if(length(subdata)>nrow(data)){
      newdata<-matrix(rowMeans(as.matrix(subdata)),ncol=1)
      colnames(newdata)<-overlapAnnot[i,"model_id"]
    }else{
      newdata<-matrix(subdata,ncol=1)
      colnames(newdata)<-overlapAnnot[i,"model_id"]
    }
    
    if(i==1){
      outdata<-newdata
    }else{
      outdata<-cbind(outdata,newdata)
    }
  }
  rownames(outdata)<-rownames(data)
  allover<-c(overlapAnnot[,c("model_id")],overlapAnnot[,c("BROAD_ID")])
  notmerge<-data[,!colnames(data)%in%unlist(allover)]
  finaldata<-cbind(notmerge[rownames(outdata),],outdata)
  return(finaldata)
}

avgOverlapBroadAnnot<-function(data,overlapAnnot,divergentCL=NULL){
  if(!is.null(divergentCL)){
    singledropSanger<-overlapAnnot[overlapAnnot$stripped_cell_line_name%in%rownames(divergentCL)[which(divergentCL[,3]==2)],"Sanger_Model_ID"]
    singledropBroad<-overlapAnnot[overlapAnnot$stripped_cell_line_name%in%rownames(divergentCL)[which(divergentCL[,3]==1)],"DepMap_ID"]
    dropcl<-c(singledropBroad,singledropSanger)
    data<-data[,!colnames(data)%in%dropcl]
    overlapAnnot<-overlapAnnot[!overlapAnnot$stripped_cell_line_name%in%rownames(divergentCL),]
    
  }
  
  for(i in 1:nrow(overlapAnnot)){
    subdata<-data[,colnames(data)%in%unlist(overlapAnnot[i,c("Sanger_Model_ID","DepMap_ID")])]
   
    if(is.matrix(subdata)){
      newdata<-matrix(rowMeans(subdata,na.rm=T),ncol=1)
      colnames(newdata)<-overlapAnnot[i,"Sanger_Model_ID"]
    }else{
      newdata<-matrix(subdata,ncol=1)
      colnames(newdata)<-overlapAnnot[i,"Sanger_Model_ID"]
    }
    if(i==1){
      outdata<-newdata
    }else{
      outdata<-cbind(outdata,newdata)
    }
  }
  rownames(outdata)<-rownames(data)
  allover<-c(overlapAnnot[,c("Sanger_Model_ID")],overlapAnnot[,c("DepMap_ID")])
  notmerge<-data[,!colnames(data)%in%unlist(allover)]
  finaldata<-cbind(notmerge[rownames(outdata),],outdata)
  return(finaldata)
}

updateRownames<-function(inputdata,Map){

  newnames<-Map[rownames(inputdata)]
  names(newnames)<-rownames(inputdata)
  newnames[is.na(newnames)]<-names(newnames)[is.na(newnames)]
  newnames<-newnames[newnames%in%names(table(newnames))[table(newnames)==1]]
  inputdata<-inputdata[names(newnames),]
  rownames(inputdata)<-newnames
  return(inputdata)
}

assembleDepletionMatrix<-function(cellLines,BFths,inputFolder='../../RESULTS/BAGEL-R_output/',BFvalues=FALSE,BFmatrix=NULL){
  
  if(!is.null(inputFolder)){
    print(paste(inputFolder,cellLines[1],'_GeneLevelBF.rdata',sep=''))
    if(file.exists(paste(inputFolder,cellLines[1],'_GeneLevelBF.rdata',sep=''))){
      load(paste(inputFolder,cellLines[1],'_GeneLevelBF.rdata',sep=''))
      genes<-rownames(GeneLevelBFs)
    }  
    
    print(head(genes))
    resMat<-foreach(i=1:length(cellLines),.combine=cbind)%dopar%{
      print(head(genes))
      
      if(file.exists(paste(inputFolder,cellLines[i],'_GeneLevelBF.rdata',sep=''))){
        load(paste(inputFolder,cellLines[i],'_GeneLevelBF.rdata',sep=''))
        if(BFvalues){
          tmp<-matrix(NA,length(genes),1,
                      dimnames = list(genes,cellLines[i]))
          tmp[rownames(GeneLevelBFs),1]<-GeneLevelBFs[,1]
        }else{
          tmp<-matrix(NA,length(genes),1,
                      dimnames = list(genes,cellLines[i]))
          tmp[rownames(GeneLevelBFs),1]<-GeneLevelBFs[,1]-BFths[cellLines[i]]
          
        }
      }
      
      #if(i==1){
      # if(!exists("resMat")){
      #   resMat<-tmp
      # }else{
      #   if(nrow(tmp)<nrow(resMat)){
      #       missingG<-setdiff(rownames(resMat),rownames(tmp))
      #       toAdd<-matrix(NA,length(missingG),dimnames = list(missingG,NULL))
      #       tmp<-rbind(tmp,toAdd)
      #       tmp<-matrix(tmp[rownames(resMat),],nrow(resMat),dimnames = list(rownames(resMat),cellLines[i]))
      #   }
      #   resMat<-cbind(resMat,tmp)
      # }
      
      
    }
    colnames(resMat)<-cellLines
  }else{
    if(!is.null(BFmatrix)){
      resMat<-foreach(i=1:length(cellLines),.combine=cbind)%dopar%{
        tmp<-matrix(BFmatrix[,cellLines[i]]-BFths[cellLines[i]],nrow(BFmatrix),1,
                    dimnames = list(rownames(BFmatrix),cellLines[i]))
      }
    }
  }
  resMat<-resMat+0
  return(resMat)
}

FindPower<-function(inputset,integratedset,fdr=0.05){
  tissues<-intersect(names(inputset),names(integratedset))
  corTest<-list()
  criteria<-list()
    for(i in 1:length(tissues)){
      i1<-inputset[[tissues[i]]]
      i2<-integratedset[[tissues[i]]]
      bothtest<-intersect(rownames(i1),rownames(i2))
      #print(head(bothtest))
      corTest[[i]]<-cor.test(i1[bothtest,"delta"],i2[bothtest,"delta"])
      i1<-i1[bothtest,]
      i2<-i2[bothtest,]
      i1$fdr<-p.adjust(i1$p,method="fdr")
      i2$fdr<-p.adjust(i2$p,method="fdr")
      f1<-i1[sign(i1$delta)==sign(i2$delta),]
      f2<-i2[sign(i1$delta)==sign(i2$delta),]
      signInt<-f2[f2$fdr<fdr&f1$fdr>fdr,]
      signInt$nposInd<-f1[f2$fdr<fdr&f1$fdr>fdr,"npos"]
      signInt$nnegInd<-f1[f2$fdr<fdr&f1$fdr>fdr,"nneg"]
      signInt$delta<-f1[f2$fdr<fdr&f1$fdr>fdr,"delta"]
      if(nrow(signInt)>0){
      signInt$tissue<-tissues[i]}
      criteria[[i]]<-signInt
    }
  names(corTest)<-tissues
  names(criteria)<-tissues
  return(list(CorTest=corTest,criteria=criteria))
}

FindSignif<-function(inputset,integratedset,fdr=0.05){
  tissues<-intersect(names(inputset),names(integratedset))
  corTest<-list()
  criteria<-list()
  for(i in 1:length(tissues)){
    i1<-inputset[[tissues[i]]]
    i2<-integratedset[[tissues[i]]]
    bothtest<-intersect(rownames(i1),rownames(i2))
    #print(head(bothtest))
    corTest[[i]]<-cor.test(i1[bothtest,"delta"],i2[bothtest,"delta"])
    i1<-i1[bothtest,]
    i2<-i2[bothtest,]
    i1$fdr<-p.adjust(i1$p,method="fdr")
    i2$fdr<-p.adjust(i2$p,method="fdr")
    
    signInt<-i2[i2$fdr<fdr&i1$fdr>fdr,]
    signInt$nposInd<-i1[i2$fdr<fdr&i1$fdr>fdr,"npos"]
    signInt$nnegInd<-i1[i2$fdr<fdr&i1$fdr>fdr,"nneg"]
    signInt$delta<-i1[i2$fdr<fdr&i1$fdr>fdr,"delta"]
      if(nrow(signInt)>0){
        signInt$tissue<-tissues[i]}
    criteria[[i]]<-signInt
    }
  names(corTest)<-tissues
  names(criteria)<-tissues
  return(list(CorTest=corTest,criteria=criteria))
}

plotAssociationExamplesCP<-function(IndividData,IntegratedData,sel,MoBEM,cmp,AnovaData,pancan=FALSE) 
{
  gene<-AnovaData[sel,"GENE"]
  feat<-AnovaData[sel,"CFE"]
  if(pancan){
    cosmic_indivd<-intersect(colnames(IndividData),cmp[,"COSMIC_ID"])
    cosmic_int<-intersect(colnames(IntegratedData),cmp[,"COSMIC_ID"])
  }else{
    tissue<-AnovaData[sel,"tissue"]
    cosmic_indivd<-intersect(colnames(IndividData),cmp[cmp$tissue==tissue,"COSMIC_ID"])
    cosmic_int<-intersect(colnames(IntegratedData),cmp[cmp$tissue==tissue,"COSMIC_ID"])
  }
  if(length(cosmic_indivd)==0){
    colnames(IndividData)<-cmp[match(colnames(IndividData),cmp$model_id),"COSMIC_ID"]
    if(pancan){
      cosmic_indivd<-intersect(colnames(IndividData),cmp[,"COSMIC_ID"])
    }else{
      
      cosmic_indivd<-intersect(colnames(IndividData),cmp[cmp$tissue==tissue,"COSMIC_ID"])
    }
  }
  if(length(cosmic_int)==0){
    colnames(IntegratedData)<-cmp[match(colnames(IntegratedData),cmp$model_id),"COSMIC_ID"]
    if(pancan){
      cosmic_int<-intersect(colnames(IntegratedData),cmp[,"COSMIC_ID"])
    }else{
      cosmic_int<-intersect(colnames(IntegratedData),cmp[cmp$tissue==tissue,"COSMIC_ID"])
    }
  }
  cosmic_indivd<-intersect(cosmic_indivd,colnames(MoBEM))
  cosmic_int<-intersect(cosmic_int,colnames(MoBEM))
  cosmic_indivd<-cosmic_indivd[!cosmic_indivd==""]
  cosmic_int<-cosmic_int[!cosmic_int==""]
  essS <- IndividData[gene, cosmic_indivd]
  essB <- IntegratedData[gene, cosmic_int]
  pattern <- MoBEM[feat, as.character(cosmic_indivd)]
  pattern2 <- MoBEM[feat, as.character(cosmic_int)]
  r1<-range(essS)
  r2<-range(essB)
  ranges<-c(min(r1[1],r2[1]),max(r1[2],r2[2]))
  cols <- rep("darkgray", length(essS))
  cols[which(pattern == 1)] <- "darkgreen"
  beeswarm(essS ~ pattern, corral = "wrap", bg = c(makeTransparent("gray"), 
                                                   makeTransparent("darkgreen")), pch = 21, col = c("gray", 
                                                                                                    "darkgreen"), cex = 1.5, ylim = range(essS), las = 2, 
           labels = c("absent", "present"), xlab = feat, ylab = "",axes=F)
  par(new = TRUE)

  boxplot(essS ~ pattern, col = NA, ylim = ranges, frame.plot = FALSE, 
          xaxt = "n",  outline = FALSE,tcl=0.5,tck=-0.01)
  cols <- rep("darkgray", length(essB))
  cols[which(pattern2 == 1)] <- "darkgreen"
  beeswarm(essB ~ pattern2, corral = "wrap", bg = c(makeTransparent("gray"), 
                                                   makeTransparent("darkgreen")), pch = 21, col = c("gray", 
                                                                                                    "darkgreen"), cex = 1.5, ylim = range(essB), las = 2, 
           labels = c("absent", "present"), xlab = feat, ylab = "",axes=F)
  par(new = TRUE)
  boxplot(essB ~ pattern2, col = NA, ylim = ranges, frame.plot = FALSE, 
          xaxt = "n", outline = FALSE,tcl=0.5,tck=-0.01)
}
ADAM2.panessprofileC<-function (depMat, display = TRUE, main_suffix = "fitness genes in at least 1 cell line", 
          xlab = "n. dependent cell lines") 
{
  depMat <- depMat[which(rowSums(depMat,na.rm=T) > 0), ]
  panessprof <- rep(0, ncol(depMat))
  names(panessprof) <- as.character(1:ncol(depMat))
  paness <- summary(as.factor(rowSums(depMat,na.rm=T)), maxsum = length(unique(as.factor(rowSums(depMat,na.rm=T)))))
  panessprof[as.character(names(paness))] <- paness
  CUMsums <- rev(cumsum(rev(panessprof)))
  names(CUMsums) <- paste(">=", names(CUMsums), sep = "")
  if (display) {
    par(mfrow = c(2, 1))
    par(mar = c(6, 4, 4, 1))
    main = paste(nrow(depMat), main_suffix)
    barplot(panessprof, ylab = "n.genes", xlab = xlab, cex.axis = 0.8, 
            cex.names = 0.8, las = 2, main = main)
    barplot(CUMsums, ylab = "n.genes", xlab = xlab, cex.axis = 0.8, 
            cex.names = 0.6, las = 2, main = "Cumulative sums")
  }
  return(list(panessprof = panessprof, CUMsums = CUMsums))
}
ADAM2.generateNullModelC<-function (depMat, ntrials = 1000, display = TRUE) 
{
  set.seed(100812)
  depMat <- depMat[which(rowSums(depMat,na.rm=T) > 0), ]
  nullProf <- matrix(NA, ntrials, ncol(depMat), dimnames = list(1:ntrials, 
                                                                1:ncol(depMat)))
  nullCumSUM <- matrix(NA, ntrials, ncol(depMat), dimnames = list(1:ntrials, 
                                                                  paste("≥", 1:ncol(depMat), sep = "")))
  print("Generating null model...")
  pb <- txtProgressBar(min = 1, max = ntrials, style = 3)
  for (i in 1:ntrials) {
    setTxtProgressBar(pb, i)
    rMat <- ADAM2.randomisedepMat(depMat)
    Ret <- ADAM2.panessprofileC(rMat, display = FALSE)
    nullProf[i, ] <- Ret$panessprof
    nullCumSUM[i, ] <- Ret$CUMsums
  }
  Sys.sleep(1)
  close(pb)
  print("")
  print("Done")
  if (display) {
    par(mfrow = c(2, 1))
    main = c(paste(ntrials, " randomised essentiality profiles of\n", 
                   nrow(depMat), " genes across ", ncol(depMat), " cell lines", 
                   sep = ""))
    boxplot(nullProf, las = 2, xlab = "n. cell lines", ylab = "genes depleted in n cell lines", 
            main = main)
    colnames(nullCumSUM) <- paste(">=", 1:ncol(nullCumSUM))
    boxplot(log10(nullCumSUM + 1), las = 2, main = "Cumulative sums", 
            xlab = "n. cell lines", ylab = "log10 [number of genes + 1]", 
            cex.axis = 0.8)
  }
  return(list(nullProf = nullProf, nullCumSUM = nullCumSUM))
}
ADAMwrapper<-function(BinaryMat){
  #changed to ADAM2.panessprofileC locally function to handle NA values
  pprofile<-ADAM2.panessprofileC(depMat=BinaryMat)
  
  # Generate a set of random profiles of number of genes depleted for a number of cell lines and corresponding
  # cumulative sums by perturbing observed data.
  #changed to ADAM2.generateNullModelC locally function to handle NA values
  nullmodel<-ADAM2.generateNullModelC(depMat=BinaryMat,ntrials = 1000)
  
  
  # Calculate log10 odd ratios of observed/expected profiles of cumulative number of fitness genes in fixed number of cell lines
  # Observed values are from the ADAM.panessprofile function and expected are the average of random set from ADAM.generateNullModle
  EO<-ADAM2.empiricalOdds(observedCumSum = pprofile$CUMsums,simulatedCumSum =nullmodel$nullCumSUM )
  
  # Calculate True positive rates for fitness genes in at least n cell lines in the observed dependency matrix,
  # with positive cases from a reference set of essential genes
  TPR<-ADAM2.truePositiveRate(BinaryMat,curated_BAGEL_essential)
  
  
  # Calculate minimum number of cell lines a gene needs to be a fitness gene in order to be considered
  # as a core-fitness gene
  crossoverpoint<-ADAM2.tradeoffEO_TPR(EO,TPR$TPR,test_set_name = 'curated BAGEL essential')
  
  #coreFitnessGenes is the list of genes predicted as core-fitness by AdAM.
  coreFitnessGenes<-rownames(BinaryMat)[rowSums(BinaryMat,na.rm=T)>=crossoverpoint]
  return(coreFitnessGenes)
}

venn.diagramCP<-
function (x, filename, height = 300, width = 300, resolution = 300, 
          imagetype = "tiff", units = "px", compression = "lzw", na = "stop", 
          main = NULL, sub = NULL, main.pos = c(0.5, 1.05), main.fontface = "plain", 
          main.fontfamily = "serif", main.col = "black", main.cex = 1, 
          main.just = c(0.5, 1), sub.pos = c(0.5, 1.05), sub.fontface = "plain", 
          sub.fontfamily = "serif", sub.col = "black", sub.cex = 1, 
          sub.just = c(0.5, 1), category.names = names(x), force.unique = TRUE, 
          print.mode = "raw", sigdigs = 3, direct.area = FALSE, area.vector = 0, 
          hyper.test = FALSE, total.population = NULL, lower.tail = TRUE, 
          ...) 
{
  time.string = gsub(":", "-", gsub(" ", "_", as.character(Sys.time())))
  if (!is.null(filename)) {
    flog.appender(appender.file(paste0(filename, ".", time.string, 
                                       ".log")), name = "VennDiagramLogger")
  }
  else {
    flog.appender(appender.file(paste0("VennDiagram", time.string, 
                                       ".log")), name = "VennDiagramLogger")
  }
  out.list = as.list(sys.call())
  out.list[[1]] <- NULL
  out.string = capture.output(out.list)
  flog.info(out.string, name = "VennDiagramLogger")
  if (direct.area) {
    if (1 == length(area.vector)) {
      list.names <- category.names
      if (is.null(list.names)) {
        list.names <- ""
      }
      grob.list <- VennDiagram::draw.single.venn(area = area.vector[1], 
                                                 category = list.names, ind = FALSE, ...)
    }
    if (3 == length(area.vector)) {
      grob.list <- VennDiagram::draw.pairwise.venn(area1 = area.vector[1], 
                                                   area2 = area.vector[2], cross.area = area.vector[3], 
                                                   category = category.names, ind = FALSE, print.mode = print.mode, 
                                                   sigdigs = sigdigs, ...)
    }
    if (7 == length(area.vector)) {
      grob.list <- VennDiagram::draw.triple.venn(area1 = 0, 
                                                 area2 = 0, area3 = 0, n12 = 0, n23 = 0, n13 = 0, 
                                                 n123 = 0, category = category.names, ind = FALSE, 
                                                 list.order = 1:3, print.mode = print.mode, sigdigs = sigdigs, 
                                                 area.vector = area.vector, direct.area = TRUE, 
                                                 ...)
    }
    if (15 == length(area.vector)) {
      grob.list <- VennDiagram::draw.quad.venn(area1 = 0, 
                                               area2 = 0, area3 = 0, area4 = 0, n12 = 0, n13 = 0, 
                                               n14 = 0, n23 = 0, n24 = 0, n34 = 0, n123 = 0, 
                                               n124 = 0, n134 = 0, n234 = 0, n1234 = 0, category = category.names, 
                                               ind = FALSE, print.mode = print.mode, sigdigs = sigdigs, 
                                               area.vector = area.vector, direct.area = TRUE, 
                                               ...)
    }
    if (31 == length(area.vector)) {
      grob.list <- VennDiagram::draw.quintuple.venn(area1 = 0, 
                                                    area2 = 0, area3 = 0, area4 = 0, area5 = 0, n12 = 0, 
                                                    n13 = 0, n14 = 0, n15 = 0, n23 = 0, n24 = 0, 
                                                    n25 = 0, n34 = 0, n35 = 0, n45 = 0, n123 = 0, 
                                                    n124 = 0, n125 = 0, n134 = 0, n135 = 0, n145 = 0, 
                                                    n234 = 0, n235 = 0, n245 = 0, n345 = 0, n1234 = 0, 
                                                    n1235 = 0, n1245 = 0, n1345 = 0, n2345 = 0, n12345 = 0, 
                                                    category = category.names, ind = FALSE, print.mode = print.mode, 
                                                    sigdigs = sigdigs, area.vector = area.vector, 
                                                    direct.area = TRUE, ...)
    }
  }
  else {
    if (force.unique) {
      for (i in 1:length(x)) {
        x[[i]] <- unique(x[[i]])
      }
    }
    if ("none" == na) {
      x <- x
    }
    else if ("stop" == na) {
      for (i in 1:length(x)) {
        if (any(is.na(x[[i]]))) {
          flog.error("NAs in dataset", call. = FALSE, 
                     name = "VennDiagramLogger")
          stop("NAs in dataset", call. = FALSE)
        }
      }
    }
    else if ("remove" == na) {
      for (i in 1:length(x)) {
        x[[i]] <- x[[i]][!is.na(x[[i]])]
      }
    }
    else {
      flog.error("Invalid na option: valid options are \"none\", \"stop\", and \"remove\"", 
                 name = "VennDiagramLogger")
      stop("Invalid na option: valid options are \"none\", \"stop\", and \"remove\"")
    }
    if (0 == length(x) | length(x) > 5) {
      flog.error("Incorrect number of elements.", call. = FALSE, 
                 name = "VennDiagramLogger")
      stop("Incorrect number of elements.", call. = FALSE)
    }
    if (1 == length(x)) {
      list.names <- category.names
      if (is.null(list.names)) {
        list.names <- ""
      }
      grob.list <- VennDiagram::draw.single.venn(area = length(x[[1]]), 
                                                 category = list.names, ind = FALSE, ...)
    }
    else if (2 == length(x)) {
      grob.list <- VennDiagram::draw.pairwise.venn(area1 = length(x[[1]]), 
                                                   area2 = length(x[[2]]), cross.area = length(intersect(x[[1]], 
                                                                                                         x[[2]])), category = category.names, ind = FALSE, 
                                                   print.mode = print.mode, sigdigs = sigdigs, ...)
    }
    else if (3 == length(x)) {
      A <- x[[1]]
      B <- x[[2]]
      C <- x[[3]]
      list.names <- category.names
      nab <- intersect(A, B)
      nbc <- intersect(B, C)
      nac <- intersect(A, C)
      nabc <- intersect(nab, C)
      grob.list <- VennDiagram::draw.triple.venn(area1 = length(A), 
                                                 area2 = length(B), area3 = length(C), n12 = length(nab), 
                                                 n23 = length(nbc), n13 = length(nac), n123 = length(nabc), 
                                                 category = list.names, ind = FALSE, list.order = 1:3, 
                                                 print.mode = print.mode, sigdigs = sigdigs, ...)
    }
    else if (4 == length(x)) {
      A <- x[[1]]
      B <- x[[2]]
      C <- x[[3]]
      D <- x[[4]]
      list.names <- category.names
      n12 <- intersect(A, B)
      n13 <- intersect(A, C)
      n14 <- intersect(A, D)
      n23 <- intersect(B, C)
      n24 <- intersect(B, D)
      n34 <- intersect(C, D)
      n123 <- intersect(n12, C)
      n124 <- intersect(n12, D)
      n134 <- intersect(n13, D)
      n234 <- intersect(n23, D)
      n1234 <- intersect(n123, D)
      grob.list <- VennDiagram::draw.quad.venn(area1 = length(A), 
                                               area2 = length(B), area3 = length(C), area4 = length(D), 
                                               n12 = length(n12), n13 = length(n13), n14 = length(n14), 
                                               n23 = length(n23), n24 = length(n24), n34 = length(n34), 
                                               n123 = length(n123), n124 = length(n124), n134 = length(n134), 
                                               n234 = length(n234), n1234 = length(n1234), category = list.names, 
                                               ind = FALSE, print.mode = print.mode, sigdigs = sigdigs, 
                                               ...)
    }
    else if (5 == length(x)) {
      A <- x[[1]]
      B <- x[[2]]
      C <- x[[3]]
      D <- x[[4]]
      E <- x[[5]]
      list.names <- category.names
      n12 <- intersect(A, B)
      n13 <- intersect(A, C)
      n14 <- intersect(A, D)
      n15 <- intersect(A, E)
      n23 <- intersect(B, C)
      n24 <- intersect(B, D)
      n25 <- intersect(B, E)
      n34 <- intersect(C, D)
      n35 <- intersect(C, E)
      n45 <- intersect(D, E)
      n123 <- intersect(n12, C)
      n124 <- intersect(n12, D)
      n125 <- intersect(n12, E)
      n134 <- intersect(n13, D)
      n135 <- intersect(n13, E)
      n145 <- intersect(n14, E)
      n234 <- intersect(n23, D)
      n235 <- intersect(n23, E)
      n245 <- intersect(n24, E)
      n345 <- intersect(n34, E)
      n1234 <- intersect(n123, D)
      n1235 <- intersect(n123, E)
      n1245 <- intersect(n124, E)
      n1345 <- intersect(n134, E)
      n2345 <- intersect(n234, E)
      n12345 <- intersect(n1234, E)
      grob.list <- VennDiagram::draw.quintuple.venn(area1 = length(A), 
                                                    area2 = length(B), area3 = length(C), area4 = length(D), 
                                                    area5 = length(E), n12 = length(n12), n13 = length(n13), 
                                                    n14 = length(n14), n15 = length(n15), n23 = length(n23), 
                                                    n24 = length(n24), n25 = length(n25), n34 = length(n34), 
                                                    n35 = length(n35), n45 = length(n45), n123 = length(n123), 
                                                    n124 = length(n124), n125 = length(n125), n134 = length(n134), 
                                                    n135 = length(n135), n145 = length(n145), n234 = length(n234), 
                                                    n235 = length(n235), n245 = length(n245), n345 = length(n345), 
                                                    n1234 = length(n1234), n1235 = length(n1235), 
                                                    n1245 = length(n1245), n1345 = length(n1345), 
                                                    n2345 = length(n2345), n12345 = length(n12345), 
                                                    category = list.names, ind = FALSE, print.mode = print.mode, 
                                                    sigdigs = sigdigs, ...)
    }
    else {
      flog.error("Invalid size of input object", name = "VennDiagramLogger")
      stop("Invalid size of input object")
    }
  }


      pdf(file = filename, height = height, width = width)
  
    print(grid.draw(grob.list))
    dev.off()
    
  
  return(grob.list)
}

PowerPlot<-function(inputset,integratedset,fdr=0.05){
  tissues<-intersect(names(inputset),names(integratedset))
  values<-list()
  valuenames<-c()
  j=1
  for(i in 1:length(tissues)){

    i1<-inputset[[tissues[i]]]
    i2<-integratedset[[tissues[i]]]
    bothtest<-intersect(rownames(i1),rownames(i2))

    i1<-i1[bothtest,]
    i2<-i2[bothtest,]
    i1$fdr<-p.adjust(i1$p,method="fdr")
    i2$fdr<-p.adjust(i2$p,method="fdr")
    sel<-(which(i2$fdr<fdr&(i2$fdr<i1$fdr)))
    ii1<-i1[sel,]
    ii2<-i2[sel,]

    if(length(sel)>0){
      deltaIntS<-Sig_All[intersect(retained,rownames(Sig_Sanger)),"effect_size"]
      deltaS<-Sig_Sanger[intersect(retained,rownames(Sig_Sanger)),"effect_size"]
      deltaB<-Sig_Broad[intersect(retained,rownames(Sig_Broad)),"effect_size"]
      esdata<-data.frame(Integrated=c(deltaIntS,deltaIntS),Individual=c(deltaS,deltaB),site=c(rep("Sanger",length(deltaS)),rep("Broad",length(deltaB))))
      
      esplot<-ggplot(esdata,aes(x=Integrated,y=Individual,colour=site))+geom_point(size=3)+theme_bw()+ylab("Individual effect sizes")+xlab("Integrated effect sizes")+geom_abline(slope=1,intercept=0)+ theme(legend.position = c(0.1, 0.9))
      pdf(paste0(dir.Results,"/EffectSizes_IntvIndiv_CCR.pdf"))
      print(esplot)
      dev.off()
      plot(log(ii1$fdr-ii2$fdr),ii1$delta-ii2$delta,main=tissues[i])

      temp<-cbind(ii2,log(ii1$fdr-ii2$fdr),ii1$delta-ii2$delta,(ii1$delta-ii2$delta)/ii1$delta,ii1$fdr,tissues[i])
      rownames(temp)<-bothtest[sel]
      colnames(temp)<-c(colnames(ii2),"diff_fdr","diff_delta","percent_diff_delta","orig_fdr","tissue")
      values[[j]]<-temp
      valuenames<-c(valuenames,tissues[i])
      #do tests:
      #wilcox and return median values
      j=j+1
    }

  }

  names(values)<-valuenames
  return(values)
}

negativePvalue<-function(FCmat,nonessentials){
  noness<-intersect(rownames(FCmat),nonessentials)
  nonessMat<-FCmat[noness,]
  meanvecs<-colMeans(nonessMat)
  sdvecs<-sqrt(colVars(nonessMat))
  pvalmat<-sapply(1:ncol(FCmat),function(x) pnorm(FCmat[,x],mean=meanvecs[x],sd=sdvecs[x]))
  return(pvalmat)
  
}

GSEAfunction<-function(rankedlist,geneset){
  guse<-intersect(rankedlist,geneset)
  Inset<-(rankedlist%in%guse)*1
  sizeset<-length(guse)
  N_r<-sizeset
  PosScores<-Inset/N_r
  NegScores<-(!rankedlist%in%guse)*1
  N_Nh<-length(rankedlist)-sizeset
  NegScores<-NegScores/N_Nh
  Allscores<-PosScores-NegScores
  RunningSum<-cumsum(Allscores)
  maxdev<-RunningSum[which.max(abs(RunningSum))]
  return(list(ESscore=maxdev,RunningSum=RunningSum))

}


fdrCurve<-function(ranklist,controlset,featurepoint=1000){
  npoints<-length(ranklist)
  ranklist<-sort(ranklist,decreasing = FALSE)
  checkset<-seq(from=1,to=npoints,length.out = featurepoint)
  fprs<-sapply(checkset,function(x) fpr(names(ranklist[1:x]),controlset))

  return(unlist(fprs))
}

fdrCurveP<-function(ranklist,controlset,featurepoint=1000){
  npoints<-length(ranklist)
  ranklist<-sort(ranklist,decreasing = FALSE)
  checkset<-seq(from=1,to=npoints,length.out = featurepoint)

  propfp<-sapply(checkset,function(x) fr(names(ranklist[1:x]),controlset))
  return(unlist(propfp))
}

fpr<-function(testgenes,controlgenes){
  sum(controlgenes%in%testgenes)/length(testgenes)
}

fr<-function(testgenes,controlgenes){
  sum(controlgenes%in%testgenes)/length(controlgenes)
}
#https://github.com/cancerdatasci/ceres/blob/master/R/scale_to_essentials.R
scale_to_essentials <- function(ge_fit){
  
  
  essential_indices <- which(row.names(ge_fit) %in% ceres::hart_essentials[["Gene"]])
  nonessential_indices <- which(row.names(ge_fit) %in% ceres::hart_nonessentials[["Gene"]])
  
  scaled_ge_fit <- ge_fit %>%
    apply(2, function(x){
      
      (x - median(x[nonessential_indices], na.rm=T)) %>%
        divide_by(median(x[nonessential_indices], na.rm=T) - median(x[essential_indices], na.rm=T))
      
    })
  
  return(scaled_ge_fit)
  
}

FeatureTissue<-function(feature,MoBEM,CL1,CL2,annot){
  outDF<-NULL
  for(i in 1:length(feature)){
    
    set1<-MoBEM[feature[i],intersect(CL1,colnames(MoBEM))]
    set1<-names(set1)[set1==1]
    set2<-MoBEM[feature[i],intersect(CL2,colnames(MoBEM))]
    set2<-names(set2)[set2==1]
    tissue1<-unique(cmp[cmp$COSMIC_ID%in%set1,"tissue"])
    tissue2<-unique(cmp[cmp$COSMIC_ID%in%set2,"tissue"])
    lc<-identical(tissue1,tissue2)
    temp<-data.frame(feature=feature[i],tissue1=paste(tissue1,collapse=","),tissue2=paste(tissue2,collapse = ","),noFeat1=length(set1),noFeat2=length(set2),lineagecheck=lc)
    outDF<-rbind(outDF,temp)
  }
  return(outDF)
}


FindNewSignif<-function(inputset,integratedset,fdr=0.05){
  tissues<-intersect(names(inputset),names(integratedset))
  corTest<-list()
  criteria<-list()
  for(i in 1:length(tissues)){
    i1<-inputset[[tissues[i]]]
    i2<-integratedset[[tissues[i]]]
    newtest<-setdiff(rownames(i2),rownames(i1))
    #print(head(bothtest))


    i2$fdr<-p.adjust(i2$p,method="fdr")
    intNew<-i2[newtest,]
    
    signInt<-intNew[intNew$fdr<fdr,]

    if(nrow(signInt)>0){
      signInt$tissue<-tissues[i]}
    criteria[[i]]<-signInt
  }

  names(criteria)<-tissues
  return(criteria=criteria)
}


plotSingleAssociation<-function(Data,sel,MoBEM,cmp,AnovaData) 
{
  gene<-AnovaData[sel,"GENE"]
  feat<-AnovaData[sel,"CFE"]
  tissue<-AnovaData[sel,"tissue"]
  cosmic_indivd<-intersect(colnames(Data),cmp[cmp$tissue==tissue,"COSMIC_ID"])

  cosmic_indivd<-intersect(cosmic_indivd,colnames(MoBEM))

  cosmic_indivd<-cosmic_indivd[!cosmic_indivd==""]

  essS <- Data[gene, cosmic_indivd]

  pattern <- MoBEM[feat, as.character(cosmic_indivd)]

  
  cols <- rep("darkgray", length(essS))
  cols[which(pattern == 1)] <- "darkgreen"
  beeswarm(essS ~ pattern, corral = "wrap", bg = c(makeTransparent("gray"), 
                                                   makeTransparent("darkgreen")), pch = 21, col = c("gray", 
                                                                                                    "darkgreen"), cex = 1.5, ylim = range(essS), las = 2, 
           labels = c("absent", "present"), xlab = feat, ylab = "",axes=F)
  par(new = TRUE)
  
  boxplot(essS ~ pattern, col = NA, ylim = range(essS), frame.plot = FALSE, 
          xaxt = "n",  outline = FALSE,tcl=0.5,tck=-0.01)
 
}


SingleFeatureAnova<-function (feature, dataset,fname) {
  usecl<-intersect(names(feature),colnames(dataset))
  feature<-feature[usecl]
  if(length(usecl)>5){
    c1<-sum(feature)>2
    c2<-sum(feature==0)>2

    if(c1&c2){
      res<-NULL
      subdataset <- dataset[, usecl]
      cfePattern <- feature
        subSetWT <- subdataset[, cfePattern == 0]
        subSetALT <- subdataset[, cfePattern == 1]
        effect_size <- unlist(lapply(seq_len(nrow(subSetALT)), function(i) {
          cohens_d(subSetALT[i, ], subSetWT[i, ])
        }))
        delta <- unlist(lapply(seq_len(nrow(subSetALT)), function(i) {
          mean(subSetALT[i, ]) - mean(subSetWT[i, ])
        }))
        pval <- unlist(lapply(seq_len(nrow(subSetALT)), function(i) {
          RES <- t.test(subSetALT[i, ], subSetWT[i, ], var.equal = TRUE)
          RES$p.value
        }))
        npos=sum(cfePattern)
        nneg=length(cfePattern)-npos
        fdr=p.adjust(pval,method="fdr")
        temp<- data.frame(CFE = fname, GENE = rownames(dataset), delta = delta, 
                          effect_size = effect_size, p = pval, npos=npos,nneg=nneg,fdr=fdr ,stringsAsFactors = FALSE)
        res<-rbind(res,temp)
        return(res)
      }else{
        #no features to run analysis with
        return(NULL)
      }
  }else{
    #not enough cell lines with molecular data to run analysis
    return(NULL)
  }
}

GeneCor<-function(expressiondata,depdata){
  genes<-intersect(rownames(expressiondata),rownames(depdata))
  expressiondata<-expressiondata[genes,]
  depdata<-depdata[genes,]
  ngenes<-nrow(expressiondata)
  gcor<-c()
  gcorp<-c()
  for(i in 1:ngenes){
    temp<-cor.test(expressiondata[i,],depdata[i,])
    gcor[i]<-temp$estimate
    gcorp[i]<-temp$p.value
  }
  names(gcor)<-genes
  names(gcorp)<-genes
  keep<-names(gcor)[gcor<0]
  gcor<-gcor[keep]
  gcorp<-gcorp[keep]
  mna<-!is.na(gcor)
  gcor<-gcor[mna]
  gcorp<-gcorp[mna]
  fdr<-p.adjust(gcorp,method="fdr")
  #mna<-!is.na(fdr)
  #gcor[mna]
  #gcorp[mna]
  #fdr<-fdr[mna]
  return(data.frame(estimate=gcor,pvalues=gcorp,fdr=fdr))
}

BinaryVenn<-function(BinaryCCR,BinaryCERES,BinaryCCRJ,CCRcorrected,CEREScorrected,CCRJcorrected,dir.Results,output){
  DepletedonceCCRJ<-rownames(CCRJ_corrected)[rowSums(BinaryCCRJ)!=0]
  DepletedonceCCR<-rownames(CCR_corrected)[rowSums(BinaryCCR)!=0]
  DepletedonceCERES<-rownames(CERES_corrected)[rowSums(BinaryCERES)!=0]
  
  dCCRJ<-length(DepletedonceCCRJ)
  dCCR<-length(DepletedonceCCR)
  dCERES<-length(DepletedonceCERES)
  
  cat(paste("CCR number genes depleted:",output,dCCR),file=paste0(dir.Results,"/Figure4Log.txt"),sep="\n",append=TRUE)
  cat(paste("CERES number genes depleted:",output,dCERES),file=paste0(dir.Results,"/Figure4Log.txt"),sep="\n",append=TRUE)
  cat(paste("CCR_JACKS number genes depleted:",output,dCCRJ),file=paste0(dir.Results,"/Figure4Log.txt"),sep="\n",append=TRUE)
  
  allBinary<-BinaryCCR+BinaryCCRJ+BinaryCERES
  all<-sum(allBinary==3,na.rm = T)
  sum(BinaryCCR,na.rm = T)
  sum(BinaryCCRJ,na.rm = T)
  sum(BinaryCERES,na.rm = T)
  CCR_CERES<-BinaryCCR+BinaryCERES
  CCR_CCRJ<-BinaryCCR+BinaryCCRJ
  CERES_CCRJ<-BinaryCERES+BinaryCCRJ
  c1_3<-sum(CCR_CERES==2,na.rm = T)-sum(allBinary==3,na.rm = T)
  c1_2<-sum(CCR_CCRJ==2,na.rm = T)-sum(allBinary==3,na.rm = T)
  c2_3<-sum(CERES_CCRJ==2,na.rm = T)-sum(allBinary==3,na.rm = T)
  c1<-sum(BinaryCCR4,na.rm = T)-(sum(CCR_CERES==2,na.rm = T)-sum(allBinary==3,na.rm = T))-(sum(CCR_CCRJ==2,na.rm = T)-sum(allBinary==3,na.rm = T))-sum(allBinary==3,na.rm = T)
  c2<-sum(BinaryCCRJ4,na.rm = T)-(sum(CCR_CCRJ==2,na.rm = T)-sum(allBinary==3,na.rm = T))-(sum(CERES_CCRJ==2,na.rm = T)-sum(allBinary==3,na.rm = T))-sum(allBinary==3,na.rm = T)
  c3<-sum(BinaryCERES4,na.rm = T)-(sum(CCR_CERES==2,na.rm = T)-sum(allBinary==3,na.rm = T))-(sum(CERES_CCRJ==2,na.rm = T)-sum(allBinary==3,na.rm = T))-sum(allBinary==3,na.rm = T)
  
  venn<-euler(c("A&B&C"=all,A=c1,B=c2,C=c3,"A&B"=c1_2,"A&C"=c1_3,
                "B&C"=c2_3))
  return(venn)
}

DepletionPlots<-function(BinaryCCR,BinaryCERES,BinaryCCRJ,dir.Results,output){
  Depletedccr<-sum(rowSums(BinaryCCR,na.rm = T)>0,na.rm = T)
  Depletedceres<-sum(rowSums(BinaryCERES,na.rm = T)>0,na.rm = T)
  DepletedccrJ<-sum(rowSums(BinaryCCRJ,na.rm = T)>0,na.rm = T)

  
  mCCR<-median(colSums(BinaryCCR,na.rm = T))
  mCERES<-median(colSums(BinaryCERES,na.rm = T))
  mCCRJ<-median(colSums(BinaryCCRJ,na.rm = T))
  
  cat(paste("CCR median number depletions per CL:",output,mCCR),file=paste0(dir.Results,"/Figure4Log.txt"),sep="\n",append=TRUE)
  cat(paste("CERES median number depletions per CL:",output,mCERES),file=paste0(dir.Results,"/Figure4Log.txt"),sep="\n",append=TRUE)
  cat(paste("CCR_JACKS median number depletions per CL:",output,mCCRJ),file=paste0(dir.Results,"/Figure4Log.txt"),sep="\n",append=TRUE)
  
  
  plotnt1<-rbind(data.frame(x="median",v=mCCR,method="CCR"),data.frame(x="median",v=mCERES,method="CERES"),data.frame(x="median",v=mCCRJ,method="CCRJ"))
  plotnt2<-rbind(data.frame(x="Dep1",v=dCCR,method="CCR"),data.frame(x="Dep1",v=dCERES,method="CERES"),data.frame(x="Dep1",v=dCCRJ,method="CCRJ"))
  
  normdepplot<-ggplot(plotnt1,aes(x=x,y=v,fill=method))+geom_bar(stat="identity",position="dodge")+scale_fill_manual(values=PipelineColours)+theme(axis.text.x=element_text(angle=45))+ylab("Number Genes")+xlab("Norm LRT")+theme_bw()+theme(legend.position = c(0.8,0.8),legend.background = element_blank())
  normdepplot2<-ggplot(plotnt2,aes(x=x,y=v,fill=method))+geom_bar(stat="identity",position="dodge")+scale_fill_manual(values=PipelineColours)+theme(axis.text.x=element_text(angle=45))+ylab("Number Genes")+xlab("Norm LRT")+theme_bw()+theme(legend.position = c(0.8,0.8),legend.background = element_blank())
  return(list(plot1=normdepplot,plot2=normdepplot2))
}

ProportionDepletions<-function(BinaryCCR,BinaryCERES,BinaryCCRJ,output){
  
  bothvCCR<-sum(BinaryCCR[,colnames(BinaryCERES)]*BinaryCERES,na.rm = T)/sum(BinaryCCR[,colnames(BinaryCERES)],na.rm = T)
  bothvCERES<-sum(BinaryCCR[,colnames(BinaryCERES)]*BinaryCERES,na.rm = T)/sum(BinaryCERES,na.rm = T)
  ccrJvCCR<-sum(BinaryCCR[,colnames(BinaryCCRJ)]*BinaryCCRJ,na.rm = T)/sum(BinaryCCR[,colnames(BinaryCCRJ)],na.rm = T)
  ccrvCCRJ<-sum(BinaryCCR4[,colnames(BinaryCCRJ)]*BinaryCCRJ,na.rm = T)/sum(BinaryCCRJ,na.rm = T)
  ccrJvCERES<-sum(BinaryCCRJ[,colnames(BinaryCERES)]*BinaryCERES,na.rm = T)/sum(BinaryCERES,na.rm = T)
  ceresvCCRJ<-sum(BinaryCCRJ[,colnames(BinaryCERES)]*BinaryCERES,na.rm = T)/sum(BinaryCCRJ[,colnames(BinaryCERES)],na.rm = T)
  
  
  
  cat(paste("Prop CCR found in CERES:",output,bothvCCR),file=paste0(dir.Results,"/Figure4Log.txt"),sep="\n",append=TRUE)
  cat(paste("Prop CERES found in CCR:",output,bothvCERES),file=paste0(dir.Results,"/Figure4Log.txt"),sep="\n",append=TRUE)
  cat(paste("Prop CCR_JACKS found in CCR:",output,ccrvCCRJ),file=paste0(dir.Results,"/Figure4Log.txt"),sep="\n",append=TRUE)
  cat(paste("Prop CCR found in CCR_JACKS:",output,ccrJvCCR),file=paste0(dir.Results,"/Figure4Log.txt"),sep="\n",append=TRUE)
  cat(paste("Prop CERES found in CCR_JACKS:",output,ccrJvCERES),file=paste0(dir.Results,"/Figure4Log.txt"),sep="\n",append=TRUE)
  cat(paste("Prop CCR_JACKS found in CERES:",output,ccrJvCCR),file=paste0(dir.Results,"/Figure4Log.txt"),sep="\n",append=TRUE)
  
  
}

AovFactors<-function(InputDF){
  AovRes<-aov(data~Batch+PreProc)
  
}
cohens_d <- function(x, y) {
  lx <- length(x)- 1
  ly <- length(y)- 1
  
  md  <- abs(mean(x) - mean(y))        ## mean difference (numerator)
  csd <- lx * var(x) + ly * var(y)
  csd <- csd/(lx + ly)
  csd <- sqrt(csd)                     ## common sd computation
  
  cd  <- md/csd                        ## cohen's d
  
  return(cd)
}
