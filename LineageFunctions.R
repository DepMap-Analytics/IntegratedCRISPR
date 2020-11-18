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
# Parameters
global <- list(
  n_genes = 'all', # set to 'all' to use all protein coding genes found in both datasets 
  umap_n_neighbors = 10, # num nearest neighbors used to create UMAP plot
  umap_min_dist = 0.5, # min distance used to create UMAP plot
  mnn_k_CL = 5, # number of nearest neighbors of tumors in the cell line data
  mnn_k_tumor = 50, # number of nearest neighbors of cell lines in the tumor data
  top_DE_genes_per = 1000, # differentially expressed genes with a rank better than this is in the cell line or tumor data
  # are used to identify mutual nearest neighbors in the MNN alignment step
  remove_cPCA_dims = c(1,2,3,4), # which cPCA dimensions to regress out of the data 
  distance_metric = 'euclidean', # distance metric used for the UMAP projection
  mod_clust_res = 5, # resolution parameter used for clustering the data
  mnn_ndist = 3, # ndist parameter used for MNN
  n_PC_dims = 70, # number of PCs to use for dimensionality reduction
  reduction.use = 'umap', # 2D projection used for plotting
  fast_cPCA = 10 # to run fast cPCA (approximate the cPCA eigenvectors instead of calculating all) set this to a value >= 4
)

create_Seurat_object <- function(exp_mat, ann, type = NULL) {
  #exp_mat should be genes x cell lines
  seu_obj <- Seurat::CreateSeuratObject(exp_mat,
                                        min.cells = 0,
                                        min.features = 0,
                                        meta.data = ann %>%
                                          magrittr::set_rownames(ann$sampleID))
  if(!is.null(type)) {
    seu_obj@meta.data$type <- type
  }
  # mean center the data, important for PCA
  seu_obj <- Seurat::ScaleData(seu_obj, features = rownames(Seurat::GetAssayData(seu_obj)), do.scale = F)
  
  seu_obj %<>% Seurat::RunPCA(assay='RNA',
                              features = rownames(Seurat::GetAssayData(seu_obj)),
                              npcs = global$n_PC_dims, verbose = F)
  
  seu_obj %<>% Seurat::RunUMAP(assay = 'RNA', dims = 1:global$n_PC_dims,
                               reduction = 'pca',
                               n.neighbors = global$umap_n_neighbors,
                               min.dist =  global$umap_min_dist,
                               metric = global$distance_metric, verbose=F)
  umapcoords<-FetchData(seu_obj,vars=c("UMAP_1","UMAP_2"))
  return(umapcoords)
}
LineageUMAP<-function(){
  CIoutput<-tCI[,TCGAbreast$sampleID]
  #assume have data want to use in CIoutput
  umapData<-create_Seurat_object(CIoutput,CIannot[CIannot$sampleID%in%TCGAbreast$sampleID,])
  
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
global <- list(
  n_genes = 'all', # set to 'all' to use all protein coding genes found in both datasets 
  umap_n_neighbors = 10, # num nearest neighbors used to create UMAP plot
  umap_min_dist = 0.5, # min distance used to create UMAP plot
  mnn_k_CL = 5, # number of nearest neighbors of tumors in the cell line data
  mnn_k_tumor = 50, # number of nearest neighbors of cell lines in the tumor data
  top_DE_genes_per = 1000, # differentially expressed genes with a rank better than this is in the cell line or tumor data
  # are used to identify mutual nearest neighbors in the MNN alignment step
  remove_cPCA_dims = c(1,2,3,4), # which cPCA dimensions to regress out of the data 
  distance_metric = 'euclidean', # distance metric used for the UMAP projection
  mod_clust_res = 5, # resolution parameter used for clustering the data
  mnn_ndist = 3, # ndist parameter used for MNN
  n_PC_dims = 50, # number of PCs to use for dimensionality reduction
  reduction.use = 'umap', # 2D projection used for plotting
  fast_cPCA = 20 # to run fast cPCA (approximate the cPCA eigenvectors instead of calculating all) set this to a value >= 4
)
create_Seurat_object <- function(exp_mat, ann, type = NULL) {
  #exp_mat should be genes x cell lines
  seu_obj <- Seurat::CreateSeuratObject(exp_mat,
                                        min.cells = 0,
                                        min.features = 0,
                                        meta.data = ann %>%
                                          magrittr::set_rownames(ann$model_id))
  if(!is.null(type)) {
    seu_obj@meta.data$type <- type
  }
  # mean center the data, important for PCA
  seu_obj <- Seurat::ScaleData(seu_obj, features = rownames(Seurat::GetAssayData(seu_obj)), do.scale = F)
  
  seu_obj %<>% Seurat::RunPCA(assay='RNA',
                              features = rownames(Seurat::GetAssayData(seu_obj)),
                              npcs = global$n_PC_dims, verbose = F)
  
  seu_obj %<>% Seurat::RunUMAP(assay = 'RNA', dims = 1:global$n_PC_dims,
                               reduction = 'pca',
                               n.neighbors = global$umap_n_neighbors,
                               min.dist =  global$umap_min_dist,
                               metric = global$distance_metric, verbose=F)
  umapcoords<-FetchData(seu_obj,vars=c("UMAP_1","UMAP_2"))
  return(umapcoords)
}

