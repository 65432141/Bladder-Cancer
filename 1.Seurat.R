  setwd("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2")
  library(dplyr)
  library(Seurat)
  library(SeuratObject)
  library(patchwork)
  library(ggplot2)
  library(data.table)
  library(harmony)
  library(stringr)
  library(RhpcBLASctl)
  library(utils)
  library(ggalluvial)
  library(tidyr)
  library(RColorBrewer)

{ 
  dir_name=list.files("sample_17样本/")
  dir_name
  T0_file<-dir_name[grep("T0",dir_name)]
  T1_file<-dir_name[grep("T1",dir_name)]
  Ta_file<-dir_name[grep("Ta",dir_name)]
  T2_file<-dir_name[grep("T2",dir_name)]
  T3_file<-dir_name[grep("T3",dir_name)]
  T4_file<-dir_name[grep("T4",dir_name)]
  
  T0_scRNAlist <- list()
    for(i in 1:length(T0_file)){
    counts <- Read10X(data.dir = paste("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4/sample_17样本/",T0_file[i], sep = ""))
    T0_scRNAlist[[i]] <- CreateSeuratObject(counts, project = T0_file[i],min.cells = 3, min.features = 300)
    print(i)
    }
    length(T0_scRNAlist)
  T0_scRNAlist_2 <- merge(T0_scRNAlist[[1]], T0_scRNAlist[2:length(T0_scRNAlist)])

  T1_scRNAlist <- list()
    for(i in 1:length(T1_file)){
    counts <- Read10X(data.dir = paste("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4/sample_17样本/",T1_file[i], sep = ""))
    T1_scRNAlist[[i]] <- CreateSeuratObject(counts, project = T1_file[i],min.cells = 3, min.features = 300)
    print(i)
    }
    length(T1_scRNAlist)
  T1_scRNAlist_2 <- merge(T1_scRNAlist[[1]], T1_scRNAlist[2:length(T1_scRNAlist)])
   
  Ta_scRNAlist <- list()
    for(i in 1:length(Ta_file)){
    counts <- Read10X(data.dir = paste("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/sample_17样本/",Ta_file[i], sep = ""))
    Ta_scRNAlist <- CreateSeuratObject(counts, project = Ta_file[i],min.cells = 3, min.features = 300)
    print(i)
    }
    class(Ta_scRNAlist)
    Ta_scRNAlist_2 <- Ta_scRNAlist
  
  T2_scRNAlist <- list()
    for(i in 1:length(T2_file)){
    counts <- Read10X(data.dir = paste("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4/sample_17样本/",T2_file[i], sep = ""))
    T2_scRNAlist[[i]] <- CreateSeuratObject(counts, project = T2_file[i],min.cells = 3, min.features = 300)
    print(i)
    }
    length(T2_scRNAlist)
  T2_scRNAlist_2 <- merge(T2_scRNAlist[[1]], T2_scRNAlist[2:length(T2_scRNAlist)])
   
   T3_scRNAlist <- list()
    for(i in 1:length(T3_file)){
    counts <- Read10X(data.dir = paste("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4/sample_17样本/",T3_file[i], sep = ""))
    T3_scRNAlist[[i]] <- CreateSeuratObject(counts, project = T3_file[i],min.cells = 3, min.features = 300)
    print(i)
    }
  T3_scRNAlist_2 <- merge(T3_scRNAlist[[1]], T3_scRNAlist[2:length(T3_scRNAlist)])

  T4_scRNAlist <- list()
    for(i in 1:length(T4_file)){
    counts <- Read10X(data.dir = paste("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4/sample_17样本/",T4_file[i], sep = ""))
    T4_scRNAlist <- CreateSeuratObject(counts, project = T4_file[i],min.cells = 3, min.features = 300)
    print(i)
    }
    
}
{ 
  dir_name=list.files("GSE135337_8样本/GSE135337_8样本")
  dir_name
  Ta_GSE135337<-dir_name[grep("Ta",dir_name)]
  T1_GSE135337<-dir_name[grep("T1",dir_name)]
  T2_GSE135337<-dir_name[grep("T2",dir_name)]
  T3_GSE135337<-dir_name[grep("T3",dir_name)]
  path="GSE135337_8样本/GSE135337_8样本/"

  samples=Ta_GSE135337
  samples
  dir <- file.path(path,samples)
  names(dir) <- samples
  options(scipen = 200)
  Ta_GSE135337list <- list()
  for(i in 1:length(dir)){
    counts=fread(dir[i],sep="\t",header=T,check.names=F,quote = F)
    counts=as.data.frame(counts)
    counts<-counts[!duplicated(counts[,2]),]
    rownames(counts) <- counts[,2]
    counts=counts[,-(1:2)]
    Ta_GSE135337list[[i]] <- CreateSeuratObject(counts = counts,project = "Ta_GSE135337", min.cells=3, min.features=200, names.delim = "_")
    print(i)
  }
  length(Ta_GSE135337list)
  Ta_GSE135337_2 <-merge(Ta_GSE135337list[[1]], 
                          y = c(Ta_GSE135337list[2:length(Ta_GSE135337list)]),
                          project = "Ta_GSE135337list")
  dim(Ta_GSE135337_2)
  table(Ta_GSE135337_2@meta.data$orig.ident)

   samples=T1_GSE135337
  samples
  dir <- file.path(path,samples)
  names(dir) <- samples
  options(scipen = 200)
  T1_GSE135337list <- list()
  for(i in 1:length(dir)){
    counts=fread(dir[i],sep="\t",header=T,check.names=F,quote = F)
    counts=as.data.frame(counts)
    counts<-counts[!duplicated(counts[,2]),]
    rownames(counts) <- counts[,2]
    counts=counts[,-(1:2)]
    T1_GSE135337list[[i]] <- CreateSeuratObject(counts = counts,project = "T1_GSE135337", min.cells=3, min.features=200, names.delim = "_")
    print(i)
  }
  length(T1_GSE135337list)
  T1_GSE135337_2 <-merge(T1_GSE135337list[[1]], 
                          y = c(T1_GSE135337list[2:length(T1_GSE135337list)]),
                          project = "T1_GSE135337list")
  dim(T1_GSE135337_2)
   
   samples=T2_GSE135337
  samples
  dir <- file.path(path,samples)
  names(dir) <- samples
  options(scipen = 200)
  T2_GSE135337list <- list()
  for(i in 1:length(dir)){
    counts=fread(dir[i],sep="\t",header=T,check.names=F,quote = F)
    counts=as.data.frame(counts)
    counts<-counts[!duplicated(counts[,2]),]
    rownames(counts) <- counts[,2]
    counts=counts[,-(1:2)]
    T2_GSE135337list <- CreateSeuratObject(counts = counts,project = "T2_GSE135337", min.cells=3, min.features=200, names.delim = "_")
    print(i)
  }
  class(T2_GSE135337list)
    dim(T2_GSE135337list)

  samples=T3_GSE135337
  samples
  dir <- file.path(path,samples)
  names(dir) <- samples
  options(scipen = 200)
  T3_GSE135337list <- list()
  for(i in 1:length(dir)){
    counts=fread(dir[i],sep="\t",header=T,check.names=F,quote = F)
    counts=as.data.frame(counts)
    counts<-counts[!duplicated(counts[,2]),]
    rownames(counts) <- counts[,2]
    counts=counts[,-(1:2)]
    T3_GSE135337list[[i]] <- CreateSeuratObject(counts = counts,project = "T3_GSE135337", min.cells=3, min.features=200, names.delim = "_")
    print(i)
  }
  length(T3_GSE135337list)
  T3_GSE135337_2 <-merge(T3_GSE135337list[[1]], 
                          y = c(T1_GSE135337list[2:length(T3_GSE135337list)]),
                          project = "T3_GSE135337list")
  dim(T3_GSE135337_2)



}

{

    Ta <- merge(Ta_GSE135337_2, 
                      y = c(Ta_scRNAlist_2), 
                      add.cell.ids = c("Ta_GSE135337_2","Ta_scRNAlist_2"), 
                      project = "Ta")
    rownames(Ta)
    unique(sapply(X = strsplit(colnames(Ta), split = "_"), FUN = "[", 1))
    table(Ta$orig.ident)
    saveRDS(Ta, file = "Ta_object_v4.rds")

    T1 <- merge(T1_GSE135337_2, 
                      y = c(T1_scRNAlist_2), 
                      add.cell.ids = c("T1_GSE135337_2","T1_scRNAlist_2"), 
                      project = "T1")
    rownames(T1)
    unique(sapply(X = strsplit(colnames(T1), split = "_"), FUN = "[", 1))
    table(T1$orig.ident)
    rownames(TaT1)
    colnames(TaT1)
    saveRDS(T1, file = "T1_object_v4.rds")

    T2 <- merge(T2_GSE135337list, 
                      y = T2_scRNAlist_2, 
                      add.cell.ids = c("T2_GSE135337list","T2_scRNAlist_2"), 
                      project = "T2")
    rownames(T2)
    unique(sapply(X = strsplit(colnames(T2), split = "_"), FUN = "[", 1))
    table(T2$orig.ident)
    saveRDS(T2, file = "T2_object_v4.rds")
    
    T3 <- merge(T3_GSE135337_2, 
                      y = T3_scRNAlist_2, 
                      add.cell.ids = c("T3_GSE135337_2","T3_scRNAlist_2"), 
                      project = "T3")
    rownames(T3)
    unique(sapply(X = strsplit(colnames(T3), split = "_"), FUN = "[", 1))
    table(T3$orig.ident)
    saveRDS(T3, file = "T3_object_v4.rds")
    T3<-readRDS("T3_object_v4.rds")  25957 32326
    
    T4<-T4_scRNAlist
    saveRDS(T4, file = "T4_object_v4.rds")
    rownames(T4)
    T4<-readRDS("T4_object_v4.rds")  22183  6464
}

{TaT1T2T3T4 <- merge(Ta, 
                      y = c(T1,T2,T3,T4), 
                      add.cell.ids = c("Ta","T1","T2","T3","T4"), 
                      project = "TaT1T2T3T4")

    relation_data<-TaT1T2T3T4@meta.data
    relation_data$group<-""
    relation_data[grep("T4",rownames(TaT1T2T3T4@meta.data)),dim(relation_data)[2]]="T4"
    relation_data[grep("T3",rownames(TaT1T2T3T4@meta.data)),dim(relation_data)[2]]="T3"
    relation_data[grep("T2",rownames(TaT1T2T3T4@meta.data)),dim(relation_data)[2]]="T2"
    relation_data[grep("T1",rownames(TaT1T2T3T4@meta.data)),dim(relation_data)[2]]="T1"
    relation_data[grep("Ta",rownames(TaT1T2T3T4@meta.data)),dim(relation_data)[2]]="Ta"
    TaT1T2T3T4@meta.data<-relation_data
    table(TaT1T2T3T4@meta.data$group) 
    saveRDS(TaT1T2T3T4, file = "TaT1T2T3T4.rds")
    TaT1T2T3T4<-readRDS("TaT1T2T3T4.rds")
}

{
    TaT1T2T3T4<-readRDS(file="TaT1T2T3T4.rds")
    TaT1T2T3T4 <- NormalizeData(TaT1T2T3T4, normalization.method = "LogNormalize", scale.factor = 10000)
    TaT1T2T3T4 <- FindVariableFeatures(TaT1T2T3T4,selection.method = "vst", nfeatures = 2000)
    TaT1T2T3T4 <- ScaleData(TaT1T2T3T4, features = rownames(TaT1T2T3T4))
    TaT1T2T3T4 <- RunPCA(TaT1T2T3T4, verbose = FALSE)
    rownames(TaT1T2T3T4)
  TaT1T2T3T4_harmony <- RunHarmony(TaT1T2T3T4, group.by.vars = "group")

        TaT1T2T3T4_harmony[["percent.mt"]] <- PercentageFeatureSet(TaT1T2T3T4_harmony, pattern = "^MT-")
        TaT1T2T3T4_harmony <- subset(TaT1T2T3T4_harmony, subset = nFeature_RNA > 500 & percent.mt < 20)  
        TaT1T2T3T4_harmony[["percent.mt"]] <- PercentageFeatureSet(TaT1T2T3T4_harmony, pattern = "^MT-") 
        TaT1T2T3T4_harmony$group <- factor(TaT1T2T3T4_harmony$group, levels = c("Ta", "T1", "T2", "T3", "T4"))
        pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/all_MT.pdf", width = 10, height = 8)
        p<-VlnPlot(TaT1T2T3T4_harmony, group.by = "group",features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
        print(p)
        dev.off()

  TaT1T2T3T4_harmony <- FindNeighbors(TaT1T2T3T4_harmony, dims = 1:30)
        TaT1T2T3T4_harmony <- FindClusters(TaT1T2T3T4_harmony, resolution = 0.5)
        TaT1T2T3T4_umap <- RunUMAP(TaT1T2T3T4_harmony, dims = 1:30)
        saveRDS(TaT1T2T3T4_umap, file = "TaT1T2T3T4_umap.rds")
        TaT1T2T3T4_umap<-readRDS("TaT1T2T3T4_umap.rds")

  p3 <- DimPlot(TaT1T2T3T4_umap,group.by = "group",reduction = "umap",combine = TRUE,label.size = 2)+
              theme(axis.title.x = element_text(size = 30),
              axis.title.y = element_text(size = 30),
              legend.text = element_text(size = 30),
              plot.title = element_text(size = 30))
  ggsave(plot=p3,filename="/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/降维/TaT1T2T3T4_umap_revise.pdf", width = 10, height = 10)


  p3 <- DimPlot(TaT1T2T3T4_umap,reduction = "umap",combine = TRUE,label.size = 2)+
              theme(axis.title.x = element_text(size = 30),
              axis.title.y = element_text(size = 30),
              legend.text = element_text(size = 20),
              plot.title = element_text(size = 30))
  ggsave(plot=p3,filename="/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/降维/TaT1T2T3T4_umap_30_revise.pdf", width = 10, height = 10)
}