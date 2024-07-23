#Fibroblast
{
  setwd("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/Fibroblast")
  TaT1T2T3T4_umap_手动<-readRDS("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/注释/TaT1T2T3T4_umap_手动.rds")
  sub_Fibroblast = TaT1T2T3T4_umap_手动[,TaT1T2T3T4_umap_手动$celltype  %in% "Fibroblast"]
  dim(sub_Fibroblast)#29959  1796 
  table(sub_Fibroblast@meta.data$group)
  sub_Fibroblast <- NormalizeData(sub_Fibroblast, normalization.method = "LogNormalize", scale.factor = 10000)
          sub_Fibroblast <- FindVariableFeatures(sub_Fibroblast,selection.method = "vst", nfeatures = 2000)
          sub_Fibroblast <- ScaleData(sub_Fibroblast, features = rownames(sub_Fibroblast ))
          sub_Fibroblast <- RunPCA(sub_Fibroblast, verbose = FALSE)
          sub_Fibroblast <- RunHarmony(sub_Fibroblast , group.by.vars = "group")
          sub_Fibroblast <- FindNeighbors(sub_Fibroblast, dims = 1:30)
          sub_Fibroblast <- FindClusters(sub_Fibroblast, resolution = 0.5)
          sub_Fibroblast_UMAP <- RunUMAP(sub_Fibroblast, dims = 1:30)
  saveRDS(sub_Fibroblast_UMAP, file = "/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/Fibroblast/sub_Fibroblast_UMAP.rds")  
  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/Fibroblast/sub_Fibroblast_UMAP.pdf",width=10,height=10,)
            p3 <- DimPlot(sub_Fibroblast_UMAP,group.by = "group",reduction = "umap",combine = FALSE,label.size = 2)
            print(p3)
            dev.off()
  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/Fibroblast/sub_Fibroblast_UMAP_11.pdf",width=10,height=10)
            p3 <- DimPlot(sub_Fibroblast_UMAP,reduction = "umap",combine = FALSE,label = TRUE,label.size = 2 )
            print(p3)
            dev.off()
  available_assays <- Assays(sub_NK_UMAP)
  dim(sub_NK_UMAP)
  sub_NK_UMAP<-JoinLayers(sub_NK_UMAP)
  markers_Fibroblast<- FindAllMarkers(object = sub_Fibroblast_UMAP, test.use="wilcox" ,
                            only.pos = TRUE,
                            logfc.threshold = 0.25)
  saveRDS(markers_Fibroblast, file = "markers_Fibroblast.rds") 
  write.csv(markers_Fibroblast, file = "/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/NK/markers_Fibroblast.csv", row.names = TRUE)  
 

  t<-table(Idents(sub_Fibroblast_UMAP), sub_Fibroblast_UMAP$group)
  tt<-prop.table(table(Idents(sub_Fibroblast_UMAP), sub_Fibroblast_UMAP$group))
  write.csv(t, file = "/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/Fibroblast/Fiber各簇占比.csv", row.names = TRUE)
  write.csv(tt, file = "/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/Fibroblast/Fiber各簇占比%.csv", row.names = TRUE)

  cell.prop<-as.data.frame(prop.table(table(sub_Fibroblast_UMAP@meta.data$seurat_clusters, sub_Fibroblast_UMAP$group)))
    head(cell.prop)
    dim(cell.prop)
    colnames(cell.prop)<-c("clusters","group","proportion")
    pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/Fibroblast/cluster堆积图.pdf",width=6,height=10)
    p<-ggplot(cell.prop,aes(group,proportion,fill=clusters))+
    geom_bar(stat="identity",position="fill")+
    ggtitle("")+
    theme_bw()+
    theme(axis.ticks.length=unit(0.5,'cm'))+
    guides(fill=guide_legend(title=NULL))
    print(p)
    dev.off()   

   
  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/Fibroblast/dCAF.pdf",width=30,height=15) 
  p4 <- DotPlot(sub_Fibroblast_UMAP, features = genes_to_check,
              assay='RNA' ,group.by = 'seurat_clusters' ) + coord_flip()+ggtitle("")
  print(p4)
  dev.off()
  genes_to_check = c("FCGR3A","CX3CR1","CEP78","CD47","NCR3","GNG2","CXCR4","IFNG","NEU1","FKBP4","IL2RB","XCL1","XCL2","KLRB1","PLCG2","CD160","CCL4","REL","SYTL3","TGFB1","KRT81","KRT86","ITGA1","STMN1","MKI67")
  genes_to_check = c("TUBA1B", "MKI67","E2F","G2M")
                    # ## mCAF   MMP11  COL1A2 4
                    # 'CSPG4',#iCAF PLA2G2A  2
                    # 'ENPP3','FCER1A','FCGR1A','FCGR2A', 'FCGR2B'#vCAF NOTCH3
                    # #tCAF PDPN MME TMEM158  NDRG1 VEGFA
                    # #ifnCAF IL32 CXCL9 CXCL10
                    # #apCAF HLA-DRA HLA-DRB1 CD74
                    # nature9   Cancer-associated fibroblast classification in single-cell and spatial proteomics data
                    # mCAF 'MMP11',"COL1A2","COL10A1","COL11A1","COL8A1","COL12A1","COL3A1","COL5A2" 0 4 
                    # iCAF 'IL6',"PLA2G2A","CFD","C3","CD34" 2 
                    # vCAF "COL18A1","NOTCH3","MCAM","CD146","RGS5" 3 
                    # tCAF 'PDPN',"MME","TMEM158","NDRG1","ENO1","HSPH1","HSP90AA1","IX"  6 11
                    # ifnCAF 'IL32',"CXCL9","CXCL10","CXCL11","IDO1"," IL2","IL6"  cancel
                    # apCAF "HLA-B", "HLA-DRB1", "HLA-C", "HLA-A", "MICA", "HLA-DQB1", "HLA-DQA1", "HLA-DPB1", "TAP2",
                    #  "HLA-H", "MICB", "HLA-DPA1", "HLA-DRB3", "HLA-K", "HLA-G", "HLA-DMB", "HLA-E", "HLA-L", "HLA-DRB4", 
                    #  "TAP1", "HLA-DRA", "HLA-DMA", "HLA-DOB", "HLA-P", "HLA-DRB5", "HLA-F", "HLA-DOA", "HLA-J", "HLA-V" 5 10 11
                    # rCAF CCL21 CCL19 cancel
                    # dCAF TUBA1B MKI67 9

  celltype=data.frame(ClusterID=0:11,
                                    celltype='unkown')

                celltype[celltype$ClusterID %in% c(0,4),2]="mCAF"
                celltype[celltype$ClusterID %in% c(2),2]="iCAF"
                celltype[celltype$ClusterID %in% c(3),2]="vCAF"
                celltype[celltype$ClusterID %in% c(5,10,11),2]="apCAF"
                celltype[celltype$ClusterID %in% c(6),2]="tCAF"
                celltype[celltype$ClusterID %in% c(1,7,8,9),2]="Unidentified_CAF"

  table(celltype)
                head(celltype)
                colnames(celltype)
                table(celltype$celltype)
                sce.in=sub_Fibroblast_UMAP
                sce.in@meta.data$手动celltype = "NA"
                for(i in 1:nrow(celltype)){
                  sce.in@meta.data[which(sce.in@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
                table(sce.in@meta.data$celltype)
                sub_Fibroblast_UMAP_2=sce.in
                saveRDS(sub_Fibroblast_UMAP_2,file="/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/Fibroblast/sub_Fibroblast_UMAP_2.rds")
                pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/Fibroblast/sub_Fibroblast_UMAP_2_A_new.pdf",width=15,height=15) 
                p.dim.cell=DimPlot(sub_Fibroblast_UMAP_2, reduction = "umap", group.by = "celltype",label = T,pt.size = 1,label.size = 7)+
                                    theme(legend.text = element_text(size = 18),legend.title = element_text(size = 50)+)
                print(p.dim.cell)
                dev.off()
                pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/Fibroblast/sub_Fibroblast_UMAP_2_B_new.pdf",width=15,height=15) 
                p.dim.cell=DimPlot(sub_Fibroblast_UMAP_2, reduction = "umap", group.by = "group",label = T,pt.size = 1,label.size = 7) +
                                    theme(legend.text = element_text(size = 18),legend.title = element_text(size = 50))
                print(p.dim.cell)
                dev.off()

                pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/Fibroblast/sub_Fibroblast_UMAP_2_A_new_revise.pdf",width=10,height=10) 
                p.dim.cell=DimPlot(sub_Fibroblast_UMAP_2, reduction = "umap", group.by = "celltype",label = T,pt.size = 1,label.size = 7)+
                                    theme(axis.title.x = element_text(size = 30),
                                          axis.title.y = element_text(size = 30),
                                          legend.text = element_text(size = 20),
                                          plot.title = element_text(size = 30),
                                          legend.title = element_text(size = 50))
                print(p.dim.cell)
                dev.off()
                pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/Fibroblast/sub_Fibroblast_UMAP_2_B_new_revise.pdf",width=10,height=10) 
                p.dim.cell=DimPlot(sub_Fibroblast_UMAP_2, reduction = "umap", group.by = "group",label = T,pt.size = 1,label.size = 7) +
                                    theme(axis.title.x = element_text(size = 30),
                                          axis.title.y = element_text(size = 30),
                                          legend.text = element_text(size = 20),
                                          plot.title = element_text(size = 30),
                                          legend.title = element_text(size = 50))
                print(p.dim.cell)
                dev.off()

  library(dplyr)
  anotation_uamp<-readRDS("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/注释/TaT1T2T3T4_umap_手动.rds")
  sub_Fibroblast<-readRDS("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/Fibroblast/sub_Fibroblast_UMAP.rds")
  anno_meta<-read.csv("单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/Fibroblast/成纤维细胞簇.csv",header=F)
  anotation_uamp@meta.data$sub_class<-anotation_uamp@meta.data$celltype
  sub_Fibroblast@meta.data %>% head(.) %>% rownames(.)
  anotation_uamp@meta.data %>% head(.) %>% rownames(.)
}

{
  matrix_meta<-c()
  for(i in 1:dim(anno_meta)[1] )
  {
    tem_barcode<-sub_Fibroblast@meta.data[which(anno_meta[i,1]==sub_Fibroblast@meta.data$seurat_clusters),] %>% rownames(.)
    anotation_uamp@meta.data[match(tem_barcode,rownames(anotation_uamp@meta.data)),match("sub_class",colnames(anotation_uamp@meta.data))]<-anno_meta[i,2]
  }
  unique(anotation_uamp@meta.data$sub_class)
  TaT1T2T3T4_umap_手动2<-anotation_uamp
  saveRDS(TaT1T2T3T4_umap_手动2,file = "/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/注释/TaT1T2T3T4_umap_手动2")
}

{
  TaT1T2T3T4_cafs<-readRDS(file="/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/注释/TaT1T2T3T4_umap_手动2")
  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/注释/手动注释_umap_cafs.pdf",width=15,height=15) 
  p.dim.cell=DimPlot(TaT1T2T3T4_cafs, reduction = "umap", group.by = "sub_class",label = T,pt.size = 1) 
  print(p.dim.cell)
  dev.off()
  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/注释/手动注释2_umap_cell.pdf",width=15,height=15) 
  p.dim.cell=DimPlot(TaT1T2T3T4_umap_手动, reduction = "umap", group.by = "group",label = T,pt.size = 1) 
  print(p.dim.cell)
  dev.off()
}

{
    sub_Fibroblast_UMAP_2<- readRDS("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/Fibroblast/sub_Fibroblast_UMAP_2.rds")
    tab_num<-table(sub_Fibroblast_UMAP_2@meta.data$celltype, sub_Fibroblast_UMAP_2$group)
    tab_num_percent <- sweep(tab_num, 2, colSums(tab_num), FUN = "/")
    cell.prop<-as.data.frame(tab_num_percent)
    head(cell.prop)
    dim(cell.prop)
    colnames(cell.prop)<-c("clusters","group","proportion")
    pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/Fibroblast/celltype堆积图.pdf",width=6,height=10)
    p<-ggplot(cell.prop,aes(group,proportion,fill=clusters))+
    geom_bar(stat="identity",position="fill")+
    ggtitle("")+
    theme_bw()+
    theme(axis.ticks.length=unit(0.5,'cm'))+
    guides(fill=guide_legend(title=NULL))
    print(p)
    dev.off()  

    library(ggalluvial)
    colnames(cell.prop)<-c("celltype","group","proportion")
    pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/Fibroblast/celltype堆积图_revise.pdf",width=6,height=6)
    p<-ggplot(cell.prop,aes(group,proportion,fill=celltype,stratum = celltype, alluvium = celltype),lable=TRUE)+
    ggtitle("")+
    theme_bw()+
    theme(axis.ticks.length=unit(0.5,'cm'),
          axis.text.x = element_text(size = 13))+
    guides(fill=guide_legend(title="Cell Type"))+
    geom_stratum(width = 0.5, color='white')+
    geom_alluvium(alpha = 0.5,
                width = 0.5,
                curve_type = "linear")
    print(p)
    dev.off()   


}
{ 
  TaT1T2T3T4_umap_手动2<-readRDS("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/注释/TaT1T2T3T4_umap_手动2.rds")
  sub_Ta = TaT1T2T3T4_umap_手动2[,TaT1T2T3T4_umap_手动2$group  %in% "Ta"]
  sub_T1 = TaT1T2T3T4_umap_手动2[,TaT1T2T3T4_umap_手动2$group  %in% "T1"]
  sub_T2 = TaT1T2T3T4_umap_手动2[,TaT1T2T3T4_umap_手动2$group  %in% "T2"]
  sub_T3 = TaT1T2T3T4_umap_手动2[,TaT1T2T3T4_umap_手动2$group  %in% "T3"]
  sub_T4 = TaT1T2T3T4_umap_手动2[,TaT1T2T3T4_umap_手动2$group  %in% "T4"]

  dim(sub_Ta) 
  sub_Ta<-JoinLayers(sub_Ta)
  sub_Ta[["RNA"]] <- as(sub_Ta[["RNA"]], "Assay")
  chat_cell<-sub_Ta
  data.input <- sub_Ta[["RNA"]]@data # normalized data matrix
  labels <- Idents(sub_Ta)
  meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
  cellchat <- createCellChat(object = chat_cell, group.by = "sub_class", assay = "RNA")
  load("CellChatDB.human.rda")
  CellChatDB <- CellChatDB.human
  showDatabaseCategory(CellChatDB)
  dplyr::glimpse(CellChatDB$interaction)
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  future::plan("multiprocess", workers = 1) 
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  load("PPI.human.rda")
  cellchat <- projectData(cellchat, PPI.human)
  cellchat@data.project[1:4,1:4]
  cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE) 
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  groupSize <- as.numeric(table(cellchat@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/caf_Ta_1.pdf",width=15 ,height=10)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, vertex.label.cex=2,
                  weight.scale = T, label.edge= F, title.name = "Number of interactions")
  dev.off()
  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/caf_Ta_2.pdf",width=18,height=12)
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, vertex.label.cex=3,
                  weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()
  
  dim(sub_T1) 
  sub_T1<-JoinLayers(sub_T1)
  sub_T1[["RNA"]] <- as(sub_T1[["RNA"]], "Assay")
  chat_cell<-sub_T1
  data.input <- sub_T1[["RNA"]]@data # normalized data matrix
  labels <- Idents(sub_T1)
  meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
  cellchat <- createCellChat(object = chat_cell, group.by = "sub_class", assay = "RNA")
  load("CellChatDB.human.rda")
  CellChatDB <- CellChatDB.human
  showDatabaseCategory(CellChatDB)
  dplyr::glimpse(CellChatDB$interaction)
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  future::plan("multiprocess", workers = 1) 
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  load("PPI.human.rda")
  cellchat <- projectData(cellchat, PPI.human)
  cellchat@data.project[1:4,1:4]
  cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE) 
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  groupSize <- as.numeric(table(cellchat@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/caf_T1_1.pdf",width=15 ,height=10)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, vertex.label.cex=2,
                  weight.scale = T, label.edge= F, title.name = "Number of interactions")
  dev.off()
  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/caf_T1_2.pdf",width=18,height=12)
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, vertex.label.cex=2,
                  weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()
  
dim(sub_T2) 
  sub_T2<-JoinLayers(sub_T2)
  sub_T2[["RNA"]] <- as(sub_T2[["RNA"]], "Assay")
  chat_cell<-sub_T2
  data.input <- sub_T2[["RNA"]]@data # normalized data matrix
  labels <- Idents(sub_T2)
  meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
  cellchat <- createCellChat(object = chat_cell, group.by = "sub_class", assay = "RNA")
  load("CellChatDB.human.rda")
  CellChatDB <- CellChatDB.human
  showDatabaseCategory(CellChatDB)
  dplyr::glimpse(CellChatDB$interaction)
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  future::plan("multiprocess", workers = 1) 
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  load("PPI.human.rda")
  cellchat <- projectData(cellchat, PPI.human)
  cellchat@data.project[1:4,1:4]
  cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE) 
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  groupSize <- as.numeric(table(cellchat@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/caf_T2_1.pdf",width=15 ,height=10)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, vertex.label.cex=2,
                  weight.scale = T, label.edge= F, title.name = "Number of interactions")
  dev.off()
  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/caf_T2_2.pdf",width=18,height=12)
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, vertex.label.cex=2,
                  weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()

dim(sub_T3) 
  sub_T3<-JoinLayers(sub_T3)
  sub_T3[["RNA"]] <- as(sub_T3[["RNA"]], "Assay")
  chat_cell<-sub_T3
  data.input <- sub_T3[["RNA"]]@data # normalized data matrix
  labels <- Idents(sub_T3)
  meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
  cellchat <- createCellChat(object = chat_cell, group.by = "sub_class", assay = "RNA")
  load("CellChatDB.human.rda")
  CellChatDB <- CellChatDB.human
  showDatabaseCategory(CellChatDB)
  dplyr::glimpse(CellChatDB$interaction)
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  future::plan("multiprocess", workers = 1) 
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  load("PPI.human.rda")
  cellchat <- projectData(cellchat, PPI.human)
  cellchat@data.project[1:4,1:4]
  cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE) 
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  groupSize <- as.numeric(table(cellchat@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/caf_T3_1.pdf",width=15 ,height=10)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, vertex.label.cex=2,
                  weight.scale = T, label.edge= F, title.name = "Number of interactions")
  dev.off()
  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/caf_T3_2.pdf",width=18,height=12)
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, vertex.label.cex=2,
                  weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()


dim(sub_T4) 
  sub_T4<-JoinLayers(sub_T4)
  sub_T4[["RNA"]] <- as(sub_T4[["RNA"]], "Assay")
  chat_cell<-sub_T4
  data.input <- sub_T4[["RNA"]]@data # normalized data matrix
  labels <- Idents(sub_T4)
  meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
  cellchat <- createCellChat(object = chat_cell, group.by = "sub_class", assay = "RNA")
  load("CellChatDB.human.rda")
  CellChatDB <- CellChatDB.human
  showDatabaseCategory(CellChatDB)
  dplyr::glimpse(CellChatDB$interaction)
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  future::plan("multiprocess", workers = 1) 
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  load("PPI.human.rda")
  cellchat <- projectData(cellchat, PPI.human)
  cellchat@data.project[1:4,1:4]
  cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE) 
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)  groupSize <- as.numeric(table(cellchat@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/caf_T4_1.pdf",width=15 ,height=10)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, vertex.label.cex=2,
                  weight.scale = T, label.edge= F, title.name = "Number of interactions")
  dev.off()
  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/caf_T4_2.pdf",width=18,height=12)
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, vertex.label.cex=2,
                  weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()

  #CAFs的cellchat
{
    TaT1T2T3T4_umap_手动2<-readRDS(file="/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/注释/TaT1T2T3T4_umap_手动2")
    TaT1T2T3T4_umap_手动2<-JoinLayers(TaT1T2T3T4_umap_手动2)
    TaT1T2T3T4_umap_手动2[["RNA"]] <- as(TaT1T2T3T4_umap_手动2[["RNA"]], "Assay")

    chat_cell<-TaT1T2T3T4_umap_手动2

    data.input <- TaT1T2T3T4_umap_手动2[["RNA"]]@data # normalized data matrix
    labels <- Idents(TaT1T2T3T4_umap_手动2)
    meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
    cellchat <- createCellChat(object = chat_cell, group.by = "sub_class", assay = "RNA")
    load("CellChatDB.human.rda")
    CellChatDB <- CellChatDB.human
    showDatabaseCategory(CellChatDB)
    dplyr::glimpse(CellChatDB$interaction)
    CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
    cellchat@DB <- CellChatDB.use
    cellchat <- subsetData(cellchat)
    future::plan("multiprocess", workers = 1) 
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    load("PPI.human.rda")
    cellchat <- projectData(cellchat, PPI.human)
    cellchat@data.project[1:4,1:4]
    cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE) 
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    saveRDS(cellchat,file = "/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/注释/cellchat2_cafs.rds")

    cellchat<-readRDS("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/cellchat2.rds")
    df.net <- subsetCommunication(cellchat)
    write.csv(df.net, "/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/df.all.csv")
    df.net1 <- subsetCommunication(cellchat,slot.name = "netP")
    levels(cellchat@idents)
    df.net2 <- subsetCommunication(cellchat, sources.use = c("Epithelial cell"), targets.use = c("Fibroblast" ,"Tcell")) 
    df.net3 <- subsetCommunication(cellchat, signaling = c("CCL", "TGFb"))
    write.csv(df.net1, "/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/dfnet1.all.csv")
    write.csv(df.net2, "/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/dfnet2.all.csv")
    write.csv(df.net3, "/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/dfnet3.all.csv")
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    groupSize <- as.numeric(table(cellchat@idents))
    par(mfrow = c(1,2), xpd=TRUE)

    pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/F1.pdf",width=12,height=8)
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                    weight.scale = T, label.edge= F, title.name = "Number of interactions")
    dev.off()

    pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/F2.pdf",width=12,height=8)
    netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                    weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
    dev.off()

    pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/F3.pdf",width=12,height=8)
    mat <- cellchat@net$weight
    par(mfrow = c(3,5), xpd=TRUE)
    for (i in 1:nrow(mat)) {
      mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
      mat2[i, ] <- mat[i, ]
      netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
    }
    dev.off()
    cellchat@netP$pathways
    pathways.show <- c("HH")  
    levels(cellchat@idents)   
    #[1] "Epi"        "Myeloid"    "Fibroblast" "T"          "Endo"       "un"     
    vertex.receiver = c(1,4,7,8,13,14,15) 
    pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/F5_hh.pdf",width=12,height=8)
    netVisual_aggregate(cellchat, signaling = pathways.show,  
                        vertex.receiver = vertex.receiver,layout = "hierarchy")
    dev.off()

}
}