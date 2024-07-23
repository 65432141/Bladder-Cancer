{
  setwd("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/")
  library(CellChat)
  library(tidyverse)
  library(Seurat)

  TaT1T2T3T4_umap_手动1<-readRDS("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/注释/TaT1T2T3T4_umap_手动.rds")
  TaT1T2T3T4_umap_手动1<-JoinLayers(TaT1T2T3T4_umap_手动1)
  TaT1T2T3T4_umap_手动1[["RNA"]] <- as(TaT1T2T3T4_umap_手动1[["RNA"]], "Assay")

  chat_cell<-TaT1T2T3T4_umap_手动1

  data.input <- TaT1T2T3T4_umap_手动1[["RNA"]]@data # normalized data matrix
  labels <- Idents(TaT1T2T3T4_umap_手动1)
  meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
  cellchat <- createCellChat(object = chat_cell, group.by = "celltype", assay = "RNA")
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
  cellchat<-readRDS("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/cellchat")
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

  
  saveRDS(cellchat,"/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/cellchat.rds")
  groupSize <- as.numeric(table(cellchat@idents))
  par(mfrow = c(1,2), xpd=TRUE)

  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/1.pdf",width=12,height=8)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                  weight.scale = T, label.edge= F, title.name = "Number of interactions")
  dev.off()

  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/2.pdf",width=12,height=8)
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                  weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()

  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/3.pdf",width=12,height=8)
  mat <- cellchat@net$weight
  par(mfrow = c(2,5), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
  dev.off()

  cellchat@netP$pathways

  pathways.show <- c("WNT")  
  levels(cellchat@idents)   
  #[1] "Epi"        "Myeloid"    "Fibroblast" "T"          "Endo"       "un"     
  vertex.receiver = c(1,2,4,5) 
  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/4.pdf",width=12,height=8)
  netVisual_aggregate(cellchat, signaling = pathways.show,  
                      vertex.receiver = vertex.receiver,layout = "hierarchy")
  dev.off()                    


  
  TaT1T2T3T4_umap_手动1<-readRDS("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/注释/TaT1T2T3T4_umap_手动.rds")
  sub_Ta = TaT1T2T3T4_umap_手动1[,TaT1T2T3T4_umap_手动1$group  %in% "Ta"]
  sub_T1 = TaT1T2T3T4_umap_手动1[,TaT1T2T3T4_umap_手动1$group  %in% "T1"]
  sub_T2 = TaT1T2T3T4_umap_手动1[,TaT1T2T3T4_umap_手动1$group  %in% "T2"]
  sub_T3 = TaT1T2T3T4_umap_手动1[,TaT1T2T3T4_umap_手动1$group  %in% "T3"]
  sub_T4 = TaT1T2T3T4_umap_手动1[,TaT1T2T3T4_umap_手动1$group  %in% "T4"]

  dim(sub_Ta) 
  sub_Ta<-JoinLayers(sub_Ta)
  sub_Ta[["RNA"]] <- as(sub_Ta[["RNA"]], "Assay")
  chat_cell<-sub_Ta
  data.input <- sub_Ta[["RNA"]]@data # normalized data matrix
  labels <- Idents(sub_Ta)
  meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
  cellchat <- createCellChat(object = chat_cell, group.by = "celltype", assay = "RNA")
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

  cellchat<-readRDS("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/cellchat")
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
  saveRDS(cellchat,"/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/cellchat_Ta.rds")
  groupSize <- as.numeric(table(cellchat@idents))
  par(mfrow = c(1,2), xpd=TRUE)

  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/Ta_1.pdf",width=15 ,height=10)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, vertex.label.cex=2,
                  weight.scale = T, label.edge= F, title.name = "Number of interactions")
  dev.off()

  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/Ta_2.pdf",width=18,height=12)
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, vertex.label.cex=3,
                  weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()
  
  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/Ta_3.pdf",width=12,height=8)
  mat <- cellchat@net$weight
  par(mfrow = c(2,5), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
  dev.off()

  cellchat@netP$pathways
  pathways.show <- c("WNT")  
  levels(cellchat@idents)   
  #[1] "Epi"        "Myeloid"    "Fibroblast" "T"          "Endo"       "un"     
  vertex.receiver = c(1,2,4,5) 
  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/Ta_4.pdf",width=12,height=8)
  netVisual_aggregate(cellchat, signaling = pathways.show,  
                      vertex.receiver = vertex.receiver,layout = "hierarchy")
  dev.off()
  

  dim(sub_T4) 
  sub_Ta<-JoinLayers(sub_T4)
  sub_T4[["RNA"]] <- as(sub_T4[["RNA"]], "Assay")
  chat_cell<-sub_T4
  data.input <- sub_T4[["RNA"]]@data # normalized data matrix
  labels <- Idents(sub_T4)
  meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
  cellchat <- createCellChat(object = chat_cell, group.by = "celltype", assay = "RNA")
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

  cellchat<-readRDS("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/cellchat")
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
  saveRDS(cellchat,"/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/cellchat_T4.rds")
  groupSize <- as.numeric(table(cellchat@idents))
  par(mfrow = c(1,2), xpd=TRUE)

  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/T4_1.pdf",width=15 ,height=10)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize,vertex.label.cex=2, 
                  weight.scale = T, label.edge= F, title.name = "Number of interactions")
  dev.off()

  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/T4_2.pdf",width=18,height=12)
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, vertex.label.cex=3,
                  weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()
  
  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/T4_3.pdf",width=12,height=8)
  mat <- cellchat@net$weight
  par(mfrow = c(2,5), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
  dev.off()

dim(sub_T1) 
  sub_T1<-JoinLayers(sub_T1)
  sub_T1[["RNA"]] <- as(sub_T1[["RNA"]], "Assay")
  chat_cell<-sub_T1
  data.input <- sub_T1[["RNA"]]@data # normalized data matrix
  labels <- Idents(sub_T1)
  meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
  cellchat <- createCellChat(object = chat_cell, group.by = "celltype", assay = "RNA")
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

  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/T1_1.pdf",width=15 ,height=10)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, vertex.label.cex=2, 
                  weight.scale = T, label.edge= F, title.name = "Number of interactions")
  dev.off()

  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/T1_2.pdf",width=18,height=12)
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, vertex.label.cex=3, 
                  weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()

  dim(sub_T2) 
  sub_T2<-JoinLayers(sub_T2)
  sub_T2[["RNA"]] <- as(sub_T2[["RNA"]], "Assay")
  chat_cell<-sub_T2
  data.input <- sub_T2[["RNA"]]@data # normalized data matrix
  labels <- Idents(sub_T2)
  meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
  cellchat <- createCellChat(object = chat_cell, group.by = "celltype", assay = "RNA")
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

  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/T2_1.pdf",width=15 ,height=10)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, vertex.label.cex=2, 
                  weight.scale = T, label.edge= F, title.name = "Number of interactions")
  dev.off()

  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/T2_2.pdf",width=18,height=12)
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, vertex.label.cex=3, 
                  weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()

  dim(sub_T3) 
  sub_T3<-JoinLayers(sub_T3)
  sub_T3[["RNA"]] <- as(sub_T3[["RNA"]], "Assay")
  chat_cell<-sub_T3
  data.input <- sub_T3[["RNA"]]@data # normalized data matrix
  labels <- Idents(sub_T3)
  meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
  cellchat <- createCellChat(object = chat_cell, group.by = "celltype", assay = "RNA")
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

  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/T3_1.pdf",width=20 ,height=20)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                  weight.scale = T, label.edge= F, title.name = "Number of interactions")
  dev.off()

  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/cellchat/T3_2.pdf",width=12,height=8)
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                  weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()
}