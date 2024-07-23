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
  library(SingleR)
  library(devtools)
  
  T1T2T3T4_umap<-readRDS(file="T1T2T3T4_umap.rds")
  available_assays <- Assays(T1T2T3T4_umap)
  dim(T1T2T3T4_umap)
  T1T2T3T4_umap<-JoinLayers(T1T2T3T4_umap)
  markers_T1T2T3T4 <- FindAllMarkers(object = T1T2T3T4_umap, test.use="wilcox" ,
                            only.pos = TRUE,
                            logfc.threshold = 0.25)
  saveRDS(markers_T1T2T3T4, file = "markers_T1T2T3T4.rds") 
  write.csv(markers_T1T2T3T4, file = "/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4/markers.csv", row.names = TRUE)                         
{
  class(TaT1T2T3T4_umap@meta.data)
  colnames(TaT1T2T3T4_umap@meta.data)
          assays(TaT1T2T3T4_umap)
  meta=TaT1T2T3T4_umap@meta.data 
  TaT1T2T3T4_umap<-JoinLayers(TaT1T2T3T4_umap)
  SingleR <- GetAssayData(TaT1T2T3T4_umap, layer="data")
  load(file="BlueprintEncode_bpe.se_human.RData")
  load(file="HumanPrimaryCellAtlas_hpca.se_human.RData")
  hesc <- SingleR(test = SingleR, ref = hpca.se, labels = hpca.se$label.main) 
  table(hesc$labels,meta$seurat_clusters) 


  TaT1T2T3T4_umap@meta.data$labels <- hesc$labels
  saveRDS(TaT1T2T3T4_umap, file = "TaT1T2T3T4_umap_自动注释.rds")

  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/注释/自动注释.pdf",width=15,height=15) 
  p<-DimPlot(TaT1T2T3T4_umap, reduction = "umap", group.by = "labels",label = T,label.size = 5) 
  print(p)
  dev.off()
  
  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/注释/自动注释2.pdf",width=10,height=10)        
        p4 <- DimPlot(TaT1T2T3T4_umap,reduction = "umap",combine = FALSE,label.size = 2)
        print(p4)
        dev.off()
}

{
  t<-table((Idents(TaT1T2T3T4_umap)), TaT1T2T3T4_umap$labels)
  write.csv(t, file = "/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/各细胞亚群对应 簇.csv", row.names = TRUE)

  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R//T1T2T3T4_2/注释/TVSNKcell2.pdf",width=30,height=15) 
  p4 <- DotPlot(TaT1T2T3T4_umap, features = genes_to_check,
              assay='RNA' ,group.by = 'seurat_clusters' ) + coord_flip()+ggtitle("")
  print(p4)
  dev.off()
  genes_to_check = c("FCGR3A","CX3CR1","CEP78","CD47","NCR3","GNG2","CXCR4","IFNG","NEU1","FKBP4","IL2RB","XCL1","XCL2","KLRB1","PLCG2","CD160","CCL4","REL","SYTL3","TGFB1","KRT81","KRT86","ITGA1","STMN1","MKI67")
  genes_to_check = c('PTPRC','CD3D','CD3E','CD4','CD8A',## Tcells "CD3D", "CD3E", "CD8A"
                    'CSPG4',#Pericyte
                   'ENPP3','FCER1A','FCGR1A','FCGR2A', 'FCGR2B',#Mast cell "TPSAB1","FCER1A","MS4A2","KIT"
                   'AQP3', 'MUC5AC', 'MUC5B', 'PIGR', 'SCGB1A1',#Airway secretory cell
                   'THBD', 'CD1C','CLEC4C','CD83','HLA.DQA2','HLA.DQA1',#Dendritic cell
                   'KRT19','ITGAE', 'VCAM1', 'IL6R', 'ANPEP', 'CD24', 'CDH1',#Epithelial cell  KRT20,"KRT20","SHOX2","KRT14","KRT19" 归癌细胞
                   'MUC1','ABCA3', 'LPCAT1', 'NAPSA', 'SFTPB', 'SFTPC', 'SLC34A2',#lung Epithelial cell ("EPCAM","EpCAM","KRT18","KRT19","NMP22","CFHR1","CFHR3")
                   'ACTA2', 'PDGFRA', 'PDGFRB','THY1',#Fibroblast "COL1A1","FGF7","MME","FAP"
                   'CD19', 'CD79A', 'MS4A1', # B cellsc  ("CD19", "CD79A", "MS4A1","IGHM","MS4A3","IGHD")
                   'TAGLN2','CD5',#Activated B cell
                   'CD27', 'CD38','LY9', 'LAIR1', 'ICAM1','KIT',  # plasma 
                   'NKG7','GNLY',#NK
                   'CD8B', 'ZNF683','FCGR3A', 'FCGR3B', 'NCAM1', 'KLRB1',#NKT
                   'CXCR5','CCR7',#memory T cell
                   'CD6','IL7R', 'IL2RA', 'IKZF2',#Treg
                   'GP1BA','SELL','IFNG','CXCR3','IL17A','IL4','GATA3',#T helper cell
                   'CD33', 'ENTPD1',#MDSC
                   'S100A8','S100A9','S100A12',#Neutrophil
                   'CD68',  'CD163','MRC1','MSR1','CXCL10','CCL18', ## Macrophage (belong to monocyte)("CD68","CD163","MRC1","MSR1","CXCL10","CCL18")（"CD14","CD68")
                   'PECAM1','VWF','MCAM','CD34','ESM1', ## Endothelial  "PECAM1","VWF","CD34","ESM1"
                   'ALDH1A1','KRT18', 'PROM1',## Cancer stem cell
                   'HHIP','SFTPC','SFTPA','SFTPC','LAMP3',
                   'MDK','SFTPB')

                celltype=data.frame(ClusterID=0:30,
                                    celltype='unkown')

                celltype[celltype$ClusterID %in% c(16,21,27),2]="Bcell"
                celltype[celltype$ClusterID %in% c(13,30),2]="Endothelial cell"
                celltype[celltype$ClusterID %in% c(14),2]="Fibroblast"
                celltype[celltype$ClusterID %in% c(4,9,12,19,28),2]="Tcell" 
                celltype[celltype$ClusterID %in% c(5,10),2]="monocytes/Macrophages"
                celltype[celltype$ClusterID %in% c(0,1,3,6,7,8,11,15,20,22,23,24,25),2]="Epithelial cell"
                celltype[celltype$ClusterID %in% c(17),2]="Neurons"
                celltype[celltype$ClusterID %in% c(26),2]="CMP"
                celltype[celltype$ClusterID %in% c(18,29),2]="Tissue_stem_cells"
                celltype[celltype$ClusterID %in% c(2),2]="NK"

                
                table(celltype)
                head(celltype)
                colnames(celltype)
                table(celltype$celltype)
                sce.in=TaT1T2T3T4_umap
                sce.in@meta.data$手动celltype = "NA"
                for(i in 1:nrow(celltype)){
                  sce.in@meta.data[which(sce.in@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
                table(sce.in@meta.data$celltype)
                TaT1T2T3T4_umap_手动=sce.in
                TaT1T2T3T4_umap_手动<-readRDS("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/注释/TaT1T2T3T4_umap_手动.rds")
                write.csv(data,file="/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4/降维/T1T2T3T4手动注释_matadata.csv")
                pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/注释/手动注释_umap_cell_new.pdf",width=15,height=15) 
                p.dim.cell=DimPlot(TaT1T2T3T4_umap_手动, reduction = "umap", group.by = "celltype",label = T,pt.size = 1,label.size=5) +
                                    theme(legend.text = element_text(size = 18),legend.title = element_text(size = 50))
                print(p.dim.cell)
                dev.off()
              p3 <- DimPlot(TaT1T2T3T4_umap_手动, reduction = "umap", group.by = "celltype",label = T,pt.size = 1,label.size=5)+
              theme(axis.title.x = element_text(size = 30),
              axis.title.y = element_text(size = 30),
              legend.text = element_text(size = 15),
              plot.title = element_text(size = 30))
              ggsave(plot=p3,filename="/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/注释/手动注释_umap_cell_new_revise.pdf", width = 12, height = 12)
              
              TaT1T2T3T4_umap_手动$celltype<-factor(x=TaT1T2T3T4_umap_手动$celltype,levels=c("Epithelial cell","Tcell","monocytes/Macrophages","NK","Bcell","Endothelial cell",
                                                    "Fibroblast","Tissue_stem_cells","Neurons","CMP"))

              p3<-DoHeatmap(TaT1T2T3T4_umap_手动,
                        features = c("KRT19","CD3D","LYZ","SYTL3","CD79A","VWF","DCN","MYL9","HES6","TPSAB1"),
                        group.by="celltype",
                        assay="RNA",
                        group.colors = c("#C77CFF" , "#7CAE00" , "#00BFC4" ,"#F8766D" ,"#AB82FF", "#90EE90" , "#00CD00","#008B8B","#FFA500","#3B5DFF"),
                        size=2
                        )+scale_fill_gradientn(colors = c("white","grey","firebrick3"))
              ggsave(plot=p3,filename="/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/注释/doheatmap_revise.pdf", width = 10, height = 10)

              TaT1T2T3T4_umap_手动
                pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/注释/手动注释2_umap_cell.pdf",width=15,height=15) 
                p.dim.cell=DimPlot(TaT1T2T3T4_umap_手动, reduction = "umap", group.by = "group",label = T,pt.size = 1) 
                print(p.dim.cell)
                dev.off()
                saveRDS(TaT1T2T3T4_umap_手动,file="/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/注释/TaT1T2T3T4_umap_手动.rds")

    cell.prop<-as.data.frame(prop.table(table(TaT1T2T3T4_umap_手动$celltype, TaT1T2T3T4_umap_手动$group)))
    pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/堆积图/group_手动堆积图2.pdf",width=6,height=10)
    p<-ggplot(cell.prop,aes(group,proportion,fill=celltype,stratum = celltype, alluvium = celltype))+
    ggtitle("")+
    theme_bw()+
    theme(axis.ticks.length=unit(0.5,'cm'))+
    guides(fill=guide_legend(title=NULL))+
    geom_stratum(width = 0.5, color='white')+
    geom_alluvium(alpha = 0.5,
                width = 0.5,
                curve_type = "linear")
    print(p)
    dev.off()   
}