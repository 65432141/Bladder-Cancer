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
  
  TaT1T2T3T4_umap_手动<-readRDS("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/注释/TaT1T2T3T4_umap_手动.rds")
  sub_cancer_手动 = TaT1T2T3T4_umap_手动[,TaT1T2T3T4_umap_手动$celltype  %in% "Epithelial cell"]
  sub_cancer_手动 <- NormalizeData(sub_cancer_手动, normalization.method = "LogNormalize", scale.factor = 10000)
          sub_cancer_手动 <- FindVariableFeatures(sub_cancer_手动,selection.method = "vst", nfeatures = 2000)
          sub_cancer_手动 <- ScaleData(sub_cancer_手动, features = rownames(sub_cancer_手动 ))
          sub_cancer_手动 <- RunPCA(sub_cancer_手动, verbose = FALSE)
          sub_cancer_手动 <- RunHarmony(sub_cancer_手动 , group.by.vars = "group")
          sub_cancer_手动 <- FindNeighbors(sub_cancer_手动, dims = 1:30)
          sub_cancer_手动 <- FindClusters(sub_cancer_手动, resolution = 0.5)
          sub_cancer_umap手动 <- RunUMAP(sub_cancer_手动, dims = 1:30)
  p3 <- DimPlot(sub_cancer_umap手动,group.by = "group",reduction = "umap",combine = TRUE,label.size = 2)+
              theme(axis.title.x = element_text(size = 30),
              axis.title.y = element_text(size = 30),
              legend.text = element_text(size = 30),
              plot.title = element_text(size = 30))
  ggsave(plot=p3,filename="/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/癌细胞亚群/sub_cancer_umap手动_revise.pdf", width = 10, height = 10)
  p3 <- DimPlot(sub_cancer_umap手动,reduction = "umap",combine = TRUE,label.size = 5)+
              theme(axis.title.x = element_text(size = 30),
              axis.title.y = element_text(size = 30),
              legend.text = element_text(size = 20),
              plot.title = element_text(size = 30))
  ggsave(plot=p3,filename="/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/癌细胞亚群/sub_cancer_umap_26手动_revise.pdf", width = 10, height = 10)

  T1T2T3T4_umap<-readRDS(file="T1T2T3T4_umap.rds")
  markers_T1T2T3T4 <- FindAllMarkers(object = T1T2T3T4_umap, test.use="wilcox" ,
                            only.pos = TRUE,
                            logfc.threshold = 0.25)
  saveRDS(markers_T1T2T3T4, file = "markers_T1T2T3T4.rds") 
  write.csv(markers_T1T2T3T4, file = "/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4/markers.csv", row.names = TRUE)                         

{
  library(gplots)
  library(ggplot2)
  library(monocle)
  library(dplyr)  
  setwd("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/拟时序")
  pd <- new('AnnotatedDataFrame', data = sub_cancer_umap手动@meta.data)
  fData <- data.frame(gene_short_name = row.names(sub_cancer_umap手动1), row.names = row.names(sub_cancer_umap手动1))
  fd <- new('AnnotatedDataFrame', data = fData)
  monocle_cds <- newCellDataSet(sub_cancer_umap手动1,
                                phenoData = pd,
                                featureData = fd,
                                lowerDetectionLimit = 0.5,
                                expressionFamily = negbinomial.size())
  monocle_cds <- estimateSizeFactors(monocle_cds) 
  monocle_cds <- estimateDispersions(monocle_cds)
  monocle_cds=detectGenes(monocle_cds,min_expr = 3)
  saveRDS(monocle_cds,file="/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/拟时序/subcancer_monocle_cds.rds")
  monocle_cds<-readRDS("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/拟时序/subcancer_monocle_cds.rds")
  
  print(head(fData(monocle_cds)))
  expressed_genes <- row.names(subset(fData(monocle_cds), num_cells_expressed >= 10))
  print(head(pData(monocle_cds)))

  diff_test_res <- differentialGeneTest(monocle_cds[expressed_genes,],
                                        fullModelFormulaStr = "~ seurat_clusters")
  head(diff_test_res)                         
  ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
  saveRDS(ordering_genes,file="/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/拟时序/ordering_genes.rds")
  ordering_genes<- readRDS("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/拟时序/ordering_genes.rds")
  monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)
  monocle_cds <- reduceDimension(monocle_cds, max_components = 2,
                              method = 'DDRTree')
  monocle_cds <- orderCells(monocle_cds)
  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/拟时序//拟时序_clusters.pdf",width=30,height=30)
  p3<-plot_cell_trajectory(monocle_cds, color_by = "seurat_clusters")
  print(p3)
  dev.off()

  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/拟时序/拟时序stage.pdf",width=30,height=30)
  p3<-plot_cell_trajectory(monocle_cds,color_by = "State")
  print(p3)
  dev.off()

  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/拟时序/拟时序Pseudotime.pdf",width=30,height=30)
  p3<-plot_cell_trajectory(monocle_cds,color_by = "Pseudotime")
  print(p3)
  dev.off()


  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/拟时序/拟时序GROUP.pdf",width=30,height=30)
  p3<-plot_cell_trajectory(monocle_cds,color_by = "group")
  print(p3)
  dev.off()

  expressed_genes=row.names(subset(fData(monocle_cds),num_cells_expressed>=10))
  pseudotime_de <- differentialGeneTest(monocle_cds[expressed_genes,],
                                        fullModelFormulaStr = "~sm.ns(Pseudotime)")
  pseudotime_de <- pseudotime_de[order(pseudotime_de$qval), ]
  states_de <- differentialGeneTest(monocle_cds[expressed_genes,],
                                    fullModelFormulaStr = "~State")
  states_de <- states_de[order(states_de$qval), ]
  pseudotime_de<-readRDS("pseudotime_de.rds")
  time_gene<-pseudotime_de%>%pull(gene_short_name)%>%as.character()
  time_gene100<-top_n(pseudotime_de,n=100,desc(qval))%>%pull(gene_short_name)%>%as.character()
  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/拟时序/pseudotime_de.pdf",width=10,height=10)
  p3<-plot_pseudotime_heatmap(monocle_cds[time_gene100],num_clusters=4,show_rownames=T,return_heatmap=T)
  print(p3)
  dev.off()

  pseudotime_de %>% arrange(qval) %>% head() %>% select(gene_short_name) -> my_pseudotime_gene
  my_pseudotime_gene=my_pseudotime_gene[,1]
  my_pseudotime_gene
  pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/拟时序/TOP502.pdf",width=10,height=10)
  p3<- plot_pseudotime_heatmap(monocle_cds[gene_to_cluster,] add_annotation_col = ac,show_rownames = TRUE,return_heatmap = TRUE)
  print(p3)
  dev.off()

}

{
  setwd("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/拟时序/")
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(data.table)
  library(harmony)
  library(stringr)
  library(RhpcBLASctl)
  library(utils)
  library(CytoTRACE)
  library(RhpcBLASctl)
  library(utils)
  library(CytoTRACE)
    dat <- GetAssayData(sub_cancer_umap手动,assay = "RNA",layer = "counts")
    dat <- as.data.frame(dat)
    results <- CytoTRACE(dat, ncores = 8, subsamplesize = 1000)
    pheno <-  as.character(sub_cancer_umap手动$seurat_clusters)  
    names(pheno) <- colnames(sub_cancer_umap手动)
    saveRDS(results,"/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/拟时序/cytotrace_result_subcancer手动.RDS")
    results<-readRDS("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/拟时序/cytotrace_result_subcancer手动.RDS")
    pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/拟时序/cytotrace_dat2.pdf",width=20,height=20)
    p3 <- plotCytoGenes(results, numOfGenes = 10)
    print(p3)
    dev.off()
}  

{ 
    setwd("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/CIBERSORT/")
    library(GSVA) 
    library(tinyarray)
    library(tidyverse)
    library(FactoMineR)
    library(ggpubr)
    library(tidyr)
    geneSet <- read.csv("ss.csv",header = F,sep = ",")
    colnames(geneSet)
    rownames(geneSet)
    class(geneSet)
    geneSet=t(geneSet)
    colnames(geneSet)=geneSet[1,]
    geneSet=geneSet[-1,]
    a <- geneSet
    a <- a[1:nrow(a),]
    set <- colnames(a)
    l <- list()
    for (i in set) 
    { 
      x <- as.character(a[,i]) 
      x <- x[nchar(x)!=0] 
      x <- as.character(x) 
      l[[i]] <-x 
    }
    TaT1T2T3T4_umap_手动<-readRDS("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/注释/TaT1T2T3T4_umap_手动.rds")
    dat <- GetAssayData(TaT1T2T3T4_umap_手动,assay = "RNA",layer = "counts")
    re <- gsva(dat, l, method="ssgsea",
              mx.diff=FALSE, verbose=FALSE)
    sub_cancer_umap手动1<-readRDS("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/癌细胞亚群/sub_cancer_umap手动.rds")
    dat1 <- GetAssayData(sub_cancer_umap手动1,assay = "RNA",layer = "counts")
    re1 <- gsva(dat1, l, method="ssgsea",
              mx.diff=FALSE, verbose=FALSE)
     saveRDS(re1,"/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/CIBERSORT/re1.RDS")
     write.csv(re1,file ="/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/CIBERSORT/re1.csv")
    Stage<-matrix("",1,dim(re1)[2])
    Stage[1,grep("Ta",colnames(re1))]<-'Ta'
    Stage[1,grep("T1",colnames(re1))]<-'T1'
    Stage[1,grep("T2",colnames(re1))]<-'T2'
    Stage[1,grep("T3",colnames(re1))]<-'T3'
    Stage[1,grep("T4",colnames(re1))]<-'T4'
    Stage<-as.character(Stage)
    Stage<-factor(Stage,levels=unique(Stage))

    pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/CIBERSORT/ssGSEA_boxplot2.pdf",width=12,height=8)
    draw_boxplot_nopoint(re1,Stage)
    dev.off()

    pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/CIBERSORT/ssGSEA_veen.pdf",width=12,height=8)
    draw_pca(re,Stage)
    dev.off()
}

{
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(msigdbr)
    library(GSVA)
    library(pheatmap)
    library(ggplot2)
    setwd("/data5/shihong/fengmh/data/CWJ/BC/first_text/单细胞功能/机器学习补图")
    sub_umap<-readRDS("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/癌细胞亚群/sub_cancer_umap手动.rds")
    sub_umap[["RNA"]] <- as(sub_umap[["RNA"]], "Assay")
    expr_data <- as.matrix(sub_umap@assays$RNA@data)
    file_path<-dir("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/癌细胞亚群/T5_markers")
    stage<-sapply(file_path,function(x){
        a<-unlist(strsplit(x,split="_"))[1]
    })
    for(i in 1:length(stage)){
        markers_Fibroblast<-read.csv(paste0("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/癌细胞亚群/T5_markers/",file_path[i]))
        rownames(markers_Fibroblast)<-markers_Fibroblast[,1]
        markers<-markers_Fibroblast[,-1]
        exp_data<-expr_data[match(rownames(markers),rownames(expr_data)),]
        up <-rownames(markers[intersect(which(markers [,1]<0.05),which(markers [,2]>=0.25)),])
        down <-rownames(markers[intersect(which(markers [,1]<0.05),which(markers [,2]<=(-0.25))),])
        res1_total <- c(up,down)
        gs = bitr(res1_total, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
        head(gs)
        ego.bp = enrichGO(gene=gs$ENTREZID, 
                            OrgDb = org.Hs.eg.db,
                            ont= "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff= 0.05,
                            qvalueCutoff= 1,
                            readable= TRUE
                            )                    
        head(ego.bp)
        write.csv(ego.bp, file =paste0("差异分析/cancer_GO_",stage[i],".csv"))
        pdf(paste0("差异分析/cancer_GO_",stage[i],"_revise.pdf"), height =15, width = 25)
        print(barplot(ego.bp, showCategory=25,title=paste0("Cancer GO Term of ",stage[i]))+theme(axis.text.y = element_text(size = 20)))
        dev.off()


        kk1 <- enrichKEGG(gene = gs$ENTREZID, 
                            organism = "hsa",
                            keyType = 'ENTREZID',
                            pAdjustMethod = "BH", 
                            pvalueCutoff = 1 ,
                            use_internal_data = T
                        )
        write.csv(kk1, file =paste0("差异分析/cancer_KEGG_",stage[i],".csv"))
        pdf(paste0("差异分析/cancer_KEGG_",stage[i],"_revise.pdf"), height =15, width = 25)
        print(barplot(kk1, showCategory=25,title=paste0("Cancer KEGG Term of ",stage[i]))+theme(axis.text.y = element_text(size = 20)))
        dev.off()
        print(i)
    }

    setwd("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/hot_gene")
    setwd("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/hot_gene")
    library(clusterProfiler)
    library(enrichplot)
    copper<-read.csv("hot_gene_gmtfile.csv")
    file_path<-"/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/癌细胞亚群/T5_markers"
    file<-list.files(file_path)
    groups<-c("Ta","T1","T2","T3","T4")
    for(i in 1:length(file)){
        markers<-read.csv(paste0(file_path,"/",file[i]))
        Copper_geneList= markers$avg_log2FC
        names(Copper_geneList)= toupper(markers[,1])
        Copper_geneList=sort(Copper_geneList,decreasing = T)
        Copper_egmt <- GSEA(Copper_geneList, TERM2GENE=copper, 
                minGSSize = 1,
                pvalueCutoff = 1,
                verbose=FALSE)
        write.csv(Copper_egmt@result,file=paste0("GSEA/",groups[i],"/",groups[i],"_total_GSEA.csv"))
        pdf(paste0("GSEA/",groups[i],"/",groups[i],"_total_GSEA.pdf"),width=10,height=8)
        p<-gseaplot2(Copper_egmt,1,color="red",pvalue_table = T,title="",base_size=10,ES_geom="line") 
        print(p)
        dev.off()
        terms<-unique(copper[,1])
        for(j in 1:length(terms)){
            single_term<-copper[which(copper[,1]==terms[j]),]
            if(any(!is.na(match(single_term[,2],names(Copper_geneList))))){
                Copper_egmt_single <- GSEA(Copper_geneList, TERM2GENE=single_term, 
                        minGSSize = 1,
                        pvalueCutoff = 1,
                        verbose=FALSE)
                write.csv(Copper_egmt_single@result,file=paste0("GSEA/",groups[i],"/",groups[i],"_",terms[j],"_GSEA.csv"))
                pdf(paste0("GSEA/",groups[i],"/",groups[i],"_",terms[j],"_GSEA.pdf"),width=10,height=8)
                p<-gseaplot2(Copper_egmt_single,1,color="red",pvalue_table = T,title="",base_size=10,ES_geom="line") 
                print(p)
                dev.off()
            }
            print(j)
        }
        print(i)
    }
}




library(estimate)
library(dplyr)
library(tidyr)
library(tibble)
library(ggpubr)
library(RColorBrewer)
library(car)
setwd("/data5/shihong/fengmh/data/CWJ/BC/first_text/单细胞功能")

T_sce<-readRDS("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/癌细胞亚群/sub_cancer_umap手动.rds")
T_sce[["RNA"]] <- as(T_sce[["RNA"]], "Assay")
sub_cancer_tpm=as.data.frame(T_sce[["RNA"]]@counts) %>% apply(2,function(x){x/sum(x) * 10000})
estimate <- function(dat,pro){
  input.f=paste0(pro,'_estimate_input.txt')
  output.f=paste0(pro,'_estimate_gene.gct')
  output.ds=paste0(pro,'_estimate_score.gct')
  write.table(dat,file = input.f,sep = '\t',quote = F)
  filterCommonGenes(input.f=input.f,
                    output.f=output.f ,
                    id="GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds=output.ds,
                platform="illumina") 
  scores=read.table(output.ds,skip = 2,header = T,check.names = F)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
 library(stringr)
  rownames(scores)=str_replace_all(rownames(scores),'[.]','-') 
  write.csv(scores,file="Stromal_Immune_ESTIMATE.Score.csv")
  return(scores)
}
pro='T_'
estimate(sub_cancer_tpm,pro)


scores<-read.csv(file="/data5/shihong/fengmh/data/CWJ/BC/first_text/单细胞功能/estimate/Stromal_Immune_ESTIMATE.Score.csv")

a <- scores 
rownames(a)<-a[,1]
a<-a[,-1]
group <- substr(rownames(a), 1, 2)
unique(group)
a$group <- group
aa<-a
a <- a %>% rownames_to_column("sample")
b <- gather(a,key=category,value = score,-c(group,sample))
b1=as.data.frame(lapply(b$score,as.numeric)) %>% t() %>% as.data.frame()
b$score <- b1$V1

{
  head(b)
  dim(b)
  compare_group<-list()
  k<-1
  for(j in 1:length(unique(b$category))){
    for(i in 1:length(unique(b$group)))
    {
      compare_group[[k]]<-b[intersect(which(b$group==unique(b$group)[i]),which(b$category==unique(b$category)[j])),]
      print(k)
      k<-k+1
    }
  }
  
t_test_results<-matrix("",30,6)
colnames(t_test_results)<-c("category","group1","group2","F.p.value","t.p.value","t'.p.value")
k<-1
  for(i in 1:3){
    for( j in 1:10)
    {
      pairs <- combn(((5*i)-4):(5*i), 2)
      var_results<-var.test(compare_group[[pairs[1,j]]]$score, compare_group[[pairs[2,j]]]$score)
      t_test_results[k,4]<-var_results$statistic
      if(var_results$p.value<0.05 && !(var_results$conf.int[1]<=1 && var_results$conf.int[2] >=1)){
        text_result<-t.test(compare_group[[pairs[1,j]]]$score,compare_group[[pairs[2,j]]]$score,var.equal=FALSE)
        t_test_results[k,6]<-c(text_result$p.value)
      }else if(var_results$p.value>0.05 && (var_results$conf.int[1]<=1 && var_results$conf.int[2] >=1)){
        text_result<-t.test(compare_group[[pairs[1,j]]]$score,compare_group[[pairs[2,j]]]$score,var.equal=TRUE)
        t_test_results[k,5]<-c(text_result$p.value)
      }else{
        t_test_results[k,5:6]<-c("","")
      }
      t_test_results[k,1:3]<-c(compare_group[[pairs[1,j]]][1,3],compare_group[[pairs[1,j]]][1,2],compare_group[[pairs[2,j]]][1,2])
      print(k)
      k<-k+1
    }
  }
save(t_test_results,file="estimate/t_test_result.Rdata")
write.csv(t_test_results,file="estimate/t_test_result.csv")
}

{
  pdf("estimate/scores_new.pdf",width=10)
  plot <- ggboxplot(b, x = "group", y = "score",
            color = "group", palette = "jco",
            facet.by = "category", short.panel.labs = FALSE)
  print(plot)
  dev.off()
}


{
    setwd("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/CIBERSORT")
    library(estimate)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggpubr)
    library(RColorBrewer)
    re1_csv<-read.csv("re1.csv")
    dim(re1_csv)
    rownames(re1_csv)<-re1_csv[,1]
    re1_csv<-re1_csv[,-1]
    re1_data<-t(re1_csv)
    a <- as.data.frame(re1_data) 
    rownames(a)
    group <- substr(rownames(a), 1, 2)
    unique(group)
    a$group <- group
    aa<-a
    a <- a %>% rownames_to_column("sample")
    b <- gather(a,key=celltype,value = expression,-c(group,sample))

    compare_means(expression ~ celltype, data = b,ref.group = ".all.",method="t.test")
    my_comparisons<-list(c("Ta","T1"),c("T1","T2"),c("T2","T3"),c("T3","T4"))
    pdf("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/CIBERSORT/cibersort_new.pdf",width=20)
    plot <- ggboxplot(b, x = "celltype", y = "expression",
                color = "group", palette = "jco",bxp.errorbar.width=0.4)+ 
                stat_compare_means(label = "p.signif", method = "t.test",
                                    ref.group = ".all.")+      
                font("xy.text", size = 4,  face = "bold")
    print(plot)
    dev.off()
    celltypes<-unique(b$celltype)
    groups<-unique(b$group)
    
    {
        library(ggplot2)
        library(scales)
        b$group<-factor(b$group,levels=c("Ta","T1","T2","T3","T4"))
        pdf("ssGSEA_boxplot_new_order.pdf",width=15)
            p<-ggplot(data=b, aes(x=celltype, y=expression)) + 
            geom_boxplot(aes(fill = group), width = 0.8,outlier.size=0.2) +
            theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=10),
                                plot.title = element_text(hjust = 0.5)) +
            labs(title="ssGSEA Score of T Stage")
        print(p)
        dev.off()



    celltype_data_list<-list()
    for(i in 1:length(groups)){
        celltype_data_list[[i]]<-b[which(b$group==groups[i]),]
        compare_means(expression ~ celltype, data = celltype_data_list[[i]],ref.group = ".all.",method="t.test")
        if(grepl("",celltypes[i])){
        cellname<-gsub(" ","_",groups[i])
        }else{
        cellname<-groups[i]
        }
        pdf(paste0("celltype_ssGSEA/",cellname,".pdf"),width=20)
        plot <- ggboxplot(celltype_data_list[[i]], x = "celltype", y = "expression",
                color = "celltype", palette = c("#FF0000", "#FF4500", "#FF8C00",
                "#FFD700", "#ADFF2F", "#7FFF00", "#00FF00", "#00FA9A", "#00FFFF",
                    "#1E90FF", "#0000FF", "#8A2BE2", "#9400D3", "#9932CC", "#8B008B",
                    "#800080", "#8B0000", "#A52A2A", "#B22222", "#DC143C", "#FF0000",
                    "#FF4500", "#FF8C00", "#FFD700", "#ADFF2F", "#7FFF00", "#00FF00", "#800080"),
                    bxp.errorbar.width=0.4)+ 
                stat_compare_means(label = "p.signif", method = "t.test",
                                    ref.group = ".all.")+      
                font("xy.text", size = 4,  face = "bold")
        print(plot)
        dev.off()
        print(i)
    }
    } 
}


