# Lym subtypes

#1. library packages
library(dplyr)
library(Seurat)
library(ggsci)
library(viridis)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(monocle)

data<-readRDS("T_subType.RDS")

#2. datasets dimension reduction and re-cluster

data <- RunPCA(object = data,verbose = FALSE)
ElbowPlot(data,ndims = 50) 
dims=1:20
data <- FindNeighbors(data,dim=dims)
data <- FindClusters(data,resolution = 0.8)
data <- RunUMAP (data,reduction="pca", dims = dims)
DimPlot(data,label = T,pt.size = 1)

#3. marker gene expression

gene<-c(
"Cd3d","Cd3e","Cd3g",
"Cd4","Cd8a","Cd8b1", # T cell markers
"Ncr1", # NK cell markers
"Tcf7","Ccr7","Lef1", # T Naive markers
"S100a4","Cd40lg","Cxcr6","Cxcr3", # T Memory markers
"Nkg7","Fasl","Klrd1","Klrg1","Ccl5","Gzma","Gzmb","Gzmk", # T effector markers
"Mki67","Stmn1", # Proliferation markers
"Trdc" # gdT cell markers
)

DefaultAssay(data)<-"SCT"
FeaturePlot(data,features=gene,col = c("grey","red"))

#4. subtype annotation
ann <- read.csv("cluster_subtype.csv",header=T,row.names=1)
data$subType <- ann[as.character(data$seurat_clusters),]
data$subType <- factor(data$subType,levels=c("CD4_Tn",
                                             "CD4_Tm",
											 "CD4_Tpro",
											 "CD8_Tn",
											 "CD8_Tem",
											 "CD8_TEMRA",
											 "CD8_Tpro",
											 "γδT",
											 "NK"))

color<-pal_d3()(10)[c(10,6,8,9,5,3,7,1,2,4)]

DimPlot(data,group.by="subType")+scale_color_manual(values=color)

#5. draw heatmap plot

dat<-AverageExpression(data)
pheatmap(dat$SCT[gene,],scale="row",cluster_rows=F,cluster_cols=F,angle_col=45,border_color=NA,cellwidth=20,cellheight=20,color=viridis(100),gaps_col=c(3,7,8),gaps_row=c(7,10,14,22,24))

#6. save metadata and Lym.RDS

meta <- data@meta.data
wriet.csv(meta,"Lym_meta.csv")
saveRDS(data,"Lym_anno.RDS")

#7. draw the proportions of each group in different subtypes

dtype<-as.data.frame(table(data@meta.data[,c('subType','sample_type')]))
p1 <- ggplot(subset(dtype,sample_type=="Control"),aes(x=sample_type,y=Freq,fill=subType))+
		geom_bar(stat = "identity",position = 'fill')+
		theme_classic()+
		scale_fill_manual(values=color)+
		coord_polar(theta="y")+
		ggtitle("Control (n=184)")
p2 <- ggplot(subset(dtype,sample_type=="Model"),aes(x=sample_type,y=Freq,fill=subType))+
		geom_bar(stat = "identity",position = 'fill')+
		theme_classic()+
		scale_fill_manual(values=color)+
		coord_polar(theta="y")+
		ggtitle("Model (n=2277)")
p3 <- ggplot(subset(dtype,sample_type=="ART"),aes(x=sample_type,y=Freq,fill=subType))+
		geom_bar(stat = "identity",position = 'fill')+
		theme_classic()+
		scale_fill_manual(values=color)+
		coord_polar(theta="y")+
		ggtitle("ART")
p1|p2|p3

#8. Lym subtypes different expressed gene (DEGs) 

Idents(data)<-"subType"
subType = unique(data$subType)
dat = data.frame(gene="",logFC=0,adj.P=0,Celltype="",State="",Com="")
for (i in subType) {
 subdata=subset(data,idents=i)
 Idents(subdata)<-"sample_type"
 marker<-FindMarkers(subdata,ident.1="Model",ident.2="Control",min.pct=0.25,logfc.threshold=0.25)
 submarker = data.frame(gene=row.names(marker),logFC=marker$avg_log2FC,adj.P=marker$p_val_adj,Celltype=i,State=ifelse(marker$p_val_adj<0.05,ifelse(marker$avg_log2FC>0.25,"Up",ifelse(marker$avg_log2FC< -0.25,"Down","No")),"No"),Com="ModvsCon") # Model vs Control
 dat=rbind(dat,submarker)
 
 marker<-FindMarkers(subdata,ident.1="ART",ident.2="Control",min.pct=0.25,logfc.threshold=0.25)
 submarker = data.frame(gene=row.names(marker),logFC=marker$avg_log2FC,adj.P=marker$p_val_adj,Celltype=i,State=ifelse(marker$p_val_adj<0.05,ifelse(marker$avg_log2FC>0.25,"Up",ifelse(marker$avg_log2FC< -0.25,"Down","No")),"No"),Com="ARTvsCon") # ART vs Control
 dat=rbind(dat,submarker)
}
dat<-dat[-1,]
dat<-subset(dat,State != "No")
write.csv(dat,"Lym-subType-Type-marker.csv",quote=F,row.names=F)

#9. Cd8+ T cells trajectory analysis

mat<-subset(data,subtype %in% c("CD8_Tn", "CD8_Tem", "CD8_TEMRA", "CD8_Tpro"))
DefaultAssay(mat) <- "RNA"
expr_matrix <- as(as.matrix(mat@assays$RNA@counts), 'sparseMatrix')
p_data <- mat@meta.data
f_data <- data.frame(gene_short_name = row.names(mat),row.names = row.names(mat))
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix,phenoData = pd,featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1)
tmp <- cds
deg <- c(unique(subset(dat,Celltype %in% c("CD8_Tn", "CD8_Tem", "CD8_TEMRA", "CD8_Tpro") & State != "Stable")[,'gene']))
cds <- setOrderingFilter(cds, deg)
table(cds@featureData@data[["use_for_ordering"]])

cds <- reduceDimension(cds, max_components = 2,method = 'DDRTree')
cds <- orderCells(cds)

saveRDS(cds,"CD8_union_monocle.RDS")

color<-pal_d3()(10)[c(10,6,8,9,5,3,7,1,2,4)][c(4:7)]

plot_cell_trajectory(cds,color_by="Pseudotime")
plot_cell_trajectory(cds,color_by="State")+scale_color_jama()
plot_cell_trajectory(cds,color_by="sample_type",cell_size=1)+scale_color_manual(values=c("#3B4992","#EE0000","#008B45"))
plot_cell_trajectory(cds,color_by="subtype",cell_size=1)+scale_color_manual(values=color)


diff_test_res <- differentialGeneTest(cds[deg,], fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
plot <- plot_pseudotime_heatmap(cds[sig_gene_names,],
                num_clusters = 4,
                cores = 1,
                show_rownames = F,
                return_heatmap=TRUE)
OrderGene<-as.data.frame(cutree(plot$tree_row,k=4))
colnames(OrderGene)<-"Cluster"
c1<-row.names(subset(OrderGene,Cluster==1))
c2<-row.names(subset(OrderGene,Cluster==2))
c3<-row.names(subset(OrderGene,Cluster==3))
c4<-row.names(subset(OrderGene,Cluster==4))

eg<-bitr(c1,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mm.eg.db")
genelist<-eg$ENTREZID
genelist<-unique(genelist)
c1_bp<-enrichGO(genelist,OrgDb = org.Mm.eg.db,ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,keyType = "ENTREZID",readable=T)

eg<-bitr(c2,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mm.eg.db")
genelist<-eg$ENTREZID
genelist<-unique(genelist)
c2_bp<-enrichGO(genelist,OrgDb = org.Mm.eg.db,ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,keyType = "ENTREZID",readable=T)

eg<-bitr(c3,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mm.eg.db")
genelist<-eg$ENTREZID
genelist<-unique(genelist)
c3_bp<-enrichGO(genelist,OrgDb = org.Mm.eg.db,ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,keyType = "ENTREZID",readable=T)

eg<-bitr(c4,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mm.eg.db")
genelist<-eg$ENTREZID
genelist<-unique(genelist)
c4_bp<-enrichGO(genelist,OrgDb = org.Mm.eg.db,ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,keyType = "ENTREZID",readable=T)


gene <- 
c(
"Ifng",
"Gzma",
"Klrd1",
"Gzmb"
)

color<-pal_d3()(10)[c(10,6,8,9,5,3,7,1,2,4)][c(4:7)]
plot_genes_in_pseudotime(cds[gene,],ncol=2,color_by="subType")+scale_color_manual(values=color)+theme(text=element_text(size=18))


