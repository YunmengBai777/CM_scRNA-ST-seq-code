# 1.scRNA Rawdata process

# 1.1 library packages

library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(ggsci)
options(future.globals.maxSize = 60000 * 1024^2)

# 1.2 load 10X genomics datasets

setwd("/public/Count/scRNA/brain_cm/1.cluster/test")

samplename = c(
"/BC1/outs/filtered_feature_bc_matrix", 
"/BC2/outs/filtered_feature_bc_matrix", 
"/BC3/outs/filtered_feature_bc_matrix", 
"/BM1/outs/filtered_feature_bc_matrix",
"/BM2/outs/filtered_feature_bc_matrix",     
"/BM3/outs/filtered_feature_bc_matrix", 
"/BA1/outs/filtered_feature_bc_matrix",
"/BA2/outs/filtered_feature_bc_matrix",
"/BA3/outs/filtered_feature_bc_matrix")

names(samplename) = c('BC1','BC2','BC3','BM1','BM2','BM3','BA1','BA2','BA3')
print(samplename)
proj <- list()

print("1.2 load 10X genomics datasets done!")

# 1.3 Perform SCTransform on multiple samples

for(i in 1:length(samplename)){
  print(names(samplename[i]))	
  counts <- Read10X(data.dir = samplename[i])
  newproj <- CreateSeuratObject(counts = counts, min.cells = 3, project = names(samplename[i]))
  print(newproj) # print datasets
  newproj$sample <- names(samplename[i])
  newproj[["percent.mt"]] <- PercentageFeatureSet(object = newproj, pattern = "mt-")
  newproj[["percent.hba"]] <- PercentageFeatureSet(object = newproj, pattern = "Hba-")
  newproj[["percent.hbb"]] <- PercentageFeatureSet(object = newproj, pattern = "Hbb-")
  #saveRDS(newproj,file=paste0(names(samplename[i]),".rds"))
  newproj <- subset(x = newproj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 500 & nCount_RNA < 30000 & percent.mt < 15)
  newproj<-SCTransform(newproj,return.only.var.genes = FALSE,assay = "RNA",verbose = FALSE)
  proj[[names(samplename[i])]] <- newproj
}

print("1.3 Perform SCTransform on multiple samples done!")

# 1.4 Draw a violin plot of brain_CM_scRNA_before_filter

# pdf("brain_CM_scRNA_before_filter.pdf",height = 5,width = 15)
# scc_merge <- merge(PC1,c(PC2,PC3,PM1,PM2,PM3,PA1,PA2,PA3,
#                          LC1,LC2,LC3,LM1,LM2,LM3,LA1,LA2,LA3,
#                          SC1,SC2,SC3,SM1,SM2,SM3,SA1,SA2,SA3))
# VlnPlot(scc_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hb"), ncol = 5,pt.size = 0) 
# dev.off()

# 1.5 Use CCA to integrate multiple samples

object_list = proj
selfeatures <- SelectIntegrationFeatures(object.list = object_list, nfeatures = 3000)
scc.list <- PrepSCTIntegration(object.list = object_list, anchor.features = selfeatures, verbose = FALSE)
scc.anchors <- FindIntegrationAnchors(object.list = scc.list, normalization.method = "SCT",anchor.features = selfeatures, verbose = FALSE)
scc_integrated <- IntegrateData(anchorset = scc.anchors, normalization.method = "SCT",verbose = FALSE)
remove(object_list,scc.anchors,scc.list)
print(scc_integrated) 
print(head(scc_integrated@meta.data))

# 1.6 Draw a violin plot of brain_CM_scRNA_after_filter

pdf("brain_CM_scRNA_after_filter.pdf",height = 5,width = 10)
VlnPlot(scc_integrated, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.hba","percent.hbb"), ncol = 5,pt.size = 0)
dev.off()

saveRDS(scc_integrated,file="brain_CM_scRNA_integrated_raw.rds")


# 1.7 Perform dimension reduction and UMAP visulations

scc_integrated <- RunPCA(object = scc_integrated,verbose = FALSE)
scc_integrated <- FindNeighbors(scc_integrated,dim=1:20)
scc_integrated <- FindClusters(scc_integrated,resolution = 0.8)
scc_integrated <- RunUMAP (scc_integrated,reduction="pca", dims = 1:20)
print(scc_integrated)


# 1.8 Add meta data information

scc_integrated$data_type = "scRNA"

sample_type <- scc_integrated$sample
sample_type <- as.character(sample_type)
sample_type[scc_integrated$sample %in% c('BC1','BC2','BC3')] <- 'Control'
sample_type[scc_integrated$sample %in% c('BM1','BM2','BM3')] <- 'Model'
sample_type[scc_integrated$sample %in% c('BA1','BA2','BA3')] <- 'ART'
table(sample_type)

sample_type <- factor(sample_type,levels = c("Control","Model","ART"))
scc_integrated$sample_type <- sample_type

#  1.9 saveRDS and plot the Umap
saveRDS(scc_integrated,"brain_CM_scRNA_integrated_Umap.rds")

brain_CM_scRNA_integrated_Umap <- DimPlot(scc_integrated,label = TRUE,reduction = "umap")
pdf("brain_CM_scRNA_integrated_Umap.pdf",height = 10,width = 10)
print(brain_CM_scRNA_integrated_Umap)
dev.off()
print("done!")


#  1.10 cell type annotation

scc_integrated <- readRDS("brain_CM_scRNA_integrated_Umap.rds")

rds@active.assay <- "SCT"

gene <- c("Cx3cr1","Aldoc","Ccdc153","Kcnj8","Sox10","Sox11",
          "Cd79a","Cldn5","Mbp","Ttr","Cd33","Cd209a","Dcn",
          "Cd68","Alas2","Slc47a1","Plac8","Pf4","S100a9")

pdf("feature_plot.pdf",width = 15,height = 15)
FeaturePlot(rds,features=gene)
dev.off()

pdf("umap_split.pdf",width = 10,height = 5)
DimPlot(rds, split.by ="sample_type" ,pt.size=0.25,label = TRUE)
dev.off()

scc_integrated2 <- RenameIdents(scc_integrated, `0` = "MG", `1` = "MG", `2` = "EC", `3` = "MG", `4` = "MG",
 `5` = "ASC", `6` = "MG", `7` = "CPC", `8`="ASC",`9` = "T", `10`="EC",`11` = "MNC", 
 `12` = "NEUT",`13` = "EC",`14` = "OLG",`15` = "T",`16` = "OLG",`17` = "EPC",`18` = "EC",`19` = "Hbb_VC",
 `20` = "EC",`21` = "ASC",`22` = "MNC",`23` = "Hbb_VC",`24` = "ASC",`25` = "VSMC",`26` = "EPC",`27` = "Neu",
 `28` = "ABC",`29` = "MAC",`30` = "MG",`31` = "EC",`32` = "ASC",`33` = "B",`34` = "unknow",`35` = "Fibroblast",
 `36` = "MAC",`37` = "T")

pdf("umap_named.pdf",width = 10,height = 10)
DimPlot(rds2, pt.size=1,label = TRUE)
dev.off()

pdf("umap_split.pdf",width = 16,height = 8)
DimPlot(rds2, split.by ="sample_type" ,pt.size=0.25,label = TRUE)
dev.off()

saveRDS(rds,"CM_scRNA_names.rds")

# 1.11 cell type and stastics

scc_integrated <-readRDS("scc_integrated_anno.rds")
type <- as.character(scc_integrated@active.ident)
# type$cell_type <- type$`scc_integrated@active.ident`
type[type == "Neu"] = "Neuron"
type[type == "OLG"] = "Oligo"
type[type == "ASC"] = "Astro"
type[type == "EC"] = "Endo"
type[type == "Fibroblast"] = "Fibro"
type[type == "MNC"] = "Mono"
type[type == "MG"] = "Micro"
type[type == "MAC"] = "Macro"
type[type == "NEUT"] = "Neutro"
type[type == "Hbb_VC"] = "RBC"
type[type == "T"] = "T/NK cell"
type[type == "B"] = "B cell"
type[type == "unknow"] = "Unknow"
table(type)
cell_type <- factor(type,levels = c("Neuron","Oligo","Astro","Endo","VSMC","RBC","Fibro","ABC","EPC","CPC",
                                  "Mono","Micro","Macro","Neutro","T/NK cell","B cell","Unknow"))
table(cell_type)

scc_integrated$cell_type <- cell_type
scc_integrated@active.ident <- scc_integrated$cell_type

DimPlot(scc_integrated)

DimPlot(scc_integrated,label = T,label.size = 5,cols = c("Neuron" = '#E95C59',
                                                             "Oligo" = '#53A85F',
                                                             "Astro" = '#F1BB72',
                                                             "Endo" = '#F3B1A0',
                                                             "VSMC" = '#D6E7A3',
                                                             "RBC" = '#57C3F3',
                                                             "Fibro" = '#E63863',
                                                             "Endo" = '#E4C755',
                                                             "ABC" = '#E59CC4',
                                                             "EPC" = '#AB3282',
                                                             "CPC" = '#23452F',
                                                             "Mono" = '#BD956A',
                                                             "Micro" = 'orange',
                                                             "Macro" = '#58A4C3',
                                                             "Neutro" = "#00BFC4",
                                                             "T/NK cell" = 'red',
                                                             "B cell" = '#8C549C',
                                                             "Unknow" = "black"
                                                             ),group.by = "cell_type") + 
  theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) # + labs(title =  paste0("Cell annotation"))


DimPlot(scc_integrated,label = F,label.size = 5,cols = c("Neuron" = '#E95C59',
                                                         "Oligo" = '#53A85F',
                                                         "Astro" = '#F1BB72',
                                                         "Endo" = '#F3B1A0',
                                                         "VSMC" = '#D6E7A3',
                                                         "RBC" = '#57C3F3',
                                                         "Fibro" = '#E63863',
                                                         "Endo" = '#E4C755',
                                                         "ABC" = '#E59CC4',
                                                         "EPC" = '#AB3282',
                                                         "CPC" = '#23452F',
                                                         "Mono" = '#BD956A',
                                                         "Micro" = 'orange',
                                                         "Macro" = '#58A4C3',
                                                         "Neutro" = "#00BFC4",
                                                         "T/NK cell" = 'red',
                                                         "B cell" = '#8C549C',
                                                         "Unknow" = "black"
),group.by = "cell_type",split.by = "sample_type") + 
  theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) # + labs(title =  paste0("Cell annotation"))

saveRDS(scc_integrated,file="scc_integrated_anno.rds")

# 1.12 save metadata info

Metadata <- scc_integrated@meta.data # 91455
write.csv(Metadata,"CM_scRNA_metadata.csv")

# 1.13.Analyzing the proportions of cell types in differ group

table(scc_integrated$sample_type)
# Control   Model     ART 
# 21684   41602   28169

table(scc_integrated$sample)
# BA1   BA2   BA3   BC1   BC2   BC3   BM1   BM2   BM3 
# 9168  9939  9062  8773  6824  6087 10522 21487  9593

table(scc_integrated$cell_type)

dfsam <- as.data.frame(table(scc_integrated$sample_type,scc_integrated$cell_type))

control <- dfsam[dfsam$Var1 == "Control",]
myLabel = as.vector(control$Var2)
myLabel = paste(myLabel, "(", round(control$Freq / sum(control$Freq) * 100, 2), "%)"
                , sep = "")  

p1 <- ggplot(control, aes(x = "", y = Freq,fill = Var2),colour = cols ) + 
  geom_bar(stat = "identity")+
  # coord_polar(theta = "y") + 
  theme_bw() + 
  labs(x = "", y = "")+
  theme(axis.ticks = element_blank())+
  theme(axis.text.x = element_blank()) + 
  theme(panel.grid=element_blank()) +    
  theme(panel.border=element_blank()) + 
  theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + 
  scale_fill_manual(values = c("Neuron" = '#E95C59',
                               "Oligo" = '#53A85F',
                               "Astro" = '#F1BB72',
                               "Endo" = '#F3B1A0',
                               "VSMC" = '#D6E7A3',
                               "RBC" = '#57C3F3',
                               "Fibro" = '#E63863',
                               "ABC" = '#E59CC4',
                               "EPC" = '#AB3282',
                               "CPC" = '#23452F',
                               "Mono" = '#BD956A',
                               "Micro" = 'orange',
                               "Macro" = '#58A4C3',
                               "Neutro" = "#00BFC4",
                               "T/NK cell" = 'red',
                               "B cell" = '#8C549C',
                               "Unknow" = "black"), 
                    name = "",labels = myLabel)

p1


model <- dfsam[dfsam$Var1 == "Model",]
myLabel = as.vector(model$Var2)
myLabel = paste(myLabel, "(", round(model$Freq / sum(model$Freq) * 100, 2), "%)"
                , sep = "") 

p2 <- ggplot(model, aes(x = "", y = Freq,fill = Var2),colour = cols ) + 
  geom_bar(stat = "identity")+
  # coord_polar(theta = "y") + 
  theme_bw() + 
  labs(x = "", y = "")+
  theme(axis.ticks = element_blank())+
  theme(axis.text.x = element_blank()) + 
  theme(panel.grid=element_blank()) +   
  theme(panel.border=element_blank()) + 
  theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + 
  scale_fill_manual(values = c("Neuron" = '#E95C59',
                               "Oligo" = '#53A85F',
                               "Astro" = '#F1BB72',
                               "Endo" = '#F3B1A0',
                               "VSMC" = '#D6E7A3',
                               "RBC" = '#57C3F3',
                               "Fibro" = '#E63863',
                               "ABC" = '#E59CC4',
                               "EPC" = '#AB3282',
                               "CPC" = '#23452F',
                               "Mono" = '#BD956A',
                               "Micro" = 'orange',
                               "Macro" = '#58A4C3',
                               "Neutro" = "#00BFC4",
                               "T/NK cell" = 'red',
                               "B cell" = '#8C549C',
                               "Unknow" = "black"), 
                    name = "",labels = myLabel)
p2

treatment <- dfsam[dfsam$Var1 == "ART",]
myLabel = as.vector(treatment$Var2)
myLabel = paste(myLabel, "(", round(treatment$Freq / sum(treatment$Freq) * 100, 2), "%)"
                , sep = "") 

p3 <- ggplot(treatment, aes(x = "", y = Freq,fill = Var2)) + 
  geom_bar(stat = "identity")+
  # coord_polar(theta = "y") + 
  theme_bw() + 
  labs(x = "", y = "")+
  theme(axis.ticks = element_blank())+
  theme(axis.text.x = element_blank()) + 
  theme(panel.grid=element_blank()) + 
  theme(panel.border=element_blank()) + 
  theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + 
  scale_fill_manual(values = c("Neuron" = '#E95C59',
                               "Oligo" = '#53A85F',
                               "Astro" = '#F1BB72',
                               "Endo" = '#F3B1A0',
                               "VSMC" = '#D6E7A3',
                               "RBC" = '#57C3F3',
                               "Fibro" = '#E63863',
                               "ABC" = '#E59CC4',
                               "EPC" = '#AB3282',
                               "CPC" = '#23452F',
                               "Mono" = '#BD956A',
                               "Micro" = 'orange',
                               "Macro" = '#58A4C3',
                               "Neutro" = "#00BFC4",
                               "T/NK cell" = 'red',
                               "B cell" = '#8C549C',
                               "Unknow" = "black"), 
                    name = "",labels = myLabel)
p3

p1+p2+p3 #5*15

# 3.2 perform DEG analysis

scRNA_markers = FindAllMarkers(scc_integrated,only.pos =T,logfc.threshold = 0.25)

# scRNA_marker <- read.csv("scRNA_markers.csv")
# scRNA_markers_top50 <- scRNA_markers %>% group_by(cluster) %>% top_n(n =50, wt = avg_log2FC)

scRNA_Markers_num <- as.data.frame(table(scRNA_marker_all$cluster))
colnames(scRNA_Markers_num) <- c("cell_type","Markers_num")
head(scRNA_Markers_num)
table(scRNA_Markers_num$cell_type) 

scRNA_Markers_num$cell_type <- factor(scRNA_Markers_num$cell_type,levels = c("Neuron","Oligo","Astro","Endo","VSMC","RBC","Fibro","ABC","EPC","CPC",
                                    "Mono","Micro","Macro","Neutro","T/NK cell","B cell","Unknow"))

theme <- theme(panel.background = element_blank(), 
               panel.grid.major.x = element_line(colour = "black"), 
               axis.line.x = element_line(colour = "black"),
               axis.title.y = element_blank())

ggplot(data = scRNA_Markers_num) + geom_point(aes(x = cell_type, y = Markers_num,color = cell_type), size = 5) +
  geom_segment(aes(x = cell_type, y = 0,xend = cell_type, yend = Markers_num)) + 
  coord_flip() + theme + ylab("Cell markers number") + 
  scale_color_manual(values = c("Neuron" = '#E95C59',
                               "Oligo" = '#53A85F',
                               "Astro" = '#F1BB72',
                               "Endo" = '#F3B1A0',
                               "VSMC" = '#D6E7A3',
                               "RBC" = '#57C3F3',
                               "Fibro" = '#E63863',
                               "ABC" = '#E59CC4',
                               "EPC" = '#AB3282',
                               "CPC" = '#23452F',
                               "Mono" = '#BD956A',
                               "Micro" = 'orange',
                               "Macro" = '#58A4C3',
                               "Neutro" = "#00BFC4",
                               "T/NK cell" = 'red',
                               "B cell" = '#8C549C',
                               "Unknow" = "black"), 
                    name = "") + scale_y_continuous(limits =c(0,1800,500) ,expand = c(0,0))

# 1.14 select DEGs of each cell type

scRNA_marker <- read.csv("scRNA_markers.csv")
table(scRNA_marker$cluster)
scRNA_marker$cluster

Unknow_genes <- scRNA_marker[scRNA_marker$cluster %in% "Unknow",]$gene
ABC_genes <- scRNA_marker[scRNA_marker$cluster %in% "ABC",]$gene
Astro_genes <- scRNA_marker[scRNA_marker$cluster %in% "Astro",]$gene
B_cell_genes <- scRNA_marker[scRNA_marker$cluster %in% "B cell",]$gene
CPC_genes <- scRNA_marker[scRNA_marker$cluster %in% "CPC",]$gene
Endo_genes <- scRNA_marker[scRNA_marker$cluster %in% "Endo",]$gene
EPC_genes <- scRNA_marker[scRNA_marker$cluster %in% "EPC",]$gene
Fibro_genes <- scRNA_marker[scRNA_marker$cluster %in% "Fibro",]$gene
Macro_genes <- scRNA_marker[scRNA_marker$cluster %in% "Macro",]$gene
Micro_genes <- scRNA_marker[scRNA_marker$cluster %in% "Micro",]$gene
Mono_genes <- scRNA_marker[scRNA_marker$cluster %in% "Mono",]$gene
Neuron_genes <- scRNA_marker[scRNA_marker$cluster %in% "Neuron",]$gene
Neutro_genes <- scRNA_marker[scRNA_marker$cluster %in% "Neutro",]$gene
Oligo_genes <- scRNA_marker[scRNA_marker$cluster %in% "Oligo",]$gene
RBC_genes <- scRNA_marker[scRNA_marker$cluster %in% "RBC",]$gene
T_NK_cell_genes <- scRNA_marker[scRNA_marker$cluster %in% "T/NK cell",]$gene
VSMC_genes <- scRNA_marker[scRNA_marker$cluster %in% "VSMC",]$gene

# 1.15 perform GO enrichment of  DEGs of each cell type
Unknow_genes_GO <- enrichGO(gene = Unknow_genes, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP",pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
ABC_genes_GO <- enrichGO(gene = ABC_genes, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP",pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
Astro_genes_GO <- enrichGO(gene = Astro_genes, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP",pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
B_cell_genes_GO <- enrichGO(gene = B_cell_genes, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP",pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
CPC_genes_GO <- enrichGO(gene = CPC_genes, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
Endo_genes_GO <- enrichGO(gene = Endo_genes, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP",pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
EPC_genes_GO <- enrichGO(gene = EPC_genes, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP",pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)   
Fibro_genes_GO <- enrichGO(gene = Fibro_genes, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP",pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
Macro_genes_GO <- enrichGO(gene = Macro_genes, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP",pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
Micro_genes_GO <- enrichGO(gene = Micro_genes, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP",pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
Mono_genes_GO <- enrichGO(gene = Mono_genes, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP",pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
Neuron_genes_GO <- enrichGO(gene = Neuron_genes, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP",pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
Neutro_genes_GO <- enrichGO(gene = Neutro_genes, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP",pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)  
Oligo_genes_GO <- enrichGO(gene = Oligo_genes, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP",pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
RBC_genes_GO <- enrichGO(gene = RBC_genes, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
T_NK_cell_genes_GO <- enrichGO(gene = T_NK_cell_genes, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP",pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
VSMC_genes_GO <- enrichGO(gene = VSMC_genes, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP",pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)

# 1.16 select top20 GO pathways of each cell type

Unknow_genes_GO_top20 <- Unknow_genes_GO@result[c(1:20),c("Description","pvalue")]
B_cell_genes_GO_top20 <- B_cell_genes_GO@result[c(1:20),c("Description","pvalue")]
T_NK_cell_genes_GO_top20 <- T_NK_cell_genes_GO@result[c(1:20),c("Description","pvalue")]
Neutro_genes_GO_top20<- Neutro_genes_GO@result[c(1:20),c("Description","pvalue")]
Macro_genes_GO_top20<- Macro_genes_GO@result[c(1:20),c("Description","pvalue")] 
Micro_genes_GO_top20<- Micro_genes_GO@result[c(1:20),c("Description","pvalue")]
Mono_genes_GO_top20<- Mono_genes_GO@result[c(1:20),c("Description","pvalue")]
CPC_genes_GO_top20<- CPC_genes_GO@result[c(1:20),c("Description","pvalue")]
EPC_genes_GO_top20<- EPC_genes_GO@result[c(1:20),c("Description","pvalue")]
ABC_genes_GO_top20 <- ABC_genes_GO@result[c(1:20),c("Description","pvalue")]
Fibro_genes_GO_top20<- Fibro_genes_GO@result[c(1:20),c("Description","pvalue")]
RBC_genes_GO_top20<- RBC_genes_GO@result[c(1:20),c("Description","pvalue")]
VSMC_genes_GO_top20<- VSMC_genes_GO@result[c(1:20),c("Description","pvalue")]
Endo_genes_GO_top20<- Endo_genes_GO@result[c(1:20),c("Description","pvalue")]
Astro_genes_GO_top20<- Astro_genes_GO@result[c(1:20),c("Description","pvalue")]
Oligo_genes_GO_top20<- Oligo_genes_GO@result[c(1:20),c("Description","pvalue")]
Neuron_genes_GO_top20 <- Neuron_genes_GO@result[c(1:30),c("Description","pvalue")]

Unknow_genes_GO_top20$subtype = "Unknow"
B_cell_genes_GO_top20$subtype = "B cell"
T_NK_cell_genes_GO_top20$subtype = "T/NK_cell"
Neutro_genes_GO_top20$subtype = "Neutro"
Macro_genes_GO_top20$subtype = "Macro"
Micro_genes_GO_top20$subtype = "Micro"
Mono_genes_GO_top20$subtype = "Mono"
CPC_genes_GO_top20$subtype = "CPC"
EPC_genes_GO_top20$subtype = "EPC"
ABC_genes_GO_top20$subtype = "ABC"
Fibro_genes_GO_top20$subtype = "Fibro"
RBC_genes_GO_top20$subtype = "RBC"
VSMC_genes_GO_top20$subtype = "VSMC"
Endo_genes_GO_top20$subtype = "Endo"
Astro_genes_GO_top20$subtype = "Astro"
Oligo_genes_GO_top20$subtype = "Oligo"
Neuron_genes_GO_top20$subtype = "Neuron"

# 1.17 merge the GO pathways

All_GO_top20 <- rbind(
  Unknow_genes_GO_top20,B_cell_genes_GO_top20,T_NK_cell_genes_GO_top20,
  Neutro_genes_GO_top20,Macro_genes_GO_top20,Micro_genes_GO_top20,
  Mono_genes_GO_top20,CPC_genes_GO_top20,EPC_genes_GO_top20,
  ABC_genes_GO_top20,Fibro_genes_GO_top20,RBC_genes_GO_top20,
  VSMC_genes_GO_top20,Endo_genes_GO_top20,Astro_genes_GO_top20,
  Oligo_genes_GO_top20,Neuron_genes_GO_top20)
write.csv(All_GO_top20,"All_GO_top20.CSV")

All_GO_top$"-log10Pvalue" <- -log10(All_GO_top$pvalue)
All_GO_top$item <- row.names(All_GO_top)

# 1.18 GO pathways visualization

All_GO_top$subtype <- factor(All_GO_top$subtype,levels = rev(c("Unknow","B cell","T/NK_cell","Neutro","Macro","Micro","Mono","CPC","EPC",     
                                         "ABC","Fibro","RBC","VSMC","Endo","Astro","Oligo","Neuron")))

ggbarplot(All_GO_top, x="subtype", y="-log10Pvalue", 
          fill = "subtype", color = "white",
          label = All_GO_top$Description,
          palette =  c("Neuron" = '#E95C59',
                       "Oligo" = '#53A85F',
                       "Astro" = '#F1BB72',
                       "Endo" = '#F3B1A0',
                       "VSMC" = '#D6E7A3',
                       "RBC" = '#57C3F3',
                       "Fibro" = '#E63863',
                       "ABC" = '#E59CC4',
                       "EPC" = '#AB3282',
                       "CPC" = '#23452F',
                       "Mono" = '#BD956A',
                       "Micro" = 'orange',
                       "Macro" = '#58A4C3',
                       "Neutro" = "#00BFC4",
                       "T/NK cell" = 'red',
                       "B cell" = '#8C549C',
                       "Unknow" = "black"),
          order = rev(c("Unknow","B cell","T/NK_cell","Neutro","Macro","Micro","Mono","CPC","EPC",     
                    "ABC","Fibro","RBC","VSMC","Endo","Astro","Oligo","Neuron")),
          x.text.angle=0, 
          xlab = NULL) + coord_flip()
