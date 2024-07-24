# 4.Seurat FindTransferAnchors

# 4.1 Loading packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
options(future.globals.maxSize = 40000 * 1024^2)

# 4.2 Loading ST dataset
scRNA <- readRDS("CM_ST_dataset.rds")

# DimPlot(scRNA)
table(scRNA$cell_type)

ST = readRDS("CM_scRNA-seq_dataset.rds")
DimPlot(ST,reduction = "tsne")

scRNA_reference <- SCTransform(scRNA, ncells = 3000, verbose = FALSE)
scRNA_reference <- RunPCA(scRNA_reference,verbose = FALSE)
scRNA_reference <- RunUMAP(scRNA_reference,dims = 1:30)

# 4.3 FindTransferAnchors
anchors <- FindTransferAnchors(reference = scRNA_reference, query = ST, normalization.method = "SCT",dims = 1:30, reference.reduction = "pca")
predictions.assay <- TransferData(anchorset = anchors, refdata = scRNA$cell_type, 
                                  prediction.assay = TRUE)#, weight.reduction = scRNA_reference[["pca"]], dims = 1:30)
ST[["predictions"]] <- predictions.assay
DefaultAssay(ST) <- "predictions"

# 4.4 FeaturePlot visualization

FeaturePlot(ST,c("Neuron", "Oligo","Astro","Endo",
                 "VSMC","RBC","Fibro","ABC"),reduction = "tsne",pt.size = 0.5,ncol = 4,cols = c("grey","red"))

FeaturePlot(ST,c("EPC","CPC","Mono","Micro",
                 "Macro","Neutro","T/NK cell","B cell"),reduction = "tsne",pt.size = 0.5,ncol = 4,cols = c("grey","red"))


# 4.5 SpatialFeaturePlot visualization

# (1)BC1_ST_A
pdf("/data/rawdata/chenjiayun/09_Chen_CM_scRNA_ST/2.ST_data/4.TransferAnchors/BC1_ST_A.pdf",height = 50,width = 10)
SpatialFeaturePlot(ST, features = c("Neuron", "Oligo","Astro","Endo",
                                    "VSMC","RBC","Fibro","ABC",
                                    "EPC","CPC","Mono","Micro",
                                    "Macro","Neutro","T/NK cell","B cell"), 
                   pt.size.factor = 1.6, ncol = 2, crop = TRUE,images = c("anterior1"))
dev.off()

# (2)BC1_ST_P
pdf("/data/rawdata/chenjiayun/09_Chen_CM_scRNA_ST/2.ST_data/4.TransferAnchors/BC1_ST_P.pdf",height = 50,width = 10)
SpatialFeaturePlot(ST, features = c("Neuron", "Oligo","Astro","Endo",
                                    "VSMC","RBC","Fibro","ABC",
                                    "EPC","CPC","Mono","Micro",
                                    "Macro","Neutro","T/NK cell","B cell"), 
                   pt.size.factor = 1.6, ncol = 2, crop = TRUE,images = c("posterior1"))
dev.off()

# (3)BC2_ST_A
pdf("/data/rawdata/chenjiayun/09_Chen_CM_scRNA_ST/2.ST_data/4.TransferAnchors/BC2_ST_A.pdf",height = 50,width = 10)
SpatialFeaturePlot(ST, features = c("Neuron", "Oligo","Astro","Endo",
                                    "VSMC","RBC","Fibro","ABC",
                                    "EPC","CPC","Mono","Micro",
                                    "Macro","Neutro","T/NK cell","B cell"), 
                   pt.size.factor = 1.6, ncol = 2, crop = TRUE,images = c("anterior2"))
dev.off()

# (4)BC2_ST_P
pdf("/data/rawdata/chenjiayun/09_Chen_CM_scRNA_ST/2.ST_data/4.TransferAnchors/BC2_ST_P.pdf",height = 50,width = 10)
SpatialFeaturePlot(ST, features = c("Neuron", "Oligo","Astro","Endo",
                                    "VSMC","RBC","Fibro","ABC",
                                    "EPC","CPC","Mono","Micro",
                                    "Macro","Neutro","T/NK cell","B cell"), 
                   pt.size.factor = 1.6, ncol = 2, crop = TRUE,images = c("posterior2"))
dev.off()

# (5)BM1_ST_A
pdf("/data/rawdata/chenjiayun/09_Chen_CM_scRNA_ST/2.ST_data/4.TransferAnchors/BM1_ST_A.pdf",height = 50,width = 10)
SpatialFeaturePlot(ST, features = c("Neuron", "Oligo","Astro","Endo",
                                    "VSMC","RBC","Fibro","ABC",
                                    "EPC","CPC","Mono","Micro",
                                    "Macro","Neutro","T/NK cell","B cell"), 
                   pt.size.factor = 1.6, ncol = 2, crop = TRUE,images = c("BM1_ST_A"))
dev.off()

# (6)BM1_ST_P
pdf("/data/rawdata/chenjiayun/09_Chen_CM_scRNA_ST/2.ST_data/4.TransferAnchors/BM1_ST_P.pdf",height = 50,width = 10)
SpatialFeaturePlot(ST, features = c("Neuron", "Oligo","Astro","Endo",
                                    "VSMC","RBC","Fibro","ABC",
                                    "EPC","CPC","Mono","Micro",
                                    "Macro","Neutro","T/NK cell","B cell"), 
                   pt.size.factor = 1.6, ncol = 2, crop = TRUE,images = c("BM1_ST_P"))
dev.off()

# (7)BM2_ST_A
pdf("/data/rawdata/chenjiayun/09_Chen_CM_scRNA_ST/2.ST_data/4.TransferAnchors/BM2_ST_A.pdf",height = 50,width = 10)
SpatialFeaturePlot(ST, features = c("Neuron", "Oligo","Astro","Endo",
                                    "VSMC","RBC","Fibro","ABC",
                                    "EPC","CPC","Mono","Micro",
                                    "Macro","Neutro","T/NK cell","B cell"), 
                   pt.size.factor = 1.6, ncol = 2, crop = TRUE,images = c("BM2_ST_A"))
dev.off()

# (8)BM2_ST_P
pdf("/data/rawdata/chenjiayun/09_Chen_CM_scRNA_ST/2.ST_data/4.TransferAnchors/BM2_ST_P.pdf",height = 50,width = 10)
SpatialFeaturePlot(ST, features = c("Neuron", "Oligo","Astro","Endo",
                                    "VSMC","RBC","Fibro","ABC",
                                    "EPC","CPC","Mono","Micro",
                                    "Macro","Neutro","T/NK cell","B cell"), 
                   pt.size.factor = 1.6, ncol = 2, crop = TRUE,images = c("BM2_ST_P"))
dev.off()

# (9)BA1_ST_A
pdf("/data/rawdata/chenjiayun/09_Chen_CM_scRNA_ST/2.ST_data/4.TransferAnchors/BA1_ST_A.pdf",height = 50,width = 10)
SpatialFeaturePlot(ST, features = c("Neuron", "Oligo","Astro","Endo",
                                    "VSMC","RBA","Fibro","ABA",
                                    "EPC","CPC","Mono","Micro",
                                    "Macro","Neutro","T/NK cell","B cell"), pt.size.factor = 1.6, ncol = 2, crop = TRUE,images = c("BA1_ST_A"))
dev.off()

# (10)BA1_ST_P
pdf("/data/rawdata/chenjiayun/09_Chen_CM_scRNA_ST/2.ST_data/4.TransferAnchors/BA1_ST_P.pdf",height = 50,width = 10)
SpatialFeaturePlot(ST, features = c("Neuron", "Oligo","Astro","Endo",
                                    "VSMC","RBA","Fibro","ABA",
                                    "EPC","CPC","Mono","Micro",
                                    "Macro","Neutro","T/NK cell","B cell"), 
                   pt.size.factor = 1.6, ncol = 2, crop = TRUE,images = c("BA1_ST_P"))
dev.off()

# (11)BA1_ST_A
pdf("/data/rawdata/chenjiayun/09_Chen_CM_scRNA_ST/2.ST_data/4.TransferAnchors/BA2_ST_A.pdf",height = 50,width = 10)
SpatialFeaturePlot(ST, features = c("Neuron", "Oligo","Astro","Endo",
                                    "VSMC","RBA","Fibro","ABA",
                                    "EPC","CPC","Mono","Micro",
                                    "Macro","Neutro","T/NK cell","B cell"), 
                   pt.size.factor = 1.6, ncol = 2, crop = TRUE,images = c("BA2_ST_A"))
dev.off()

# (12)BA1_ST_P
pdf("/data/rawdata/chenjiayun/09_Chen_CM_scRNA_ST/2.ST_data/4.TransferAnchors/BA2_ST_P.pdf",height = 50,width = 10)
SpatialFeaturePlot(ST, features = c("Neuron", "Oligo","Astro","Endo",
                                    "VSMC","RBA","Fibro","ABA",
                                    "EPC","CPC","Mono","Micro",
                                    "Macro","Neutro","T/NK cell","B cell"), 
                   pt.size.factor = 1.6, ncol = 2, crop = TRUE,images = c("BA2_ST_P"))
dev.off()


# 4.6 For specific cell type

pdf("/data/rawdata/chenjiayun/09_Chen_CM_scRNA_ST/2.ST_data/4.TransferAnchors/ST_T.pdf",height = 20,width = 30)
SpatialFeaturePlot(ST, features = c("T/NK cell"), pt.size.factor = 1.6, ncol = 4, crop = TRUE,images = c("anterior1","posterior1","anterior2","posterior2"))/
  SpatialFeaturePlot(ST, features = c("T/NK cell"), pt.size.factor = 1.6, ncol = 4, crop = TRUE,images = c("BM1_ST_A","BM1_ST_P","BM2_ST_A","BM2_ST_P"))/
  SpatialFeaturePlot(ST, features = c("T/NK cell"), pt.size.factor = 1.6, ncol = 4, crop = TRUE,images = c("BA1_ST_A","BA1_ST_P","BA2_ST_A","BA2_ST_P"))
dev.off()


# 4.5 Dotplot visualization

AverageExpression_value <-  AverageExpression(ST,assays = "predictions", 
                                              features = c("Neuron","Oligo","Astro","Endo","VSMC","RBC","Fibro","ABC","EPC",
                                                           "CPC","Mono","Micro","Macro","Neutro","T/NK cell","B cell"),
                                              return.seurat = FALSE, group.by = c("SingleR.labels"))


DotPlot(ST,group.by = "SingleR.labels",features = rev(c("Neuron","Oligo","Astro","Endo","VSMC","RBC","Fibro","ABC","EPC",
                                                        "CPC","Mono","Micro","Macro","Neutro","T/NK cell","B cell"))) +
  coord_flip() + xlab(NULL) + ylab(NULL) + scale_color_gradientn(colours = viridis(20), guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"), name = "Module socre") +
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1) +
  annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) # 5*8



