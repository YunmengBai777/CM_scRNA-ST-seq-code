# 1.Process and integration

# 1.1 Library packages

library(Seurat)
library(ggplot2)
library(RCurl)
library(cowplot)
library(dplyr)
library(hdf5r)
options(future.globals.maxSize = 10000 * 1024^2)  # set allowed size to 2K MiB

# 1.2 Load datasets

BC1_ST_A <- readRDS("brain1.rds")
names(BC1_ST_A@images) = "BC1_ST_A"

BC1_ST_P <- readRDS("brain2.rds")
names(BC1_ST_P@images) = "BC1_ST_P"

BC2_ST_A <- readRDS("brain3.rds")
names(BC2_ST_P@images) = "BC2_ST_P"
BC2_ST_P <- readRDS("brain4.rds")

BA1_ST_A <- Load10X_Spatial(
    data.dir = "BA1_ST_A/outs/",
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "BA1_ST_A",
    filter.matrix = TRUE,
    to.upper = FALSE)

BA1_ST_P <- Load10X_Spatial(
    data.dir = "BA1_ST_P/outs/",
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "BA1_ST_P",
    filter.matrix = TRUE,
    to.upper = FALSE)

BA2_ST_A <- Load10X_Spatial(
    data.dir = "BA2_ST_A/outs/",
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "BA2_ST_A",
    filter.matrix = TRUE,
    to.upper = FALSE)

BA2_ST_P <- Load10X_Spatial(
    data.dir = "BA2_ST_P/outs/",
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "BA2_ST_P",
    filter.matrix = TRUE,
    to.upper = FALSE)

BM1_ST_A <- Load10X_Spatial(
    data.dir = "BM1_ST_A/outs/",
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "BM1_ST_A",
    filter.matrix = TRUE,
    to.upper = FALSE)

BM1_ST_P <- Load10X_Spatial(
    data.dir = "BM1_ST_P/outs/",
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "BM1_ST_P",
    filter.matrix = TRUE,
    to.upper = FALSE)

BM2_ST_A <- Load10X_Spatial(
    data.dir = "BM2_ST_A/outs/",
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "BM2_ST_A",
    filter.matrix = TRUE,
    to.upper = FALSE)

BM2_ST_P <- Load10X_Spatial(
    data.dir = "BM2_ST_P/outs/",
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "BM2_ST_P",
    filter.matrix = TRUE,
    to.upper = FALSE)

# 1.3 Add metadata

BC1_ST_A$type <- "Control"
BC1_ST_P$type <- "Control"
BC2_ST_A$type <- "Control"
BC2_ST_P$type <- "Control"
BM1_ST_A$type <- "Model"
BM1_ST_P$type <- "Model"
BM2_ST_A$type <- "Model"
BM2_ST_P$type <- "Model"
BA1_ST_A$type <- "ART"
BA1_ST_P$type <- "ART"
BA2_ST_A$type <- "ART"
BA2_ST_P$type <- "ART"

BC1_ST_A$sample <- "BC1"
BC1_ST_P$sample <- "BC1"
BC2_ST_A$sample <- "BC2"
BC2_ST_P$sample <- "BC2"
BM1_ST_A$sample <- "BM1"
BM1_ST_P$sample <- "BM1"
BM2_ST_A$sample <- "BM2"
BM2_ST_P$sample <- "BM2"
BA1_ST_A$sample <- "BA1"
BA1_ST_P$sample <- "BA1"
BA2_ST_A$sample <- "BA2"
BA2_ST_P$sample <- "BA2"

BC1_ST_A$data_type <- "ST"
BC1_ST_P$data_type <- "ST"
BC2_ST_A$data_type <- "ST"
BC2_ST_P$data_type <- "ST"
BM1_ST_A$data_type <- "ST"
BM1_ST_P$data_type <- "ST"
BM2_ST_A$data_type <- "ST"
BM2_ST_P$data_type <- "ST"
BA1_ST_A$data_type <- "ST"
BA1_ST_P$data_type <- "ST"
BA2_ST_A$data_type <- "ST"
BA2_ST_P$data_type <- "ST"

BC1_ST_A$region <- "anterior"
BC1_ST_P$region <- "posterior"
BC2_ST_A$region <- "anterior"
BC2_ST_P$region <- "posterior"
BM1_ST_A$region <- "anterior"
BM1_ST_P$region <- "posterior"
BM2_ST_A$region <- "anterior"
BM2_ST_P$region <- "posterior"
BA1_ST_A$region <- "anterior"
BA1_ST_P$region <- "posterior"
BA2_ST_A$region <- "anterior"
BA2_ST_P$region <- "posterior"

BC1_ST_A$slice <- "BC1_A"
BC1_ST_P$slice <- "BC1_P"
BC2_ST_A$slice <- "BC2_A"
BC2_ST_P$slice <- "BC2_P"
BM1_ST_A$slice <- "BM1_A"
BM1_ST_P$slice <- "BM1_P"
BM2_ST_A$slice <- "BM2_A"
BM2_ST_P$slice <- "BM2_P"
BA1_ST_A$slice <- "BA1_A"
BA1_ST_P$slice <- "BA1_P"
BA2_ST_A$slice <- "BA2_A"
BA2_ST_P$slice <- "BA2_P"

# 1.4 Evaluate MT percent

BC1_ST_A <- PercentageFeatureSet(BC1_ST_A, "^mt-", col.name = "percent_mito")
BC1_ST_P <- PercentageFeatureSet(BC1_ST_P, "^mt-", col.name = "percent_mito")
BC2_ST_A <- PercentageFeatureSet(BC2_ST_A, "^mt-", col.name = "percent_mito")
BC2_ST_P <- PercentageFeatureSet(BC2_ST_P, "^mt-", col.name = "percent_mito")
BM1_ST_A <- PercentageFeatureSet(BM1_ST_A, "^mt-", col.name = "percent_mito")
BM1_ST_P <- PercentageFeatureSet(BM1_ST_P, "^mt-", col.name = "percent_mito")
BM2_ST_A <- PercentageFeatureSet(BM2_ST_A, "^mt-", col.name = "percent_mito")
BM2_ST_P <- PercentageFeatureSet(BM2_ST_P, "^mt-", col.name = "percent_mito")
BA1_ST_A <- PercentageFeatureSet(BA1_ST_A, "^mt-", col.name = "percent_mito")
BA1_ST_P <- PercentageFeatureSet(BA1_ST_P, "^mt-", col.name = "percent_mito")
BA2_ST_A <- PercentageFeatureSet(BA2_ST_A, "^mt-", col.name = "percent_mito")
BA2_ST_P <- PercentageFeatureSet(BA2_ST_P, "^mt-", col.name = "percent_mito")

# 1.5 Merge datasets

brain.merge <- merge(BA1_ST_A,c(BA1_ST_P,BA2_ST_A,BA2_ST_P,BC1_ST_A,BC1_ST_P,
BC2_ST_A,BC2_ST_P,BM1_ST_A,BM1_ST_P,BM2_ST_A,BM2_ST_P))
brain.integrated$type  <- factor(brain.integrated$type, levels = c("Control","Model","ART"))
brain.merge$slice = factor(brain.merge$slice, levels = c("BC1_A", "BC1_P","BC2_A", "BC2_P","BM1_A","BM1_P",
"BM2_A","BM2_P","BA1_A","BA1_P","BA2_A","BA2_P"))

# 1.6 Datasets QC

pdf("brain_CM_ST_before_filter.pdf",height = 5,width = 12)
VlnPlot(brain.merge, features = c("nFeature_Spatial","nCount_Spatial","percent_mito"), pt.size = 0,group.by = "slice")
dev.off()

brain.merge.filter <- subset(brain.merge, subset = nFeature_Spatial > 200 & nFeature_Spatial < 8000 & nCount_Spatial > 1000 & nCount_Spatial < 60000 & percent_mito < 30)
brain.merge.filter

pdf("brain_CM_ST_afer_filter.pdf",height = 5,width = 12)
VlnPlot(brain.merge.filter, features = c("nFeature_Spatial","nCount_Spatial","percent_mito"), pt.size = 0,group.by = "slice")
dev.off()

# 1.7 Integrated datasets 

st.list = list(BA1_ST_A = BA1_ST_A,     
BA1_ST_P = BA1_ST_P,         
BA2_ST_A = BA2_ST_A,
BA2_ST_P = BA2_ST_P,         
BC1_ST_A = BC1_ST_A,
BC1_ST_P = BC1_ST_P,
BC2_ST_A = BC2_ST_A,         
BC2_ST_P = BC2_ST_P,
BM1_ST_A = BM1_ST_A,
BM1_ST_P = BM1_ST_P,
BM2_ST_A = BM2_ST_A,
BM2_ST_P = BM2_ST_P)
st.list = lapply(st.list, SCTransform, assay = "Spatial", method = "poisson")

st.features = SelectIntegrationFeatures(st.list, nfeatures = 3000, verbose = FALSE)
st.list <- PrepSCTIntegration(object.list = st.list, anchor.features = st.features, verbose = FALSE)
int.anchors <- FindIntegrationAnchors(object.list = st.list, normalization.method = "SCT", verbose = FALSE, anchor.features = st.features)
brain.integrated <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT", verbose = FALSE)

# 1.8 Dimension reduction 

brain.integrated <- RunPCA(object = brain.integrated,verbose = FALSE)

pdf("ElbowPlot.pdf",height = 8,width = 8)
ElbowPlot(brain.integrated,ndims=50)
dev.off()

brain.integrated <- FindNeighbors(brain.integrated,dim=1:40)
brain.integrated <- FindClusters(brain.integrated,resolution = 0.8)
brain.integrated <- RunUMAP(brain.integrated,reduction="pca", dims = 1:40)
brain.integrated <- RunTSNE(brain.integrated,reduction="pca", dims = 1:40)

# 1.9 Umap and Tsne visulization

pdf("brain_CM_ST_umap_cluster.pdf",height = 8 ,width = 8)
DimPlot(brain.integrated,reduction = "umap",group.by="integrated_snn_res.0.8",label = F,pt.size = 1) + ggtitle("Cluster")
dev.off()

pdf("brain_CM_ST_umap_region.pdf",height = 8,width = 8)
DimPlot(brain.integrated,reduction = "umap" ,group.by="region",label = F,pt.size = 1) + ggtitle("Region")
dev.off()

pdf("brain_CM_ST_umap_sample.pdf",height = 8,width = 8)
DimPlot(brain.integrated,reduction = "umap" ,group.by="slice",label = F,pt.size = 1) + ggtitle("Sample")
dev.off()

pdf("brain_CM_ST_umap_type.pdf",height = 8,width = 8)
DimPlot(brain.integrated,reduction = "umap" ,group.by="type",label = F,pt.size = 1,cols = c("#3B4992","#EE0000","#008B45")) + ggtitle("Type")
dev.off()

pdf("brain_CM_ST_tsne_cluster.pdf",height = 8 ,width = 8)
DimPlot(brain.integrated,reduction = "tsne",group.by="integrated_snn_res.0.8",label = F,pt.size = 1) + ggtitle("Cluster")
dev.off()

pdf("brain_CM_ST_tsne_region.pdf",height = 8,width = 8)
DimPlot(brain.integrated,reduction = "tsne" ,group.by="region",label = F,pt.size = 1) + ggtitle("Region")
dev.off()

pdf("brain_CM_ST_tsne_sample.pdf",height = 8,width = 8)
DimPlot(brain.integrated,reduction = "tsne" ,group.by="slice",label = F,pt.size = 1) + ggtitle("Sample")
dev.off()

pdf("brain_CM_ST_tsne_type.pdf",height = 8,width = 8)
DimPlot(brain.integrated,reduction = "tsne" ,group.by="type",label = F,pt.size = 1,cols = c("#3B4992","#EE0000","#008B45")) + ggtitle("Type")
dev.off()

# 1.10 Save datasets

saveRDS(brain.integrated,"CM_ST_dataset.rds")