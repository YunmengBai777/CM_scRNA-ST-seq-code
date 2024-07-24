# 2.SingleR annotation

# 2.1 Loading packages
library(Seurat)
library(ggplot2)
library(RCurl)
library(cowplot)
library(dplyr)
library(SingleR)
library(ggsci)
library(viridis)

# 2.2 Loading dataset

ST <- readRDS("CM_ST_dataset.rds")
ST
DimPlot(ST,label = T,label.size = 5)

# 2.3 SingleR annotation

ref.se <- MouseRNAseqData()
combin.data <- ST
expdata <- combin.data[["Spatial"]]@data
anno.cell.main <- SingleR(test=expdata,ref = ref.se,labels = ref.se$label.main)  
table(anno.cell.main$labels)
# Adipocytes      Astrocytes Endothelial cells  Epithelial cells       Fibroblasts 
#     2             141                99               127                54 
# Macrophages    Microglia         Monocytes         Neurons  Oligodendrocytes 
#    3               1                19             26762             10461

# 2.4 Dimplot visualization

DimPlot(combin.data,label = F,group.by = "SingleR.labels",
        label.size = 5,reduction = "tsne",pt.size = 1)

DimPlot(combin.data,label = F,group.by = "SingleR.labels",
        label.size = 5,reduction = "tsne",
        split.by = "type",pt.size = 1)

# 2.5 Change spot labels

table(combin.data$SingleR.labels)
type <- as.character(combin.data$SingleR.labels)
type[type %in% c("Monocytes","Macrophages","Microglia","Adipocytes")] = "Immune"
type[type %in% c("Astrocytes")] = "Astro"
type[type %in% c("Endothelial cells")] = "Endo"
type[type %in% c("Epithelial cells")] = "CPC"
type[type %in% c("Fibroblasts")] = "Fibro"
type[type %in% c("Oligodendrocytes")] = "Oligo"
type[type %in% c("Neurons")] = "Neuron"

table(type)
# Astro    CPC   Endo  Fibro Immune Neuron  Oligo 
#  141    127     99     54     25  26762  10461

combin.data$SingleR.labels <- factor(type,levels = c("Neuron","Oligo","Astro","CPC","Endo","Fibro","Immune"))

table(combin.data$SingleR.labels)
# Neuron  Oligo  Astro    CPC   Endo  Fibro Immune 
# 26762  10461    141    127     99     54     25

combin.data@active.ident <- combin.data$SingleR.labels
combin.data$type <- factor(combin.data$type,levels = c("Control","Model","ART"))

ST <- combin.data

table(ST$type,ST$SingleR.labels)
table(ST$slice,ST$SingleR.labels)

# 2.6 Save dataset and metadata

write.csv(ST@meta.data,"CM_ST_metadata.csv")
saveRDS(ST,"1.Data/brain.integrated.anno.rds")

