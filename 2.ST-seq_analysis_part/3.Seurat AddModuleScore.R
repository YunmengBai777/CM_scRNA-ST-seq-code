# 3 Seurat_AddModuleScore 

# 3.1 Loading packages
library(Seurat)
library(ggplot2)
library(RCurl)
library(cowplot)
library(dplyr)
library(ggsci)
library(viridis)

# 3.2 Loading ST dataset
ST <- readRDS("CM_ST_dataset.rds")
ST

# 3.3 Making cell type maker geneset as module from scRNA-seq dataset

scRNA_marker <- read.csv("scRNA_markers.csv")
table(scRNA_marker$cluster)

Neuron <- list(scRNA_marker[scRNA_marker$cluster %in% "Neuron",]$gene)
Oligo <- list(scRNA_marker[scRNA_marker$cluster %in% "Oligo",]$gene)
Astro <- list(scRNA_marker[scRNA_marker$cluster %in% "Astro",]$gene)
Endo <- list(scRNA_marker[scRNA_marker$cluster %in% "Endo",]$gene)
Mural <- list(scRNA_marker[scRNA_marker$cluster %in% "VSMC",]$gene)
RBC <- list(scRNA_marker[scRNA_marker$cluster %in% "RBC",]$gene)
Fibro <- list(scRNA_marker[scRNA_marker$cluster %in% "Fibro",]$gene)
ABC <- list(scRNA_marker[scRNA_marker$cluster %in% "ABC",]$gene)
EPC <- list(scRNA_marker[scRNA_marker$cluster %in% "EPC",]$gene)
CPC <- list(scRNA_marker[scRNA_marker$cluster %in% "CPC",]$gene)
Mono <- list(scRNA_marker[scRNA_marker$cluster %in% "Mono",]$gene)
Micro <- list(scRNA_marker[scRNA_marker$cluster %in% "Micro",]$gene)
Macro <- list(scRNA_marker[scRNA_marker$cluster %in% "Macro",]$gene)
Neutro <- list(scRNA_marker[scRNA_marker$cluster %in% "Neutro",]$gene)
T_NK_cell <- list(scRNA_marker[scRNA_marker$cluster %in% "T/NK cell",]$gene)
B_cell <- list(scRNA_marker[scRNA_marker$cluster %in% "B cell",]$gene)

# 3.4 Add modulesocre to ST dataset

ST <- AddModuleScore(object = ST,features = c(Neuron,Oligo,Astro,Endo,Mural,RBC,Fibro,ABC,EPC,CPC,
                                              Mono,Micro,Macro,Neutro,T_NK_cell,B_cell),
                      assay="SCT",name=c("Neuron","Oligo","Astro","Endo","Mural","RBC","Fibro","ABC","EPC","CPC",
                                         "Mono","Micro","Macro","Neutro","T_NK_cell","B_cell"))
colnames(head(ST@meta.data))

meta <- ST@meta.data[,c(14,16:31)]

AverageExpression_value <-  AverageExpression(ST,# assays = "SCT", 
                                              features = c("Neuron1","Oligo2","Astro3","Endo4","Mural5","RBC6","Fibro7","ABC8","EPC9",
                                                           "CPC10","Mono11","Micro12","Macro13","Neutro14","T_NK_cell15","B_cell16"),
                                              return.seurat = FALSE, group.by = c("SingleR.labels"))

DotPlot(ST,group.by = "SingleR.labels",features = rev(c("Neuron1","Oligo2","Astro3","EPC9",
                                                          "CPC10","Endo4","Mural5","RBC6","Fibro7","ABC8","Mono11","Micro12","Macro13","Neutro14","T_NK_cell15","B_cell16"))) +
  coord_flip() + xlab(NULL) + ylab(NULL) + scale_color_gradientn(colours = viridis(20), guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"), name = "Module socre") +
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1) +
  annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) # 5*8



# 3.5 FeaturePlot visualization

FeaturePlot(ST,features = c("Neuron1","Oligo2","Astro3","CPC10",
                            "Endo4","Fibro7","Micro12","T_NK_cell15"),
            cols=c("white","red"),order=T,pt.size=0.01,reduction="tsne",ncol = 4) # 18*8

FeaturePlot(ST,features = c("Mural5","RBC6","ABC8","EPC9",
                            "Mono11","Macro13","Neutro14","B_cell16"),
            cols=c("white","red"),order=T,pt.size=0.01,reduction="tsne",ncol = 4) # 18*8
