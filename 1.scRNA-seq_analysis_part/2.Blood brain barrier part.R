# 2.BBB part

# 2.1 library packages

library(Seurat)
library(ggplot2)
library(ggsci)
library(viridis)
library(dplyr)
options(stringsAsFactors = F)

# 2.2 load scRNA datasets

scRNA <- readRDS("CM_scRNA_dataset.rds")
DimPlot(scRNA)
data <- subset(scRNA,subset = cell_type %in% c("Endo","VSMC","Astro"))
DimPlot(data)
data

# 2.3 reduction

data <- SCTransform(data,return.only.var.genes = FALSE,assay = "RNA",verbose = FALSE)
data <- RunPCA(object = data,verbose = FALSE)
ElbowPlot(data)
data <- FindNeighbors(data,dim=1:20)
data <- FindClusters(data,resolution = 0.8)
data <- RunUMAP(data,reduction="pca", dims = 1:20)

DimPlot(data,label = T)
DimPlot(data,group.by = "cell_type")

DoHeatmap(data,c("Cldn5","Kcnj8","Aldoc","Acta2"),group.by = "cell_type")

# 2.3.1 Endo

Endo <- subset(data,SCT_snn_res.0.8 %in% c(2:7,9:12,14:15,20,21,24))
DimPlot(Endo,group.by = "SCT_snn_res.0.8",split.by = "sample_type",label = T)

DoHeatmap(Endo,c("Cldn5","Kdr","Tagln","Myh11","Acta2","Crip1","Myl9"))

"Bmx","Efnb2","Vegfc","Sema3g","Gkn3"  # Endo_Art,
"Mfsd2a","Tfrc","Slc7a5","Slc16a1" # Endo_Cap
"Nr2f2","Vwf","Slc38a5","Ackr1","Lrg1","Vwf","Ch25h","Vcam1","Icam1" # Endo_Ven,
"Ccl3","Ccl4","Ccl5" # # Endo_Inflam
"Hba-a1","Hba-a2" # Endo_Ap,
"Acta2","Myl9","Myh11","Tagln" # aSMC, 
"Acta2","Abcc9","Pdgfrb","Vtn", # vSMC,
"Kcnj8" # Peri

# 2.3.2 Astrocyte

Astro <- subset(data,SCT_snn_res.0.8 %in% c(0,1,8,17,18,22,23))
DimPlot(Astro,group.by = "SCT_snn_res.0.8",split.by = "sample_type",label = T)

"Aldoc","Aqp4","Aldh1l1","Slc1a2","S100b","Gfap","Apoe","Clu"
"Lcn2","Ifitm3","Timp1","B2M","Vim","Gfap","Psmb8","Bst2"
"Apoe","Vegfa","Mmp9","Edn1","Ccl2","Cxcl10"
"Mfge8","Igfbp2","Nnat","Dbi","Agt","Nupr1","Dbi","Prss56","Slc43a3",
"Timp1","Thbs4","Igfbp5","Fxyd6","Cd9","Itih3","Sparc","Nkx6-2"

Astro_DEG <- FindAllMarkers(Astro,only.pos = T)

Astro_DEG_top10 <- Astro_DEG %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC)
Astro_DEG_top10$gene
DoHeatmap(Astro,Astro_DEG_top10$gene)

AverageExpression_value <-  AverageExpression(Endo, assays = "SCT", 
                                              features = c("Cldn5","Flt1","Slc2a1","Pecam1","Vegfc","Arl15","Igfbp3", # Arteriole
                                                           "Atp10a", "Syne1","Abcb1a","Npipb6a","Cmtm8","Angpt2","Itm2a","Inpp5d", # Capillary
                                                           "Tshz2","Ackr1","Anxa2","Pdgfrb","Adgrg6","Aff3" # Venule
                                              ),return.seurat = FALSE, group.by = "ident")

pheatmap::pheatmap(AverageExpression_value$SCT,
                   cluster_rows = F,cluster_cols = T,
                   border_color = "black",
                   angle_col = 45,
                   # gaps_col = c(1,5,6,7),
                   gaps_row = c(7,14),
                   cellwidth = 20,cellheight = 10,
                   scale = "row",viridis(100))

data2 <- subset(data,subset = SCT_snn_res.0.8 ==  "13",invert =T)
DimPlot(data2)

# 2.4 Rename the cluster

new.cluster.ids <- c("Astro_C1","Astro_C2","Endo_Cap","Endo_Cap","Endo_Cap","Endo_Inflam","Endo_Inflam","Endo_Cap","Astro_C2","Endo_Art","Endo_Inflam","Endo_Art","vSMC","Endo_Ven","Endo_Art","Pericyte","Astro_C3","Astro_C4","aSMC","vSMC","Endo_Ven","Astro_C2","Astro_C5","Endo_Ap")

names(new.cluster.ids) <- levels(data2)
data2 <- RenameIdents(data2, new.cluster.ids)

table(data2@active.ident)
# Astro_C1    Astro_C2    Endo_Cap Endo_Inflam    Endo_Art        aSMC    Endo_Ven    Pericyte    Astro_C3 
# 3527        4516        8177        4229        2663        1335        1041         597         483 
# Astro_C4        vSMC    Astro_C5     Endo_Ap 
# 478         451         217         212 


levels(data2) <- c("Astro_C1","Astro_C2","Astro_C3","Astro_C4","Astro_C5",
                  "Endo_Art","Endo_Cap","Endo_Ven","Endo_Inflam","Endo_Ap",
                  "aSMC","vSMC","Pericyte")

DimPlot(data2,pt.size = 1,reduction = "umap",label = T,label.size = 5) #+ scale_color_lancet()

saveRDS(data2,"7.BBB/BBB_anno.rds")


theme <- theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               # panel.border = element_blank(),
               axis.title = element_blank(),  
               axis.text = element_blank(), 
               axis.ticks = element_blank(),
               panel.background = element_rect(fill = 'white'),
               plot.background=element_rect(fill="white"))

cell_type_cols <- c("#FF34B3","#BC8F8F","#20B2AA","#FF3030","#FFA500",
                    "#FFC1C1","#FF6A6A","#CD5C5C", "#AB82FF","#90EE90",
                    "#00CD00","#008B8B","#6495ED") 

DimPlot(data2,pt.size = 1,reduction = "umap",label = T,label.size = 5,cols = cell_type_cols) + theme +  
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1) +
  annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)

# heatmap
AverageExpression_value <-  AverageExpression(data2, assays = "SCT", 
                                              features = c("Aldoc","Aldh1l1","Slc1a2","S100b","Aqp4",
                                                           "Agt","Itih3","Nnat","Sparc","Igfbp5",  # C0
                                                           "Mfge8","Ptn", # C1,C8
                                                           "Cx3cr1","Csf1r", # C17
                                                           "Gria1","Hopx",  # C18
                                                           "Gfap","Meg3",
                                                           "Cldn5","Kdr","Pecam1",
                                                           "Bmx","Efnb2","Vegfc","Sema3g","Gkn3",
                                                           "Mfsd2a","Tfrc","Slc7a5","Slc16a1",
                                                           "Vwf","Slc38a5","Ackr1","Lrg1","Vwf","Ch25h","Vcam1","Icam1",
                                                           "Ccl3","Ccl4","Ccl5",
                                                           "Hba-a1","Hba-a2",
                                                           "Acta2","Myl9","Myh11","Tagln", # aSMC, C12,20
                                                           "Acta2","Abcc9","Pdgfrb","Vtn", # vSMC, C19
                                                           "Kcnj8" # Peri, C16
                                              ),return.seurat = FALSE, group.by = "ident")

pheatmap::pheatmap(AverageExpression_value$SCT,
                   cluster_rows = F,cluster_cols = F,
                   border_color = "black",
                   angle_col = 45,
                   gaps_col = c(5,10),
                   gaps_row = c(18,42),
                   cellwidth = 20,cellheight = 10,
                   scale = "row",viridis(100))

data2$cell_type_sub <- data2@active.ident

Astro <- subset(data2,subset = cell_type_sub %in%  c("Astro_C1","Astro_C2","Astro_C3","Astro_C4","Astro_C5"))

VlnPlot(Astro,c("Apoe","Vegfa","Edn1","C1qc"),slot = "scale.data",
        group.by = "cell_type_sub",pt.size = 0,sort = F,cols = c("Astro_C1"="#FF34B3",
                                                                 "Astro_C2"="#BC8F8F",
                                                                 "Astro_C3"="#20B2AA",
                                                                 "Astro_C4"="#FF3030",
                                                                 "Astro_C5"="#FFA500"),ncol = 2) # 6*6 
table(Astro$sample_type,Astro$cell_type_sub)

Endo <- subset(data2,subset = cell_type_sub %in%  c("Endo_Art","Endo_Cap","Endo_Ven","Endo_Inflam","Endo_Ap"))

VlnPlot(Endo,c("Edn1","Vcam1","Icam1","Hif1a","Slc2a1","Ldha","Cxcl10","Ccl4","Ccl5"
               # "Procr","Vegfa"
               ),slot = "scale.data",
        group.by = "cell_type_sub",pt.size = 0,sort = F,cols = c("#3B4992","#EE0000","#008B45"),
        split.by = "sample_type",stack = T,flip = T) + xlab("")  # 6*6 

cols = c("Endo_Art"="#FFC1C1",
         "Endo_Cap"="#FF6A6A",
         "Endo_Ven"="#CD5C5C",
         "Endo_Inflam"="#AB82FF",
         "Endo_Ap"="#90EE90")

table(data2$sample_type)

# 2.5 cell proportion

dfsam <- as.data.frame(table(data2$sample_type,data2@active.ident))

control <- dfsam[dfsam$Var1 == "Control",]
myLabel = as.vector(control$Var2)
myLabel = paste(myLabel, "(", round(control$Freq / sum(control$Freq) * 100, 2), "%)"
                , sep = "")  

p1 <- ggplot(control, aes(x = "", y = Freq,fill = Var2),colour = cols ) + 
  geom_bar(stat = "identity")+
  coord_polar(theta = "y") +
  theme_bw() + 
  labs(x = "", y = "")+
  theme(axis.ticks = element_blank())+
  theme(axis.text.x = element_blank()) + 
  theme(panel.grid=element_blank()) +    
  theme(panel.border=element_blank()) + 
  theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + ggtitle("Control(n=8081)") +
  scale_fill_manual(values = cell_type_cols,name = "",labels = myLabel)

p1

model <- dfsam[dfsam$Var1 == "Model",]
myLabel = as.vector(model$Var2)
myLabel = paste(myLabel, "(", round(model$Freq / sum(model$Freq) * 100, 2), "%)"
                , sep = "") 

p2 <- ggplot(model, aes(x = "", y = Freq,fill = Var2),colour = cols ) + 
  geom_bar(stat = "identity")+
  coord_polar(theta = "y") +
  theme_bw() + 
  labs(x = "", y = "")+
  theme(axis.ticks = element_blank())+
  theme(axis.text.x = element_blank()) + 
  theme(panel.grid=element_blank()) +    
  theme(panel.border=element_blank()) + 
  theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + ggtitle("Model(n=11742)") +
  scale_fill_manual(values = cell_type_cols,name = "",labels = myLabel)
p2

treatment <- dfsam[dfsam$Var1 == "ART",]
myLabel = as.vector(treatment$Var2)
myLabel = paste(myLabel, "(", round(treatment$Freq / sum(treatment$Freq) * 100, 2), "%)"
                , sep = "") 

p3 <- ggplot(treatment, aes(x = "", y = Freq,fill = Var2),colour = cols ) + 
  geom_bar(stat = "identity")+
  coord_polar(theta = "y") +
  theme_bw() + 
  labs(x = "", y = "")+
  theme(axis.ticks = element_blank())+
  theme(axis.text.x = element_blank()) + 
  theme(panel.grid=element_blank()) +    
  theme(panel.border=element_blank()) + 
  theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + ggtitle("ART(n=8103)") +
  scale_fill_manual(values = cell_type_cols,name = "",labels = myLabel)
p3

p1 + p2 + p3 #15*5


# 2.6 Endo markers

Endo <- subset(BBB,subset = subtype %in% c("Endo_Art","Endo_Cap","Endo_Ven","Endo_Inflam","Endo_Ap"))

table(Endo$subtype)

VlnPlot(Endo,c("Cldn5","Edn1","Vcam1","Icam1",
               "Cxcl10","Ccl3","Ccl4","Ccl5",
               "Hba-a1","Hba-a2"),
        group.by = "subtype",split.by = "sample_type",
        pt.size = 0,stack = T,
        cols = c("#3B4992","#EE0000","#008B45"),flip = T) + xlab("")


library(clusterProfiler)
library(org.Mm.eg.db)

# 2.7.DEGs and GO enrichment

Endo@active.ident <- Endo$sample_type

Endo_MvC_DEG <- FindMarkers(Endo,ident.1 = "Model",ident.2 = "Control",
                              min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
Endo_MvC_DEG$change = ifelse(Endo_MvC_DEG$p_val_adj < 0.05 & abs(Endo_MvC_DEG$avg_log2FC) >= 0.25, 
                               ifelse(Endo_MvC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Endo_MvC_DEG$change)
# Down Stable     Up 
# 510      3    364 

Endo_MvC_DEG$type <- "Endo_MvC"

Endo_AvC_DEG <- FindMarkers(Endo,ident.1 = "ART",ident.2 = "Control",
                              min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
Endo_AvC_DEG$change = ifelse(Endo_AvC_DEG$p_val_adj < 0.05 & abs(Endo_AvC_DEG$avg_log2FC) >= 0.25, 
                               ifelse(Endo_AvC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Endo_AvC_DEG$change)

Endo_AvC_DEG$type <- "Endo_AvC"

intersect(rownames(Endo_MvC_DEG[Endo_MvC_DEG$change == "Up",]),
          rownames(Endo_AvC_DEG[Endo_AvC_DEG$change == "Up",]))

intersect(rownames(Endo_MvC_DEG[Endo_MvC_DEG$change == "Down",]),
          rownames(Endo_AvC_DEG[Endo_AvC_DEG$change == "Down",]))

# Merge
ALL_DEG <- rbind(Endo_MvC_DEG,Endo_AvC_DEG)
ALL_DEG$type  <- factor(ALL_DEG$type,levels = c("Endo_MvC","Endo_AvC"))

ALL_DEG_sub <- ALL_DEG[ALL_DEG$change != "Stable",]

table(ALL_DEG$type)
table(ALL_DEG$label)

ggplot() + geom_point(ALL_DEG_sub, mapping=aes(x=type, y= avg_log2FC, color=change),
                      width=0.4,size=-log10(p_val)) +
  theme_test()    +
  theme(axis.text.x=element_text(size=10,angle=30,face ='bold'),
        axis.text.y=element_text(size=10,face ='bold'),
        axis.title.x=element_text(size = 10),
        axis.title.y=element_text(size = 14)) + 
  xlab(NULL) + ylab('avg_log2FC')  + 
  ylim(-1,1)+
  geom_hline(yintercept = 0,lty=1.5,lwd=2,alpha=0.5)+
  scale_colour_manual(values = c("#6B9AC7","#E24D36")) +
  guides(color=guide_legend(override.aes = list(size=6))) + #+ coord_flip() # 5*6
  geom_label_repel(ALL_DEG_sub, mapping=aes(x=type, y= avg_log2FC,label = label),
                   size = 4,
                   box.padding = unit(2, "lines"),
                   point.padding = unit(0.1, "lines"), 
                   segment.color = "black", 
                   show.legend = FALSE,
                   max.overlaps = Inf)

# 2.8. GO enrichment

## up-regularate pathways
Endo_MvC_DEG_up <- rownames(Endo_MvC_DEG[Endo_MvC_DEG$change == "Up",])

Endo_MvC_DEG_up_GO <- enrichGO(gene = Endo_MvC_DEG_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                         ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
barplot(Endo_MvC_DEG_up_GO,shwocha)

Endo_MvC_DEG_up_GO@result[c(1:20),c("Description","pvalue","Count")]

Endo_MvC_DEG_up_GO_5 <- Endo_MvC_DEG_up_GO@result[c(1,2,6,8,13),c("Description","pvalue","Count")]
Endo_MvC_DEG_up_GO_5$log10Pvalue <- -log10(Endo_MvC_DEG_up_GO_5$pvalue)
Endo_MvC_DEG_up_GO_5$subtype <- "Endo_MvC_DEG_up"
Endo_MvC_DEG_up_GO_5

Endo_AvC_DEG_up <- rownames(Endo_AvC_DEG[Endo_AvC_DEG$change == "Up",])

Endo_AvC_DEG_up_GO <- enrichGO(gene = Endo_AvC_DEG_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                                 ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
dotplot(Endo_AvC_DEG_up_GO)

Endo_AvC_DEG_up_GO_5 <- Endo_AvC_DEG_up_GO@result[c(1,17,21,26,27),c("Description","pvalue","Count")]
Endo_AvC_DEG_up_GO_5$log10Pvalue <- -log10(Endo_AvC_DEG_up_GO_5$pvalue)
Endo_AvC_DEG_up_GO_5$subtype <- "Endo_AvC_DEG_up"
Endo_AvC_DEG_up_GO_5

Go_Endo_DEG_up <- rbind(Endo_MvC_DEG_up_GO_5,Endo_AvC_DEG_up_GO_5)
Go_Endo_DEG_up$item <- row.names(Go_Endo_DEG_up)
Go_Endo_DEG_up$subtype <- factor(Go_Endo_DEG_up$subtype,levels = c("Endo_MvC_DEG_up","Endo_AvC_DEG_up"))

p1 <- ggplot(Go_Endo_DEG_up,aes(x=subtype,y=Description)) +
  geom_point(aes(size=Count,color=log10Pvalue)) + 
  scale_color_gradient(low="black",high = "red") + 
  theme_bw() + ylab("")+ xlab("")+
  theme(axis.text.y=element_text(size=12,colour="black"), 
        axis.title.y=element_text(size = 12,face="bold",colour="black"), 
        axis.text.x=element_text(size=12,face="bold",colour="black",vjust = 0.5, hjust = 0.4,angle = 45), 
        axis.title.x=element_text(size = 12,face="bold",colour="black"),
        legend.title=element_text(size=12,face="bold",colour="black"),
        # legend.text=element_text(size=12,face="bold",colour="black"), 
        title = element_text(size = 13,face="bold",colour="black"),
        # panel.border = element_blank(),axis.line = element_line(colour = "black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) # 5*7
p1

## down-regularate pathways
Endo_MvC_DEG_down <- rownames(Endo_MvC_DEG[Endo_MvC_DEG$change == "Down",])

Endo_MvC_DEG_down_GO <- enrichGO(gene = Endo_MvC_DEG_down, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                                   ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
dotplot(Endo_MvC_DEG_down_GO)

Endo_AvC_DEG_down <- rownames(Endo_AvC_DEG[Endo_AvC_DEG$change == "Down",])

Endo_AvC_DEG_down_GO <- enrichGO(gene = Endo_AvC_DEG_down, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                                   ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
dotplot(Endo_AvC_DEG_down_GO)

Endo_MvC_DEG_down_GO@result[c(1:10),c("Description")]
Endo_MvC_DEG_down_GO_5 <- Endo_MvC_DEG_down_GO@result[c(1:5),c("Description","pvalue","Count")]
Endo_MvC_DEG_down_GO_5$log10Pvalue <- -log10(Endo_MvC_DEG_down_GO_5$pvalue)
Endo_MvC_DEG_down_GO_5$subtype <- "Endo_MvC_DEG_down"
Endo_MvC_DEG_down_GO_5

Endo_AvC_DEG_down_GO@result[c(1:10),c("Description")] 
Endo_AvC_DEG_down_GO_5 <- Endo_AvC_DEG_down_GO@result[c(1:5),c("Description","pvalue","Count")]
Endo_AvC_DEG_down_GO_5$log10Pvalue <- -log10(Endo_AvC_DEG_down_GO_5$pvalue)
Endo_AvC_DEG_down_GO_5$subtype <- "Endo_AvC_DEG_down"
Endo_AvC_DEG_down_GO_5

Go_Endo_DEG_down <- rbind(Endo_MvC_DEG_down_GO_5,Endo_AvC_DEG_down_GO_5)
Go_Endo_DEG_down$item <- row.names(Go_Endo_DEG_down)
Go_Endo_DEG_down$subtype <- factor(Go_Endo_DEG_down$subtype,levels = c("Endo_MvC_DEG_down","Endo_AvC_DEG_down"))

p2 <- ggplot(Go_Endo_DEG_down,aes(x=subtype,y=Description)) +
  geom_point(aes(size=Count,color=log10Pvalue)) + 
  scale_color_gradient(low="black",high = "blue") + 
  theme_bw() + ylab("")+ xlab("")+
  theme(axis.text.y=element_text(size=12,colour="black"), 
        axis.title.y=element_text(size = 12,face="bold",colour="black"), 
        axis.text.x=element_text(size=12,face="bold",colour="black",vjust = 0.5, hjust = 0.4,angle = 45), 
        axis.title.x=element_text(size = 12,face="bold",colour="black"),
        legend.title=element_text(size=12,face="bold",colour="black"),
        # legend.text=element_text(size=12,face="bold",colour="black"), 
        title = element_text(size = 13,face="bold",colour="black"),
        # panel.border = element_blank(),axis.line = element_line(colour = "black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) # 5*7
p2

p1 / p2