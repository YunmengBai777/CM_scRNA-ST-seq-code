# 6.Neuron type

# 6.1 library packages

library(Seurat)
library(ggplot2)
library(ggsci)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)
options(stringsAsFactors = F)

# 6.2 load scRNA datasets

scc_integrated <- readRDS("CM_scRNA_dataset.rds")
Neuron <- subset(scc_integrated,idents = c("Neuron"))
remove(scc_integrated)

Neuron <- SCTransform(Neuron,return.only.var.genes = FALSE,assay = "RNA",verbose = FALSE)
Neuron <- RunPCA(object = Neuron,verbose = FALSE)
Neuron <- FindNeighbors(Neuron,dim=1:20)
Neuron <- FindClusters(Neuron,resolution = 0.8)
Neuron <- RunUMAP(Neuron,reduction="pca", dims = 1:20)
# Neuron <- RunTSNE(Neuron,reduction="pca", dims = 1:20)

Neuron_markers = FindAllMarkers(Neuron,only.pos =T,logfc.threshold = 0.25)
Neuron_markers_top10 <- Neuron_markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC)

DoHeatmap(Neuron_sub,c("Sox11","Nrp1","Ntn1","Snap25",
                       "Dcx","Tubb3","Tbr1",
                       "Rbfox3","Map2","Dlg4", 
                       # "Mki67",
                       # "Ascl1", #"Pax6",
                       # "Eno2","Uchl1","Vsnl1","S100b", 
                       "Camk2a", "Cbln2", "Ldb2", "Slc17a6","Slc17a7","Calb2", #  excitatory markers
                       "Gad1","Gad2","Lhfpl3", "Pcdh15","Slc32a1" #  inhibitory markers
))

# 6.3 Rename the cluster

new.cluster.ids <- c("Neuron_IM","Neuron_IM","Neuron_Pro","Neuron_IP","Neuron_IM","Neuron_M_ex","Neuron_M_in","Neuron_M_ex","Neuron_M_ex")

names(new.cluster.ids) <- levels(Neuron)
Neuron <- RenameIdents(Neuron, new.cluster.ids)

table(Neuron@active.ident)
# Neuron_IM  Neuron_Pro   Neuron_IP Neuron_M_ex Neuron_M_in 
# 589         169         138         235         122

levels(Neuron) <- c("Neuron_Pro","Neuron_IP","Neuron_IM","Neuron_M_ex","Neuron_M_in")
DimPlot(Neuron,pt.size = 1,reduction = "tsne",label = T,label.size = 5) + scale_color_lancet()

saveRDS(Neuron,"Neuron/Neuron_anno.rds")

# 6.4 AverageExpression value

AverageExpression_value <-  AverageExpression(Neuron, assays = "SCT", 
                                              features = c(
                                                "Sox11","Ntn1","Nrxn3", # Neuron markers 
                                                "Mki67","Ascl1", #"Pax6", 
                                                "Dcx","Tubb3",#"Tbr1","Neurod1", # immature markers
                                                "Rbfox3","Dlg4", # mature markers
                                                "Camk2a","Slc17a6","Slc17a7", #  excitatory markers
                                                "Gad1","Gad2","Slc32a1", #  inhibitory markers
                                                "Eno2","Uchl1"#,"Vsnl1","S100b" # injury markers
                                              ),return.seurat = FALSE, group.by = "ident")

pheatmap::pheatmap(AverageExpression_value$SCT,
                   cluster_rows = F,cluster_cols = F,
                   border_color = "black",
                   angle_col = 45,
                   # gaps_col = c(1,5,6,7),
                   gaps_row = c(3,4,5,7,9,12,15,17),
                   cellwidth = 20,cellheight = 20,
                   scale = "row",viridis(100))


# 6.5 cell proportion

dfsam <- as.data.frame(table(Neuron$sample_type,Neuron@active.ident))

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
  theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + ggtitle("Control (n=648)") +
  scale_fill_lancet(name = "",labels = myLabel)

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
  theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + ggtitle("Model (n=300)") +
  scale_fill_lancet(name = "",labels = myLabel)
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
  theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + ggtitle("ART (n=305)") +
  scale_fill_lancet(name = "",labels = myLabel)
p3

p1 / p2 / p3 #5*15

table(Neuron$sample_type)


# 6.6 Brain injury markers comparision

VlnPlot(subset(subset(Neuron,idents = c("Neuron_IM","Neuron_M_ex")),Eno2>0),
        "Eno2",group.by = "cell_subtype",split.by = "sample_type",
        pt.size = 0,y.max = 2,cols = c("#3B4992","#EE0000","#008B45"))/
  VlnPlot(subset(subset(Neuron,idents = c("Neuron_IM","Neuron_M_ex")),Uchl1>0),
        "Uchl1",group.by = "cell_subtype",split.by = "sample_type",
        pt.size = 0,cols = c("#3B4992","#EE0000","#008B45"))


Neuron@active.ident <- Neuron$sample_type
VlnPlot(subset(Neuron,Eno2>0),
        "Eno2",#group.by = "cell_subtype",split.by = "sample_type",
        pt.size = 0,y.max = 2,cols = c("#3B4992","#EE0000","#008B45")) /
  VlnPlot(subset(Neuron,Uchl1>0),
          "Uchl1",#group.by = "cell_subtype",split.by = "sample_type",
          pt.size = 0,cols = c("#3B4992","#EE0000","#008B45"))

# 6.7 DEGs analysis

Neuron@active.ident <- Neuron$sample_type

Neuron_MvC_DEG <- FindMarkers(Neuron,ident.1 = "Model",ident.2 = "Control",
                              min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
Neuron_MvC_DEG$change = ifelse(Neuron_MvC_DEG$p_val_adj < 0.05 & abs(Neuron_MvC_DEG$avg_log2FC) >= 0.25, 
                               ifelse(Neuron_MvC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Neuron_MvC_DEG$change)
# Down Stable     Up 
# 39     78     21

Neuron_MvC_DEG$type <- "Neuron_MvC"

Neuron_AvC_DEG <- FindMarkers(Neuron,ident.1 = "ART",ident.2 = "Control",
                              min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
Neuron_AvC_DEG$change = ifelse(Neuron_AvC_DEG$p_val_adj < 0.05 & abs(Neuron_AvC_DEG$avg_log2FC) >= 0.25, 
                               ifelse(Neuron_AvC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Neuron_AvC_DEG$change)
# Down Stable     Up 
# 24     83     56

Neuron_AvC_DEG$type <- "Neuron_AvC"

intersect(rownames(Neuron_MvC_DEG[Neuron_MvC_DEG$change == "Up",]),
          rownames(Neuron_AvC_DEG[Neuron_AvC_DEG$change == "Up",]))
# "Cirbp" "Ccl5"  "H2-K1" "Rbm3"  "B2m"   "H2-D1"

intersect(rownames(Neuron_MvC_DEG[Neuron_MvC_DEG$change == "Down",]),
          rownames(Neuron_AvC_DEG[Neuron_AvC_DEG$change == "Down",]))
# "Gm42418" "Ttr"     "Plp1"    "Lars2"   "mt-Atp8" "mt-Nd4l" "Cst3"    "Enpp2"   "Hexb"! 
# "Mbp"     "Dbi"!

Neuron_MvC_DEG$label = ifelse(rownames(Neuron_MvC_DEG) %in% c("Ccl4","Ccl5","H2-K1","B2m","H2-D1",
                                                              "Hexb","Plp1","Dbi","Mbp"), 
                              rownames(Neuron_MvC_DEG),'')

Neuron_AvC_DEG$label = ifelse(rownames(Neuron_AvC_DEG) %in% c("Ccl5","H2-K1","B2m","H2-D1",
                                                              "Hexb","Plp1","Dbi","Mbp"), 
                              rownames(Neuron_AvC_DEG),'')

## Merge
ALL_DEG <- rbind(Neuron_MvC_DEG,Neuron_AvC_DEG)
ALL_DEG$type  <- factor(ALL_DEG$type,levels = c("Neuron_MvC","Neuron_AvC"))

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

# 6.8 GO enrichment

## up-regularate pathways
Neuron_MvC_DEG_up <- rownames(Neuron_MvC_DEG[Neuron_MvC_DEG$change == "Up",])

Neuron_MvC_DEG_up_GO <- enrichGO(gene = Neuron_MvC_DEG_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                          ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
barplot(Neuron_MvC_DEG_up_GO)

Neuron_MvC_DEG_up_GO_5 <- Neuron_MvC_DEG_up_GO@result[c(1,2,6,8,13),c("Description","pvalue","Count")]
Neuron_MvC_DEG_up_GO_5$log10Pvalue <- -log10(Neuron_MvC_DEG_up_GO_5$pvalue)
Neuron_MvC_DEG_up_GO_5$subtype <- "Neuron_MvC_DEG_up"
Neuron_MvC_DEG_up_GO_5

Neuron_AvC_DEG_up <- rownames(Neuron_AvC_DEG[Neuron_AvC_DEG$change == "Up",])

Neuron_AvC_DEG_up_GO <- enrichGO(gene = Neuron_AvC_DEG_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                                 ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
dotplot(Neuron_AvC_DEG_up_GO)

Neuron_AvC_DEG_up_GO_5 <- Neuron_AvC_DEG_up_GO@result[c(1,17,21,26,27),c("Description","pvalue","Count")]
Neuron_AvC_DEG_up_GO_5$log10Pvalue <- -log10(Neuron_AvC_DEG_up_GO_5$pvalue)
Neuron_AvC_DEG_up_GO_5$subtype <- "Neuron_AvC_DEG_up"
Neuron_AvC_DEG_up_GO_5

Go_Neuron_DEG_up <- rbind(Neuron_MvC_DEG_up_GO_5,Neuron_AvC_DEG_up_GO_5)
Go_Neuron_DEG_up$item <- row.names(Go_Neuron_DEG_up)
Go_Neuron_DEG_up$subtype <- factor(Go_Neuron_DEG_up$subtype,levels = c("Neuron_MvC_DEG_up","Neuron_AvC_DEG_up"))

p1 <- ggplot(Go_Neuron_DEG_up,aes(x=subtype,y=Description)) +
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
Neuron_MvC_DEG_down <- rownames(Neuron_MvC_DEG[Neuron_MvC_DEG$change == "Down",])

Neuron_MvC_DEG_down_GO <- enrichGO(gene = Neuron_MvC_DEG_down, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                                   ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
dotplot(Neuron_MvC_DEG_down_GO)

Neuron_AvC_DEG_down <- rownames(Neuron_AvC_DEG[Neuron_AvC_DEG$change == "Down",])

Neuron_AvC_DEG_down_GO <- enrichGO(gene = Neuron_AvC_DEG_down, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                                   ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
dotplot(Neuron_AvC_DEG_down_GO)

Neuron_MvC_DEG_down_GO@result[c(1:10),c("Description")]
Neuron_MvC_DEG_down_GO_5 <- Neuron_MvC_DEG_down_GO@result[c(1:5),c("Description","pvalue","Count")]
Neuron_MvC_DEG_down_GO_5$log10Pvalue <- -log10(Neuron_MvC_DEG_down_GO_5$pvalue)
Neuron_MvC_DEG_down_GO_5$subtype <- "Neuron_MvC_DEG_down"
Neuron_MvC_DEG_down_GO_5

Neuron_AvC_DEG_down_GO@result[c(1:10),c("Description")] 
Neuron_AvC_DEG_down_GO_5 <- Neuron_AvC_DEG_down_GO@result[c(1:5),c("Description","pvalue","Count")]
Neuron_AvC_DEG_down_GO_5$log10Pvalue <- -log10(Neuron_AvC_DEG_down_GO_5$pvalue)
Neuron_AvC_DEG_down_GO_5$subtype <- "Neuron_AvC_DEG_down"
Neuron_AvC_DEG_down_GO_5

Go_Neuron_DEG_down <- rbind(Neuron_MvC_DEG_down_GO_5,Neuron_AvC_DEG_down_GO_5)
Go_Neuron_DEG_down$item <- row.names(Go_Neuron_DEG_down)
Go_Neuron_DEG_down$subtype <- factor(Go_Neuron_DEG_down$subtype,levels = c("Neuron_MvC_DEG_down","Neuron_AvC_DEG_down"))

p2 <- ggplot(Go_Neuron_DEG_down,aes(x=subtype,y=Description)) +
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



