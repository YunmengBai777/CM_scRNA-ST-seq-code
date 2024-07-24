# Myeloid subtypes

#1. library packages
library(dplyr)
library(Seurat)
library(ggsci)
library(viridis)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(pheatmap)

data<-readRDS("CM_scRNA_names.rds")
cell<- c("Micro","Macro","Mono","Neutro")
data<-subset(data,cell_type %in% as.character(cell))

#2. datasets dimension reduction and re-cluster
data <- RunPCA(data, verbose = FALSE)
ElbowPlot(data,ndims=50)
dims=1:40
data <- FindNeighbors(data, dims = dims, verbose = FALSE)
data <- FindClusters(data, resolution = 0.8, verbose = FALSE)
data <- RunUMAP(data, dims = dims, verbose = FALSE)
data <- RunTSNE(data, dims = dims, verbose = FALSE)

MajorColor <- c("#58a4c3","#ffa500","#bd9569","#00bfc4")
DimPlot(data,group.by="cell_type",label=TRUE,label.size=8)+scale_color_manual(values=MajorColor)+NoLegend()

#3. marker gene expression

gene<-c(
"Cd68","Aif1","Cx3cr1","Plac8","S100a9",
"Tmem119","P2ry12",#Home MG
"H2-Oa","H2-M3",#Ap MG
"Ccl3","Ccl4",#inflam MG
"Jun","Fos",#IEG MG
"Snrnp70","Luc7l2",#RNA splicing MG
'Mrc1',"Cd163",#home Macro
"H2-Aa","H2-Ab1",#Ap Macro
"Ccr2","Ly6c2",#classical Mono
"Ear2","Ace",#non-classical Mono
"Ltf","Ngp",##Home Neutro
"Il1b","Ly6g",#Activate Neutro
"Cxcl12","Igfbp7"#Chemo Neutro
)

DefaultAssay(data)<-"SCT"
DotPlot(data,features=rev(gene),cols=c("grey","red"))+coord_flip()+theme(axis.text.x=element_text(angle=45,hjust=1))+labs(x="",y="")

#4. subtype annotation

ann <- read.csv("cluster_subtype.csv",header=T,row.names=1)
data$subtype <- ann[as.character(data$seurat_clusters),]
data$subtype <- factor(data$subtype,levels=c("Micro_Home",
											 "Micro_Ap",
											 "Micro_Inflam",
											 "Micro_IEG",
											 "Micro_RNASplic",
											 "Macro_Home",
											 "Macro_Ap",
											 "Mono_classical",
											 "Mono_non-classical",
											 "Neutro_Home",
											 "Neutro_Activated",
											 "Neutro_Chemo"))
color <- c("#00A1D5FF","#B24745FF","#79AF97FF","#6A6599FF","#80796BFF","#1F77B4FF","#FF7F0EFF","#2CA02CFF","#9467BDFF","#DF8F44FF","#8C564BFF","#374E55FF")

#5. draw heatmap plot

dat <- AverageExpression(data)
pheatmap(dat$SCT[gene,],scale="row",cluster_cols=F,cluster_rows=F,border_color=NA,angle_col=45,color=viridis(100),gaps_col=c(5,7,9),gaps_row=c(5,15,19,23))

#6. save metadata and Lym.RDS

meta <- data@meta.data
wriet.csv(meta,"Mye_meta.csv")
saveRDS(data,"Mye_anno.RDS")

#7. draw the proportions of each group in different subtypes

dtype<-as.data.frame(table(data@meta.data[,c('subtype','sample_type')]))
p1 <- ggplot(subset(dtype,sample_type=="Control"),aes(x=sample_type,y=Freq,fill=subtype))+
		geom_bar(stat = "identity",position = 'fill')+
		theme_classic()+
		scale_fill_manual(values=color)+
		coord_polar(theta="y")+
		ggtitle("Control (n=8137)")
p2 <- ggplot(subset(dtype,sample_type=="Model"),aes(x=sample_type,y=Freq,fill=subtype))+
		geom_bar(stat = "identity",position = 'fill')+
		theme_classic()+
		scale_fill_manual(values=color)+
		coord_polar(theta="y")+
		ggtitle("Model (n=16575)")
p3 <- ggplot(subset(dtype,sample_type=="ART"),aes(x=sample_type,y=Freq,fill=subtype))+
		geom_bar(stat = "identity",position = 'fill')+
		theme_classic()+
		scale_fill_manual(values=color)+
		coord_polar(theta="y")+
		ggtitle("ART (n=12388)")
p1|p2|p3


#8. draw heatmap plot of key genes

gene <- c(
"Il1a","Ccl3","Ccl4","Csf1","Ifnb1",
"Tnf","Nfkb1","Icam1","Nfkbid","Il6","Ccl12",
"Ifng","Ccl2",
"H2-Aa","H2-Ab1",
"H2-K1","H2-D1",
"Il12b","Cox5b","C3","Ccl22",
"Il10","Isg15","Il1b"
)

data$Group <- paste(data$subtype,data$sample_type,sep="_")
dat <- AverageExpression(data,group.by="Group")
dat <- dat$SCT
a <- dat[,colnames(as.data.frame(dat))[c(17,18,16,2,3,1,23,24,22,26,27,25,29,30,28)]]
pheatmap(a[gene,],scale="row",cluster_cols=F,color=viridis(100),angle_col=45,border_color="black",gaps_col=c(3,6,9,12),cellwidth=20,cellheight=15,cluster_rows=F)


#9. DEG of Micro_Inflam cells in different groups and GO enrichment

Micro_Inflam <- subset(data,subtype=="Micro_Inflam")
Idents(Micro_Inflam) <- "sample_type"

MVC <- FindMarkers(PrepSCTFindMarkers(Micro_Inflam),ident.1 = "Model",ident.2 = "Control",
                              logfc.threshold = 0,only.pos = F,min.diff.pct = 0)
MVC$type = ifelse(MVC$p_val_adj < 0.05 & abs(MVC$avg_log2FC) >= 0.25, 
                             ifelse(MVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
MVC$type <- factor(MVC$type,levels=c("Up","Stable","Down"))
MVC$Com="ModvsCon"
MVC$Gene <- row.names(MVC)
gene <- c(subset(MVC,type=="Up")[1:5,'Gene'],subset(MVC,type=="Down")[1:5,'Gene'])
MVC$Label <- ""
MVC[gene,"Label"]=gene
p1 <- ggplot(MVC, aes(x= -log10(p_val_adj),y=avg_log2FC))+geom_point(aes(color=type))+geom_hline(yintercept=c(-0.25,0.25),linetype = 'dashed')+geom_vline(xintercept= -log10(0.05),linetype = 'dashed')+ggtitle("Micro_Inflam_DEGs\n(ModvsCon)")+theme_bw()+geom_text_repel(aes(label=Label))


AVC <- FindMarkers(PrepSCTFindMarkers(Micro_Inflam),ident.1 = "ART",ident.2 = "Control",
                              logfc.threshold = 0,only.pos = F,min.diff.pct = 0)
AVC$type = ifelse(AVC$p_val_adj < 0.05 & abs(AVC$avg_log2FC) >= 0.25, 
                             ifelse(AVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
AVC$type <- factor(AVC$type,levels=c("Up","Stable","Down"))
AVC$Com="ARTvsCon"
AVC$Gene <- row.names(AVC)
gene <- c(subset(AVC,type=="Up")[1:5,'Gene'],subset(AVC,type=="Down")[1:5,'Gene'])
AVC$Label <- ""
AVC[gene,"Label"]=gene
p2 <- ggplot(AVC, aes(x= -log10(p_val_adj),y=avg_log2FC))+geom_point(aes(color=type))+geom_hline(yintercept=c(-0.25,0.25),linetype = 'dashed')+geom_vline(xintercept= -log10(0.05),linetype = 'dashed')+ggtitle("Micro_Inflam_DEGs\n(ARTvsCon)")+theme_bw()+geom_text_repel(aes(label=Label))

pdf("Micro_Inflam_DEG.pdf",width=8,height=4)
p1|p2
dev.off()

dat <- rbind(MVC,AVC)
write.csv(dat,"Micro_Inflam_DEG.csv")
a <- subset(dat,type != "Stable")
a$Group <- paste(a$Com,a$type,sep="_")

library(clusterProfiler)
ids=bitr(a$Gene,'SYMBOL','ENTREZID','org.Mm.eg.db')
a =merge(a,ids,by.x='Gene',by.y='SYMBOL')
gcSample=split(a$ENTREZID, a$Group)
bp <- compareCluster(gcSample,fun="enrichGO", OrgDb="org.Mm.eg.db", ont= "BP")

cnetplot(bp)

enrichmentNetwork(bp@result)

#10. Ap Score distribution
 
antigen_processing_and_presentation_of_peptide_antigen <- unique(c("Ctsl","Ctss","H2-Aa","H2-Ab1","Slc11a1","B2m","Bag6","H2-D1","H2-K1","H2-L","H2-Q7","H2-Q10","H2-T23","Marchf1","Marchf8","Mr1","Abcb9","Azgp1","Calr","Cd74","Ctse","Erap1","Fcer1g","Fcgr1","Fcgr2b","Fcgr3","Gm7030","Gm8909","Gm11127","H2-DMa","H2-DMb1","H2-DMb2","H2-Ea","H2-Eb1","H2-Eb2","H2-M1","H2-M2","H2-M3","H2-M5","H2-M9","H2-M10.1","H2-M10.2","H2-M10.3","H2-M10.4","H2-M10.5","H2-M10.6","H2-M11","H2-Oa","H2-Ob","H2-Q1","H2-Q2","H2-Q4","H2-Q6","H2-Q8","H2-Q9","H2-T3","H2-T22","H2-T24","Hfe","Ide","Ifi30","Mfsd6","Pdia3","Pikfyve","Pycard","Tap1","Tap2","Tapbp","Tapbpl","Traf6","Trem2","Unc93b1"))

data <- AddModuleScore(object = data,features = list(antigen_processing_and_presentation_of_peptide_antigen),assay="SCT",name='antigen_processing_and_presentation_of_peptide_antigen')

pdf("sc_Ap.pdf")
VlnPlot(data,"antigen_processing_and_presentation_of_peptide_antigen1",
        group.by = "sample_type",pt.size = 0,sort = F,cols = c("#3B4992","#EE0000","#008B45")) + 
  ggtitle("antigen processing and presentation of \n peptide antigen module score in ScRNA dataset") + xlab("")
dev.off()

STdata <- readRDS("STdata.RDS")
STdata <- AddModuleScore(object = STdata,features = list(antigen_processing_and_presentation_of_peptide_antigen),assay="SCT",name='antigen_processing_and_presentation_of_peptide_antigen')
VlnPlot(STdata,"antigen_processing_and_presentation_of_peptide_antigen1",
        group.by = "sample_type",pt.size = 0,sort = F,cols = c("#3B4992","#EE0000","#008B45")) + 
  ggtitle("antigen processing and presentation of \n peptide antigen module score in ST dataset") + xlab("")
