# 5. library packages
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)
library(ggsci)
library(VennDiagram)

# 5.1 loading scRNA-seq data
data<-readRDS("CM_scRNA_dataset.rds")#scRNA dataset
meta <- data@meta.data

Control_meta <- subset(meta,sample_type=="Control")
Model_meta <- subset(meta,sample_type=="Model")
ART_meta <- subset(meta,sample_type=="ART")

Control_data<-as.matrix(data@assays$SCT@data[,as.character(row.names(Control_meta))])
Model_data<-as.matrix(data@assays$SCT@data[,as.character(row.names(Model_meta))])
ART_data<-as.matrix(data@assays$SCT@data[,as.character(row.names(ART_meta))])

color = c(
"Endo" = '#F3B1A0',
"Macro" = '#58A4C3',
"Micro" = 'orange',
"Mono" = '#BD956A',
"Neuron" = '#E95C59',
"T/NK cell" = 'red')

# 5.2 cell-cell communication in each group

cellchat <- createCellChat(object = Control_data, meta = Control_meta, group.by = "cell_type")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 40) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat) # project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat) # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10) # Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat) # Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",color.use=color)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",color.use=color)
saveRDS(cellchat, "ScRNA_Control_cellchat.RDS")

cellchat <- createCellChat(object = Model_data, meta = Model_meta, group.by = "cell_type")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 40) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat) # project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat) # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10) #Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat) #Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",color.use=color)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",color.use=color)
saveRDS(cellchat, "ScRNA_Model_cellchat.RDS")

cellchat <- createCellChat(object = ART_data, meta = ART_meta, group.by = "cell_type")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 40) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat) # project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat) # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10) #Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat) #Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",color.use=color)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",color.use=color)
saveRDS(cellchat, "ScRNA_ART_cellchat.RDS")

# 5.3 Merge cell-cell communication in scRNA-seq dataset

control<-readRDS( "Control_cellchat.RDS")
model<-readRDS( "Model_cellchat.RDS") 
art<-readRDS( "ART_cellchat.RDS")

object.list <- list(Control=control,Model=model,ART=art)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
netVisual_heatmap(cellchat,color.use=color)

netVisual_heatmap(cellchat,comparison=c("Control","Model"),title.name = "Different number of interactions of MVC\nin scRNA-seq dataset",font.size=16,font.size.title = 18)+ggtitle("ModvsCon")

netVisual_heatmap(cellchat,comparison=c("Control","ART"),title.name = "Different number of interactions of AVC\nin scRNA-seq dataset",font.size=16,font.size.title = 18)+ggtitle("ARTvsCon")

# 5.4 Different enriched pathways in scRNA-seq and ST-seq datasets

SC_MvsC_Path=c("MHC-I","LCK","IFN-II","PARs","ANNEXIN","CD86","ncWNT","ICOS","CALCR","TWEAK","VCAM","THY1","ITGAL-ITGB2","LAIR1")
SC_AvsC_Path=c("MHC-I","LCK","CSF","CX3C","CADM","CD86","PARs","IFN-II","ncWNT","SN","ANNEXIN","WNT","ICOS","APRIL","FASLG","CALCR","IL1","COMPLEMENT","ITGAL-ITGB2","CCL","ICAM","VCAM","SEMA4","LAIR1","SEMA7")

rankNet(sc_cellchat, mode = "comparison", stacked = T,do.stat =TRUE,comparison = c(1,2,3),color.use=c("#3B4992","#EE0000","#008B45" ),signaling =unique(c(SC_MvsC_Path,SC_AvcC_Path)))+coord_flip()

# 5.5 different pathway overlap across different groups in scRNA-seq and ST-seq datastes

venn.diagram(x=list(SC_MvsC=SC_MvsC_Path,SC_AvsC=SC_AvsC_Path,ST_MvsC=ST_MvsC_Path,ST_AvsC=ST_AvsC_Path),filename=NULL,fill=pal_jama()(4),cex = 1.5, cat.col = 'black', cat.cex = 1.5, cat.fontface = "bold", margin = 0.05,rotation.degree = 0)
grid.draw(plot)

# 5.6 ligand-receptor pairs expression across different groups in scRNA-seq dataset 

Path <- c("MHC-I","ITGAL-ITGB2","IFN-II","ICAM","CCL","VCAM")
pairLR.use <- extractEnrichedLR(sc_cellchat, signaling = Path)
pairLR.use  <- pairLR.use[c(8,19:24,103:106,115,116,123:125),,drop=FALSE]
netVisual_bubble(cellchat_sub, comparison = c(1, 2, 3),  angle.x = 45, remove.isolate = F, pairLR.use = pairLR.use, color.text = c("#3B4992","#EE0000","#008B45" ))






