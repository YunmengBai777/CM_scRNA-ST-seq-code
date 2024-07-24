The spatiotemporal transcriptional profiling of murine brain during cerebral malaria progression and artesunate treatment
---
This repository contains code from data analysis of scRNA-seq and ST-seq split into 2 parts:

- scRNA-seq analysis
- ST-seq analysis

Citation
---
The code in this repository pertains to publication:
**The spatiotemporal transcriptional profiling of murine brain during cerebral malaria progression and artesunate treatment**
Jiayun Chen, Yunmeng Bai, Xueling He, Lina Chen, Yin Kwan Wong, Chen Wang, Peng Gao, Guangqing Chen, Liting Xu, Jichao Sun, Chengchao Xu,Jigang Wang

Data
---
The raw single cell RNA sequencing and spatial transcriptome sequencing data files were deposited in the Genome Sequence Archive (GSA) under accession number **CRA007721** and **CRA007982**, respectively. 

Code description
---
scRNA-seq_analysis_part

**1.scRNA dataset process and annotation.R** | Code contains functions for preprocessing and annotation steps of scRNA-seq dataset: Quality control and SCTransform of each sample; Multiple samples integration; Dimension reduction and UMAP visulations; Cell type annotation and stastics; Cell type specific DEG analysis; GO enrichment of DEGs in each cell type.

**2.Blood brain barrier part.R** | This part investigates the transcription programs of blood brain barrier (BBB) related cells: Subtype identification; UMPA visulations; cellular proportions comparison; DEGs analysis and GO enrichment; Module scores evaluation and correlation calculation.

**3.Myeloid cell part.R** | This part investigates the transcription programs of myeloid cells:  Subtype identification; UMPA visulations; cellular proportions comparison; DEGs analysis and GO enrichment; Module scores evaluation.

**4.Lymphocytes part.R** | This part investigates the transcription programs of T and NK cells: Subtype identification; UMPA visulations; cellular proportions comparison; Pseudotime trajectory inference; GO enrichment.

**5.Cell-Cell communication.R** | Code contains functions for cell-cell communication analysis: Generate cellchat objects; Identify different enriched pathways; Visulations of key ligand-receptor pairs by *netVisual_bubble* function.

**6.Neuron part.R** | This part investigates the transcription programs of neurons: Subtype identification; UMPA visulations; cellular proportions comparison; DEG analysis; GSEA analysis; Module scores evaluation.

ST-seq_analysis_part

**1.Process and integration.R** | Code contains functions for preprocessing and annotation steps of ST-seq dataset: Quality control and SCTransform of each sample; Multiple samples integration; Dimension reduction and UMAP visulations.

**2.SingleR annotation.R** | This part applied SingleR methods to assign each spot's identity based on the MouseRNAseqData reference dataset.

**3.Seurat AddModuleScore.R** | This part evaluated cellular module scores of spots by *AddModuleScore* function based on cluster-specific markers of scRNA-seq dataset.

**4.Seurat FindTransferAnchors.R** | This part calcaulated prediction scores of spots by *TransferAnchors* function based on scRNA-seq dataset. 

**5.Spot type visualization.R** | Code shows the visualization of spot types by *DotPlot* and *SpatialDimPlot* functions.

**6.SPATA.R** | Code shows the enrichment scores of key pathways at different zones of brain slides via the SPATA analysis.

**7.cell2location.py** | Code shows the spatially resolved fine-grained cell types by integrating ST-seq and scRNA-seq reference of cell type:  Estimate reference cell type signatures via negative binomial regression model; Find shared genes and subset both ST-seq and scRNA-seq reference signatures; Create and train the model; Visualize cell abundance in spatial coordinates.

