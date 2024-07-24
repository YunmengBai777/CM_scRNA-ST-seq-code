# 6.1 Loading package

library(SPATA)
library(hdf5r)
library(magrittr)
library(ggplot2)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(Seurat)
library(gsva) 
library(homologene)
library(ggsci)

# 6.2 Loading dataset

BM1_ST_A <-  initiateSpataObject_10X(input_paths = c("1.Data/1.ST_data/BM1_ST_A/"), sample_names = c("BM1_ST_A"))
BM1_ST_P <-  initiateSpataObject_10X(input_paths = c("1.Data/1.ST_data/BM1_ST_P/"), sample_names = c("BM1_ST_P"))
BM2_ST_A <-  initiateSpataObject_10X(input_paths = c("1.Data/1.ST_data/BM2_ST_A/"), sample_names = c("BM2_ST_A"))
BM2_ST_P <-  initiateSpataObject_10X(input_paths = c("1.Data/1.ST_data/BM2_ST_P/"), sample_names = c("BM2_ST_P"))
BA1_ST_A <-  initiateSpataObject_10X(input_paths = c("1.Data/1.ST_data/BA1_ST_A/"), sample_names = c("BA1_ST_A"))
BA1_ST_P <-  initiateSpataObject_10X(input_paths = c("1.Data/1.ST_data/BA1_ST_P/"), sample_names = c("BA1_ST_P"))
BA2_ST_A <-  initiateSpataObject_10X(input_paths = c("1.Data/1.ST_data/BA2_ST_A/"), sample_names = c("BA2_ST_A"))
BA2_ST_P <-  initiateSpataObject_10X(input_paths = c("1.Data/1.ST_data/BA2_ST_P/"), sample_names = c("BA2_ST_P"))

# 6.3 Gene ID mapping

BA1_ST_A_genename <- BA1_ST_A@data@counts@Dimnames[[1]]
BA1_ST_A_genename_m2h <- homologene(BA1_ST_A_genename, inTax = 10090, outTax = 9606)
BA1_ST_A@data@counts <- BA1_ST_A@data@counts[BA1_ST_A_genename_m2h$`10090`,]
BA1_ST_A@data@norm_exp <- BA1_ST_A@data@norm_exp[BA1_ST_A_genename_m2h$`10090`,]
table(row.names(BA1_ST_A@data@counts) == rownames(BA1_ST_A@data@norm_exp))
row.names(BA1_ST_A@data@counts) <- BA1_ST_A_genename_m2h$`9606`
row.names(BA1_ST_A@data@norm_exp) <- BA1_ST_A_genename_m2h$`9606`

BA1_ST_P_genename <- BA1_ST_P@data@counts@Dimnames[[1]]
BA1_ST_P_genename_m2h <- homologene(BA1_ST_P_genename, inTax = 10090, outTax = 9606)
BA1_ST_P@data@counts <- BA1_ST_P@data@counts[BA1_ST_P_genename_m2h$`10090`,]
BA1_ST_P@data@norm_exp <- BA1_ST_P@data@norm_exp[BA1_ST_P_genename_m2h$`10090`,]
table(row.names(BA1_ST_P@data@counts) == rownames(BA1_ST_P@data@norm_exp))
row.names(BA1_ST_P@data@counts) <- BA1_ST_P_genename_m2h$`9606`
row.names(BA1_ST_P@data@norm_exp) <- BA1_ST_P_genename_m2h$`9606`

BA2_ST_A_genename <- BA2_ST_A@data@counts@Dimnames[[1]]
BA2_ST_A_genename_m2h <- homologene(BA2_ST_A_genename, inTax = 10090, outTax = 9606)
BA2_ST_A@data@counts <- BA2_ST_A@data@counts[BA2_ST_A_genename_m2h$`10090`,]
BA2_ST_A@data@norm_exp <- BA2_ST_A@data@norm_exp[BA2_ST_A_genename_m2h$`10090`,]
table(row.names(BA2_ST_A@data@counts) == rownames(BA2_ST_A@data@norm_exp))
row.names(BA2_ST_A@data@counts) <- BA2_ST_A_genename_m2h$`9606`
row.names(BA2_ST_A@data@norm_exp) <- BA2_ST_A_genename_m2h$`9606`

BA2_ST_P_genename <- BA2_ST_P@data@counts@Dimnames[[1]]
BA2_ST_P_genename_m2h <- homologene(BA2_ST_P_genename, inTax = 10090, outTax = 9606)
BA2_ST_P@data@counts <- BA2_ST_P@data@counts[BA2_ST_P_genename_m2h$`10090`,]
BA2_ST_P@data@norm_exp <- BA2_ST_P@data@norm_exp[BA2_ST_P_genename_m2h$`10090`,]
table(row.names(BA2_ST_P@data@counts) == rownames(BA2_ST_P@data@norm_exp))
row.names(BA2_ST_P@data@counts) <- BA2_ST_P_genename_m2h$`9606`
row.names(BA2_ST_P@data@norm_exp) <- BA2_ST_P_genename_m2h$`9606`

BM1_ST_A_genename <- BM1_ST_A@data@counts@Dimnames[[1]]
BM1_ST_A_genename_m2h <- homologene(BM1_ST_A_genename, inTax = 10090, outTax = 9606)
BM1_ST_A@data@counts <- BM1_ST_A@data@counts[BM1_ST_A_genename_m2h$`10090`,]
BM1_ST_A@data@norm_exp <- BM1_ST_A@data@norm_exp[BM1_ST_A_genename_m2h$`10090`,]
table(row.names(BM1_ST_A@data@counts) == rownames(BM1_ST_A@data@norm_exp))
row.names(BM1_ST_A@data@counts) <- BM1_ST_A_genename_m2h$`9606`
row.names(BM1_ST_A@data@norm_exp) <- BM1_ST_A_genename_m2h$`9606`

BM1_ST_P_genename <- BM1_ST_P@data@counts@Dimnames[[1]]
BM1_ST_P_genename_m2h <- homologene(BM1_ST_P_genename, inTax = 10090, outTax = 9606)
BM1_ST_P@data@counts <- BM1_ST_P@data@counts[BM1_ST_P_genename_m2h$`10090`,]
BM1_ST_P@data@norm_exp <- BM1_ST_P@data@norm_exp[BM1_ST_P_genename_m2h$`10090`,]
table(row.names(BM1_ST_P@data@counts) == rownames(BM1_ST_P@data@norm_exp))
row.names(BM1_ST_P@data@counts) <- BM1_ST_P_genename_m2h$`9606`
row.names(BM1_ST_P@data@norm_exp) <- BM1_ST_P_genename_m2h$`9606`

BM2_ST_A_genename <- BM2_ST_A@data@counts@Dimnames[[1]]
BM2_ST_A_genename_m2h <- homologene(BM2_ST_A_genename, inTax = 10090, outTax = 9606)
BM2_ST_A@data@counts <- BM2_ST_A@data@counts[BM2_ST_A_genename_m2h$`10090`,]
BM2_ST_A@data@norm_exp <- BM2_ST_A@data@norm_exp[BM2_ST_A_genename_m2h$`10090`,]
table(row.names(BM2_ST_A@data@counts) == rownames(BM2_ST_A@data@norm_exp))
row.names(BM2_ST_A@data@counts) <- BM2_ST_A_genename_m2h$`9606`
row.names(BM2_ST_A@data@norm_exp) <- BM2_ST_A_genename_m2h$`9606`

BM2_ST_P_genename <- BM2_ST_P@data@counts@Dimnames[[1]]
BM2_ST_P_genename_m2h <- homologene(BM2_ST_P_genename, inTax = 10090, outTax = 9606)
BM2_ST_P@data@counts <- BM2_ST_P@data@counts[BM2_ST_P_genename_m2h$`10090`,]
BM2_ST_P@data@norm_exp <- BM2_ST_P@data@norm_exp[BM2_ST_P_genename_m2h$`10090`,]
table(row.names(BM2_ST_P@data@counts) == rownames(BM2_ST_P@data@norm_exp))
row.names(BM2_ST_P@data@counts) <- BM2_ST_P_genename_m2h$`9606`
row.names(BM2_ST_P@data@norm_exp) <- BM2_ST_P_genename_m2h$`9606`

# 6.4 PlotSurface in SPATA

BM1_ST_A_plot <- 
  (plotSurface(object = BM1_ST_A,of_sample = "BM1_ST_A",color_to = "AIF1",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = T,method_gs = "mean") + NoLegend() + 
     plotSurface(object = BM1_ST_P,of_sample = "BM1_ST_P",color_to = "AIF1",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = F,method_gs = "mean") + NoLegend() + ggtitle("BM1")) / 
  (plotSurface(object = BM1_ST_A,of_sample = "BM1_ST_A",color_to = "HLA-DQA1",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = T,method_gs = "mean") + NoLegend() + 
     plotSurface(object = BM1_ST_P,of_sample = "BM1_ST_P",color_to = "HLA-DQA1",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = F,method_gs = "mean") + NoLegend())  /
  (plotSurface(object = BM1_ST_A,of_sample = "BM1_ST_A",color_to = "BP.GO_MYELOID_LEUKOCYTE_ACTIVATION",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = T,method_gs = "mean") + NoLegend() + 
     plotSurface(object = BM1_ST_P,of_sample = "BM1_ST_P",color_to = "BP.GO_MYELOID_LEUKOCYTE_ACTIVATION",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = F,method_gs = "mean")+ NoLegend() ) / 
  (plotSurface(object = BM1_ST_A,of_sample = "BM1_ST_A",color_to = "BP.GO_ANTIGEN_PROCESSING_AND_PRESENTATION",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = T,method_gs = "mean") + NoLegend() + 
     plotSurface(object = BM1_ST_P,of_sample = "BM1_ST_P",color_to = "BP.GO_ANTIGEN_PROCESSING_AND_PRESENTATION",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = F,method_gs = "mean")+ NoLegend()) 

BA1_ST_A_plot <- 
  (plotSurface(object = BA1_ST_A,of_sample = "BA1_ST_A",color_to = "AIF1",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = T,method_gs = "mean") + NoLegend() + 
     plotSurface(object = BA1_ST_P,of_sample = "BA1_ST_P",color_to = "AIF1",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = F,method_gs = "mean")  + ggtitle("BA1")) / 
  (plotSurface(object = BA1_ST_A,of_sample = "BA1_ST_A",color_to = "HLA-DQA1",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = F,method_gs = "mean") + NoLegend() + 
     plotSurface(object = BA1_ST_P,of_sample = "BA1_ST_P",color_to = "HLA-DQA1",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = F,method_gs = "mean"))  /
  (plotSurface(object = BA1_ST_A,of_sample = "BA1_ST_A",color_to = "BP.GO_MYELOID_LEUKOCYTE_ACTIVATION",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = F,method_gs = "mean") + NoLegend() + 
     plotSurface(object = BA1_ST_P,of_sample = "BA1_ST_P",color_to = "BP.GO_MYELOID_LEUKOCYTE_ACTIVATION",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = F,method_gs = "mean")) / 
  (plotSurface(object = BA1_ST_A,of_sample = "BA1_ST_A",color_to = "BP.GO_ANTIGEN_PROCESSING_AND_PRESENTATION",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = F,method_gs = "mean") + NoLegend() + 
     plotSurface(object = BA1_ST_P,of_sample = "BA1_ST_P",color_to = "BP.GO_ANTIGEN_PROCESSING_AND_PRESENTATION",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = F,method_gs = "mean")) 

BM1_ST_A_plot | BA1_ST_A_plot # 8*8

BM2_ST_A_plot <- 
  (plotSurface(object = BM2_ST_A,of_sample = "BM2_ST_A",color_to = "AIF1",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = T,method_gs = "mean") + NoLegend() + 
     plotSurface(object = BM2_ST_P,of_sample = "BM2_ST_P",color_to = "AIF1",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = F,method_gs = "mean") + NoLegend() + ggtitle("BM2")) / 
  (plotSurface(object = BM2_ST_A,of_sample = "BM2_ST_A",color_to = "HLA-DQA1",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = T,method_gs = "mean") + NoLegend() + 
     plotSurface(object = BM2_ST_P,of_sample = "BM2_ST_P",color_to = "HLA-DQA1",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = F,method_gs = "mean") + NoLegend())  /
  (plotSurface(object = BM2_ST_A,of_sample = "BM2_ST_A",color_to = "BP.GO_MYELOID_LEUKOCYTE_ACTIVATION",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = T,method_gs = "mean") + NoLegend() + 
     plotSurface(object = BM2_ST_P,of_sample = "BM2_ST_P",color_to = "BP.GO_MYELOID_LEUKOCYTE_ACTIVATION",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = F,method_gs = "mean")+ NoLegend() ) / 
  (plotSurface(object = BM2_ST_A,of_sample = "BM2_ST_A",color_to = "BP.GO_ANTIGEN_PROCESSING_AND_PRESENTATION",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = T,method_gs = "mean") + NoLegend() + 
     plotSurface(object = BM2_ST_P,of_sample = "BM2_ST_P",color_to = "BP.GO_ANTIGEN_PROCESSING_AND_PRESENTATION",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = F,method_gs = "mean")+ NoLegend()) 

BA2_ST_A_plot <- 
  (plotSurface(object = BA2_ST_A,of_sample = "BA2_ST_A",color_to = "AIF1",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = F,method_gs = "mean") + NoLegend() + 
     plotSurface(object = BA2_ST_P,of_sample = "BA2_ST_P",color_to = "AIF1",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = F,method_gs = "mean")  + ggtitle("BA2")) / 
  (plotSurface(object = BA2_ST_A,of_sample = "BA2_ST_A",color_to = "HLA-DQA1",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = F,method_gs = "mean") + NoLegend() + 
     plotSurface(object = BA2_ST_P,of_sample = "BA2_ST_P",color_to = "HLA-DQA1",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = F,method_gs = "mean"))  /
  (plotSurface(object = BA2_ST_A,of_sample = "BA2_ST_A",color_to = "BP.GO_MYELOID_LEUKOCYTE_ACTIVATION",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = F,method_gs = "mean") + NoLegend() + 
     plotSurface(object = BA2_ST_P,of_sample = "BA2_ST_P",color_to = "BP.GO_MYELOID_LEUKOCYTE_ACTIVATION",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = F,method_gs = "mean")) / 
  (plotSurface(object = BA2_ST_A,of_sample = "BA2_ST_A",color_to = "BP.GO_ANTIGEN_PROCESSING_AND_PRESENTATION",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = F,method_gs = "mean") + NoLegend() + 
     plotSurface(object = BA2_ST_P,of_sample = "BA2_ST_P",color_to = "BP.GO_ANTIGEN_PROCESSING_AND_PRESENTATION",pt_size = 0.1,pt_clrsp = "magma",smooth = TRUE,smooth_span = 0.01,display_title = F,method_gs = "mean")) 

BM2_ST_A_plot | BA2_ST_A_plot # 8*8

# 6.5 plotSurfaceComparison in SPATA

((plotSurfaceComparison(object = BM1_ST_A,of_sample = "BM1_ST_A",variables = "BP.GO_ANTIGEN_PROCESSING_AND_PRESENTATION_VIA_MHC_CLASS_IB",
                        pt_size =  0.01,pt_clrsp = "magma"  ,method_gs = "gsva") + 
    plotSurfaceComparison(object = BM1_ST_P,of_sample = "BM1_ST_P",variables = "BP.GO_ANTIGEN_PROCESSING_AND_PRESENTATION_VIA_MHC_CLASS_IB",
                          pt_size =  0.01,pt_clrsp = "magma"  ,method_gs = "gsva"))) / 
  ((plotSurfaceComparison(object = BA2_ST_A,of_sample = "BA2_ST_A",variables = "BP.GO_ANTIGEN_PROCESSING_AND_PRESENTATION_VIA_MHC_CLASS_IB",
                          pt_size =  0.01,pt_clrsp = "magma"  ,method_gs = "mean") + 
      plotSurfaceComparison(object = BA2_ST_P,of_sample = "BA2_ST_P",variables = "BP.GO_ANTIGEN_PROCESSING_AND_PRESENTATION_VIA_MHC_CLASS_IB",
                            pt_size =  0.01,pt_clrsp = "magma"  ,method_gs = "mean"))) /
  ((plotSurfaceComparison(object = BM2_ST_A,of_sample = "BM2_ST_A",variables = "BP.GO_ANTIGEN_PROCESSING_AND_PRESENTATION_VIA_MHC_CLASS_IB",
                          pt_size =  0.01,pt_clrsp = "magma"  ,method_gs = "gsva") + 
      plotSurfaceComparison(object = BM2_ST_P,of_sample = "BM2_ST_P",variables = "BP.GO_ANTIGEN_PROCESSING_AND_PRESENTATION_VIA_MHC_CLASS_IB",
                            pt_size =  0.01,pt_clrsp = "magma"  ,method_gs = "gsva"))) / 
  ((plotSurfaceComparison(object = BA1_ST_A,of_sample = "BA1_ST_A",variables = "BP.GO_ANTIGEN_PROCESSING_AND_PRESENTATION_VIA_MHC_CLASS_IB",
                          pt_size =  0.01,pt_clrsp = "magma"  ,method_gs = "mean") + 
      plotSurfaceComparison(object = BA1_ST_P,of_sample = "BA1_ST_P",variables = "BP.GO_ANTIGEN_PROCESSING_AND_PRESENTATION_VIA_MHC_CLASS_IB",
                            pt_size =  0.01,pt_clrsp = "magma"  ,method_gs = "mean"))) 
