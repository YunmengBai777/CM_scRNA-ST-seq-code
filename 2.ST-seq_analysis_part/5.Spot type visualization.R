# 5.Spot type visualization

# 5.1 Markers gene (Dotplot)

ST$level <- factor(ST$SingleR.labels,levels = c("Neuron","Oligo","Astro","CPC","Endo","Fibro","Immune"))

DotPlot(ST,features = rev(c("Nrgn","Nrn1", # Neu
                        "Mbp","Olig1", # Oligo
                        "Aldoc","Gfap",#Astro
                        "Ttr","Foxj1", # Epi
                        "Cldn5", # Endo
                        "Dcn", # Fibro
                        "Cd68","Aif1",#"Cx3cr1","Tmem119",
                        "Cd3d","Nkg7" # lymph
)),group.by = "level")  +
  xlab(NULL) + ylab(NULL) + scale_color_gradientn(colours = viridis(20), 
  guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"), name = "Module socre") +
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1) +
  annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) + coord_flip() # 5*8


# 5.2 SpatialDimPlot

BC1_ST_A <- (SpatialDimPlot(ST,images = c("anterior2"),cols = cols,alpha = 0,stroke = 0) + NoLegend() + 
               ggtitle("BC1_ST_A") + theme(panel.border = element_rect(fill = NA,colour = "black"))) /
  (SpatialDimPlot(ST,images = c("anterior2"),cols = cols,alpha = 1,stroke = 0) + NoLegend() + 
     theme(panel.border = element_rect(fill = NA,colour = "black")))
BC1_ST_A

BC1_ST_P <- (SpatialDimPlot(ST,images = c("posterior2"),cols = cols,alpha = 0,stroke = 0) + NoLegend() + 
               ggtitle("BC1_ST_P") + theme(panel.border = element_rect(fill = NA,colour = "black"))) /
  (SpatialDimPlot(ST,images = c("posterior2"),cols = cols,alpha = 1,stroke = 0) + NoLegend() + 
     theme(panel.border = element_rect(fill = NA,colour = "black")))
BC1_ST_P

BC1_ST_A | BC1_ST_P

BM1_ST_A <- (SpatialDimPlot(ST,images = c("BM1_ST_A"),cols = cols,alpha = 0,stroke = 0) + NoLegend() + 
               ggtitle("BM1_ST_A") + theme(panel.border = element_rect(fill = NA,colour = "black"))) /
  (SpatialDimPlot(ST,images = c("BM1_ST_A"),cols = cols,alpha = 1,stroke = 0) + NoLegend() + 
     theme(panel.border = element_rect(fill = NA,colour = "black")))
BM1_ST_A

BM1_ST_P <- (SpatialDimPlot(ST,images = c("BM1_ST_P"),cols = cols,alpha = 0,stroke = 0) + NoLegend() + 
               ggtitle("BM1_ST_P") + theme(panel.border = element_rect(fill = NA,colour = "black"))) /
  (SpatialDimPlot(ST,images = c("BM1_ST_P"),cols = cols,alpha = 1,stroke = 0) + NoLegend() + 
     theme(panel.border = element_rect(fill = NA,colour = "black")))
BM1_ST_P

BM1_ST_A | BM1_ST_P

BA1_ST_A <- (SpatialDimPlot(ST,images = c("BA1_ST_A"),cols = cols,alpha = 0,stroke = 0) + NoLegend() + 
               ggtitle("BA1_ST_A") + theme(panel.border = element_rect(fill = NA,colour = "black"))) /
  (SpatialDimPlot(ST,images = c("BA1_ST_A"),cols = cols,alpha = 1,stroke = 0) + NoLegend() + 
     theme(panel.border = element_rect(fill = NA,colour = "black")))
BA1_ST_A

BA1_ST_P <- (SpatialDimPlot(ST,images = c("BA1_ST_P"),cols = cols,alpha = 0,stroke = 0) + NoLegend() + 
               ggtitle("BA1_ST_P") + theme(panel.border = element_rect(fill = NA,colour = "black"))) /
  (SpatialDimPlot(ST,images = c("BA1_ST_P"),cols = cols,alpha = 1,stroke = 0) + #NoLegend() + 
     theme(panel.border = element_rect(fill = NA,colour = "black")))
BA1_ST_P

BA1_ST_A | BA1_ST_P

BC2_ST_A <- (SpatialDimPlot(ST,images = c("anterior2"),cols = cols,alpha = 0,stroke = 0) + NoLegend() + 
  ggtitle("BC2_ST_A") + theme(panel.border = element_rect(fill = NA,colour = "black"))) /
  (SpatialDimPlot(ST,images = c("anterior2"),cols = cols,alpha = 1,stroke = 0) + NoLegend() + 
 theme(panel.border = element_rect(fill = NA,colour = "black")))
BC2_ST_A

BC2_ST_P <- (SpatialDimPlot(ST,images = c("posterior2"),cols = cols,alpha = 0,stroke = 0) + NoLegend() + 
    ggtitle("BC2_ST_P") + theme(panel.border = element_rect(fill = NA,colour = "black"))) /
  (SpatialDimPlot(ST,images = c("posterior2"),cols = cols,alpha = 1,stroke = 0) + NoLegend() + 
     theme(panel.border = element_rect(fill = NA,colour = "black")))
BC2_ST_P

BC2_ST_A | BC2_ST_P

BM2_ST_A <- (SpatialDimPlot(ST,images = c("BM2_ST_A"),cols = cols,alpha = 0,stroke = 0) + NoLegend() + 
               ggtitle("BM2_ST_A") + theme(panel.border = element_rect(fill = NA,colour = "black"))) /
  (SpatialDimPlot(ST,images = c("BM2_ST_A"),cols = cols,alpha = 1,stroke = 0) + NoLegend() + 
     theme(panel.border = element_rect(fill = NA,colour = "black")))
BM2_ST_A

BM2_ST_P <- (SpatialDimPlot(ST,images = c("BM2_ST_P"),cols = cols,alpha = 0,stroke = 0) + NoLegend() + 
               ggtitle("BM2_ST_P") + theme(panel.border = element_rect(fill = NA,colour = "black"))) /
  (SpatialDimPlot(ST,images = c("BM2_ST_P"),cols = cols,alpha = 1,stroke = 0) + NoLegend() + 
     theme(panel.border = element_rect(fill = NA,colour = "black")))
BM2_ST_P

BM2_ST_A | BM2_ST_P

BA2_ST_A <- (SpatialDimPlot(ST,images = c("BA2_ST_A"),cols = cols,alpha = 0,stroke = 0) + NoLegend() + 
               ggtitle("BA2_ST_A") + theme(panel.border = element_rect(fill = NA,colour = "black"))) /
  (SpatialDimPlot(ST,images = c("BA2_ST_A"),cols = cols,alpha = 1,stroke = 0) + NoLegend() + 
     theme(panel.border = element_rect(fill = NA,colour = "black")))
BA2_ST_A

BA2_ST_P <- (SpatialDimPlot(ST,images = c("BA2_ST_P"),cols = cols,alpha = 0,stroke = 0) + NoLegend() + 
               ggtitle("BA2_ST_P") + theme(panel.border = element_rect(fill = NA,colour = "black"))) /
  (SpatialDimPlot(ST,images = c("BA2_ST_P"),cols = cols,alpha = 1,stroke = 0) + #NoLegend() + 
     theme(panel.border = element_rect(fill = NA,colour = "black")))
BA2_ST_P

BA2_ST_A | BA2_ST_P

# (SpatialDimPlot(ST,images = c("anterior2"),alpha = 1,stroke = 0,group.by = "integrated_snn_res.0.8") + NoLegend() + 
#    theme(panel.border = element_rect(fill = NA,colour = "black"))) / 
# (SpatialDimPlot(ST,images = c("posterior2"),alpha = 1,stroke = 0,group.by = "integrated_snn_res.0.8") + NoLegend() + 
#    theme(panel.border = element_rect(fill = NA,colour = "black"))) / 

# 5.3 Markers gene (Dotplot)