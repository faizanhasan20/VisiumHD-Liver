## Load Libraries
library(here)
library(Seurat)

library(reshape2)
library(dplyr)
library(purrr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

library(sp)
library(scales)
library(viridis)
library(colorspace)

source("scripts/00_pretty_plots.R")


localdir <- "/media/redgar/Seagate Portable Drive/visiumHD_liver/reseq/MacParland_Diana__C107_D1/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8))

# Prepare spot coordinates
coords <- as.data.frame(cbind(object@meta.data, object@images$slice1$centroids@coords))
coords$bin<-rownames(coords)

gene_exp<-FetchData(object, vars=c("SLCO1B3","LGALS4"))
gene_exp$cell<-rownames(gene_exp)
coords_exp<-cbind(coords, gene_exp)


rm(object)
gc()





#######
## layer iamge and bins in a plot
#######
coords_exp<-coords_exp[order(coords_exp$SLCO1B3),]
coords_exp$SLCO1B3[which(coords_exp$SLCO1B3==0)]<-NA

# # Plot image and overlay points
# ggplot() +
#   scale_fill_identity() +
#   geom_point(data = coords_exp, aes(x = y, y = -x, color=SLCO1B3),  size = 0.3) +
#   coord_fixed() + scale_color_viridis()


bin_size <- 13.4543/4
unit_per_um <- bin_size / 2

# Choose scale bar size in μm
scale_length_um <- 1000
scale_length_plot <- scale_length_um * unit_per_um

x_start <- max(coords_exp$y) - (scale_length_um*1.8)
y_start <- (-max(coords_exp$x))+380


SLCO1B3_exp_plot<-ggplot() +
  # bin values
  geom_point(data = coords_exp, aes(x = y, y = -x, color=SLCO1B3),  size = 0.15) +
  # Scale bar
  annotate("segment",x = x_start,xend = x_start + scale_length_plot,y = y_start,yend = y_start,
           color = "black",linewidth = 0.75) +
  annotate("text",x = x_start + scale_length_plot / 2,
           y = y_start - 150,label = paste0(scale_length_um, " µm"),
           color = "black",size = 3.5,hjust = 0.5) +
  # Black border around image
  geom_rect(aes(xmin = min(coords_exp$y), xmax = max(coords_exp$y), 
                                    ymin = -min(coords_exp$x), ymax = -max(coords_exp$x)),
            fill = NA,color = "black",linewidth = 1) +
  coord_fixed() +  scale_color_viridis(na.value = "grey90", name="Counts") +
  theme_void()
save_plts(SLCO1B3_exp_plot, "SLCO1B3_exp_plot", w=7, h=8)



coords_exp<-coords_exp[order(coords_exp$LGALS4),]
coords_exp$LGALS4[which(coords_exp$LGALS4==0)]<-NA

LGALS4_exp_plot<-ggplot() +
  # bin values
  geom_point(data = coords_exp, aes(x = y, y = -x, color=LGALS4),  size = 0.15) +
  # Scale bar
  annotate("segment",x = x_start,xend = x_start + scale_length_plot,y = y_start,yend = y_start,
           color = "black",linewidth = 0.75) +
  annotate("text",x = x_start + scale_length_plot / 2,
           y = y_start - 150,label = paste0(scale_length_um, " µm"),
           color = "black",size = 3.5,hjust = 0.5) +
  # Black border around image
  geom_rect(aes(xmin = min(coords_exp$y), xmax = max(coords_exp$y), 
                ymin = -min(coords_exp$x), ymax = -max(coords_exp$x)),
            fill = NA,color = "black",linewidth = 1) +
  coord_fixed() +  scale_color_viridis(na.value = "grey90", name="Counts") +
  theme_void()
save_plts(LGALS4_exp_plot, "LGALS4_exp_plot", w=7, h=8)

