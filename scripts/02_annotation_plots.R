library(shiny)
library(ggplot2)
library(scales)
library(dplyr)
library(viridis)
library(cowplot)
library(here)
library(Seurat)
library(tiff)
library(reshape2)
library(ggforce)


source(here("scripts/00_pretty_plots.R"))
#######################
## C107 Annotation Plots
#######################

annotation_C107_8um<-read.csv(here("data/annotation_C107_8um_with_banksy.csv"))


########################
## cell counts
########################

cell_counts<-annotation_C107_8um %>% 
  group_by(final_anno) %>% 
  summarise(count=length(unique(X))) %>% 
  mutate(countT= sum(count)) %>%
  group_by(final_anno, add=TRUE) %>%
  mutate(per=100*count/countT)

cell_counts$final_anno<-factor(cell_counts$final_anno, levels=c(
  "layer 1 hepatocyte","layer 2 hepatocyte","layer 3 hepatocyte","layer 4 hepatocyte",
  "layer 5 hepatocyte","layer 6 hepatocyte","layer 7 hepatocyte","layer 8 hepatocyte",
  "layer 9 hepatocyte",
  
  "T cell", "CD4-positive, alpha-beta T cell","CD8-positive, alpha-beta T cell",
  "hepatic pit cell",  "natural killer cell","neutrophil", "neutrophil2",
  "mast cell", "mature B cell" ,"plasma cell", "plasmacytoid dendritic cell",
  "Kupffer cell", "macrophage", "monocyte", "conventional dendritic cell",
  "erythrocyte",
  
  "endothelial cell of artery","endothelial cell of pericentral hepatic sinusoid",
  "endothelial cell of periportal hepatic sinusoid","vein endothelial cell" ,
  
  "fibroblast","hepatic stellate cell", "intrahepatic cholangiocyte",                     
  
  "Low_UMI", "Unannotated","unknown" ))

C107_count<-ggplot(cell_counts[which(cell_counts$final_anno!="Low_UMI"),], aes(final_anno, count, fill=final_anno))+geom_bar(stat = "identity", color="black")+
  theme_bw()+fillscale_cellType+xlab("")+ylab("8um Bin Count") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
C107_count
save_plts(C107_count, "C107_counts", w=6,h=8)


annotation_C107_8um[which(annotation_C107_8um$final_anno=="Unannotated"),]
annotation_C107_8um[which(annotation_C107_8um$final_anno=="unknown"),]

#################
## plot spatially
#################
localdir <- "/media/redgar/Seagate Portable Drive/visiumHD_liver/reseq/MacParland_Diana__C107_D1/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8))

# Prepare spot coordinates
coords <- as.data.frame(cbind(object@meta.data, object@images$slice1$centroids@coords))
coords$bin<-rownames(coords)
rm(object)
gc()

plt_annotation_spatial<-merge(coords,annotation_C107_8um, by.x="bin", by.y="bin")

bin_size <- 13.4543
half_bin <- bin_size / 2
unit_per_um <- bin_size / 8

# Choose scale bar size in μm
scale_length_um <- 1000
scale_length_plot <- scale_length_um * unit_per_um
x_start <- (-max(plt_annotation_spatial$x))
y_start <- max(plt_annotation_spatial$y)

get_leg = function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

plot_bins<-ggplot() +
  geom_point(data = plt_annotation_spatial,
             aes(x = y, y = -x, color = final_anno),
             size = 0.25) +
  # Scale bar
  annotate("segment",x = y_start-scale_length_plot-100,xend = y_start-100 ,y = x_start+ 300,yend = x_start+300,
           color = "black",linewidth = 1) +
  annotate("text",x = y_start - (scale_length_plot / 2)-100,y = x_start + 200,label = paste0(scale_length_um, " µm"),
           color = "black",size = 3,hjust = 0.5) +
  # Black border around image
  geom_rect(aes(xmin = min(plt_annotation_spatial$y), xmax = max(plt_annotation_spatial$y), ymin = (-min(plt_annotation_spatial$x)), ymax = (-max(plt_annotation_spatial$x))),
            fill = NA,
            color = "black",
            linewidth = 1) +
  coord_fixed() +
  colscale_cellType+
  theme_void()


full_tissue<-plot_grid(plot_bins+theme(legend.position = "none"), get_leg(plot_bins+guides(color = guide_legend(ncol=4,override.aes = list(size = 2)))), ncol=1, rel_heights = c(3,1))
full_tissue
save_plts(full_tissue, "C107_full_tisse", w=15, h=20)


#################
## plot spatially  spot class
#################
localdir <- "/media/redgar/Seagate Portable Drive/visiumHD_liver/reseq/MacParland_Diana__C107_D1/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8))

# Prepare spot coordinates
coords <- as.data.frame(cbind(object@meta.data, object@images$slice1$centroids@coords))
coords$bin<-rownames(coords)
rm(object)
gc()

plt_annotation_spatial<-merge(coords[,c("bin","x","y")],annotation_C107_8um, by.x="bin", by.y="bin")

bin_size <- 13.4543
half_bin <- bin_size / 2
unit_per_um <- bin_size / 8

# Choose scale bar size in μm
scale_length_um <- 1000
scale_length_plot <- scale_length_um * unit_per_um
x_start <- (-max(plt_annotation_spatial$x))
y_start <- max(plt_annotation_spatial$y)

get_leg = function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

plt_annotation_spatial$spot_class<-as.factor(plt_annotation_spatial$spot_class)
levels(plt_annotation_spatial$spot_class)<-c("Doublet\n(certain)","Doublet\n(uncertain)",
                                             "Low UMI","Reject","Singlet")
plt_annotation_spatial$spot_class<-factor(plt_annotation_spatial$spot_class, levels=c("Singlet","Doublet\n(certain)","Doublet\n(uncertain)",
                                                                                      "Low UMI","Reject"))

plot_bins<-ggplot() +
  geom_point(data = plt_annotation_spatial,
             aes(x = y, y = -x, color = spot_class),
             size = 0.25) +
  # Scale bar
  annotate("segment",x = y_start-scale_length_plot-100,xend = y_start-100 ,y = x_start+ 300,yend = x_start+300,
           color = "black",linewidth = 1) +
  annotate("text",x = y_start - (scale_length_plot / 2)-100,y = x_start + 200,label = paste0(scale_length_um, " µm"),
           color = "black",size = 3,hjust = 0.5) +
  # Black border around image
  geom_rect(aes(xmin = min(plt_annotation_spatial$y), xmax = max(plt_annotation_spatial$y), ymin = (-min(plt_annotation_spatial$x)), ymax = (-max(plt_annotation_spatial$x))),
            fill = NA,
            color = "black",
            linewidth = 1) +
  coord_fixed() +scale_color_manual(values=c("#7FA8F0","#2B3D41","#4C5F6B","#FFA630","#B8336A"))+
  theme_void()


full_tissue<-plot_grid(plot_bins+theme(legend.position = "none"), get_leg(plot_bins+guides(color = guide_legend(ncol=4,override.aes = list(size = 2)))), ncol=1, rel_heights = c(3,1))
full_tissue
save_plts(full_tissue, "C107_full_tisse_spotclass", w=15, h=20)

## bar plot of spot class
cell_counts<-annotation_C107_8um %>% 
  group_by(spot_class) %>% 
  summarise(count=length(unique(X))) %>% 
  mutate(countT= sum(count)) %>%
  group_by(spot_class, add=TRUE) %>%
  mutate(per=100*count/countT)

cell_counts$spot_class<-as.factor(cell_counts$spot_class)
levels(cell_counts$spot_class)<-c("Doublet\n(certain)","Doublet\n(uncertain)",
                                  "Low UMI","Reject","Singlet")
cell_counts$spot_class<-factor(cell_counts$spot_class, levels=c("Singlet","Doublet\n(certain)","Doublet\n(uncertain)",
                                                                "Low UMI","Reject"))


C107_count_spot_class<-ggplot(cell_counts, aes(spot_class, count, fill=spot_class))+geom_bar(stat = "identity", color="black")+
  theme_bw()+scale_fill_manual(values=c("#7FA8F0","#2B3D41","#4C5F6B","#FFA630","#B8336A"))+
  xlab("")+ylab("8um Bin Count") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
C107_count_spot_class
save_plts(C107_count_spot_class, "C107_counts_spotclass", w=3,h=4)


## add zoom area
## coordiantes
xmin<-6920
xmax<-7900
ymin<-10160
ymax<-12480


plot_bins<-ggplot() +
  geom_point(data = plt_annotation_spatial,
             aes(x = y, y = -x, color = final_anno),
             size = 0.1) +
  # Scale bar
  annotate("segment",x = y_start-scale_length_plot-100,xend = y_start-100 ,y = x_start- 200,yend = x_start-200,
           color = "black",linewidth = 1) +
  annotate("text",x = y_start - (scale_length_plot / 2)-100,y = x_start -300,label = paste0(scale_length_um, " µm"),
           color = "black",size = 5,hjust = 0.5) +
  # Black border around image
  geom_rect(aes(xmin = min(plt_annotation_spatial$y), xmax = max(plt_annotation_spatial$y), ymin = (-min(plt_annotation_spatial$x)), ymax = (-max(plt_annotation_spatial$x))),
            fill = NA,
            color = "black",
            linewidth = 1) +
  # zoom area
  geom_rect(aes(xmin = ymin, xmax = ymax, ymin = -xmin, ymax = -xmax),
            fill = NA,
            color = "black",
            linewidth = 0.75) +
  coord_fixed() +
  colscale_cellType+
  theme_void()+
  theme(legend.position = "none")

plot_bins
save_plts(plot_bins, "C107_full_tisse_zoom_box", w=15, h=15)


#######################
## zoom with outlines
#######################

## coordiantes
xmin<-6920
xmax<-7900
ymin<-10160
ymax<-12480


# load highest res image
tiff_path<-here("/media/redgar/Seagate Portable Drive/visiumHD_liver/high_res_HnE/H1-NC8X8J2_D1.tif")
tiff_res <- readTIFF(tiff_path, as.is = TRUE)
dim(tiff_res)

## from seurat import above
x_range <- xmin:xmax
y_range <- ymin:ymax



# Step 3: Subset the image
# TIFFs in R are typically in [row, col, channel] format
# Note: y is rows (height), x is columns (width)
zoomed_img <- tiff_res[x_range, y_range, , drop = FALSE]  # if it's RGB
dim(zoomed_img)

# Create RGB array to raster
rgb_img <- rgb(zoomed_img[,,1],
               zoomed_img[,,2],
               zoomed_img[,,3])

# Convert to matrix with same dimensions
dim(rgb_img) <- dim(zoomed_img)[1:2]

# Step 5: Plot
plot(1, type = "n", xlim = c(0, dim(rgb_img)[2]), ylim = c(0, dim(rgb_img)[1]), xlab = "", ylab = "", axes = FALSE)
rasterImage(rgb_img, 0, 0, dim(rgb_img)[2], dim(rgb_img)[1])

rm(tiff_res)
gc()

### bin annotation data
coords<-plt_annotation_spatial[which(plt_annotation_spatial$x>xmin & plt_annotation_spatial$x<xmax & plt_annotation_spatial$y>ymin & plt_annotation_spatial$y<ymax),]

# Get image dimensions
h <- dim(rgb_img)[1]
w <- dim(rgb_img)[2]

# Convert matrix to data frame
df <- melt(rgb_img)
colnames(df) <- c("y", "x", "color")
df$y <- h - df$y + 1  # Flip y-axis


# Flip y for coords too
coords$y_centred <- (coords$y - ymin)+1
coords$x_centred <- (xmin - (coords$x)+1)+(xmax-xmin)

# Plot image and overlay points
ggplot() +
  geom_raster(data = df, aes(x = x, y = y, fill = color), show.legend = FALSE) +
  scale_fill_identity() +
  geom_point(data = coords, aes(x = y_centred, y = x_centred, color=nCount_Spatial.008um),  size = 0.3) +
  coord_fixed() + scale_color_distiller(palette = "Spectral")

## check alignemt
plot_grid(SpatialFeaturePlot(zoom_seurat, features = "nCount_Spatial.008um", pt.size.factor = 4) +
            theme(legend.position = "right"),
          ggplot() +
            geom_raster(data = df, aes(x = x, y = y, fill = color), show.legend = FALSE) +
            scale_fill_identity() +
            geom_point(data = coords, aes(x = y_centred, y = x_centred, color=nCount_Spatial.008um),  size = 0.3) +
            coord_fixed() + scale_color_distiller(palette = "Spectral"), ncol=2 , rel_widths = c(0.5,1) )

bin_size <- 13.4543
half_bin <- bin_size / 2

coords_rect <- coords |> 
  mutate(
    xmin = y_centred - half_bin,
    xmax = y_centred + half_bin,
    ymin = x_centred - half_bin,
    ymax = x_centred + half_bin
  )

unit_per_um <- bin_size / 8

# Choose scale bar size in μm
scale_length_um <- 50
scale_length_plot <- scale_length_um * unit_per_um

x_start <- 10
y_start <- -10

image_bounds <- data.frame(
  xmin = min(df$x),
  xmax = max(df$x),
  ymin = min(df$y),
  ymax = max(df$y)
)






aling_views<-plot_grid(
  ggplot() +
    # Raster image
    geom_raster(data = df, aes(x = x, y = y, fill = color), show.legend = FALSE, alpha=1) +
    scale_fill_identity() +
    # Transparent bin rectangles
    geom_rect(data = coords_rect,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              color = alpha("white", 0.3),fill = NA,linewidth = 0.1) +
    # Scale bar
    annotate("segment", x = x_start, xend = x_start + scale_length_plot, y = y_start-10, yend = y_start-10,
             color = "black", linewidth = 0.5) +
    annotate("text", x = x_start + scale_length_plot / 2, y = y_start - 30, label = paste0(scale_length_um, " µm"),
             color = "black", size = 2, hjust = 0.5) +
    # Black border around image
    geom_rect(data = image_bounds,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = NA,color = "black",linewidth = 1) +
    # PV and CV
    annotate("text", x = 470, y = 550, label = "CV", fontface=2, color = "black", size = 5, hjust = 0.5) +
    annotate("text", x = 2020, y = 490, label = "PV", fontface=2, color = "black", size = 5, hjust = 0.5) +
    geom_ellipse(aes(x0 = 460, y0 = 540, a = 280, b = 120, angle = -45 ), color="black") +
    geom_ellipse(aes(x0 = 2025, y0 = 500, a = 70, b = 70, angle = 0 ), color="black") +
    coord_fixed() +
    colscale_cellType +
    theme_void(),
  
  ggplot() +
    # Centroid points
    geom_point(data = coords,aes(x = y_centred, y = x_centred, color = final_anno),
               size = 1, shape=15) +
    # Scale bar
    annotate("segment", x = x_start, xend = x_start + scale_length_plot, y = y_start-10, yend = y_start-10,
             color = "black", linewidth = 0.5) +
    annotate("text", x = x_start + scale_length_plot / 2, y = y_start - 30, label = paste0(scale_length_um, " µm"),
             color = "black", size = 2, hjust = 0.5) +
    # Black border around image
    geom_rect(data = image_bounds,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = NA,color = "black",linewidth = 1) +
    # PV and CV
    annotate("text", x = 470, y = 550, label = "CV", fontface=2, color = "black", size = 5, hjust = 0.5) +
    annotate("text", x = 2020, y = 490, label = "PV", fontface=2, color = "black", size = 5, hjust = 0.5) +
    geom_ellipse(aes(x0 = 460, y0 = 540, a = 280, b = 120, angle = -45 ), color="black") +
    geom_ellipse(aes(x0 = 2025, y0 = 500, a = 70, b = 70, angle = 0 ), color="black") +
    coord_fixed() +
    colscale_cellType +
    theme_void()+theme(legend.position = "none"),
  ncol=1, align="v")
aling_views

save_plts(aling_views, "C107_PV_CV_tisse", w=8, h=8)



## save cell type legend
leg_ex<-ggplot() +
  geom_point(data = coords,aes(x = y, y = x, color = final_anno),
             size = 2.5, shape=15) +
  colscale_cellType +
  theme_void()

save_plts( get_leg(leg_ex+guides(color = guide_legend(ncol=4,override.aes = list(size = 6)))), "HLiCA_cell_type_legend", w=9, h=3)


## for annotation overview
localdir <- "/media/redgar/Seagate Portable Drive/visiumHD_liver/reseq/MacParland_Diana__C107_D1/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(16))

coords <- as.data.frame(cbind(object@meta.data, object@images$slice1$centroids@coords))
coords$bin<-rownames(coords)
rm(object)
gc()


min_max_df <- coords %>%  summarize(    min_x = min(x, na.rm = TRUE),    max_x = max(x, na.rm = TRUE),    min_y = min(y, na.rm = TRUE),    max_y = max(y, na.rm = TRUE) )        


UMI_number16<-ggplot() +
  geom_point(aes(y, -x),coords, size=1, color="black")+
  geom_point(data = coords,
             aes(x = y, y = -x, color = nCount_Spatial.016um),
             size = 0.1) + 
  geom_rect(data=min_max_df, aes(xmin = min_y, xmax = max_y, ymin = -min_x, ymax = -max_x),
            fill = NA,color = "black",linewidth = 1)+
  theme_void()+ theme(strip.text = element_text(size=15))+
  scale_color_gradientn(colors = grey.colors(10), transform="sqrt", name = "Number of UMI\nper 8um bin") 
UMI_number16
save_plts(UMI_number16, "grey_for_concept", w=6, h=5)
