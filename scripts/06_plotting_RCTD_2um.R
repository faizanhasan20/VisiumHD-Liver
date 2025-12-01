library(tiff)
library(here)

#remotes::install_github(repo = "satijalab/seurat", ref = "visium-hd")
#https://satijalab.org/seurat/articles/visiumhd_commands_intro
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(arrow)
library(reshape2)
library(cowplot)
library(ggforce)


source(here("scripts/00_pretty_plots.R"))


# load highest res image
tiff_path<-here("/media/redgar/Seagate Portable Drive/visiumHD_liver/high_res_HnE/H1-NC8X8J2_D1.tif")
tiff_res <- readTIFF(tiff_path, as.is = TRUE)
dim(tiff_res)


## coordiantes (tiny zoom)
      # xmin<-7400
      # xmax<-7600
      # ymin<-12115
      # ymax<-12315

## coordiantes (CV)
xmin<-7700
xmax<-8090
ymin<-7300
ymax<-7775


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




### load data from Faizan RCTD
#object <- readRDS("data/zoom_2um_C107_RCTD.rds")
object <- readRDS("data/zoomCV_2um_C107_RCTD.rds")
object$first_type<-gsub("bin ","layer ",object$first_type)

plt_hack<-cbind(object@meta.data, object@images$slice1$centroids@coords)
zoom<-plt_hack[which(plt_hack$x>xmin & plt_hack$x<xmax & plt_hack$y>ymin & plt_hack$y<ymax),]

object$cell<-colnames(object)
zoom_seurat <- subset(object, subset = cell %in% rownames(zoom))

# Prepare spot coordinates
coords <- as.data.frame(cbind(zoom_seurat@meta.data, zoom_seurat@images$slice1$centroids@coords))

# want to show gene expression
counts <- zoom_seurat@assays$Spatial$counts





#######
## layer iamge and bins in a plot
#######
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
  geom_point(data = coords, aes(x = y_centred, y = x_centred, color=nCount_Spatial.002um),  size = 2.1, shape=15) +
  coord_fixed() + scale_color_distiller(palette = "Spectral")


## calculate bin size by comparing consecutive centroid coordinates
# 13.4543/4
# coords[which(coords$y_centred>800.95),]
# 564.1025-560.7389
# 157.0807


bin_size <- 13.4543/4
half_bin <- bin_size / 2

coords_rect <- coords |> 
  mutate(
    xmin = y_centred - half_bin,
    xmax = y_centred + half_bin,
    ymin = x_centred - half_bin,
    ymax = x_centred + half_bin
  )

unit_per_um <- bin_size / 2

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


# 
# #save all necessary elements
# save(counts,df,coords, coords_rect, image_bounds,x_start,y_start,scale_length_um,  xmin,  xmax, ymin, ymax,scale_length_plot,
#      file="/media/redgar/Seagate Portable Drive/visiumHD_liver/zoom_2um_allobjects_RCTD.RData")
# 
# load("/media/redgar/Seagate Portable Drive/visiumHD_liver/zoom_2um_allobjects_RCTD.RData")
# ########## tile fill color


########################
## Gene expression plot
########################
gene<-"ALB"

coords$gene<-counts[which(rownames(counts)==gene),]

ggplot() +
  # Raster image
  geom_raster(data = df, aes(x = x, y = y, fill = color), show.legend = FALSE, alpha=1) +
  scale_fill_identity() +
  
  # Transparent bin rectangles
  geom_rect(
    data = coords_rect,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    color = alpha("white", 0.3),
    fill = NA,
    linewidth = 0.1
  ) +
  # bin values
  geom_point(data = coords,
             aes(x = y_centred, y = x_centred, color = gene),
             size = 2.1, shape=15, alpha=0.5) +
  # Scale bar
  annotate("segment",
           x = x_start,
           xend = x_start + scale_length_plot,
           y = y_start,
           yend = y_start,
           color = "black",
           linewidth = 1) +
  annotate("text",
           x = x_start + scale_length_plot / 2,
           y = y_start - 10,
           label = paste0(scale_length_um, " µm"),
           color = "black",
           size = 3,
           hjust = 0.5) +
  
  # Black border around image
  geom_rect(data = image_bounds,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA,
            color = "black",
            linewidth = 1) +
  
  coord_fixed() +
  scale_colour_viridis_c(option = "C", direction = -1) +
  theme_void()



########################
## "first_type" and "spot_class" 
########################

get_leg = function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

plot_bins<-ggplot() +
  # Raster image
  geom_raster(data = df, aes(x = x, y = y, fill = color), show.legend = FALSE, alpha=1) +
  scale_fill_identity() +
  
  # Transparent bin rectangles
  geom_rect(
    data = coords_rect,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    color = alpha("white", 0.3),
    fill = NA,
    linewidth = 0.1
  ) +
  # bin values
  geom_point(data = coords,
             aes(x = y_centred, y = x_centred, color = first_type),
             size = 2.1, shape=15, alpha=1) +
  # Scale bar
  annotate("segment",
           x = x_start,
           xend = x_start + scale_length_plot,
           y = y_start,
           yend = y_start,
           color = "black",
           linewidth = 1) +
  annotate("text",
           x = x_start + scale_length_plot / 2,
           y = y_start - 10,
           label = paste0(scale_length_um, " µm"),
           color = "black",
           size = 3,
           hjust = 0.5) +
  
  # Black border around image
  geom_rect(data = image_bounds,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA,
            color = "black",
            linewidth = 1) +
  coord_fixed() +
  colscale_cellType+
  theme_void()

plot_grid(plot_bins+theme(legend.position = "none"), get_leg(plot_bins), ncol=1, rel_heights = c(3,1))


unique(coords$first_type)[which(!(unique(coords$first_type)%in%color_possibilities_celltype))]



plot_bins<-ggplot() +
  geom_raster(data = df, aes(x = x, y = y, fill = color), show.legend = FALSE, alpha=1) +
  scale_fill_identity() +
  geom_rect(data = coords_rect,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),color = alpha("white", 0.3),fill = NA, linewidth = 0.1) +
  geom_point(data = coords[which(coords$first_type%in%c("endothelial cell of periportal hepatic sinusoid","vein endothelial cell","endothelial cell of artery","endothelial cell of pericentral hepatic sinusoid")),],
             aes(x = y_centred, y = x_centred, color = first_type),size = 1.5, shape=15, alpha=1) +
  annotate("segment",x = x_start,xend = x_start + scale_length_plot,y = y_start,yend = y_start,color = "black",linewidth = 1) +
  annotate("text",x = x_start + scale_length_plot / 2,y = y_start - 10,label = paste0(scale_length_um, " µm"),color = "black",size = 3,hjust = 0.5) +
  geom_rect(data = image_bounds,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),fill = NA,color = "black",linewidth = 1) +
  coord_fixed() +  colscale_cellType+  theme_void()

plot_grid(plot_bins+theme(legend.position = "none"), get_leg(plot_bins), ncol=1, rel_heights = c(3,1))



#######################
### 8 bin annotation data
#######################
annotation_C107_8um<-read.csv(here("data/metadata_for_C107.csv"))
annotation_C107_8um[which(annotation_C107_8um$final_anno=="Unannotated"),]
annotation_C107_8um[which(annotation_C107_8um$final_anno=="unknown"),]
annotation_C107_8um$final_anno<-gsub("bin ","layer ",annotation_C107_8um$final_anno)


localdir <- "/media/redgar/Seagate Portable Drive/visiumHD_liver/reseq/MacParland_Diana__C107_D1/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8))

# Prepare spot coordinates
coords_8um <- as.data.frame(cbind(object@meta.data, object@images$slice1$centroids@coords))
coords_8um$bin<-rownames(coords_8um)
rm(object)
gc()

plt_annotation_spatial<-merge(coords_8um,annotation_C107_8um, by.x="bin", by.y="X")

coords_8um<-plt_annotation_spatial[which(plt_annotation_spatial$x>xmin & plt_annotation_spatial$x<xmax & plt_annotation_spatial$y>ymin & plt_annotation_spatial$y<ymax),]
coords_8um$y_centred <- (coords_8um$y - ymin)+1
coords_8um$x_centred <- (xmin - (coords_8um$x)+1)+(xmax-xmin)


## save cell type legend
leg_ex<-ggplot() +
  geom_point(aes(x = c(color_possibilities_celltype,NA), y = 1, color =  c(color_possibilities_celltype,NA)),
             size = 2.5, shape=15) +
  colscale_cellType +
  theme_void()

aling_views<-plot_grid(
  ggplot() +
    # Raster image
    geom_raster(data = df, aes(x = x, y = y, fill = color), show.legend = FALSE, alpha=1) +
    scale_fill_identity() +
    # Transparent bin rectangles
    #geom_rect(data = coords_rect,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),color = alpha("white", 0.3),fill = NA,linewidth = 0.1) +
    # Scale bar
    annotate("segment", x = x_start, xend = x_start + scale_length_plot, y = y_start-10, yend = y_start-10,color = "black", linewidth = 0.5) +
    annotate("text", x = x_start + scale_length_plot / 2, y = y_start - 30, label = paste0(scale_length_um, " µm"),color = "black", size = 2, hjust = 0.5) +
    # Black border around image
    geom_rect(data = image_bounds,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),fill = NA,color = "black",linewidth = 1) +
    # CV annotation
    annotate("text", x = 305, y = 190, label = "CV", fontface=2, color = "black", size = 4, hjust = 0.5) +
    geom_ellipse(aes(x0 = 305, y0 = 187, a = 40, b = 25, angle = 45 ), color="black") +
    coord_fixed() +    colscale_cellType +    theme_void(),
  
  ggplot() +
    # Centroid points
    geom_point(data = coords_8um,aes(x = y_centred, y = x_centred, color = final_anno),size = 3, shape=15) +
    # Scale bar
    annotate("segment", x = x_start, xend = x_start + scale_length_plot, y = y_start-10, yend = y_start-10,color = "black", linewidth = 0.5) +
    annotate("text", x = x_start + scale_length_plot / 2, y = y_start - 30, label = paste0(scale_length_um, " µm"),color = "black", size = 2, hjust = 0.5) +
    # Black border around image
    geom_rect(data = image_bounds,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),fill = NA,color = "black",linewidth = 1) +
    # CV annotation
    annotate("text", x = 305, y = 190, label = "CV", fontface=2, color = "black", size = 4, hjust = 0.5) +
    geom_ellipse(aes(x0 = 305, y0 = 187, a = 40, b = 25, angle = 45 ), color="black") +
    coord_fixed() +    colscale_cellType +    theme_void()+theme(legend.position = "none"),
  
  ggplot() +
    # Centroid points
    geom_point(data = coords,aes(x = y_centred, y = x_centred, color = first_type),size = 0.3, shape=15) +
    # Scale bar
    annotate("segment", x = x_start, xend = x_start + scale_length_plot, y = y_start-10, yend = y_start-10,color = "black", linewidth = 0.5) +
    annotate("text", x = x_start + scale_length_plot / 2, y = y_start - 30, label = paste0(scale_length_um, " µm"),color = "black", size = 2.5, hjust = 0.5) +
    # Black border around image
    geom_rect(data = image_bounds,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),fill = NA,color = "black",linewidth = 1) +
    # CV annotation
    annotate("text", x = 305, y = 190, label = "CV", fontface=2, color = "black", size = 4, hjust = 0.5) +
    geom_ellipse(aes(x0 = 305, y0 = 187, a = 40, b = 25, angle = 45 ), color="black") +
    coord_fixed() +    colscale_cellType +    theme_void()+theme(legend.position = "none"),
  ncol=1, align="v")

aling_views

save_plts(plot_grid(aling_views, get_leg(leg_ex+guides(color = guide_legend(ncol=2,override.aes = list(size = 2)))),
                    ncol=2, rel_widths  = c(1,1)), "C107_CV_2um_8um_tisse", w=12, h=8)

#### big panel 3 option


# load highest res image
tiff_path<-here("/media/redgar/Seagate Portable Drive/visiumHD_liver/high_res_HnE/H1-NC8X8J2_D1.tif")
tiff_res <- readTIFF(tiff_path, as.is = TRUE)
dim(tiff_res)

## coordiantes (CV)
xmin<-7800
xmax<-8000
ymin<-7500
ymax<-7700


## from seurat import above
x_range <- xmin:xmax
y_range <- ymin:ymax


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




### load data RCTD
object <- readRDS("data/zoomCV_2um_C107_RCTD.rds")
object$first_type<-gsub("bin ","layer ",object$first_type)

plt_hack<-cbind(object@meta.data, object@images$slice1$centroids@coords)
zoom<-plt_hack[which(plt_hack$x>xmin & plt_hack$x<xmax & plt_hack$y>ymin & plt_hack$y<ymax),]

object$cell<-colnames(object)
zoom_seurat <- subset(object, subset = cell %in% rownames(zoom))

# Prepare spot coordinates
coords <- as.data.frame(cbind(zoom_seurat@meta.data, zoom_seurat@images$slice1$centroids@coords))

# want to show gene expression
counts <- zoom_seurat@assays$Spatial$counts

#######
## layer iamge and bins in a plot
#######
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

bin_size <- 13.4543/4
half_bin <- bin_size / 2

coords_rect <- coords |> 
  mutate(
    xmin = y_centred - half_bin,
    xmax = y_centred + half_bin,
    ymin = x_centred - half_bin,
    ymax = x_centred + half_bin
  )

unit_per_um <- bin_size / 2

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

#######################
### 8 bin annotation data
#######################
annotation_C107_8um<-read.csv(here("data/metadata_for_C107.csv"))
annotation_C107_8um[which(annotation_C107_8um$final_anno=="Unannotated"),]
annotation_C107_8um[which(annotation_C107_8um$final_anno=="unknown"),]
annotation_C107_8um$final_anno<-gsub("bin ","layer ",annotation_C107_8um$final_anno)


localdir <- "/media/redgar/Seagate Portable Drive/visiumHD_liver/reseq/MacParland_Diana__C107_D1/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8))

# Prepare spot coordinates
coords_8um <- as.data.frame(cbind(object@meta.data, object@images$slice1$centroids@coords))
coords_8um$bin<-rownames(coords_8um)
rm(object)
gc()

plt_annotation_spatial<-merge(coords_8um,annotation_C107_8um, by.x="bin", by.y="X")

coords_8um<-plt_annotation_spatial[which(plt_annotation_spatial$x>xmin & plt_annotation_spatial$x<xmax & plt_annotation_spatial$y>ymin & plt_annotation_spatial$y<ymax),]
coords_8um$y_centred <- (coords_8um$y - ymin)+1
coords_8um$x_centred <- (xmin - (coords_8um$x)+1)+(xmax-xmin)

## save cell type legend
leg_ex<-ggplot() +  geom_point(aes(x = c(color_possibilities_celltype,NA), y = 1, color =  c(color_possibilities_celltype,NA)),size = 2.5, shape=15) +  colscale_cellType + theme_void()




aling_views_more_zoom<-plot_grid(
  ggplot() +
    # Raster image
    geom_raster(data = df, aes(x = x, y = y, fill = color), show.legend = FALSE, alpha=1) +
    scale_fill_identity() +
    # Transparent bin rectangles
    #geom_rect(data = coords_rect,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),color = alpha("white", 0.3),fill = NA,linewidth = 0.1) +
    # Scale bar
    annotate("segment", x = x_start, xend = x_start + scale_length_plot, y = y_start-10, yend = y_start-10,color = "black", linewidth = 0.5) +
    annotate("text", x = x_start + scale_length_plot / 2, y = y_start - 30, label = paste0(scale_length_um, " µm"),color = "black", size = 2, hjust = 0.5) +
    # Black border around image
    geom_rect(data = image_bounds,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),fill = NA,color = "black",linewidth = 1) +
    # CV annotation
    annotate("text", x = 107, y = 97, label = "CV", fontface=2, color = "black", size = 5, hjust = 0.5) +
    geom_ellipse(aes(x0 = 107, y0 = 97, a = 40, b = 22, angle = 45 ), color="black") +
    coord_fixed() +    colscale_cellType +    theme_void(),
  
  ggplot() +
    # Centroid points
    geom_point(data = coords_8um,aes(x = y_centred, y = x_centred, color = final_anno),size = 3, shape=15) +
    # Scale bar
    annotate("segment", x = x_start, xend = x_start + scale_length_plot, y = y_start-10, yend = y_start-10,color = "black", linewidth = 0.5) +
    annotate("text", x = x_start + scale_length_plot / 2, y = y_start - 30, label = paste0(scale_length_um, " µm"),color = "black", size = 2, hjust = 0.5) +
    # Black border around image
    geom_rect(data = image_bounds,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),fill = NA,color = "black",linewidth = 1) +
    # CV annotation
    annotate("text", x = 107, y = 97, label = "CV", fontface=2, color = "black", size = 5, hjust = 0.5) +
    geom_ellipse(aes(x0 = 107, y0 = 97, a = 40, b = 22, angle = 45 ), color="black") +
    coord_fixed() +    colscale_cellType +    theme_void()+theme(legend.position = "none"),
  
  ggplot() +
    # Centroid points
    geom_point(data = coords,aes(x = y_centred, y = x_centred, color = first_type),size = 0.3, shape=15) +
    # Scale bar
    annotate("segment", x = x_start, xend = x_start + scale_length_plot, y = y_start-10, yend = y_start-10,color = "black", linewidth = 0.5) +
    annotate("text", x = x_start + scale_length_plot / 2, y = y_start - 30, label = paste0(scale_length_um, " µm"),color = "black", size = 2.5, hjust = 0.5) +
    # Black border around image
    geom_rect(data = image_bounds,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),fill = NA,color = "black",linewidth = 1) +
    # CV annotation
    annotate("text", x = 107, y = 97, label = "CV", fontface=2, color = "black", size = 5, hjust = 0.5) +
    geom_ellipse(aes(x0 = 107, y0 = 97, a = 40, b = 22, angle = 45 ), color="black") +
    coord_fixed() +    colscale_cellType +    theme_void()+theme(legend.position = "none"),
  ncol=1, align="v")

aling_views_more_zoom

save_plts(plot_grid(aling_views_more_zoom, get_leg(leg_ex+guides(color = guide_legend(ncol=2,override.aes = list(size = 2)))),
                    ncol=2, rel_widths  = c(1,1)), "C107_CV_2um_8um_tisse_morezoom", w=12, h=8)




##################
### more zoom
##################
tiff_path<-here("/media/redgar/Seagate Portable Drive/visiumHD_liver/high_res_HnE/H1-NC8X8J2_D1.tif")
tiff_res <- readTIFF(tiff_path, as.is = TRUE)
dim(tiff_res)

## coordiantes (CV)
xmin<-7800
xmax<-8000
ymin<-7500
ymax<-7700


## from seurat import above
x_range <- xmin:xmax
y_range <- ymin:ymax


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


object <- readRDS("data/zoomCV_2um_C107_RCTD.rds")
object$first_type<-gsub("bin ","layer ",object$first_type)

plt_hack<-cbind(object@meta.data, object@images$slice1$centroids@coords)
zoom<-plt_hack[which(plt_hack$x>xmin & plt_hack$x<xmax & plt_hack$y>ymin & plt_hack$y<ymax),]
object$cell<-colnames(object)
zoom_seurat <- subset(object, subset = cell %in% rownames(zoom))
coords <- as.data.frame(cbind(zoom_seurat@meta.data, zoom_seurat@images$slice1$centroids@coords))
counts <- zoom_seurat@assays$Spatial$counts


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

bin_size <- 13.4543/4
half_bin <- bin_size / 2

coords_rect <- coords |> 
  mutate(
    xmin = y_centred - half_bin,
    xmax = y_centred + half_bin,
    ymin = x_centred - half_bin,
    ymax = x_centred + half_bin
  )

unit_per_um <- bin_size / 2

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


annotation_C107_8um<-read.csv(here("data/metadata_for_C107.csv"))
annotation_C107_8um[which(annotation_C107_8um$final_anno=="Unannotated"),]
annotation_C107_8um[which(annotation_C107_8um$final_anno=="unknown"),]
annotation_C107_8um$final_anno<-gsub("bin ","layer ",annotation_C107_8um$final_anno)


localdir <- "/media/redgar/Seagate Portable Drive/visiumHD_liver/reseq/MacParland_Diana__C107_D1/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8))

# Prepare spot coordinates
coords_8um <- as.data.frame(cbind(object@meta.data, object@images$slice1$centroids@coords))
coords_8um$bin<-rownames(coords_8um)
rm(object)
gc()

plt_annotation_spatial<-merge(coords_8um,annotation_C107_8um, by.x="bin", by.y="X")

coords_8um<-plt_annotation_spatial[which(plt_annotation_spatial$x>xmin & plt_annotation_spatial$x<xmax & plt_annotation_spatial$y>ymin & plt_annotation_spatial$y<ymax),]
coords_8um$y_centred <- (coords_8um$y - ymin)+1
coords_8um$x_centred <- (xmin - (coords_8um$x)+1)+(xmax-xmin)

## save cell type legend
leg_ex<-ggplot() +  geom_point(aes(x = c(color_possibilities_celltype,NA), y = 1, color =  c(color_possibilities_celltype,NA)),size = 2.5, shape=15) +  colscale_cellType + theme_void()

aling_views_more_zoom<-plot_grid(
  ggplot() +
    # Raster image
    geom_raster(data = df, aes(x = x, y = y, fill = color), show.legend = FALSE, alpha=1) +
    scale_fill_identity() +
    # Transparent bin rectangles
    #geom_rect(data = coords_rect,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),color = alpha("white", 0.3),fill = NA,linewidth = 0.1) +
    # Scale bar
    annotate("segment", x = x_start, xend = x_start + scale_length_plot, y = y_start-10, yend = y_start-10,color = "black", linewidth = 0.5) +
    annotate("text", x = x_start + scale_length_plot / 2, y = y_start - 30, label = paste0(scale_length_um, " µm"),color = "black", size = 2.5, hjust = 0.5) +
    # Black border around image
    geom_rect(data = image_bounds,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),fill = NA,color = "black",linewidth = 1) +
    # CV annotation
    annotate("text", x = 107, y = 97, label = "CV", fontface=2, color = "black", size = 5, hjust = 0.5) +
    geom_ellipse(aes(x0 = 107, y0 = 97, a = 40, b = 22, angle = 45 ), color="black") +
    coord_fixed() +    colscale_cellType +    theme_void(),
  
  ggplot() +
    # Centroid points
    geom_point(data = coords_8um,aes(x = y_centred, y = x_centred, color = final_anno),size = 4.5, shape=15) +
    # Scale bar
    annotate("segment", x = x_start, xend = x_start + scale_length_plot, y = y_start-10, yend = y_start-10,color = "black", linewidth = 0.5) +
    annotate("text", x = x_start + scale_length_plot / 2, y = y_start - 30, label = paste0(scale_length_um, " µm"),color = "black", size = 2.5, hjust = 0.5) +
    # Black border around image
    geom_rect(data = image_bounds,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),fill = NA,color = "black",linewidth = 1) +
    # CV annotation
    annotate("text", x = 107, y = 97, label = "CV", fontface=2, color = "black", size = 5, hjust = 0.5) +
    geom_ellipse(aes(x0 = 107, y0 = 97, a = 40, b = 22, angle = 45 ), color="black") +
    coord_fixed() +    colscale_cellType +    theme_void()+theme(legend.position = "none"),
  
  ggplot() +
    # Centroid points
    geom_point(data = coords,aes(x = y_centred, y = x_centred, color = first_type),size = 0.8, shape=15) +
    # Scale bar
    annotate("segment", x = x_start, xend = x_start + scale_length_plot, y = y_start-10, yend = y_start-10,color = "black", linewidth = 0.5) +
    annotate("text", x = x_start + scale_length_plot / 2, y = y_start - 30, label = paste0(scale_length_um, " µm"),color = "black", size = 2.5, hjust = 0.5) +
    # Black border around image
    geom_rect(data = image_bounds,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),fill = NA,color = "black",linewidth = 1) +
    # CV annotation
    annotate("text", x = 107, y = 97, label = "CV", fontface=2, color = "black", size = 5, hjust = 0.5) +
    geom_ellipse(aes(x0 = 107, y0 = 97, a = 40, b = 22, angle = 45 ), color="black") +
    coord_fixed() +    colscale_cellType +    theme_void()+theme(legend.position = "none"),
  ncol=1, align="v")

aling_views_more_zoom

save_plts(plot_grid(aling_views_more_zoom, get_leg(leg_ex+guides(color = guide_legend(ncol=2,override.aes = list(size = 2)))),
                    ncol=2, rel_widths  = c(1,1)), "C107_CV_2um_8um_tisse_morezoom", w=12, h=8)


plot_bins<-ggplot() +
  scale_fill_identity() +
  geom_rect(data = coords_rect,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),color = alpha("white", 0.3),fill = NA, linewidth = 0.1) +
  geom_point(data = coords,
             aes(x = y_centred, y = x_centred, color = first_type),size = 1.75, shape=15, alpha=0.5) +
  geom_point(data = coords[which(coords$first_type%in%c("endothelial cell of periportal hepatic sinusoid","vein endothelial cell","endothelial cell of artery","endothelial cell of pericentral hepatic sinusoid")),],
             aes(x = y_centred, y = x_centred),color="black",size = 1.75, shape=15, alpha=1) +
  geom_point(data = coords[which(coords$first_type%in%c("endothelial cell of periportal hepatic sinusoid","vein endothelial cell","endothelial cell of artery","endothelial cell of pericentral hepatic sinusoid")),],
             aes(x = y_centred, y = x_centred, color = first_type),size = 1.15, shape=15, alpha=1) +
  annotate("segment",x = x_start,xend = x_start + scale_length_plot,y = y_start,yend = y_start,color = "black",linewidth = 1) +
  annotate("text",x = x_start + scale_length_plot / 2,y = y_start - 10,label = paste0(scale_length_um, " µm"),color = "black",size = 3,hjust = 0.5) +
  geom_rect(data = image_bounds,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),fill = NA,color = "black",linewidth = 1) +
  coord_fixed() +  colscale_cellType+  theme_void()

plot_bins

save_plts(plot_grid(plot_bins+theme(legend.position = "none"), get_leg(plot_bins), ncol=1, rel_heights = c(3,1)),
          "C107_CV_2um_8um_tisse_morezoom_endo_highlight", w=6, h=6)


plot_bins<-ggplot() +
  geom_raster(data = df, aes(x = x, y = y, fill = color), show.legend = FALSE, alpha=1) +
  scale_fill_identity() +
  geom_rect(data = coords_rect,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),color = alpha("white", 0.3),fill = NA, linewidth = 0.1) +
  geom_point(data = coords[which(coords$first_type%in%c("endothelial cell of periportal hepatic sinusoid","vein endothelial cell","endothelial cell of artery","endothelial cell of pericentral hepatic sinusoid")),],
             aes(x = y_centred, y = x_centred),color="black",size = 1.75, shape=15, alpha=1) +
  geom_point(data = coords[which(coords$first_type%in%c("endothelial cell of periportal hepatic sinusoid","vein endothelial cell","endothelial cell of artery","endothelial cell of pericentral hepatic sinusoid")),],
             aes(x = y_centred, y = x_centred, color = first_type),size = 1.15, shape=15, alpha=1) +
  annotate("segment",x = x_start,xend = x_start + scale_length_plot,y = y_start,yend = y_start,color = "black",linewidth = 1) +
  annotate("text",x = x_start + scale_length_plot / 2,y = y_start - 10,label = paste0(scale_length_um, " µm"),color = "black",size = 3,hjust = 0.5) +
  geom_rect(data = image_bounds,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),fill = NA,color = "black",linewidth = 1) +
  coord_fixed() +  colscale_cellType+  theme_void()

save_plts(plot_grid(plot_bins+theme(legend.position = "none"), get_leg(plot_bins), ncol=1, rel_heights = c(3,1)),
          "C107_CV_2um_8um_tisse_morezoom_endo_only", w=6, h=6)




aling_views_more_zoom<-plot_grid(
  ggplot() +
    # Raster image
    geom_raster(data = df, aes(x = x, y = y, fill = color), show.legend = FALSE, alpha=1) +
    scale_fill_identity() +
    annotate("segment",x = x_start,xend = x_start + scale_length_plot,y = y_start,yend = y_start,color = "black",linewidth = 1) +
    annotate("text",x = x_start + scale_length_plot / 2,y = y_start - 10,label = paste0(scale_length_um, " µm"),color = "black",size = 3,hjust = 0.5) +    geom_rect(data = image_bounds,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),fill = NA,color = "black",linewidth = 1) +
    coord_fixed() +    colscale_cellType +    theme_void(),
  
  ggplot() +
    geom_raster(data = df, aes(x = x, y = y, fill = color), show.legend = FALSE, alpha=1) +
    scale_fill_identity() +
    geom_point(data = coords_8um[which(coords_8um$final_anno%in%c("endothelial cell of periportal hepatic sinusoid","vein endothelial cell","endothelial cell of artery","endothelial cell of pericentral hepatic sinusoid")),],
               aes(x = y_centred, y = x_centred), color = "black",size = 5.5, shape=15)+
    geom_point(data = coords_8um[which(coords_8um$final_anno%in%c("endothelial cell of periportal hepatic sinusoid","vein endothelial cell","endothelial cell of artery","endothelial cell of pericentral hepatic sinusoid")),],
               aes(x = y_centred, y = x_centred, color = final_anno),size = 4.5, shape=15) +
    annotate("segment",x = x_start,xend = x_start + scale_length_plot,y = y_start,yend = y_start,color = "black",linewidth = 1) +
    annotate("text",x = x_start + scale_length_plot / 2,y = y_start - 10,label = paste0(scale_length_um, " µm"),color = "black",size = 3,hjust = 0.5) +
    geom_rect(data = image_bounds,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),fill = NA,color = "black",linewidth = 1) +
    coord_fixed() +    colscale_cellType +    theme_void()+theme(legend.position = "none"),
  
  ggplot() +
    geom_raster(data = df, aes(x = x, y = y, fill = color), show.legend = FALSE, alpha=1) +
    scale_fill_identity() +
    geom_point(data = coords[which(coords$first_type%in%c("endothelial cell of periportal hepatic sinusoid","vein endothelial cell","endothelial cell of artery","endothelial cell of pericentral hepatic sinusoid")),],
               aes(x = y_centred, y = x_centred),color="black",size = 1.75, shape=15, alpha=1) +
    geom_point(data = coords[which(coords$first_type%in%c("endothelial cell of periportal hepatic sinusoid","vein endothelial cell","endothelial cell of artery","endothelial cell of pericentral hepatic sinusoid")),],
               aes(x = y_centred, y = x_centred, color = first_type),size = 1.15, shape=15, alpha=1) +
    annotate("segment",x = x_start,xend = x_start + scale_length_plot,y = y_start,yend = y_start,color = "black",linewidth = 1) +
    annotate("text",x = x_start + scale_length_plot / 2,y = y_start - 10,label = paste0(scale_length_um, " µm"),color = "black",size = 3,hjust = 0.5) +
    geom_rect(data = image_bounds,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),fill = NA,color = "black",linewidth = 1) +
    coord_fixed() +  colscale_cellType+  theme_void(),
  ncol=1, align="v")

aling_views_more_zoom

save_plts(aling_views_more_zoom,  "C107_CV_2um_8um_tisse_morezoom_endothelial_only", w=10, h=10)



######################
## distance from vein
######################
ggplot() +
  geom_raster(data = df, aes(x = x, y = y, fill = color), show.legend = FALSE, alpha=1) +
  scale_fill_identity() +
  geom_point(data = coords[which(coords$first_type%in%c("endothelial cell of periportal hepatic sinusoid","vein endothelial cell","endothelial cell of artery","endothelial cell of pericentral hepatic sinusoid")),],
             aes(x = y_centred, y = x_centred),color="black",size = 1.75, shape=15, alpha=1) +
  geom_point(data = coords[which(coords$first_type%in%c("endothelial cell of periportal hepatic sinusoid","vein endothelial cell","endothelial cell of artery","endothelial cell of pericentral hepatic sinusoid")),],
             aes(x = y_centred, y = x_centred, color = first_type),size = 1.15, shape=15, alpha=1) +
  annotate("segment",x = x_start,xend = x_start + scale_length_plot,y = y_start,yend = y_start,color = "black",linewidth = 1) +
  annotate("text",x = x_start + scale_length_plot / 2,y = y_start - 10,label = paste0(scale_length_um, " µm"),color = "black",size = 3,hjust = 0.5) +
  geom_rect(data = image_bounds,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),fill = NA,color = "black",linewidth = 1) +
  geom_ellipse(aes(x0 = 107, y0 = 98, a = 39, b = 25, angle = 45 ), color="black") +
  coord_fixed() +  colscale_cellType+  theme_void()


# eucledian distance to vein eclipse
#Ellipse parameters
x0 <- 107
y0 <- 98
a <- 39
b <- 25
angle <- 45 * pi / 180  # convert to radians

coords_dist <- coords %>%
  filter(first_type %in% c("endothelial cell of periportal hepatic sinusoid",
                           "vein endothelial cell",
                           "endothelial cell of artery",
                           "endothelial cell of pericentral hepatic sinusoid")) %>%
  mutate(
    # translate to ellipse center
    x_t = y_centred - x0,
    y_t = x_centred - y0,
    # rotate by -angle
    x_r =  x_t * cos(-angle) - y_t * sin(-angle),
    y_r =  x_t * sin(-angle) + y_t * cos(-angle),
    # distance to ellipse boundary (approx)
    r_actual = sqrt(x_r^2 + y_r^2),
    r_ellipse = sqrt(1 / ((x_r^2 / a^2) + (y_r^2 / b^2))),
    dist_to_ellipse = abs(r_actual - r_ellipse),
    # compute closest boundary point coordinates (approx)
    scale_factor = r_ellipse / r_actual,
    x_ellipse = x0 + (x_t * scale_factor * cos(angle) - y_t * scale_factor * sin(angle)),
    y_ellipse = y0 + (x_t * scale_factor * sin(angle) + y_t * scale_factor * cos(angle))
  )

ggplot() +
  geom_raster(data = df, aes(x = x, y = y, fill = color), show.legend = FALSE, alpha = 1) +
  scale_fill_identity() +
  geom_segment(
    data = coords_dist,
    aes(x = y_centred, y = x_centred, xend = x_ellipse, yend = y_ellipse),
    linewidth = 0.5, alpha = 0.8
  ) +
  geom_point(
    data = coords_dist,
    aes(x = y_centred, y = x_centred, color = dist_to_ellipse),
    size = 1.5, shape = 15
  ) +
  geom_ellipse(aes(x0 = x0, y0 = y0, a = a, b = b, angle = 45), color = "black") +
  coord_fixed() +
  theme_void()+
  scale_color_viridis_c(option = "plasma")

ggplot(coords_dist, aes(first_type, dist_to_ellipse))+geom_violin()+geom_point()


counts <- zoom_seurat@assays$Spatial$counts
counts_endo<-as.data.frame(counts[,which(colnames(counts) %in% coords_dist$cell)])

(total_counts[rev(order(total_counts))][60:100])

names(total_counts[rev(order(total_counts))][1:60])

gene<-"CXCL2"
counts_endo_gene<-counts_endo[which(rownames(counts_endo)==gene),]
identical(colnames(counts_endo_gene), coords_dist$cell)
counts_endo_gene<-melt(counts_endo_gene)
coords_dist_gene<-merge(coords_dist, counts_endo_gene, by.x="cell", by.y="variable")

ggplot(coords_dist_gene, aes(dist_to_ellipse, value))+geom_point()+stat_smooth(method="lm")
