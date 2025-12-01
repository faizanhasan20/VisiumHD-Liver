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

source(here("scripts/00_pretty_plots.R"))


# load highest res image
tiff_path<-here("/media/redgar/Seagate Portable Drive/visiumHD_liver/high_res_HnE/H1-NC8X8J2_D1.tif")
tiff_res <- readTIFF(tiff_path, as.is = TRUE)
dim(tiff_res)


## coordiantes
xmin<-7050
xmax<-7650
ymin<-10300
ymax<-12400




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



# 8um object
annotation_C107_8um<-read.csv(here("data/metadata_for_C107.csv"))

localdir <- "/media/redgar/Seagate Portable Drive/visiumHD_liver/reseq/MacParland_Diana__C107_D1/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8))

identical(annotation_C107_8um$X, colnames(object))
rownames(annotation_C107_8um) <- annotation_C107_8um$X
annotation_C107_8um$X<-NULL

object<-AddMetaData(object, annotation_C107_8um)

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


bin_size <- 13.4543
unit_per_um <- bin_size / 8

# Choose scale bar size in μm
scale_length_um <- 100
scale_length_plot <- scale_length_um * unit_per_um

x_start <- 10
y_start <- -10

image_bounds <- data.frame(  xmin = min(df$x),  xmax = max(df$x),  ymin = min(df$y),  ymax = max(df$y))



#save all necessary elements
coords$final_anno<-gsub("bin ","layer ",coords$final_anno)

save(counts,df,coords, image_bounds,x_start,y_start,scale_length_um,  xmin,  xmax, ymin, ymax,scale_length_plot,
     file="/media/redgar/Seagate Portable Drive/visiumHD_liver/zoom_8um_allobjects_RCTD.RData", compress = "xz")

load("/media/redgar/Seagate Portable Drive/visiumHD_liver/zoom_8um_allobjects_RCTD.RData")
