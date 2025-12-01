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
xmin<-4080
xmax<-4420
ymin<-12100
ymax<-12440



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
object <- readRDS("data/subset4_20_100.rds")

plt_hack<-cbind(object@meta.data, object@images$slice1$centroids@coords)
zoom<-plt_hack[which(plt_hack$x>xmin & plt_hack$x<xmax & plt_hack$y>ymin & plt_hack$y<ymax),]

object$cell<-colnames(object)
zoom_seurat <- subset(object, subset = cell %in% rownames(zoom))

# save(zoom_seurat, file="/media/redgar/Seagate Portable Drive/visiumHD_liver/zoom_2um.RData")
# rm(object)
# gc()

SpatialFeaturePlot(zoom_seurat, features = "nCount_Spatial", pt.size.factor = 1.2) +
  theme(legend.position = "right")

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
  geom_point(data = coords, aes(x = y_centred, y = x_centred, color=nCount_Spatial),  size = 0.3) +
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



#save all necessary elements
save(counts,df,coords, coords_rect, image_bounds,x_start,y_start,scale_length_um,  xmin,  xmax, ymin, ymax,scale_length_plot,
     file="/media/redgar/Seagate Portable Drive/visiumHD_liver/zoom_2um_allobjects_RCTD.RData")
