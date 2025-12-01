
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

library(igraph)
library(tiff)
library(tidyr)



source("scripts/00_pretty_plots.R")



##################
### add image
##################
tiff_path<-here("/media/redgar/Seagate Portable Drive/visiumHD_liver/high_res_HnE/H1-NC8X8J2_D1.tif")
tiff_res <- readTIFF(tiff_path, as.is = TRUE)
dim(tiff_res)

## coordiantes (CV)
# xmin<-4000
# xmax<-4800
# ymin<-6700
# ymax<-7700

xmin<-5400
xmax<-6800
ymin<-7750
ymax<-9000


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


# Get image dimensions
h <- dim(rgb_img)[1]
w <- dim(rgb_img)[2]

# Convert matrix to data frame
df <- melt(rgb_img)
colnames(df) <- c("y", "x", "color")
df$y <- h - df$y + 1  # Flip y-axis




###################################################
#Large duct vs small duct cholangiocytes differentiation by counting cholangiocytes
###################################################
annotation_C107_8um<-read.csv(here("data/metadata_for_C107.csv"))

localdir <- "/media/redgar/Seagate Portable Drive/visiumHD_liver/reseq/MacParland_Diana__C107_D1/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8))

head(object)
object<-AddMetaData(object,annotation_C107_8um)

Idents(object)<-"final_anno"
object<-JoinLayers(object)

# Prepare spot coordinates
coords <- as.data.frame(cbind(object@meta.data, object@images$slice1$centroids@coords))
coords$bin<-rownames(coords)


# Step 1: Subset the Seurat cholangiocyte_bins to keep only cholangiocyte bins
cholangiocyte_bins <- subset(object, subset = final_anno == "intrahepatic cholangiocyte")

# Step 2: Extract spatial coordinates
chol_coords <- as.data.frame(cholangiocyte_bins@images$slice1$centroids@coords)
rownames(chol_coords)<-colnames(cholangiocyte_bins)

# Step 3: Compute pairwise Euclidean distance matrix
dist_matrix <- as.matrix(dist(chol_coords))

# Step 4: Define adjacency based on a distance threshold (adjust as needed)
distance_threshold <- 20  # Adjust based on your dataset scale
adj_matrix <- dist_matrix < distance_threshold  # Create adjacency matrix (TRUE for close neighbors)

# Ensure diagonal is FALSE (a bin is not its own neighbor)
diag(adj_matrix) <- FALSE

# Step 5: Create graph from adjacency matrix
graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")

# Step 6: Identify connected components (clusters of spatially adjacent cholangiocytes)
clusters <- components(graph)$membership

# Step 7: Assign cluster labels
cholangiocyte_bins$cholangiocyte_cluster <- paste0("cholangiocyte cluster ", clusters)

# Step 8: Add cluster labels back to the full Seurat cholangiocyte_bins
object$cholangiocyte_cluster <- NA  # Initialize column
object$cholangiocyte_cluster[colnames(cholangiocyte_bins)] <- cholangiocyte_bins$cholangiocyte_cluster


#calculate sizes 

# Step 1: Get the count of bins in each cholangiocyte cluster
bin_count <- table(cholangiocyte_bins$cholangiocyte_cluster)

# Step 2: Re-classify based on the new criteria
cholangiocyte_bins$cholangiocyte_type <- sapply(cholangiocyte_bins$cholangiocyte_cluster, function(cluster) {
  bin_num <- bin_count[as.character(cluster)]
  if (bin_num <= 2) {
    return("technical artifact")
  } else if (bin_num < 20) {
    return("small duct")
  } else if (bin_num >= 20) {
    return("large duct")
  }
})

# Step 4: Add this classification back to the full Seurat cholangiocyte_bins
object$cholangiocyte_type <- NA  # Initialize column
object$cholangiocyte_type[colnames(cholangiocyte_bins)] <- cholangiocyte_bins$cholangiocyte_type

table(object$cholangiocyte_type)

annotation_C107_8um_bile_duct<-object@meta.data
write.csv(annotation_C107_8um_bile_duct, file=here("data/metadata_for_C107_bileductsize.csv"))


# plot bile ducts annotated by size
chol_coords$bin<-rownames(chol_coords)
plt_cholangiocyte_spatial<-merge(chol_coords,cholangiocyte_bins@meta.data, by.x="bin", by.y="X")



bin_size <- 13.4543
half_bin <- bin_size / 2
unit_per_um <- bin_size / 8

# Choose scale bar size in μm
scale_length_um <- 1000
scale_length_plot <- scale_length_um * unit_per_um
x_start <- (-max(plt_cholangiocyte_spatial$x))
y_start <- max(plt_cholangiocyte_spatial$y)

get_leg = function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

plot_bins<-ggplot() +
  geom_point(data = plt_cholangiocyte_spatial,aes(x = y, y = -x, color = cholangiocyte_type),size = 0.25) +
  # Scale bar
  annotate("segment",x = y_start-scale_length_plot-100,xend = y_start-100 ,y = x_start+ 300,yend = x_start+300,color = "black",linewidth = 1) +
  annotate("text",x = y_start - (scale_length_plot / 2)-100,y = x_start + 200,label = paste0(scale_length_um, " µm"),color = "black",size = 3,hjust = 0.5) +
  # Black border around image
  geom_rect(aes(xmin = min(plt_cholangiocyte_spatial$y), xmax = max(plt_cholangiocyte_spatial$y), 
                ymin = (-min(plt_cholangiocyte_spatial$x)), ymax = (-max(plt_cholangiocyte_spatial$x))), 
            fill = NA,color = "black",linewidth = 1) +
  coord_fixed() +   theme_void()
plot_bins


rm(object)
gc()



coords_8um_all<-coords[which(coords$x>xmin & coords$x<xmax & coords$y>ymin & coords$y<ymax),]
coords_8um_all$y_centred <- (coords_8um_all$y - ymin)+1
coords_8um_all$x_centred <- (xmin - (coords_8um_all$x)+1)+(xmax-xmin)

coords_8um<-plt_cholangiocyte_spatial[which(plt_cholangiocyte_spatial$x>xmin & plt_cholangiocyte_spatial$x<xmax & plt_cholangiocyte_spatial$y>ymin & plt_cholangiocyte_spatial$y<ymax),]
coords_8um$y_centred <- (coords_8um$y - ymin)+1
coords_8um$x_centred <- (xmin - (coords_8um$x)+1)+(xmax-xmin)


bin_size <- 13.4543/4
half_bin <- bin_size / 2

coords_rect <- coords_8um |> 
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

## save cell type legend
leg_ex<-ggplot() +  geom_point(aes(x = c(color_possibilities_celltype,NA), y = 1, color =  c(color_possibilities_celltype,NA)),size = 2.5, shape=15) +  colscale_cellType + theme_void()



hulls <- coords_8um %>%
  group_by(cholangiocyte_cluster,cholangiocyte_type) %>%
  slice(chull(x_centred, y_centred))  # chull() gives convex hull indices


plot_bins<-ggplot() +
  geom_point(data = coords_8um_all,aes(x = y_centred, y = x_centred, color = final_anno),size = 2, shape=15) +
  # Scale bar
  geom_polygon(data = hulls, aes(x=y_centred, y=x_centred, group = cholangiocyte_cluster), alpha = 0, color="black") +
  annotate("segment",x = x_start,xend = x_start + scale_length_plot,y = y_start-20,yend = y_start-20,color = "black",linewidth = 1) +
  annotate("text",x = x_start + scale_length_plot / 2,y = y_start - 60,label = paste0(scale_length_um, " µm"),color = "black",size = 3,hjust = 0.5) +
  # Black border around image
  geom_rect(aes(xmin = min(coords_8um_all$y_centred), xmax = max(coords_8um_all$y_centred), ymin = (min(coords_8um_all$x_centred)), ymax = (max(coords_8um_all$x_centred))),fill = NA,color = "black",linewidth = 1) +
  # zoom area
  coord_fixed() +  colscale_cellType+  theme_void()+  theme(legend.position = "none")

ducts<-ggplot() +
  geom_raster(data = df, aes(x = x, y = y, fill = color), show.legend = FALSE, alpha=1) +
  scale_fill_identity() +
  geom_polygon(data = hulls, aes(x=y_centred, y=x_centred, color = factor(cholangiocyte_type), group = cholangiocyte_cluster), alpha = 0) +
  annotate("segment",x = x_start,xend = x_start + scale_length_plot,y = y_start-20,yend = y_start-20,color = "black",linewidth = 1) +
  annotate("text",x = x_start + scale_length_plot / 2,y = y_start - 60,label = paste0(scale_length_um, " µm"),color = "black",size = 3,hjust = 0.5) +
  geom_rect(data = image_bounds,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),fill = NA,color = "black",linewidth = 1) +
  coord_fixed() +     theme_void()+theme(legend.position = "none")+
  scale_color_manual(values=c("#ff9327","#72f507","white"))
ducts

aling_views_ducts<-plot_grid(plot_bins, ducts, ncol=1, align="v")
aling_views_ducts
save_plts(aling_views_ducts, "C107_small_large_ducts", w=8, h=8)


#############
## DEG between cholangiocyte types
#############
Idents(cholangiocyte_bins)<-"cholangiocyte_type"
cholangiocyte_bins<-JoinLayers(cholangiocyte_bins)
cholangiocyte_bins<-NormalizeData(cholangiocyte_bins)

duct_size_markers <- FindMarkers(cholangiocyte_bins, ident.1 = "large duct", ident.2="small duct")

duct_size_markers[which(duct_size_markers$p_val_adj<0.05 & duct_size_markers$avg_log2FC > 0),]
duct_size_markers[which(duct_size_markers$p_val_adj<0.05 & duct_size_markers$avg_log2FC < 0),]

rownames(duct_size_markers[which(duct_size_markers$p_val_adj<0.01),])

VlnPlot(cholangiocyte_bins, features = rownames(duct_size_markers[which(duct_size_markers$p_val_adj<0.05 & duct_size_markers$avg_log2FC > 0),]))


gene_exp<-FetchData(cholangiocyte_bins, vars=c("EGR1","LCN2"))
gene_exp$cell<-rownames(gene_exp)
gene_exp<-melt(gene_exp)

plt<-merge(gene_exp, cholangiocyte_bins@meta.data, by.x="cell", by.y="X")

DEG_duct<-ggplot(plt[which(plt$cholangiocyte_type!="technical artifact"),], aes(cholangiocyte_type, value, fill=cholangiocyte_type))+
  geom_violin(scale = "width",  trim = TRUE)+theme_bw()+facet_wrap(~variable)+ scale_fill_manual(values=c("#ff9327","#72f507","white"), name="Bile\nDuct\nType")+
  xlab("Bile Duct Type")+ylab("Gene Expression")
save_plts(DEG_duct, "DEG_duct_size", w=5, h=3)

#### bar plot of types
df_count<-as.data.frame(table(cholangiocyte_bins$cholangiocyte_type))
duct_count<-ggplot(df_count[which(df_count$Var1!="technical artifact"),], aes(Var1, Freq, fill=Var1))+
  geom_bar(stat="identity", color="black")+ scale_fill_manual(values=c("#ff9327","#72f507","white"), name="Bile\nDuct\nType")+theme_bw()+
  xlab("Bile Duct Type")+ylab("8µm Bin Count")
save_plts(duct_count, "duct_count", w=3.5, h=4)
