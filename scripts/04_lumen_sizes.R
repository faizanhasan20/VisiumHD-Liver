library(tiff)
library(here)
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(arrow)
library(reshape2)
library(cowplot)
library(ggforce)
library(viridis)
library(colorspace)
library(igraph)
library(tidyr)

library(png)
library(grid)
library(ggplot2)
library(dplyr)
library(sf)


source(here("scripts/00_pretty_plots.R"))


# Read PNG
img <- png::readPNG("/media/redgar/Seagate Portable Drive/visiumHD_liver/MacParland_Diana__C107_D1/outs/spatial/tissue_hires_image.png")

# Step 2: Extract dimensions
h <- dim(img)[1]
w <- dim(img)[2]

# Step 3: Create a dataframe of all pixel coordinates and colors
df <- expand.grid(
  y = 1:h,
  x = 1:w
)

df$R <- as.vector(img[,,1])
df$G <- as.vector(img[,,2])
df$B <- as.vector(img[,,3])

# Optional: handle alpha
if (dim(img)[3] == 4) {
  df$A <- as.vector(img[,,4])
  df <- df[df$A > 0, ]  # remove fully transparent pixels
}

# Step 4: Cluster pixels by color
set.seed(123)
k <- 5
km <- kmeans(df[, c("R", "G", "B")], centers = k)
df$cluster <- as.factor(km$cluster)

ggplot(df, aes(x = x, y = h - y, fill = cluster)) +  # h - y flips vertically
  geom_raster() +
  coord_equal() +
  theme_void() +
  scale_fill_manual(values = rgb(km$centers))


cluster_means <- df %>%
  group_by(cluster) %>%
  summarize(
    R = mean(R), G = mean(G), B = mean(B),
    n = n()
  ) %>%
  mutate(color = rgb(R, G, B))

# Identify background cluster (brightest)
background_cluster <- cluster_means %>%
  mutate(brightness = (R + G + B) / 3) %>%
  slice_max(brightness, n = 1) %>%
  pull(cluster)


df_background<-df[which(df$cluster==background_cluster),]

ggplot(df_background, aes(x = x, y = h - y)) +  
  geom_raster() +
  coord_equal() +
  theme_void() 



# Extract spatial coordinates
lumen_coords <- df_background[,c("y","x")]

# Define tile size (adjust to image scale)
tile_size <- 100

# Define tile IDs
lumen_coords <- lumen_coords %>%
  mutate(
    tile_x = floor(x / tile_size),
    tile_y = floor(y / tile_size),
    tile_id = paste(tile_x, tile_y, sep = "_")
  )

rm(df_background)
rm(df)
gc()

## only those roughly covered by visiumHD
#lumen_coords_rough_visiumHD<-lumen_coords[which(lumen_coords$x>700 & lumen_coords$x<3600 & lumen_coords$y>300 ),]
lumen_coords_filtered <- lumen_coords[
  (lumen_coords$x > 700 & lumen_coords$x < 3600 & lumen_coords$y > 600 & lumen_coords$y < 2800),  # region 4
]

keep_area<-data.frame(
  x = c(700, 3600, 3600, 3000, 3000, 2000,2000,700,700),
  y = c(2800, 2800, 600, 600, 500, 500,450,450,2800)
)


ggplot(lumen_coords_filtered, aes(x = x, y = y)) +  
  geom_raster() + 
  coord_equal() 




ggplot(lumen_coords_filtered, aes(x = x, y = h - y)) +  
  geom_raster() +
  coord_equal() +
  theme_void() 

rm(lumen_coords)
gc()

lumen_coords_tiled<-lapply(1:length(unique(lumen_coords_filtered$tile_id)), function(x){
  print(paste("tile",x, "/", length(unique(lumen_coords_filtered$tile_id))))
  
  lumen_coords_tile<-lumen_coords_filtered[which(lumen_coords_filtered$tile_id==unique(lumen_coords_filtered$tile_id)[x]),]
  # Compute pairwise Euclidean distance matrix
  dist_matrix <- as.matrix(dist(lumen_coords_tile))
  
  # Define adjacency based on a distance threshold (adjust as needed)
  distance_threshold <- 2  # Adjust based on your dataset scale
  adj_matrix <- dist_matrix < distance_threshold  # Create adjacency matrix (TRUE for close neighbors)
  
  # Ensure diagonal is FALSE (a bin is not its own neighbor)
  diag(adj_matrix) <- FALSE
  
  # Create graph from adjacency matrix
  graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  
  # Identify connected components (clusters of spatially adjacent cholangiocytes)
  clusters <- components(graph)$membership
  
  rm(dist_matrix)
  rm(adj_matrix)
  gc()
  
  # Step 7: Assign cluster labels
  lumen_coords_tile$lumen_cluster <- paste0("tile_",x,"_cluster_", clusters)
  lumen_coords_tile
})

# ggplot(lumen_coords_tile, aes(x = x, y = y, fill = factor(lumen_cluster))) +
#   geom_tile() +
#   coord_fixed() +
#   theme_void()+theme(legend.position = "none")


lumen_coords_tiled<-do.call(rbind, lumen_coords_tiled)

save(lumen_coords_tiled, file=here("data/lumen_coords_tiled.RData"))

ggplot(lumen_coords_tiled, aes(x = x, y = -y, fill = factor(lumen_cluster))) +
  geom_tile() +
  coord_fixed() +
  theme_void()+theme(legend.position = "none")




#calculate sizes 
bin_count <- as.data.frame(table(lumen_coords_tiled$lumen_cluster))
head(bin_count[rev(order(bin_count$Freq)),])

lumen_coords_tiled<-merge(lumen_coords_tiled, bin_count, by.x="lumen_cluster", by.y="Var1")

ggplot(lumen_coords_tiled, aes(x = x, y = -y, fill = Freq)) +
  geom_tile() +
  coord_fixed() 


# Step 2: Re-classify based on the new criteria
lumen_coords_tiled$lumen_type <- "sinusoid"
lumen_coords_tiled$lumen_type[which(lumen_coords_tiled$Freq>100)]<-"small vein"
lumen_coords_tiled$lumen_type[which(lumen_coords_tiled$Freq>800)]<-"large vein"

lumen_coords_tiled<-as.data.frame(lumen_coords_tiled)

ggplot(lumen_coords_tiled, aes(x = x, y = -y, fill = factor(lumen_type))) +
  geom_tile() +
  coord_fixed() +scale_fill_manual(values=c("blue","black","red"))


lumen_coords_tiled_veins<-lumen_coords_tiled[which(lumen_coords_tiled$lumen_type!="sinusoid"),]

ggplot(lumen_coords_tiled_veins, aes(x = x, y = -y, fill = factor(lumen_type))) +
  geom_tile() +
  coord_fixed() +scale_fill_manual(values=c("blue","black","red"))


#######################
## round two graph
#######################
lumen_coords2 <- lumen_coords_tiled_veins[,c("y", "x")]

# Define tile size (adjust to image scale)
tile_size <- 1000

# Define tile IDs
lumen_coords2 <- lumen_coords2 %>%
  mutate(
    tile_x = floor(x / tile_size),
    tile_y = floor(y / tile_size),
    tile_id = paste(tile_x, tile_y, sep = "_")
  )

rm(lumen_coords_tiled_veins)
rm(lumen_coords_tiled)
gc()


lumen_coords_tiled_roundtwo<-lapply(1:length(unique(lumen_coords2$tile_id)), function(x){
  print(paste("tile",x, "/", length(unique(lumen_coords2$tile_id))))
  
  lumen_coords_tile<-lumen_coords2[which(lumen_coords2$tile_id==unique(lumen_coords2$tile_id)[x]),]
  # Compute pairwise Euclidean distance matrix
  dist_matrix <- as.matrix(dist(lumen_coords_tile))
  
  # Define adjacency based on a distance threshold (adjust as needed)
  distance_threshold <- 2  # Adjust based on your dataset scale
  adj_matrix <- dist_matrix < distance_threshold  # Create adjacency matrix (TRUE for close neighbors)
  
  # Ensure diagonal is FALSE (a bin is not its own neighbor)
  diag(adj_matrix) <- FALSE
  
  # Create graph from adjacency matrix
  graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  
  # Identify connected components (clusters of spatially adjacent cholangiocytes)
  clusters <- components(graph)$membership
  
  rm(dist_matrix)
  rm(adj_matrix)
  gc()
  
  # Step 7: Assign cluster labels
  lumen_coords_tile$lumen_cluster <- paste0("tile_",x,"_cluster_", clusters)
  lumen_coords_tile
})

lumen_coords_tiled_roundtwo<-do.call(rbind, lumen_coords_tiled_roundtwo)

ggplot(lumen_coords_tiled_roundtwo, aes(x = x, y = -y, fill = factor(lumen_cluster))) +
  geom_tile() +
  coord_fixed() +
  theme_void()+theme(legend.position = "none")



#calculate sizes 
bin_count <- as.data.frame(table(lumen_coords_tiled_roundtwo$lumen_cluster))
head(bin_count[rev(order(bin_count$Freq)),])

lumen_coords_tiled_roundtwo<-merge(lumen_coords_tiled_roundtwo, bin_count, by.x="lumen_cluster", by.y="Var1")

ggplot(lumen_coords_tiled_roundtwo, aes(x = x, y = -y, fill = Freq)) +
  geom_tile() +
  coord_fixed() 



# Step 2: Re-classify based on the new criteria
lumen_coords_tiled_roundtwo$lumen_type <- "sinusoid"
lumen_coords_tiled_roundtwo$lumen_type[which(lumen_coords_tiled_roundtwo$Freq>150)]<-"small vein"
lumen_coords_tiled_roundtwo$lumen_type[which(lumen_coords_tiled_roundtwo$Freq>600)]<-"large vein"

lumen_coords_tiled_roundtwo<-as.data.frame(lumen_coords_tiled_roundtwo)

ggplot(lumen_coords_tiled_roundtwo, aes(x = x, y = -y, fill = factor(lumen_type))) +
  geom_tile() +
  coord_fixed() +scale_fill_manual(values=c("blue","black","red"))



lumen_coords_tiled_core_veins<-lumen_coords_tiled_roundtwo[which(lumen_coords_tiled_roundtwo$lumen_type!="sinusoid"),]

ggplot(lumen_coords_tiled_core_veins, aes(x = x, y = -y, fill = factor(lumen_type))) +
  geom_tile() +
  coord_fixed() +scale_fill_manual(values=c("blue","black","red"))

save(lumen_coords_tiled_core_veins, file="lumen_coords_tiled_core_veins.RData")



#################
## polygons from point to overlay onto transcription data
#################
load(here("data/lumen_coords_tiled_core_veins.RData"))

# Convert to sf points
lumen_points <- st_as_sf(lumen_coords_tiled_core_veins, coords = c("x", "y"), crs = NA)

# Create one polygon per lumen_cluster using convex hull
lumen_polygons <- lumen_points %>%
  group_by(lumen_cluster) %>%
  summarise(geometry = st_combine(geometry) |> st_convex_hull()) %>%
  ungroup()

# Check result
plot(lumen_polygons["lumen_cluster"])

#################
## overlay on VisiumHD - check nUMI and histology align
#################
localdir <- "/media/redgar/Seagate Portable Drive/visiumHD_liver/reseq/MacParland_Diana__C107_D1/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8))

# Prepare spot coordinates
coords <- as.data.frame(cbind(object@meta.data, object@images$slice1$centroids@coords))
coords$bin<-rownames(coords)
rm(object)
gc()



cols <- viridis(100)

ggplot() +
  coord_fixed() +
  geom_point(data = coords, aes(x = (y)/4.155, y = (x)/4.155, color=nCount_Spatial.008um),  size = 0.25, shape=15) + 
  geom_sf(aes(color = lumen_cluster), data=lumen_polygons, color = "black", linewidth = 0.2, alpha=0.5) +
  theme(legend.position = "none")+
  scale_color_gradientn(colors = cols, transform="sqrt", name = "Number of UMI\nper 8um bin") 
  
coords$lumen_aligned_x<-coords$y/4.159
coords$lumen_aligned_y<-coords$x/4.159


### bins in lumen
cell_points_sf <- st_as_sf(coords, coords = c("lumen_aligned_x", "lumen_aligned_y"), crs = NA)

# Spatial join: assign each point to a lumen polygon if it’s inside
points_in_lumen <- st_join(cell_points_sf, lumen_polygons, join = st_within)

bins_in_lumen<-ggplot() +
  coord_fixed() +
  geom_point(data = points_in_lumen, aes(x = y, y = -x, color=lumen_cluster),  size = 0.25, shape=15) + 
  theme(legend.position = "none")
save_plts(bins_in_lumen, "bins_in_lumen", w=10,h=10)



### bins adjacent to lumen
bin_size <- 13.4543
unit_per_um <- bin_size / 8
# adjacent will be 2 bins
(16 * unit_per_um) / 4.159

buffer_polygons <- st_buffer(lumen_polygons, dist = (16 * unit_per_um) / 4.159) # where dist_um <- dist * 4.159 / unit_per_um

plot(lumen_polygons["lumen_cluster"])
plot(buffer_polygons["lumen_cluster"])

### bins in lumen
cell_points_sf <- st_as_sf(coords, coords = c("lumen_aligned_x", "lumen_aligned_y"), crs = NA)

# Spatial join: assign each point to a lumen polygon if it’s inside
points_around_lumen <- st_join(cell_points_sf, buffer_polygons, join = st_within)

bins_around_lumen<-ggplot() +
  coord_fixed() +
  geom_point(data = points_around_lumen, aes(x = y, y = -x, color=lumen_cluster),  size = 0.25, shape=15) + 
  theme(legend.position = "none")
save_plts(bins_around_lumen, "bins_around_lumen", w=10,h=10)

neighbouring<-points_around_lumen[which(!(is.na(points_around_lumen$lumen_cluster))),]
in_lumen<-points_in_lumen[which(!(is.na(points_in_lumen$lumen_cluster))),]
neighbouring<-neighbouring[which(!(neighbouring$bin%in%in_lumen$bin)),]
neighbouring$lumen_adjacent<-neighbouring$lumen_cluster
neighbouring$lumen_cluster<-NA

points_in_lumen$lumen_adjacent<-NA
not_neighbouring<-points_in_lumen[which(!(points_in_lumen$bin%in%neighbouring$bin)),]

points_in_lumen_plus_adjacent<-rbind(not_neighbouring, neighbouring)
points_in_lumen_plus_adjacent<-points_in_lumen_plus_adjacent[match(points_in_lumen$bin,points_in_lumen_plus_adjacent$bin),]
identical(points_in_lumen_plus_adjacent$bin, points_in_lumen$bin)

head(points_in_lumen_plus_adjacent[which(!is.na(points_in_lumen_plus_adjacent$lumen_adjacent)),])

points_in_lumen_plus_adjacent<-as.data.frame(points_in_lumen_plus_adjacent)




### bins adjacent to lumen (larger to label zonation)
# expanded will be 10 bins
(80 * unit_per_um) / 4.159

buffer_polygons_lrg <- st_buffer(lumen_polygons, dist = (80 * unit_per_um) / 4.159) # where dist_um <- dist * 4.159 / unit_per_um
plot(lumen_polygons["lumen_cluster"])
plot(buffer_polygons["lumen_cluster"])
plot(buffer_polygons_lrg["lumen_cluster"])

### bins in lumen
cell_points_sf <- st_as_sf(coords, coords = c("lumen_aligned_x", "lumen_aligned_y"), crs = NA)
points_around_lumen_lrg <- st_join(cell_points_sf, buffer_polygons_lrg, join = st_within)

neighbouring_lrg<-points_around_lumen_lrg[which(!(is.na(points_around_lumen_lrg$lumen_cluster))),]
in_lumen_or_adjacent<-points_in_lumen_plus_adjacent[which(!(is.na(points_in_lumen_plus_adjacent$lumen_cluster)) | !(is.na(points_in_lumen_plus_adjacent$lumen_adjacent))),]
neighbouring_lrg<-neighbouring_lrg[which(!(neighbouring_lrg$bin%in%in_lumen_or_adjacent$bin)),]
neighbouring_lrg$lumen_expanded<-neighbouring_lrg$lumen_cluster
neighbouring_lrg$lumen_cluster<-NA
neighbouring_lrg$lumen_adjacent<-NA

points_in_lumen_plus_adjacent$lumen_expanded<-NA
not_neighbouring_lrg<-points_in_lumen_plus_adjacent[which(!(points_in_lumen_plus_adjacent$bin%in%neighbouring_lrg$bin)),]

neighbouring_lrg<-neighbouring_lrg[,match(colnames(not_neighbouring_lrg),colnames(neighbouring_lrg))]

points_in_lumen_plus_adjacentexpanded<-rbind(not_neighbouring_lrg, neighbouring_lrg)
points_in_lumen_plus_adjacentexpanded<-points_in_lumen_plus_adjacentexpanded[match(points_in_lumen$bin,points_in_lumen_plus_adjacentexpanded$bin),]
identical(points_in_lumen_plus_adjacentexpanded$bin, points_in_lumen$bin)

head(points_in_lumen_plus_adjacentexpanded[which(!is.na(points_in_lumen_plus_adjacentexpanded$lumen_expanded)),])

points_in_lumen_plus_adjacentexpanded<-as.data.frame(points_in_lumen_plus_adjacentexpanded)
save(points_in_lumen_plus_adjacentexpanded, file=here("data/points_in_lumen.RData"))


##############################
# sanity check
##############################
tiff_path<-here("/media/redgar/Seagate Portable Drive/visiumHD_liver/high_res_HnE/H1-NC8X8J2_D1.tif")
tiff_res <- readTIFF(tiff_path, as.is = TRUE)
dim(tiff_res)

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

# Get image dimensions
h <- dim(rgb_img)[1]
w <- dim(rgb_img)[2]

# Convert matrix to data frame
df <- melt(rgb_img)
colnames(df) <- c("y", "x", "color")
df$y <- h - df$y + 1  # Flip y-axis


###  nUMI for 2um in area
load(here("data/points_in_lumen.RData"))
points_in_lumen_plus_adjacentexpanded<-points_in_lumen_plus_adjacentexpanded[which(points_in_lumen_plus_adjacentexpanded$x>xmin & points_in_lumen_plus_adjacentexpanded$x<xmax & points_in_lumen_plus_adjacentexpanded$y>ymin & points_in_lumen_plus_adjacentexpanded$y<ymax),]

points_in_lumen_plus_adjacentexpanded$y_centred <- (points_in_lumen_plus_adjacentexpanded$y - ymin)+1
points_in_lumen_plus_adjacentexpanded$x_centred <- (xmin - (points_in_lumen_plus_adjacentexpanded$x)+1)+(xmax-xmin)


ggplot() +
  geom_raster(data = df, aes(x = x, y = -y, fill = color), show.legend = FALSE) +
  scale_fill_identity() +
  coord_fixed() +
  geom_point(data = points_in_lumen_plus_adjacentexpanded, aes(x = y_centred, y = -x_centred, color=lumen_cluster),  size = 2.1, shape=15) + 
  theme(legend.position = "none")

ggplot() +
  geom_raster(data = df, aes(x = x, y = -y, fill = color), show.legend = FALSE) +
  scale_fill_identity() +
  coord_fixed() +
  geom_point(data = points_in_lumen_plus_adjacentexpanded, aes(x = y_centred, y = -x_centred, color=lumen_adjacent),  size = 2.1, shape=15) + 
  theme(legend.position = "none")

bin_size <- 13.4543
half_bin <- bin_size / 2
unit_per_um <- bin_size / 8

# Choose scale bar size in μm
scale_length_um <- 74
scale_length_plot <- scale_length_um * unit_per_um
x_start <- (-max(points_in_lumen_plus_adjacentexpanded$x_centred))
y_start <- max(points_in_lumen_plus_adjacentexpanded$y_centred)

ggplot() +
  geom_raster(data = df, aes(x = x, y = -y, fill = color), show.legend = FALSE) +
  scale_fill_identity() +
  coord_fixed() +
  geom_point(data = points_in_lumen_plus_adjacentexpanded, aes(x = y_centred, y = -x_centred, color=lumen_expanded),  size = 2.1, shape=15) + 
  theme(legend.position = "none")+
  annotate("segment",x = y_start-scale_length_plot-210,xend = y_start-210 ,y = x_start+ 170,yend = x_start+170,color = "black",linewidth = 1) +
  annotate("text",x = y_start - (scale_length_plot / 2)-100,y = x_start + 200,label = paste0(scale_length_um, " µm"),color = "black",size = 3,hjust = 0.5) 
  


##############################
# sanity check II
##############################
tiff_path<-here("/media/redgar/Seagate Portable Drive/visiumHD_liver/high_res_HnE/H1-NC8X8J2_D1.tif")
tiff_res <- readTIFF(tiff_path, as.is = TRUE)
dim(tiff_res)

## coordiantes (CV)
xmin<-6920
xmax<-7900
ymin<-10160
ymax<-12480

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

# Get image dimensions
h <- dim(rgb_img)[1]
w <- dim(rgb_img)[2]

# Convert matrix to data frame
df <- melt(rgb_img)
colnames(df) <- c("y", "x", "color")
df$y <- h - df$y + 1  # Flip y-axis


###  nUMI for 2um in area
load(here("data/points_in_lumen.RData"))
points_in_lumen_plus_adjacentexpanded<-points_in_lumen_plus_adjacentexpanded[which(points_in_lumen_plus_adjacentexpanded$x>xmin & points_in_lumen_plus_adjacentexpanded$x<xmax & points_in_lumen_plus_adjacentexpanded$y>ymin & points_in_lumen_plus_adjacentexpanded$y<ymax),]

points_in_lumen_plus_adjacentexpanded$y_centred <- (points_in_lumen_plus_adjacentexpanded$y - ymin)+1
points_in_lumen_plus_adjacentexpanded$x_centred <- (xmin - (points_in_lumen_plus_adjacentexpanded$x)+1)+(xmax-xmin)


ggplot() +
  geom_raster(data = df, aes(x = x, y = -y, fill = color), show.legend = FALSE) +
  scale_fill_identity() +
  coord_fixed() +
  geom_point(data = points_in_lumen_plus_adjacentexpanded, aes(x = y_centred, y = -x_centred, color=lumen_cluster),  size = 1, shape=15) + 
  theme(legend.position = "none")

ggplot() +
  geom_raster(data = df, aes(x = x, y = -y, fill = color), show.legend = FALSE) +
  scale_fill_identity() +
  coord_fixed() +
  geom_point(data = points_in_lumen_plus_adjacentexpanded, aes(x = y_centred, y = -x_centred, color=lumen_adjacent),  size = 2.1, shape=15) + 
  theme(legend.position = "none")

ggplot() +
  geom_raster(data = df, aes(x = x, y = -y, fill = color), show.legend = FALSE) +
  scale_fill_identity() +
  coord_fixed() +
  geom_point(data = points_in_lumen_plus_adjacentexpanded, aes(x = y_centred, y = -x_centred, color=lumen_expanded),  size = 2.1, shape=15) + 
  theme(legend.position = "none")



###########################
## classify lumen by size
###########################
load(here("data/points_in_lumen.RData"))

bin_count <- as.data.frame(table(points_in_lumen_plus_adjacentexpanded$lumen_cluster))
head(bin_count[rev(order(bin_count$Freq)),])

points_in_lumen_plus_adjacentexpanded_lumen<-merge(points_in_lumen_plus_adjacentexpanded, bin_count, by.x="lumen_cluster", by.y="Var1")

ggplot(points_in_lumen_plus_adjacentexpanded_lumen, aes(x = y, y = -x, color = Freq)) +
  geom_point(size=0.5) +
  coord_fixed() 


# classify based on number of bins
points_in_lumen_plus_adjacentexpanded_lumen$lumen_type <- "small vein"
points_in_lumen_plus_adjacentexpanded_lumen$lumen_type[which(points_in_lumen_plus_adjacentexpanded_lumen$Freq>100)]<-"large vein"

points_in_lumen_plus_adjacentexpanded_lumen<-as.data.frame(points_in_lumen_plus_adjacentexpanded_lumen)

ggplot(points_in_lumen_plus_adjacentexpanded_lumen, aes(x = y, y = -x, color = factor(lumen_type))) +
  geom_point(size=0.5) +
  coord_fixed() +scale_color_manual(values=c("blue","black","red"))

adjacent_bins<-points_in_lumen_plus_adjacentexpanded[which(!is.na(points_in_lumen_plus_adjacentexpanded$lumen_adjacent)),]
adjacent_bins<-merge(adjacent_bins, bin_count, by.x="lumen_adjacent", by.y="Var1")
adjacent_bins$lumen_type <- "small vein"
adjacent_bins$lumen_type[which(adjacent_bins$Freq>100)]<-"large vein"

ggplot(adjacent_bins, aes(x = y, y = -x, color = factor(lumen_type))) +
  geom_point(size=0.5) +
  coord_fixed() +scale_color_manual(values=c("blue","black","red"))


###################
## add bin annotation
###################
annotation_C107_8um<-read.csv(here("data/metadata_for_C107.csv"))
annotation_C107_8um$final_anno<-gsub("bin ","layer ",annotation_C107_8um$final_anno)

points_in_lumen_plus_adjacentexpanded_celltype<-merge(points_in_lumen_plus_adjacentexpanded, annotation_C107_8um, by.x="bin", by.y="X")
points_in_lumen_plus_adjacentexpanded_celltype<-merge(points_in_lumen_plus_adjacentexpanded_celltype, bin_count, by.x="lumen_cluster", by.y="Var1", all.x=T)

points_in_lumen_plus_adjacentexpanded_celltype$lumen_type <- "small vein"
points_in_lumen_plus_adjacentexpanded_celltype$lumen_type[which(points_in_lumen_plus_adjacentexpanded_celltype$Freq>100)]<-"large vein"

save(points_in_lumen_plus_adjacentexpanded_celltype, file=here("data/C107_annotation_lumen_sizes.RData"))

### what is adjacent to veins?
adjacent_bins<-points_in_lumen_plus_adjacentexpanded_celltype[which(!is.na(points_in_lumen_plus_adjacentexpanded_celltype$lumen_adjacent)),]
adjacent_bins<-merge(adjacent_bins, bin_count, by.x="lumen_adjacent", by.y="Var1")
adjacent_bins$lumen_type <- "small vein"
adjacent_bins$lumen_type[which(adjacent_bins$Freq.y>100)]<-"large vein"

lumen_summary_adjacent <- adjacent_bins %>%
  group_by(lumen_type, final_anno) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(percentage = 100 * n / sum(n)) %>%
  ungroup()

ggplot(lumen_summary_adjacent, aes(final_anno, percentage, fill=lumen_type))+geom_bar(position=position_dodge(), stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


### zonation of each vein
lumen_summary_per_vein <- adjacent_bins %>%
  group_by(lumen_expanded, final_anno) %>%
  summarise(
    n = n(),
    lumen_type = first(lumen_type),  # keep one lumen_type per lumen_adjacent
    .groups = "drop_last"
  ) %>%
  mutate(percentage = 100 * n / sum(n)) %>%
  ungroup()


lumen_summary_per_vein_cholangiocytes<-lumen_summary_per_vein[which(lumen_summary_per_vein$final_anno=="intrahepatic cholangiocyte"),]

lumen_summary_per_vein_hepatocye_portal<-lumen_summary_per_vein[grep("6|7|8",lumen_summary_per_vein$final_anno),]
lumen_summary_per_vein_hepatocye_portal<-lumen_summary_per_vein_hepatocye_portal[which(lumen_summary_per_vein_hepatocye_portal$n>1),]

lumen_summary_per_vein_hepatocye_central<-lumen_summary_per_vein[grep("1|2|3",lumen_summary_per_vein$final_anno),]
lumen_summary_per_vein_hepatocye_central<-lumen_summary_per_vein_hepatocye_central[which(lumen_summary_per_vein_hepatocye_central$n>1),]

length(intersect(lumen_summary_per_vein_hepatocye_central$lumen_adjacent, lumen_summary_per_vein_hepatocye_portal$lumen_adjacent))
length(intersect(lumen_summary_per_vein_cholangiocytes$lumen_adjacent, lumen_summary_per_vein_hepatocye_portal$lumen_adjacent))
length(intersect(lumen_summary_per_vein_hepatocye_central$lumen_adjacent, lumen_summary_per_vein_cholangiocytes$lumen_adjacent))




### zonation of each vein
expanded_bins<-points_in_lumen_plus_adjacentexpanded_celltype[which(!is.na(points_in_lumen_plus_adjacentexpanded_celltype$lumen_expanded)),]
expanded_bins<-merge(expanded_bins, bin_count, by.x="lumen_expanded", by.y="Var1")
expanded_bins$lumen_type <- "small vein"
expanded_bins$lumen_type[which(expanded_bins$Freq.y>100)]<-"large vein"

lumen_summary_expanded <- expanded_bins %>%
  group_by(lumen_expanded, final_anno) %>%
  summarise(
    n = n(),
    lumen_type = first(lumen_type),  # keep one lumen_type per lumen_adjacent
    .groups = "drop_last"
  ) %>%
  mutate(percentage = 100 * n / sum(n)) %>%
  ungroup()


lumen_summary_per_vein_cholangiocytes<-lumen_summary_per_vein[which(lumen_summary_per_vein$final_anno=="intrahepatic cholangiocyte"),]
lumen_summary_per_vein_cholangiocytes<-lumen_summary_per_vein_cholangiocytes[which(lumen_summary_per_vein_cholangiocytes$percentage>2),]

lumen_summary_per_vein_hepatocye_portal<-lumen_summary_per_vein[grep("6|7|8",lumen_summary_per_vein$final_anno),]
lumen_summary_per_vein_hepatocye_portal<-lumen_summary_per_vein_hepatocye_portal[which(lumen_summary_per_vein_hepatocye_portal$percentage>10),]

lumen_summary_per_vein_hepatocye_central<-lumen_summary_per_vein[grep("1|2|3",lumen_summary_per_vein$final_anno),]
lumen_summary_per_vein_hepatocye_central<-lumen_summary_per_vein_hepatocye_central[which(lumen_summary_per_vein_hepatocye_central$percentage>10),]

length(intersect(lumen_summary_per_vein_hepatocye_central$lumen_expanded, lumen_summary_per_vein_hepatocye_portal$lumen_expanded))
length(intersect(lumen_summary_per_vein_cholangiocytes$lumen_expanded, lumen_summary_per_vein_hepatocye_portal$lumen_expanded))
length(intersect(lumen_summary_per_vein_hepatocye_central$lumen_expanded, lumen_summary_per_vein_cholangiocytes$lumen_expanded))

length(unique(lumen_summary_per_vein$lumen_expanded))
length(unique(lumen_summary_per_vein_hepatocye_central$lumen_expanded))
length(unique(lumen_summary_per_vein_hepatocye_portal$lumen_expanded))

portal<-unique(c(unique(lumen_summary_per_vein_cholangiocytes$lumen_expanded),
          unique(lumen_summary_per_vein_hepatocye_portal$lumen_expanded)[which(!(unique(lumen_summary_per_vein_hepatocye_portal$lumen_expanded)%in%unique(lumen_summary_per_vein_hepatocye_central$lumen_expanded)))]))
central<- unique(lumen_summary_per_vein_hepatocye_central$lumen_expanded)[which(!(unique(lumen_summary_per_vein_hepatocye_central$lumen_expanded)%in%unique(lumen_summary_per_vein_hepatocye_portal$lumen_expanded)))]
ambigious<-c(unique(lumen_summary_per_vein$lumen_expanded)[which(!(unique(lumen_summary_per_vein$lumen_expanded)%in%c(portal,central)))])

lumen_summary_per_vein[which(lumen_summary_per_vein$lumen_expanded==ambigious[1]),]



points_in_lumen_plus_adjacentexpanded_celltype$zonation<-NA
points_in_lumen_plus_adjacentexpanded_celltype$zonation[which(points_in_lumen_plus_adjacentexpanded_celltype$lumen_cluster%in%portal)]<-"PV"
points_in_lumen_plus_adjacentexpanded_celltype$zonation[which(points_in_lumen_plus_adjacentexpanded_celltype$lumen_cluster%in%central)]<-"CV"
points_in_lumen_plus_adjacentexpanded_celltype$zonation[which(points_in_lumen_plus_adjacentexpanded_celltype$lumen_cluster%in%ambigious)]<-"Unclear"

ggplot(points_in_lumen_plus_adjacentexpanded_celltype[which(!is.na(points_in_lumen_plus_adjacentexpanded_celltype$lumen_cluster)),], aes(x = y, y = -x, color = factor(zonation))) +
  geom_point(size=0.5) +
  coord_fixed() +scale_color_manual(values=c("#B5E48C","#184E77","grey30"))



#########################
## cell percentages in each region but lumen size
#########################

### what is directly adjacent to veins?
adjacent_bins<-points_in_lumen_plus_adjacentexpanded_celltype[which(!is.na(points_in_lumen_plus_adjacentexpanded_celltype$lumen_adjacent)),]
adjacent_bins<-merge(adjacent_bins, bin_count, by.x="lumen_adjacent", by.y="Var1")
adjacent_bins$lumen_type <- "small vein"
adjacent_bins$lumen_type[which(adjacent_bins$Freq.y>100)]<-"large vein"

lumen_summary_adjacent <- adjacent_bins %>%
  group_by(lumen_type, final_anno) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(percentage = 100 * n / sum(n)) %>%
  ungroup()

lumen_summary_adjacent$area<-"adjacent"

### what is adjacent expanded area to veins?
expanded_bins<-points_in_lumen_plus_adjacentexpanded_celltype[which(!is.na(points_in_lumen_plus_adjacentexpanded_celltype$lumen_expanded)),]
expanded_bins<-merge(expanded_bins, bin_count, by.x="lumen_expanded", by.y="Var1")
expanded_bins$lumen_type <- "small vein"
expanded_bins$lumen_type[which(expanded_bins$Freq.y>100)]<-"large vein"

lumen_summary_expanded <- expanded_bins %>%
  group_by(lumen_type, final_anno) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(percentage = 100 * n / sum(n)) %>%
  ungroup()

lumen_summary_expanded$area<-"expanded"

lumen_summary_adjacent_expanded_plot<-rbind(lumen_summary_adjacent, lumen_summary_expanded)


ggplot(lumen_summary_adjacent_expanded_plot, aes(final_anno, percentage, fill=lumen_type))+
  geom_bar(stat = "identity", position = position_dodge())+facet_wrap(~area)+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(lumen_summary_adjacent_expanded_plot, aes(lumen_type, percentage, fill=area))+
  geom_bar(stat = "identity", position = position_dodge())+
  facet_wrap(~final_anno)+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

lumen_summary_adjacent_expanded_plot$area_size<-as.factor(paste(lumen_summary_adjacent_expanded_plot$area, lumen_summary_adjacent_expanded_plot$lumen_type))
levels(lumen_summary_adjacent_expanded_plot$area_size)<-c("Large\nlumen\nadjacent","Small\nlumen\nadjacent","Large\nlumen\nproximal","Small\nlumen\nproximal")

celltype<-ggplot(lumen_summary_adjacent_expanded_plot[which(lumen_summary_adjacent_expanded_plot$final_anno%in%c("hepatic stellate cell","fibroblast")),], 
       aes(area_size, percentage, fill=area_size))+
  geom_bar(stat = "identity", position = position_dodge(), color="black")+
  facet_wrap(~final_anno)+ xlab("Proximity to lumen")+ylab("Percent of bins annotated as cell type")+
  scale_fill_manual(values=c("#184E77","#6a9942","#93b5cf","#B5E48C"))+theme_bw()+
  theme(legend.position = "none")
save_plts(celltype, "celltype_around_lumen", w=6, h=4)


#### bar plot of types
bin_count$lumen_type <- "small lumen"
bin_count$lumen_type[which(bin_count$Freq>100)]<-"large lumen"

df_count<-as.data.frame(table(bin_count$lumen_type))
df_count$Var1<-as.factor(df_count$Var1)
levels(df_count$Var1)<-c("Large","Small")

lumen_count<-ggplot(df_count, aes(Var1, Freq, fill=Var1))+
  geom_bar(stat="identity", color="black")+ scale_fill_manual(values=c("#184E77","#6a9942","white"), name="Lumen\nType")+theme_bw()+
  theme(legend.position = "none")+
  xlab("Lumen Type")+ylab("Count")
lumen_count
save_plts(lumen_count, "lumen_count", w=2, h=4)






##############################
# representative plot
##############################
tiff_path<-here("/media/redgar/Seagate Portable Drive/visiumHD_liver/high_res_HnE/H1-NC8X8J2_D1.tif")
tiff_res <- readTIFF(tiff_path, as.is = TRUE)
dim(tiff_res)

## coordiantes (CV)
xmin<-6450
xmax<-7800
ymin<-9000
ymax<-11300

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

# Get image dimensions
h <- dim(rgb_img)[1]
w <- dim(rgb_img)[2]

# Convert matrix to data frame
df <- melt(rgb_img)
colnames(df) <- c("y", "x", "color")
df$y <- h - df$y + 1  # Flip y-axis


load(here("data/points_in_lumen.RData"))
bin_count <- as.data.frame(table(points_in_lumen_plus_adjacentexpanded$lumen_cluster))
head(bin_count[rev(order(bin_count$Freq)),])

points_in_lumen_plus_adjacentexpanded<-points_in_lumen_plus_adjacentexpanded[which(points_in_lumen_plus_adjacentexpanded$x>xmin & points_in_lumen_plus_adjacentexpanded$x<xmax & points_in_lumen_plus_adjacentexpanded$y>ymin & points_in_lumen_plus_adjacentexpanded$y<ymax),]

points_in_lumen_plus_adjacentexpanded$y_centred <- (points_in_lumen_plus_adjacentexpanded$y - ymin)+1
points_in_lumen_plus_adjacentexpanded$x_centred <- (xmin - (points_in_lumen_plus_adjacentexpanded$x)+1)+(xmax-xmin)



expanded_bins_plt<-points_in_lumen_plus_adjacentexpanded[which(!is.na(points_in_lumen_plus_adjacentexpanded$lumen_expanded)),c("bin","y_centred","x_centred","lumen_expanded")]
colnames(expanded_bins_plt)[4]<-"lumen_cluster"
expanded_bins_plt$area<-"expanded"
adjacent_bins_plt<-points_in_lumen_plus_adjacentexpanded[which(!is.na(points_in_lumen_plus_adjacentexpanded$lumen_adjacent)),c("bin","y_centred","x_centred","lumen_adjacent")]
colnames(adjacent_bins_plt)[4]<-"lumen_cluster"
adjacent_bins_plt$area<-"adjacent"
plt_lumens<-rbind(adjacent_bins_plt,expanded_bins_plt)

plt_lumens<-merge(plt_lumens, bin_count, by.x="lumen_cluster", by.y="Var1")
plt_lumens$lumen_size<-"Small lumen"
plt_lumens$lumen_size[which(plt_lumens$Freq>100)]<-"Large lumen"
plt_lumens$area_size<-as.factor(paste(plt_lumens$area, plt_lumens$lumen_size))

levels(plt_lumens$area_size)<-c("Large\nlumen\nadjacent","Small\nlumen\nadjacent","Large\nlumen\nproximal","Small\nlumen\nproximal")

bin_size <- 13.4543
half_bin <- bin_size / 2
unit_per_um <- bin_size / 8

# Choose scale bar size in μm
scale_length_um <- 100
scale_length_plot <- scale_length_um * unit_per_um
x_start <- (-max(plt_lumens$x_centred))
y_start <- max(plt_lumens$y_centred)

small_lrg_lumen<-plot_grid(ggplot() +
  geom_raster(data = df, aes(x = x, y = -y, fill = color), show.legend = FALSE) +  scale_fill_identity() +  coord_fixed() +theme_void(),
  ggplot() +
    geom_raster(data = df, aes(x = x, y = -y, fill = color), show.legend = FALSE) +  scale_fill_identity() +  coord_fixed() +
    geom_point(data = plt_lumens, aes(x = y_centred, y = -x_centred, color=area_size),  size = 1, shape=15) +
    scale_color_manual(values=c("#184E77","#6a9942","#93b5cf","#B5E48C"))+theme_void()+
    annotate("segment",x = y_start-scale_length_plot-100,xend = y_start-100 ,y = x_start-100,yend = x_start-100,color = "black",linewidth = 1) +
    annotate("text",x = y_start - (scale_length_plot / 2)-100,y = x_start -150,label = paste0(scale_length_um, " µm"),color = "black",size = 3,hjust = 0.5) ,
  ncol=1, align = "v")
small_lrg_lumen

save_plts(small_lrg_lumen, "C107_small_large_lumen", w=10, h=20)





#############
## DEG between endo in lumen sizes
#############
localdir <- "/media/redgar/Seagate Portable Drive/visiumHD_liver/MacParland_Diana__C107_D1/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8))

annotation_C107_8um_additional<-read.csv(file=here("data/annotation_C107_8um_additional.csv"))
annotation_C107_8um_additional <- annotation_C107_8um_additional[match(colnames(object), annotation_C107_8um_additional$bin),]

identical(annotation_C107_8um_additional$bin, colnames(object))
rownames(annotation_C107_8um_additional) <- annotation_C107_8um_additional$X
annotation_C107_8um_additional$X<-NULL

object<-AddMetaData(object, annotation_C107_8um_additional)

endo<-subset(object, subset = final_anno%in%c("endothelial cell of periportal hepatic sinusoid",
                                              "vein endothelial cell",
                                              "endothelial cell of artery",
                                              "endothelial cell of pericentral hepatic sinusoid"))

rm(object)
gc()

endo$Lumen_size<-as.factor(endo$lumen_association)
levels(endo$Lumen_size)<-c("Large lumen","Large lumen","Small lumen","Small lumen")

Idents(endo)<-"Lumen_size"
table(Idents(endo))
endo<-JoinLayers(endo)
endo<-NormalizeData(endo)

endo_lumen_size_markers <- FindMarkers(endo, ident.1="Large lumen", ident.2="Small lumen", min.pct = 0.1)
endo_lumen_size_markers$gene<-rownames(endo_lumen_size_markers)

endo_lumen_size_markers[which(endo_lumen_size_markers$p_val_adj<0.05 & endo_lumen_size_markers$avg_log2FC > 0),]
endo_lumen_size_markers[which(endo_lumen_size_markers$p_val_adj<0.05 & endo_lumen_size_markers$avg_log2FC < 0),]

unique(endo_lumen_size_markers[which(endo_lumen_size_markers$p_val_adj<0.05 & abs(endo_lumen_size_markers$avg_log2FC) > 1),]$gene)

# negative is mainly hepatocytes
#most highly expressed and up in large
pos<-endo_lumen_size_markers[which(endo_lumen_size_markers$p_val_adj<0.05 & abs(endo_lumen_size_markers$avg_log2FC) > 1),]


endothelial_genes <- c(
  "CD34", "PECAM1", "CLDN5", "CAV1", "RAMP2", "EPAS1", "MCAM", "SOX17",
  "MMRN2", "AQP1", "TM4SF1", "TIMP3", "ADAM15", "NPNT", "MGP", "VEGFC",
  "JAG1", "JAG2", "NOTCH3", "HEY1", "FBLN5", "ELN", "FBN1", "FLNA"
)

pos_endo<-pos[which(pos$gene%in%endothelial_genes),]
pos_endo[order(pos_endo$pct.1),]



#### Pathways
source("scripts/00_GSEA_function.R")
#http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/
GO_file = here("../liver_transplant/data/Human_GO_AllPathways_noPFOCR_no_GO_iea_March_01_2025_symbol.gmt")

gene_list = endo_lumen_size_markers$avg_log2FC
names(gene_list) = endo_lumen_size_markers$gene
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]
res = GSEA(gene_list, GO_file, pval = 0.05)
## no significant pathways


VlnPlot(endo, features = rownames(endo_lumen_size_markers[which(endo_lumen_size_markers$p_val_adj<0.05 & endo_lumen_size_markers$avg_log2FC > 0),]))


gene_exp<-FetchData(endo, vars=c("PECAM1","CLDN5","CD34", "ALB"))
gene_exp$cell<-rownames(gene_exp)
gene_exp<-melt(gene_exp)

plt<-merge(gene_exp, endo@meta.data, by.x="cell", by.y="bin")

coords <- as.data.frame(cbind(endo@meta.data, endo@images$slice1$centroids@coords))
coords$bin<-rownames(coords)

plt<-merge(plt, coords[,c("bin","x","y")], by.x="cell", by.y="bin")


ggplot(plt, aes(Lumen_size, value, fill=Lumen_size))+
  geom_violin(scale = "width",  trim = TRUE)+theme_bw()+facet_wrap(~variable)+ 
  scale_fill_manual(values=c("#184E77","#6a9942"), name="Associated lumen size")+
  xlab("Associated lumen size")+ylab("Gene Expression")

ggplot(plt[which(plt$variable=="PECAM1"& !is.na(plt$Lumen_size) ),], aes(y,-x, color=value))+geom_point()
ggplot(plt[which(plt$variable=="VWF" & !is.na(plt$Lumen_size)),], aes(x,y, color=Lumen_size))+geom_point()


### expression on one region
one_geneROI<-plt[which(plt$variable=="ALB"& !is.na(plt$Lumen_size) ),]

plt_gene_ROI<-merge(plt_lumens,one_geneROI[,c("cell","value")], by.x="bin", by.y="cell")

ggplot() +
  geom_raster(data = df, aes(x = x, y = -y, fill = color), show.legend = FALSE) +  scale_fill_identity() +  coord_fixed() +
  geom_point(data = plt_lumens[which(plt_lumens$lumen_size=="Large lumen"),], aes(x = y_centred, y = -x_centred), color="#184E77",  size = 1, shape=15) +
  geom_point(data = plt_lumens[which(plt_lumens$lumen_size=="Small lumen"),], aes(x = y_centred, y = -x_centred), color="#6a9942",  size = 1, shape=15) +
  geom_point(data = plt_gene_ROI, aes(x = y_centred, y = -x_centred, color=value),  size = 1, shape=15) +
  theme_void()+
  annotate("segment",x = y_start-scale_length_plot-100,xend = y_start-100 ,y = x_start-100,yend = x_start-100,color = "black",linewidth = 1) +
  annotate("text",x = y_start - (scale_length_plot / 2)-100,y = x_start -150,label = paste0(scale_length_um, " µm"),color = "black",size = 3,hjust = 0.5) 




#############
## DEG between all bins in lumen sizes
#############
localdir <- "/media/redgar/Seagate Portable Drive/visiumHD_liver/reseq/MacParland_Diana__C107_D1/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8))

annotation_C107_8um_additional<-read.csv(file=here("data/annotation_C107_8um_additional.csv"))
annotation_C107_8um_additional <- annotation_C107_8um_additional[match(colnames(object), annotation_C107_8um_additional$bin),]

identical(annotation_C107_8um_additional$bin, colnames(object))
rownames(annotation_C107_8um_additional) <- annotation_C107_8um_additional$X
annotation_C107_8um_additional$X<-NULL

object<-AddMetaData(object, annotation_C107_8um_additional)

object$Lumen_size<-as.factor(object$lumen_association)
levels(object$Lumen_size)<-c("Large lumen","Large lumen","Small lumen","Small lumen")

Idents(object)<-"Lumen_size"
table(Idents(object))
object<-JoinLayers(object)
object<-NormalizeData(object)

lumen_size_markers <- FindMarkers(object, ident.1="Large lumen", ident.2="Small lumen", min.pct = 0.1)
lumen_size_markers$gene<-rownames(lumen_size_markers)

lumen_size_markers[which(lumen_size_markers$p_val_adj<0.05 & lumen_size_markers$avg_log2FC > 0),]
lumen_size_markers[which(lumen_size_markers$p_val_adj<0.05 & lumen_size_markers$avg_log2FC < 0),]

unique(lumen_size_markers[which(lumen_size_markers$p_val_adj<0.05 & abs(lumen_size_markers$avg_log2FC) > 1),]$gene)


lumen_size_markers[order(lumen_size_markers$pct.1),]



#### Pathways
source("scripts/00_GSEA_function.R")
#http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/
GO_file = here("../liver_transplant/data/Human_GO_AllPathways_noPFOCR_no_GO_iea_March_01_2025_symbol.gmt")

gene_list = lumen_size_markers$avg_log2FC
names(gene_list) = lumen_size_markers$gene
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]
res = GSEA(gene_list, GO_file, pval = 0.05)
## no significant pathways




gene_exp<-FetchData(object, vars=c("ORM1","ORM2","ADH4", "TDO2"))
gene_exp$cell<-rownames(gene_exp)
gene_exp<-melt(gene_exp)

plt<-merge(gene_exp, object@meta.data, by.x="cell", by.y="bin")

coords <- as.data.frame(cbind(object@meta.data, object@images$slice1$centroids@coords))
coords$bin<-rownames(coords)

plt<-merge(plt, coords[,c("bin","x","y")], by.x="cell", by.y="bin")


ggplot(plt, aes(Lumen_size, value, fill=Lumen_size))+
  geom_violin(scale = "width",  trim = TRUE)+theme_bw()+facet_wrap(~variable)+ 
  scale_fill_manual(values=c("#184E77","#6a9942"), name="Associated lumen size")+
  xlab("Associated lumen size")+ylab("Gene Expression")

ggplot(plt[which(plt$variable=="PECAM1"& !is.na(plt$Lumen_size) ),], aes(y,-x, color=value))+geom_point()
ggplot(plt[which(plt$variable=="VWF" & !is.na(plt$Lumen_size)),], aes(x,y, color=Lumen_size))+geom_point()


### expression on one region
one_geneROI<-plt[which(plt$variable=="ORM1"& !is.na(plt$Lumen_size) ),]

plt_gene_ROI<-merge(plt_lumens,one_geneROI[,c("cell","value")], by.x="bin", by.y="cell")

ggplot() +
  geom_raster(data = df, aes(x = x, y = -y, fill = color), show.legend = FALSE) +  scale_fill_identity() +  coord_fixed() +
  geom_point(data = plt_lumens[which(plt_lumens$lumen_size=="Large lumen"),], aes(x = y_centred, y = -x_centred), color="#184E77",  size = 1, shape=15) +
  geom_point(data = plt_lumens[which(plt_lumens$lumen_size=="Small lumen"),], aes(x = y_centred, y = -x_centred), color="#6a9942",  size = 1, shape=15) +
  geom_point(data = plt_gene_ROI, aes(x = y_centred, y = -x_centred, color=value),  size = 1, shape=15) +
  theme_void()+
  annotate("segment",x = y_start-scale_length_plot-100,xend = y_start-100 ,y = x_start-100,yend = x_start-100,color = "black",linewidth = 1) +
  annotate("text",x = y_start - (scale_length_plot / 2)-100,y = x_start -150,label = paste0(scale_length_um, " µm"),color = "black",size = 3,hjust = 0.5) 

