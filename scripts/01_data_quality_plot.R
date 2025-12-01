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

library(png)
library(grid)
library(sp)



source(here("scripts/00_pretty_plots.R"))



#########################
## nUMI plots
#########################
localdir <- "/media/redgar/Expansion/VisiumHD/spaceranger_out_for_upload/C113_only_spaceranger"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8))
C113 <- as.data.frame(cbind(object@meta.data, object@images$slice1$centroids@coords))
C113$bin<-rownames(C113)
C113$sample<-"C113"
rm(object)
gc()

ggplot(C113, aes(x,y, color=nCount_Spatial.008um))+geom_point(size=0.1)

polygon_coords <- data.frame(
  x = c(1300, 3000, 3000, 1500),
  y = c(200, 200, 1200, 1200)
)

poly <- Polygon(polygon_coords)
polys <- Polygons(list(poly), 1)
sp_poly <- SpatialPolygons(list(polys))

pts <- C113[, c("x", "y")]
sp_pts <- SpatialPoints(pts)
inside <- !is.na(over(sp_pts, sp_poly))

C113_only <- C113[inside, ]


localdir <- "/media/redgar/Seagate Portable Drive/visiumHD_liver/MacParland_Diana__C95-3_A1/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8))
C95 <- as.data.frame(cbind(object@meta.data, object@images$slice1$centroids@coords))
C95$bin<-rownames(C95)
C95$sample<-"C95"
rm(object)
gc()


localdir <- "/media/redgar/Seagate Portable Drive/visiumHD_liver/reseq/MacParland_Diana__C107_D1/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8))
C107 <- as.data.frame(cbind(object@meta.data, object@images$slice1$centroids@coords))
C107$bin<-rownames(C107)
C107$sample<-"C107"
rm(object)
gc()


### combined so same scale (colours)
nUMI<-rbind(C107, C95, C113_only)
nUMI$sample<-factor(nUMI$sample, levels = c("C107","C95","C113"))

nUMI_8um<-nUMI
save(nUMI_8um, file=here("data/QC_8um.RData"))

bin_size <- 13.4543
unit_per_um <- bin_size / 8

# Choose scale bar size in μm
scale_length_um <- 1000
scale_length_plot <- scale_length_um * unit_per_um
scale_length_plot_C113 <- scale_length_plot/8 # coordinate system different in C113

## boundary boxes           
min_max_df <- nUMI %>%
  group_by(sample) %>%
  summarize(
    min_x = min(x, na.rm = TRUE),
    max_x = max(x, na.rm = TRUE),
    min_y = min(y, na.rm = TRUE),
    max_y = max(y, na.rm = TRUE)
  )        
       
min_max_df$scale_bar<-c(scale_length_plot,scale_length_plot,scale_length_plot_C113)

       
cols <- viridis(100)

UMI_number<-ggplot() +
  geom_point(data = nUMI,
             aes(x = y, y = -x, color = nCount_Spatial.008um),
             size = 0.1) + facet_wrap(~sample, scales = "free", ncol=1)+
  # Black border around image
  geom_rect(data=min_max_df, aes(xmin = min_y, xmax = max_y, ymin = -min_x, ymax = -max_x),
            fill = NA,color = "black",linewidth = 1)+
  # Scale bar
  geom_segment(aes(x = max_y-scale_bar-100,xend = max_y-100 ,y = -max_x- 200,yend = -max_x-200),
           data=min_max_df, color = "black",linewidth = 1) +
  geom_text(aes(x = max_y - (scale_bar / 2)-100,y = -max_x -500,label = paste0(scale_length_um, " µm")),
            data=min_max_df, color = "black",size = 3,hjust = 0.5) +
  theme_void()+ theme(strip.text = element_text(size=15))+
  scale_color_gradientn(colors = cols, transform="sqrt", name = "Number of UMI\nper 8um bin") 
UMI_number
save_plts(UMI_number, "number_UMI_allsamples", w=6, h=15)

#flipped

UMI_number<-ggplot() +
  geom_point(data = nUMI,
             aes(x = y, y = -x, color = nCount_Spatial.008um),
             size = 0.1) + facet_wrap(~sample, scales = "free", ncol=3)+
  # Black border around image
  geom_rect(data=min_max_df, aes(xmin = min_y, xmax = max_y, ymin = -min_x, ymax = -max_x),
            fill = NA,color = "black",linewidth = 1)+
  theme_void()+ theme(strip.text = element_text(size=15))+
  scale_color_gradientn(colors = cols, transform="sqrt", name = "Number of UMI\nper 8um bin") 

UMI_number_scales<-ggplot() +
  facet_wrap(~sample, scales = "free", ncol=3)+
  # Black border around image
  geom_rect(data=min_max_df, aes(xmin = min_y, xmax = max_y, ymin = -min_x, ymax = -max_x),
            fill = NA,color = "black",linewidth = 1)+
  # Scale bar
  geom_segment(aes(x = max_y-scale_bar-100,xend = max_y-100 ,y = -max_x- 200,yend = -max_x-200),
               data=min_max_df, color = "black",linewidth = 1) +
  geom_text(aes(x = max_y - (scale_bar / 2)-100,y = -max_x -500,label = paste0(scale_length_um, " µm")),
            data=min_max_df, color = "black",size = 5,hjust = 0.5) +
  theme_void()+ theme(strip.text = element_text(size=15))+
  scale_color_gradientn(colors = cols, transform="sqrt", name = "Number of UMI\nper 8um bin") 

plot_grid(UMI_number, UMI_number_scales, ncol=1, align="v")
save_plts(plot_grid(UMI_number, UMI_number_scales, ncol=1, align="v"), "number_UMI_allsamples_option2", w=15, h=10)



######################## 
## Box plots of metrics
######################## 

### 16um
localdir <- "/media/redgar/Expansion/VisiumHD/spaceranger_out_for_upload/C113_only_spaceranger"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(16))
C113 <- as.data.frame(cbind(object@meta.data, object@images$slice1$centroids@coords))
C113$bin<-rownames(C113)
C113$sample<-"C113"
rm(object)
gc()

polygon_coords <- data.frame(
  x = c(1300, 3000, 3000, 1500),
  y = c(200, 200, 1200, 1200)
)


## confirm shape is right for this bin level
ggplot() +
  geom_point(data = C113,
             aes(x = y, y = -x, color = nCount_Spatial.016um),
             size = 0.1)+
  geom_polygon( data = polygon_coords,aes(x = y, y = -x),fill = NA,    
                color = "black",  linewidth = 0.5,inherit.aes = FALSE)

poly <- Polygon(polygon_coords)
polys <- Polygons(list(poly), 1)
sp_poly <- SpatialPolygons(list(polys))

pts <- C113[, c("x", "y")]
sp_pts <- SpatialPoints(pts)
inside <- !is.na(over(sp_pts, sp_poly))

C113_only <- C113[inside, ]

localdir <- "/media/redgar/Seagate Portable Drive/visiumHD_liver/MacParland_Diana__C95-3_A1/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(16))
C95 <- as.data.frame(cbind(object@meta.data, object@images$slice1$centroids@coords))
C95$bin<-rownames(C95)
C95$sample<-"C95"
rm(object)
gc()


localdir <- "/media/redgar/Seagate Portable Drive/visiumHD_liver/reseq/MacParland_Diana__C107_D1/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(16))
C107 <- as.data.frame(cbind(object@meta.data, object@images$slice1$centroids@coords))
C107$bin<-rownames(C107)
C107$sample<-"C107"
rm(object)
gc()

### combined so same scale (colours)
nUMI<-rbind(C107, C95, C113_only)
nUMI$sample<-factor(nUMI$sample, levels = c("C107","C95","C113"))

nUMI_16um<-nUMI
save(nUMI_16um, file=here("data/QC_16um.RData"))




### 2um
localdir <- "/media/redgar/Expansion/VisiumHD/spaceranger_out_for_upload/C113_only_spaceranger"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(2))
C113 <- as.data.frame(cbind(object@meta.data, object@images$slice1$centroids@coords))
C113$bin<-rownames(C113)
C113$sample<-"C113"
rm(object)
gc()

polygon_coords <- data.frame(
  x = c(1300, 3000, 3000, 1500),
  y = c(200, 200, 1200, 1200)
)


## confirm shape is right for this bin level
ggplot() +
  geom_point(data = C113,
             aes(x = y, y = -x, color = nCount_Spatial.002um),
             size = 0.1)+
  geom_polygon( data = polygon_coords,aes(x = y, y = -x),fill = NA,    
                color = "black",  linewidth = 0.5,inherit.aes = FALSE)

poly <- Polygon(polygon_coords)
polys <- Polygons(list(poly), 1)
sp_poly <- SpatialPolygons(list(polys))

pts <- C113[, c("x", "y")]
sp_pts <- SpatialPoints(pts)
inside <- !is.na(over(sp_pts, sp_poly))

C113_only <- C113[inside, ]
save(C113_only, file=here("data/QC_C113_only_2um.RData"))



localdir <- "/media/redgar/Seagate Portable Drive/visiumHD_liver/MacParland_Diana__C95-3_A1/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(2))
C95 <- as.data.frame(cbind(object@meta.data, object@images$slice1$centroids@coords))
C95$bin<-rownames(C95)
C95$sample<-"C95"
rm(object)
gc()
save(C95, file=here("data/QC_C95_2um.RData"))


localdir <- "/media/redgar/Seagate Portable Drive/visiumHD_liver/reseq/MacParland_Diana__C107_D1/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(2))
C107 <- as.data.frame(cbind(object@meta.data, object@images$slice1$centroids@coords))
C107$bin<-rownames(C107)
C107$sample<-"C107"
rm(object)
gc()
save(C107, file=here("data/QC_C107_2um.RData"))

load(here("data/QC_C107_2um.RData"))
load(here("data/QC_C95_2um.RData"))
load(here("data/QC_C113_only_2um.RData"))

### combined so same scale (colours)
nUMI<-rbind(C107, C95, C113_only)
nUMI$sample<-factor(nUMI$sample, levels = c("C107","C95","C113"))

nUMI_2um<-nUMI
save(nUMI_2um, file=here("data/QC_2um.RData"))




## plot all together
load(here("data/QC_2um.RData"))
load(here("data/QC_8um.RData"))
load(here("data/QC_16um.RData"))

nUMI_2um$bin<-"2um"
nUMI_8um$bin<-"8um"
nUMI_16um$bin<-"16um"

colnames(nUMI_16um)[2:3]<-c("nCount","nFeature")
colnames(nUMI_8um)[2:3]<-c("nCount","nFeature")
colnames(nUMI_2um)[2:3]<-c("nCount","nFeature")

nUMI<-rbind(nUMI_2um, nUMI_8um, nUMI_16um)


df_summary <- nUMI %>%
  group_by(bin, sample) %>%
  summarise(
    median_nCount = median(nCount, na.rm = TRUE),
    median_nFeature = median(nFeature, na.rm = TRUE),
    n_rows = n(),
    .groups = "drop"
  )

df_summary$bin<-factor(df_summary$bin, levels = rev(c("2um","8um", "16um")))
df_summary$bin<-gsub("um","",df_summary$bin)
df_summary$bin<-factor(df_summary$bin, levels=c(2,8,16))

metrcis<-plot_grid(
ggplot(df_summary, aes(bin, n_rows, fill=sample))+
  geom_bar(stat="identity", color="black")+
  facet_wrap(~sample, ncol=1)+
  scale_y_continuous(labels = label_number(accuracy = 0.1, scale_cut = cut_si(""))) +
  scale_fill_manual(values=c("#023047","#219ebc","#8ecae6"))+
  theme_bw()+  theme(strip.background = element_rect(fill="white"),
                     legend.position = "none")+
  ylab("Number Bins")+xlab("Bin size"),

ggplot(df_summary, aes(bin, median_nFeature, fill=sample))+
  geom_bar(stat="identity", color="black")+
  facet_wrap(~sample, ncol=1)+
  scale_fill_manual(values=c("#023047","#219ebc","#8ecae6"))+
  theme_bw()+  theme(strip.background = element_rect(fill="white"),
                     legend.position = "none")+
  ylab("Median Number of Genes in a Bin")+xlab("Bin size"),

ggplot(df_summary, aes(bin, median_nCount, fill=sample))+
  geom_bar(stat="identity", color="black")+
  facet_wrap(~sample, ncol=1)+
  scale_fill_manual(values=c("#023047","#219ebc","#8ecae6"))+
  theme_bw()+  theme(strip.background = element_rect(fill="white"),
        legend.position = "none")+
  ylab("Median Number of Transcripts in a Bin")+xlab("Bin size"),
ncol=3)

save_plts(metrcis, "metrics_allsamples_option1", w=6, h=6)



metrcis<-plot_grid(
  ggplot(df_summary, aes(bin, n_rows, fill=sample))+
    geom_bar(stat="identity", color="black",position="dodge")+
    scale_y_continuous(labels = label_number(accuracy = 0.1, scale_cut = cut_si(""))) +
    scale_fill_manual(values=c("#669bbc","#588157","#edafb8"))+
    theme_bw()+  theme(strip.background = element_rect(fill="white"),
                       axis.text.x = element_blank(),
                       axis.title.x = element_blank(),
                       axis.ticks.x = element_blank())+
    ylab("Number of Bins")+xlab("Bin size (µm)"),
  
  ggplot(df_summary, aes(bin, median_nFeature, fill=sample))+
    geom_bar(stat="identity", color="black",position="dodge")+
    scale_fill_manual(values=c("#669bbc","#588157","#edafb8"))+
    theme_bw()+  theme(strip.background = element_rect(fill="white"), 
                       axis.text.x = element_blank(),
                       axis.title.x = element_blank(),
                       axis.ticks.x = element_blank())+
    ylab("Genes\n(median per bin)")+xlab("Bin size (µm)"),
  
  ggplot(df_summary, aes(bin, median_nCount, fill=sample))+
    geom_bar(stat="identity", color="black",position="dodge")+
    scale_fill_manual(values=c("#669bbc","#588157","#edafb8"))+
    theme_bw()+  theme(strip.background = element_rect(fill="white"))+
    ylab("Transcripts\n(median per bin)")+xlab("Bin size (µm)"),
  ncol=1, align="v")
metrcis

save_plts(metrcis, "metrics_allsamples_option2", w=3.5, h=7)


metrcis<-plot_grid(
  ggplot(df_summary, aes(bin, n_rows, fill=sample))+
    geom_bar(stat="identity", color="black")+
    scale_y_continuous(labels = label_number(accuracy = 0.1, scale_cut = cut_si(""))) +
    scale_fill_manual(values=c("#023047","#219ebc","#8ecae6"))+
    theme_bw()+  theme(strip.background = element_rect(fill="white"),
                       legend.position = "none")+
    ylab("Number Bins")+xlab("Bin size (µm)"),
  
  ggplot(df_summary, aes(bin, median_nFeature, fill=sample))+
    geom_bar(stat="identity", color="black")+
    scale_fill_manual(values=c("#023047","#219ebc","#8ecae6"))+
    theme_bw()+  theme(strip.background = element_rect(fill="white"),
                       legend.position = "none")+
    ylab("Median Number of\nGenes in a Bin")+xlab("Bin size (µm)"),
  
  ggplot(df_summary, aes(bin, median_nCount, fill=sample))+
    geom_bar(stat="identity", color="black")+
    scale_fill_manual(values=c("#023047","#219ebc","#8ecae6"))+
    theme_bw()+  theme(strip.background = element_rect(fill="white"),
                       legend.position = "none")+
    ylab("Median Number of\nTranscripts in a Bin")+xlab("Bin size (µm)"),
  ncol=1, align="v")
metrcis

save_plts(metrcis, "metrics_allsamples_option2", w=5, h=7)
