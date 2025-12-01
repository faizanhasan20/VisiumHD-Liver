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


############# on trillium
annotation_C107_8um<-read.csv(here("/scratch/redgar25/VisiumHD/metadata_for_C107.csv"))

localdir <- "/scratch/redgar25/VisiumHD/MacParland_Diana__C107_D1/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8))

head(object)
object<-AddMetaData(object,annotation_C107_8um)
object<-JoinLayers(object)
object<-NormalizeData(object)
gc()

Idents(object) <- "zones"
zone_markers <- FindAllMarkers(object, logfc.threshold = 0, min.pct=0)
save(zone_markers, file = here("/scratch/redgar25/VisiumHD/zone_markers.RData"))




##################
load(here("data/zone_markers.RData"))

zone_markers$cluster<-as.factor(zone_markers$cluster)
levels(zone_markers$cluster)<-gsub("bin | hepatocyte","",levels(zone_markers$cluster))

zone_markers<-zone_markers[which(zone_markers$cluster!="NA"),]

plot_gene_foldchange <- function(gene_symbol, markers_df) {
  # Filter for the specified gene and exclude zone 9
  gene_data <- markers_df %>% filter(gene == gene_symbol)
  
  if (nrow(gene_data) == 0) {
    message("No data found for gene: ", gene_symbol)
    return(NULL)
  }
  
  # Ensure cluster is numeric and ordered
  gene_data$cluster <- as.numeric(as.character(gene_data$cluster))
  gene_data <- gene_data %>% arrange(cluster)
  
  # Determine the y-axis limits to span all values with padding
  y_min <- min(gene_data$avg_log2FC, na.rm = TRUE)
  y_max <- max(gene_data$avg_log2FC, na.rm = TRUE)
  y_pad <- 0.5
  
  ggplot(gene_data, aes(x = cluster, y = avg_log2FC)) +
    # Green zone above 1
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1, ymax = y_max + y_pad,
             fill = "lightgreen", alpha = 0.5) +
    # Green zone below -1
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = y_min - y_pad, ymax = -1,
             fill = "lightgreen", alpha = 0.5) +
    geom_line(color = "black", size = 1) +
    geom_point(size = 2, color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_x_continuous(breaks = sort(unique(gene_data$cluster))) +
    labs(
      x = "Zonation Zone",
      y = "log2 Fold Change") +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14)
    )
}


plot_gene_foldchange("CPS1",zone_markers)
plot_gene_foldchange("SLCO1B3",zone_markers)

plot_gene_foldchange("LGALS4",zone_markers)
zone_markers[grep("LGALS4", zone_markers$gene),]

## sig markers to plot
sig_zonated<-zone_markers[which(zone_markers$p_val_adj<0.005 & abs(zone_markers$avg_log2FC)>1),]

## remove those in signature
zone_genes <- c("CYP3A4", "ADH1B", "CYP1A2", "CYP2E1", "APOA2", "APOC1",
                       "ADH4", "ADH1A", "APOH", "AMBP", "GSTA2", "ADH1C", "SLCO1B3",
                       "AOX1", "APOA5", "DCXR", "RBP4", "OAT", "CYP2C19", "GC",
                "SERPINA1", "APOA1", "ALB", "C7", "NNMT", "HAMP", "ALDOB",
                      "ASS1", "CYP2A7", "MGP", "A2M", "FXYD2", "CCL21", "HAL",
                      "IGFBP2", "SDS", "AQP1", "CYP2A6", "FBLN1", "PTGDS")

sig_zonated_plot<-sig_zonated[which(!(sig_zonated$gene%in%zone_genes)),]

# sig across zones
table(sig_zonated_plot$gene)[which(table(sig_zonated_plot$gene)>6)]

plot_gene_foldchange("GLUL",zone_markers) # present as validating a known marker



# "new" markers to highlight
table(sig_zonated_plot$gene)[which(table(sig_zonated_plot$gene)>2)]

plot_gene_foldchange("CXCL13",zone_markers) # on protein atlas staining unclear
plot_gene_foldchange("CHI3L1",zone_markers) # not on protein atlas
plot_gene_foldchange("HSD17B13",zone_markers) # on protein atlas staining unclear
plot_gene_foldchange("LCN2",zone_markers) # on protein atlas staining unclear
plot_gene_foldchange("SDSL",zone_markers) # Maybe some staining in zonation
plot_gene_foldchange("UPP2",zone_markers) # not on protein atlas

plot_gene_foldchange("MUC6",zone_markers) # not on protein atlas

fc_plot<-plot_grid(plot_gene_foldchange("SLCO1B3",zone_markers), plot_gene_foldchange("LGALS4",zone_markers), ncol=1)
save_plts(fc_plot, "gene_fold_change", w=4, h=6)
