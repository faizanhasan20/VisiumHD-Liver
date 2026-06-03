library(shiny)
library(ggplot2)
library(viridis)
library(reshape2)
library(cowplot)

load("/media/redgar/Seagate Portable Drive/xenium_liver/shiny_minimum_data.RData")
source("/home/redgar/Documents/xenium_liver/scripts/00_pretty_plots.R")


gene<-"CYP1A2"

percentile_count <- as.numeric(0)

shiny_min_data_df_C95<-shiny_min_data_df[which(shiny_min_data_df$sample=="Adult (C95)"),]

gene_percentile_count <- quantile(shiny_min_data_df_C95[shiny_min_data_df_C95[[gene]] > 0, gene], probs = percentile_count, na.rm = TRUE)

exp_plot<-ggplot() +
  geom_point(aes(centroid_x, -centroid_y), shiny_min_data_df_C95, color = "black", size = 1.5, shape = 19) +
  geom_point(aes(centroid_x, -centroid_y), shiny_min_data_df_C95, color = "grey95", size = 0.2, shape = 19) +
  geom_point(aes(centroid_x, -centroid_y, color = (as.numeric(!!sym(gene)))), 
             shiny_min_data_df_C95[shiny_min_data_df_C95[[gene]] > gene_percentile_count, ], 
             size = 0.2, shape = 19) +theme_void()+
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), data = xenium_scale_bar_df[which(xenium_scale_bar_df$sample=="Adult (C95)"),]) +
  geom_text(aes(x = xmax - 250, y = ymin * 1.025), data = xenium_scale_bar_df[which(xenium_scale_bar_df$sample=="Adult (C95)"),], label = "0.5mm", size = 4) +
  # geom_rect(aes(xmin = min(shiny_min_data_df_C95$centroid_x), xmax = max(shiny_min_data_df_C95$centroid_x),
  #               ymin = -min(shiny_min_data_df_C95$centroid_y), ymax = -max(shiny_min_data_df_C95$centroid_y)),
  #           fill = NA,color = "black",linewidth = 1) +
  coord_fixed() +  scale_color_viridis(na.value = "grey90", name="Expression") 
exp_plot
save_plts(exp_plot, "C95_CYP1A2_xenium", w=10,h=6)
