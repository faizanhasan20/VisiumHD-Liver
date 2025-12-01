save_plts<-function(plt, name, w,h){
  ggsave(plt, file=paste(here("figures/"),name,".pdf", sep=""), w=w, h=h)
  ggsave(plt, file=paste(here("figures/jpeg/"),name,".jpeg", sep=""), w=w, h=h, bg = "white")
  ggsave(plt, file=paste(here("figures/png/"),name,".png", sep=""), w=w, h=h, bg = "white")}

save_plts_black<-function(plt, name, w,h){
  ggsave(plt, file=paste(here("figures/"),name,".pdf", sep=""), w=w, h=h, bg = "black")
  ggsave(plt, file=paste(here("figures/jpeg/"),name,".jpeg", sep=""), w=w, h=h, bg = "black")
  ggsave(plt, file=paste(here("figures/png/"),name,".png", sep=""), w=w, h=h, bg = "black")}

theme_presentation<- function(base_size = 15, base_family = "") {
  # Starts with theme_grey and then modify some parts
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      strip.background = element_blank(),
      strip.text = element_blank(),
      axis.text = element_blank(),
      axis.ticks =  element_blank(), 
      axis.title= element_blank(),
      panel.background = element_rect(fill="black"), 
      panel.border =element_blank(),  
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.margin = unit(1.0, "lines"), 
      plot.background = element_rect(fill="black"), 
      plot.title =element_text(size=15,colour="white"), 
      plot.margin = unit(c(1,  1, 1, 1), "lines"),
      legend.background=element_rect(fill='black'),
      legend.title=element_text(size=15,colour="white"),
      legend.text=element_text(size=15,colour="white"),
      legend.key = element_rect( fill = 'black'),
      legend.key.size = unit(c(1, 1), "lines"),
      axis.line = element_blank()
    )
}

th_present <- theme(axis.text=element_text(size=12),
                    axis.title=element_text(size=14),
                    strip.text = element_text(size = 12),
                    legend.text=element_text(size=12),
                    legend.title=element_text(size=14))

th <- theme(axis.text=element_text(size=12),
                    axis.title=element_text(size=14),
                    strip.text = element_text(size = 12),
                    legend.text=element_text(size=12),
                    legend.title=element_text(size=14))


## grab legened from plot
get_leg = function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}


######### 
## Cell Type Color
######### 

myColors_celltype <- c("#cccc16",
                       
                       "#c2884a","#f4a259","#a63603","#7d4f0e",
                       
                       "#0a15f2","#35cdde",
                       "#5dd45d","#89fc00", "#283618",
                       "#8e11bf","#67038f",  
                       "#d1100d","#058205",  "#155e15", 
                       "#72b01d","#04e762", 
                       "#5e1525",
                       "#ad235c","#f781b0","#ed0ce2","#f54082","#FF00A0",
                       "grey90","grey60","grey60",
                       rev(c("#184E77", "#1E6091","#1A759F", "#168AAD", "#34A0A4", "#52B69A", "#76C893", "#99D98C" ,"#B5E48C")),
                       "#184E77","#34A0A4","#B5E48C","#168AAD",
                       "black")  



color_possibilities_celltype<-c("intrahepatic cholangiocyte",
                                
                                "endothelial cell of pericentral hepatic sinusoid","endothelial cell of periportal hepatic sinusoid",
                                "vein endothelial cell","endothelial cell of artery",
                                
                                "hepatic stellate cell","fibroblast",
                                
                                "T cell","natural killer cell","hepatic pit cell",
                                "mature B cell","plasma cell", 
                                "erythrocyte","neutrophil", "neutrophil2",
                                "CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell",
                                "plasmacytoid dendritic cell",
                                "macrophage","Kupffer cell","monocyte","mast cell","conventional dendritic cell",
                                
                                "Low_UMI","Unannotated","unknown",
                                "layer 1 hepatocyte", "layer 2 hepatocyte", "layer 3 hepatocyte", "layer 4 hepatocyte", "layer 5 hepatocyte", "layer 6 hepatocyte", "layer 7 hepatocyte",
                                "layer 8 hepatocyte", "layer 9 hepatocyte",
                                #for 2um RCTD straight mapping
                                "periportal region hepatocyte","midzonal region hepatocyte","centrilobular region hepatocyte","hepatocyte",
                                "NRXN1+ HSC")  


names(myColors_celltype) <- color_possibilities_celltype
fillscale_cellType <- scale_fill_manual(name="Cell Type",
                                        values = myColors_celltype, drop = T, limits=force,na.value = "grey90")
colscale_cellType <- scale_color_manual(name="Cell Type",
                                        values = myColors_celltype, drop = T, limits=force,na.value = "grey90")






       

######### 
## highlight cell type
######### 
highlight_cell<-function(fav_celltype){
  cell_types<-unique(plot_freely$CellType)
  myColors_celltype <- c(rep("grey20", length(cell_types)))
  myColors_celltype[which(cell_types==fav_celltype)]<-"green"

  names(myColors_celltype) <- cell_types
  scale_color_manual(name="Cell Type",values = myColors_celltype, drop = T, limits=force)}


#################
## Density of a cell type
#################
cell_density<-function(bin_number,min_cell_in_bin ,subset){
  
  proportion_overall<-length(plot_freely$x[which(plot_freely$CellType==subset)])/nrow(plot_freely)
  print(paste(round(proportion_overall,3), subset, "across whole section"))
  
  # Bin width
  bin_width_x <- max(plot_freely$x)/bin_number
  bin_width_y <- max(plot_freely$y)/bin_number
  
  # Create binning grid
  x_bins <- seq(0, max(plot_freely$x), by = bin_width_x)
  y_bins <- seq(0, max(plot_freely$y), by = bin_width_y)
  bins <- list(x = x_bins, y = y_bins)
  
  # Cut points into bins
  x_cut <- cut(plot_freely$x, breaks = x_bins, include.lowest = TRUE)
  y_cut <- cut(plot_freely$y, breaks = y_bins, include.lowest = TRUE)
  
  # Create a table of bin counts
  point_counts_all <- table(x_cut, y_cut)
  point_counts_all[point_counts_all < min_cell_in_bin] <- NA
  
  # Cut points into bins  subset
  x_cut_subset <- cut(plot_freely$x[which(plot_freely$CellType==subset)], breaks = x_bins, include.lowest = TRUE)
  y_cut_subset <- cut(plot_freely$y[which(plot_freely$CellType==subset)], breaks = y_bins, include.lowest = TRUE)
  
  point_counts_subset <- table(x_cut_subset, y_cut_subset)
  
  # Display the point counts in each bin
  density_difference<-as.matrix(point_counts_subset) / as.matrix(point_counts_all)
  colnames(density_difference)<-1:bin_number
  rownames(density_difference)<-1:bin_number
  density_plt<-as.data.frame(density_difference)
  
  
  ggplot(density_plt, aes(x=x_cut_subset, y=rev(y_cut_subset), fill=Freq))+geom_tile( color="black")+
    theme(
      axis.text = element_blank(),
      axis.ticks =  element_blank(), 
      axis.title= element_blank(),
      panel.border =element_blank(),  
      legend.title=element_text(size=15,colour="black"),
      legend.text=element_text(size=15,colour="black")
    ) + 
    scale_fill_gradientn(
      colours=c("#0C6291", "#ebe8e8", "#A63446"), 
      na.value = "grey60", name = "Freq",
      values = c(0,  proportion_overall, max(density_plt$Freq, na.rm=T))/max(density_plt$Freq, na.rm=T)
    )
  
}


