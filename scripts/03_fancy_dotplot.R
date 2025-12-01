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

source("scripts/00_pretty_plots.R")



annotation_C107_8um<-read.csv(here("data/metadata_for_C107.csv"))

localdir <- "/media/redgar/Seagate Portable Drive/visiumHD_liver/reseq/MacParland_Diana__C107_D1/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8))


####################################
## Fancy dot plot
####################################


## markers
hep_genes<-c("SERPINA1", "CYP2E1","CYP1A2","HAL")
Endothelial_genes<-c( "CD34","LYVE1","TAGLN")
act_mes_genes<-c( "LAMC3","COL1A1")

MHCII_genes<-c("CD163")
mono_genes<-c( "S100A8","S100A9","LYZ")
KC_genes<-c( "MARCO","VSIG4","CD74")

immune_gene<-c("PTPRC")
T_genes<-c( "CD3D","CD8A")
NK_genes<-c( "GNLY","NKG7")

chol_genes<-c( "EPCAM")

neutro_genes<-c("FCGR3A", "DEFA3")

plasma_genes<-c( "IGHG1")
B_genes<-c( "MS4A1")

cycle_genes<-c( "HBA","HBB","MKI67",  "TOP2A")



gene_exp<-FetchData(object, vars=c(hep_genes,Endothelial_genes,act_mes_genes,chol_genes,
                                       immune_gene,
                                       T_genes,NK_genes,B_genes,plasma_genes,neutro_genes,
                                       mono_genes,MHCII_genes, KC_genes,cycle_genes))
gene_exp$cell<-rownames(gene_exp)
gene_exp<-melt(gene_exp)

rm(object)
gc()

plt<-merge(gene_exp, annotation_C107_8um, by.x="cell", by.y="X")
gc()

## summarize
scale_this <- function(x){(x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}

plt_summary<-plt %>% group_by(final_anno, variable) %>% 
  summarise(mn=mean(value), count=length(which(value>0)), percent_exp=(length(which(value>0))/length(value))*100)
plt_summary <- plt_summary %>% group_by(variable) %>%
  dplyr::mutate(scaled = scale_this(mn))
plt_summary<-as.data.frame(plt_summary)

# remove dots where 0 cell expressing marker
plt_summary<-plt_summary[(which(plt_summary$count>0)),]

plt_summary$variable<-factor(plt_summary$variable, levels=rev(c(hep_genes,Endothelial_genes,act_mes_genes,chol_genes,
                                                                immune_gene,
                                                                T_genes,NK_genes,B_genes,plasma_genes,neutro_genes,
                                                                mono_genes,MHCII_genes, KC_genes,cycle_genes)))

plt_summary$final_anno<-factor(plt_summary$final_anno, levels=c(
  "bin 1 hepatocyte", "bin 2 hepatocyte", "bin 3 hepatocyte", "bin 4 hepatocyte", "bin 5 hepatocyte", "bin 6 hepatocyte", "bin 7 hepatocyte",
  "bin 8 hepatocyte", "bin 9 hepatocyte",
  "endothelial cell of pericentral hepatic sinusoid","endothelial cell of periportal hepatic sinusoid",
  "vein endothelial cell","endothelial cell of artery",
  "fibroblast","hepatic stellate cell",
  "intrahepatic cholangiocyte",
  "CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell",
  "T cell","natural killer cell","hepatic pit cell",
  "mature B cell","plasma cell", 
  "neutrophil", "neutrophil2",
  "monocyte","macrophage","Kupffer cell","mast cell","conventional dendritic cell", "plasmacytoid dendritic cell",
  "erythrocyte", "Low_UMI","Unannotated","unknown"))

levels(plt_summary$final_anno)<-gsub("bin","layer", levels(plt_summary$final_anno))





gene_list_len<-sapply(list(hep_genes,Endothelial_genes,act_mes_genes,chol_genes,
                           immune_gene,
                           T_genes,NK_genes,B_genes,plasma_genes,neutro_genes,
                           mono_genes,MHCII_genes, KC_genes,cycle_genes), function(x) length(x))    



fancy_dotplot<-plot_grid(
  ggplot(plt_summary, aes(final_anno, variable, color=scaled, size=percent_exp))+geom_point()+
    theme_classic()+
    scale_color_continuous_sequential(palette = "Oslo", rev=F, name="Scaled\nMean\nExpression")+
    scale_size(name="Percent\nCells\nExpressing")+
    theme(axis.text.x = element_blank(),axis.title = element_blank(),axis.ticks.x = element_blank())+
    geom_hline(yintercept =c(2,3,10,12,14,18,19,20,22,25)+0.5, color="grey70")+
    geom_vline(xintercept = c(9,13,15,16,21,23,25,31,32)+0.5, color="grey70"),
  ggplot(plt_summary, aes(final_anno, y=1, fill=final_anno))+geom_tile(color="black")+
    theme_classic()+fillscale_cellType+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
          legend.position = "none",axis.line  = element_blank(),
          plot.margin = margin(t = 0,  # Top margin
                               r = 50,  # Right margin
                               b = 40,  # Bottom margin
                               l = 10)),
  ncol=1, rel_heights = c(6,2.5), align = "v", axis="lr")
fancy_dotplot

save_plts(fancy_dotplot, "C107_dot_plot_final_anno", w=9,h=12)

