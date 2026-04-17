library(data.table)
library(spacexr)
library(Seurat)

library(ggplot2)
library(cowplot)
library(here)

source("scripts/00_pretty_plots.R")



#gene map will be used to convert gene symbols to ensemble IDs for comparison to single cell
gene_map <- fread("C113_only/binned_outputs/square_008um/filtered_feature_bc_matrix/features.tsv.gz")

colnames(gene_map) <- c("ensembl_id", "gene_symbol", "feature_type")

localdir <- "C113_only"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8))

DefaultAssay(object) <- "Spatial"
SpatialFeaturePlot(object, features = "nCount_Spatial.008um")



ref <- LoadSeuratRds("b09dac61-836f-418a-ad0a-064d2bc5343d.rds")

Idents(ref) <- "cell_type"
counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$cell_type)

reference <- Reference(counts, cluster)

counts_hd <- object[["Spatial.008um"]]$count

coords <- GetTissueCoordinates(object)[, 1:2]

query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

symbol_to_ensembl <- setNames(gene_map$ensembl_id, gene_map$gene_symbol)

rownames(query@counts) <- symbol_to_ensembl[rownames(query@counts)]

RCTD <- create.RCTD(query, reference, max_cores = 40)

RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

object <- AddMetaData(object, metadata = RCTD@results$results_df)

saveRDS(object,"full_C113_RCTD.rds")













#gene map will be used to convert gene symbols to ensemble IDs for comparison to single cell
gene_map <- fread("C95/binned_outputs/square_008um/filtered_feature_bc_matrix/features.tsv.gz")

colnames(gene_map) <- c("ensembl_id", "gene_symbol", "feature_type")

localdir <- "C95"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8))

DefaultAssay(object) <- "Spatial"
#SpatialFeaturePlot(object, features = "nCount_Spatial.008um")
print(dim(object))


ref <- LoadSeuratRds("b09dac61-836f-418a-ad0a-064d2bc5343d.rds")
print(dim(ref))


Idents(ref) <- "cell_type"
counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$cell_type)

reference <- Reference(counts, cluster)

counts_hd <- object[["Spatial.008um"]]$count

coords <- GetTissueCoordinates(object)[, 1:2]

query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

symbol_to_ensembl <- setNames(gene_map$ensembl_id, gene_map$gene_symbol)

rownames(query@counts) <- symbol_to_ensembl[rownames(query@counts)]

RCTD <- create.RCTD(query, reference, max_cores = 40)

RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

object <- AddMetaData(object, metadata = RCTD@results$results_df)

saveRDS(object,"full_C95_RCTD.rds")








#############
## correlation of cell type proportions of good samples
#############
C95_object<-readRDS("full_C95_RCTD.rds")
C95_object$first_type[grep("hepatocyte",C95_object$first_type)]<-"hepatocyte"
C95_annotation_proportions<-as.data.frame(table(C95_object$first_type))
C95_annotation_proportions
C95_annotation_proportions$Freq<-(C95_annotation_proportions$Freq/sum(C95_annotation_proportions$Freq))*100
colnames(C95_annotation_proportions)<-c("CellType","C95")

C113_object<-readRDS("full_C113_RCTD.rds")
C113_object$first_type[grep("hepatocyte",C113_object$first_type)]<-"hepatocyte"
C113_annotation_proportions<-as.data.frame(table(C113_object$first_type))
C113_annotation_proportions$Freq<-(C113_annotation_proportions$Freq/sum(C113_annotation_proportions$Freq))*100
colnames(C113_annotation_proportions)<-c("CellType","C113")

C107_annotation<-read.csv("annotation_C107_8um_with_banksy.csv")
C107_annotation$final_anno[grep("hepatocyte",C107_annotation$final_anno)]<-"hepatocyte"
C107_annotation_proportions<-as.data.frame(table(C107_annotation$final_anno))
C107_annotation_proportions$Freq<-(C107_annotation_proportions$Freq/sum(C107_annotation_proportions$Freq))*100
colnames(C107_annotation_proportions)<-c("CellType","C107")

porportion_plt<-merge(C107_annotation_proportions, C95_annotation_proportions, by="CellType")
porportion_plt<-merge(porportion_plt, C113_annotation_proportions, by="CellType")


save(porportion_plt, file="proportions_all_samples.RData")


C95_object<-readRDS("full_2024_C95_RCTD.rds")
write.csv(C95_object@meta.data, file="annotation_C95_8um_RCTD.csv")

C95<-read.csv("data/annotation_C95_8um_RCTD.csv")
(nrow(C95)-length(which(is.na(C95$first_type))))/nrow(C95)

C113_object<-readRDS("full_C113_RCTD.rds")
write.csv(C113_object@meta.data, file="annotation_C113_8um_RCTD.csv")

C113<-read.csv("data/annotation_C113_8um_RCTD.csv")
(nrow(C113)-length(which(is.na(C113$first_type))))/nrow(C113)

C107_annotation<-read.csv("data/annotation_C107_8um_with_banksy.csv")
(nrow(C107_annotation)-length(which(C107_annotation$spot_class=="Low_UMI")))/nrow(C107_annotation)

table(C107_annotation$spot_class)


################
## plot correlation
################
load("data/proportions_all_samples.RData")

leg_celltype<-get_leg(ggplot(porportion_plt, aes(C107, C95))+geom_point(aes(fill=CellType),shape=21)+fillscale_cellType+theme_bw())

ct <- cor.test( porportion_plt$C107[porportion_plt$C107 > 3],  porportion_plt$C95[porportion_plt$C107 > 3])
many_C95<-ggplot(porportion_plt[which(porportion_plt$C107>3),], aes(C107, C95))+  stat_smooth(method="lm", col="grey", se=F)+
  geom_point(aes(fill=CellType),shape=21)+fillscale_cellType+theme_bw()+
  theme(legend.position = "none")+xlab("C107 (percent of bins)")+ylab("C95 (percent of bins)")+
  annotate(  "text",  x = 40,  y = 10,  label = sprintf("r = %.3f\np = %.2g", ct$estimate, ct$p.value),  hjust = 0,  size = 2.5, color="grey40")
many_C95

ct <- cor.test( porportion_plt$C107[porportion_plt$C107 < 3],  porportion_plt$C95[porportion_plt$C107 < 3])
few_C95<-ggplot(porportion_plt[which(porportion_plt$C107<3),], aes(C107, C95))+stat_smooth(method="lm", col="grey", se=F)+
  geom_point(aes(fill=CellType),shape=21)+fillscale_cellType+theme_bw()+
  theme(legend.position = "none")+xlab("C107 (percent of bins)")+ylab("C95 (percent of bins)")+
  annotate(  "text",  x = 0.3,  y = 0.12,  label = sprintf("r = %.3f\np = %.2g", ct$estimate, ct$p.value),  hjust = 0,  size = 2.5, color="grey40")
few_C95


ct <- cor.test( porportion_plt$C107[porportion_plt$C107 > 3],  porportion_plt$C113[porportion_plt$C107 > 3])
many_C113<-ggplot(porportion_plt[which(porportion_plt$C107>3),], aes(C107, C113))+  stat_smooth(method="lm", col="grey", se=F)+
  geom_point(aes(fill=CellType),shape=21)+fillscale_cellType+theme_bw()+
  theme(legend.position = "none")+xlab("C107 (percent of bins)")+ylab("C113 (percent of bins)")+
  annotate(  "text",  x = 45,  y = 8,  label = sprintf("r = %.3f\np = %.2g", ct$estimate, ct$p.value),  hjust = 0,  size = 2.5, color="grey40")
many_C113

ct <- cor.test( porportion_plt$C107[porportion_plt$C107 < 3],  porportion_plt$C113[porportion_plt$C107 < 3])
few_C113<-ggplot(porportion_plt[which(porportion_plt$C107<3),], aes(C107, C113))+stat_smooth(method="lm", col="grey", se=F)+
  geom_point(aes(fill=CellType),shape=21)+fillscale_cellType+theme_bw()+
  theme(legend.position = "none")+ xlab("C107 (percent of bins)")+ylab("C113 (percent of bins)")+
  annotate(  "text",  x = 0.3,  y = 0.15,  label = sprintf("r = %.3f\np = %.2g", ct$estimate, ct$p.value),  hjust = 0,  size = 2.5, color="grey40")
few_C113

plot_grid(few_C95, many_C95, few_C113,many_C113)

save_plts(plot_grid(few_C95, many_C95, few_C113,many_C113), "proportion_correlation", w=5,h=5)






