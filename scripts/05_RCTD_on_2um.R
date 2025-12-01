
library(data.table)
library(spacexr)
library(Seurat)


#gene map will be used to convert gene symbols to ensemble IDs for comparison to single cell
gene_map <- fread("/scratch/redgar25/VisiumHD/MacParland_Diana__C107_D1/outs/binned_outputs/square_002um/filtered_feature_bc_matrix/features.tsv.gz")

colnames(gene_map) <- c("ensembl_id", "gene_symbol", "feature_type")

localdir <- "/scratch/redgar25/VisiumHD/MacParland_Diana__C107_D1/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(2))

DefaultAssay(object) <- "Spatial"

## coordiantes (CV)
xmin<-7700
xmax<-8090
ymin<-7300
ymax<-7775

plt_hack<-cbind(object@meta.data, object@images$slice1$centroids@coords)
zoom<-plt_hack[which(plt_hack$x>xmin & plt_hack$x<xmax & plt_hack$y>ymin & plt_hack$y<ymax),]

object$cell<-colnames(object)
zoom_seurat <- subset(object, subset = cell %in% rownames(zoom))


## single-cell reference from CELLxGENE
ref <- LoadSeuratRds("/scratch/redgar25/VisiumHD/b09dac61-836f-418a-ad0a-064d2bc5343d.rds")

Idents(ref) <- "cell_type"
counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$cell_type)

reference <- Reference(counts, cluster)

counts_hd <- zoom_seurat[["Spatial.002um"]]$count

coords <- GetTissueCoordinates(zoom_seurat)[, 1:2]

query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

symbol_to_ensembl <- setNames(gene_map$ensembl_id, gene_map$gene_symbol)

rownames(query@counts) <- symbol_to_ensembl[rownames(query@counts)]

RCTD <- create.RCTD(query, reference, max_cores = 40, UMI_min = 20, UMI_min_sigma = 100)

RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

zoom_seurat <- AddMetaData(zoom_seurat, metadata = RCTD@results$results_df)

saveRDS(zoom_seurat,"/scratch/redgar25/VisiumHD/zoomCV_2um_C107_RCTD.rds")
