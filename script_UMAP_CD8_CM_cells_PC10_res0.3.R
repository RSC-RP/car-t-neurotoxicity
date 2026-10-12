# title : Script to generate UMAPs and FeaturePlots for CD8 central memory cells. 
# author: Neerja Katiyar

library(Seurat)
library(dplyr)
library(ggplot2)

seur <- readRDS(file = "./objects/PLAT_annotated_celltypes_cellbender_new_idents_chk_viewmastR.rds")

seur_CD8_CM <- subset(seur, idents = "20_CD8.CM")
sample_pre <- sub("^(.*?_.*?)_.*", "\\1", seur_CD8_CM$orig.ident)
seur_CD8_CM$indiv <- sample_pre

samples_to_remove <- c("PLAT05_S043", "PLAT05_S055")
seurat_obj <- subset(seur_CD8_CM, subset = !indiv %in% samples_to_remove)

####################################################################################
# Additional helpful columns
seurat_obj$time_point <- factor(seurat_obj$time_point, levels = c("preA", "Day7", "Day28"))

# Normalize
seurat_obj <- NormalizeData(seurat_obj)
# Find variable features
seurat_obj <- FindVariableFeatures(seurat_obj,
				selection.method = "vst",
				nfeatures = 2000)
# Scale data
seurat_obj <- ScaleData(seurat_obj)

# Run PCA
seurat_obj <- RunPCA(seurat_obj,features = VariableFeatures(seurat_obj))

pdf("./results/ElbowPlot_CD8_CM_cells.pdf")
ElbowPlot(seurat_obj)
dev.off()

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.3)

seurat_obj <- RunUMAP(seurat_obj, dims = 1:15)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
pdf("./results/DimPlot_CD8_CMcells_PC15_res0.3.pdf")
DimPlot(seurat_obj, reduction = "umap")
dev.off()

pdf("./results/DimPlot_CD8_CMcells_clusters_PC15_res0.3.pdf")
DimPlot(seurat_obj,
	reduction = "umap",
	group.by = "seurat_clusters",
	label = TRUE,
	repel = TRUE
	) +
  ggtitle("CD8 CM Clusters")
dev.off()

pdf("./results/DimPlot_CD8_CMcells_category_PC15_res0.3.pdf")
DimPlot(seurat_obj,
	reduction = "umap",
	group.by = "category",
	label = TRUE,
	repel = TRUE
	) +
  ggtitle("CD8 CM clusters with category")
dev.off()

pdf("./results/DimPlot_CD8_CMcells_indiv_PC15_res0.3.pdf")
DimPlot(seurat_obj,
	reduction = "umap",
	group.by = "indiv",
	label = TRUE,
	repel = TRUE
	) +
  ggtitle("CD8 CM clusters with individuals")
dev.off()

pdf("./results/FeaturePlot_CD8_CMcells_IL10_category_all_PC15_res0.3.pdf")
FeaturePlot(
	seurat_obj,
	features = "IL10",
	reduction = "umap",
	order = T
	)
dev.off()

pdf("./results/FeaturePlot_CD8_CMcells_IL10_indiv_all_PC15_res0.3.pdf", width=30)
FeaturePlot(
	seurat_obj,
	features = "IL10",
	reduction = "umap",
	split.by = "indiv",
	order = T
	)
dev.off()

saveRDS(seurat_obj, "./results/seur_CD8_CMcells_only_PC15_res0.3.rds")

seurat_obj <- readRDS("./results/seur_CD8_CMcells_only_PC15_res0.3.rds")

pdf("./results/FeaturePlot_CD8_CMcells_EGFR_all_PC15_res0.3.pdf")
FeaturePlot(
	seurat_obj,
	features = "EGFR",
	reduction = "umap",
	order = T
	)
dev.off()

pdf("./results/FeaturePlot_CD8_CMcells_IL10_all_PC15_res0.3.pdf")
FeaturePlot(
	seurat_obj,
	features = "IL10",
	reduction = "umap",
	order = T
	)
dev.off()


