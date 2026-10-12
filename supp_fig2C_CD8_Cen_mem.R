
library(Seurat)

seur_cm_severe_subset <- readRDS(file = "../CD8_CM_IL10_clustering_analysis/seur_cm_severe_subset_tutorial_dims10_res0.2.rds")

##############-----------####################
pdf("Featureplot_CD8_CM_IL10_dim10_res0.2.pdf")
FeaturePlot(seur_cm_severe_subset, feature = "IL10", reduction = "umap") +
	theme(plot.title = element_text(size = 16),        # Plot title
	axis.title.x = element_text(size = 16),      # X-axis title
	axis.title.y = element_text(size = 16),      # Y-axis title
	axis.text.x = element_text(size = 16),       # Y-axis labels
	axis.text.y = element_text(size = 16),       # Y-axis labels
	legend.title = element_text(size = 16),      # Legend title
	legend.text = element_text(size = 16))
dev.off()

########################3--------#############
pdf("DimPlot_seur_CD8_cm_dim10_res0.2_paper.pdf")
DimPlot(seur_cm_severe_subset, reduction = "umap") +
	theme(plot.title = element_text(size = 16),        # Plot title
	axis.title.x = element_text(size = 16),      # X-axis title
	axis.title.y = element_text(size = 16),      # Y-axis title
	axis.text.x = element_text(size = 16),       # Y-axis labels
	axis.text.y = element_text(size = 16),       # Y-axis labels
	legend.title = element_text(size = 16),      # Legend title
	legend.text = element_text(size = 16))
dev.off()

#################3-------#######################
pdf("Featureplot_CD8_CM_EGFR_dim10_res0.2.pdf")
FeaturePlot(seur_cm_severe_subset, feature = "EGFR", reduction = "umap") +
	theme(plot.title = element_text(size = 16),        # Plot title
	axis.title.x = element_text(size = 16),      # X-axis title
	axis.title.y = element_text(size = 16),      # Y-axis title
	axis.text.x = element_text(size = 16),       # Y-axis labels
	axis.text.y = element_text(size = 16),       # Y-axis labels
	legend.title = element_text(size = 16),      # Legend title
	legend.text = element_text(size = 16))
dev.off()
