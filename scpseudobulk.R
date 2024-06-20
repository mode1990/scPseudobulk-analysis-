library(Seurat)
library(ggpubr)
library(factoextra)
library(FactoMineR)

# Load data
load_data <- function(file_path) {
  return(readRDS(file_path))
}

# Normalize and prepare data
prepare_data <- function(seurat_obj, cell_types, assay_type = "RNA") {
  Idents(seurat_obj) <- seurat_obj$CellType
  seurat_obj <- subset(seurat_obj, idents = cell_types)
  DefaultAssay(seurat_obj) <- assay_type
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000, assay = assay_type)
  return(seurat_obj)
}

# Average expression to create pseudobulk
create_pseudobulk <- function(seurat_obj, group_by = 'Clone') {
  s.bulk <- AverageExpression(seurat_obj, return.seurat = TRUE, group.by = group_by)
  s.bulk <- FindVariableFeatures(s.bulk)
  s.bulk <- ScaleData(s.bulk)
  s.bulk <- RunPCA(s.bulk, features = VariableFeatures(object = s.bulk), npcs = 10)
  return(s.bulk)
}

# Save and load pseudobulk
save_pseudobulk <- function(s.bulk, file_path) {
  saveRDS(s.bulk, file_path)
}

load_pseudobulk <- function(file_path) {
  return(readRDS(file_path))
}

# Perform K-means clustering
perform_kmeans <- function(s.bulk, assay = "RNA", num_clusters = 5) {
  exp <- GetAssayData(object = s.bulk, assay = assay)
  expression_matrix <- t(as.matrix(exp))
  kmeans_result <- kmeans(expression_matrix, centers = num_clusters)
  return(kmeans_result)
}

# Plot K-means clustering result
plot_kmeans <- function(kmeans_result, expression_matrix) {
  return(fviz_cluster(kmeans_result, data = expression_matrix, stand = FALSE))
}

# Calculate and plot heatmap
plot_heatmap <- function(s.bulk, assay = "RNA") {
  distance_matrix <- dist(t(s.bulk@assays[[assay]]@scale.data), method = "euclidean")
  correlation_matrix <- 1 - as.matrix(distance_matrix)
  heatmap(correlation_matrix, 
          Colv = TRUE,   
          Rowv = TRUE,   
          col = colorRampPalette(c("blue", "white", "red"))(100),
          scale = "none",
          main = "Correlation Heatmap of RNA data")
}

# Run PCA and print variance
run_pca <- function(correlation_matrix) {
  pca_result <- prcomp(correlation_matrix, scale. = TRUE)
  pc_scores <- pca_result$x
  variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  cat("Variance explained by PC1:", variance_explained[1], "\n")
  cat("Variance explained by PC2:", variance_explained[2], "\n")
  return(list(pc_scores = pc_scores, variance_explained = variance_explained))
}
