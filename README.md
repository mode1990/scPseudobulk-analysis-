
# 🧬 scPseudobulk-analysis

**scPseudobulk-analysis** is an automated R pipeline designed for generating pseudobulk expression profiles from single-cell RNA-seq data and performing downstream unsupervised analyses using [Seurat v5](https://satijalab.org/seurat/). It includes steps such as data loading, normalization, pseudobulk aggregation, clustering, heatmap visualization, and principal component analysis (PCA).

---

## 📦 Features

- Easy integration with Seurat v5 objects (`.rds`)
- Customizable pseudobulk aggregation (e.g., by clone, sample, condition)
- Unsupervised K-means clustering
- Correlation heatmap plotting
- PCA with explained variance reporting

---

## 🗂️ Table of Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Pipeline Overview](#pipeline-overview)
  - [1. Load Data](#1-load-data)
  - [2. Prepare Data](#2-prepare-data)
  - [3. Create Pseudobulk](#3-create-pseudobulk)
  - [4. Save/Load Pseudobulk](#4-saveload-pseudobulk)
  - [5. K-means Clustering](#5-k-means-clustering)
  - [6. Visualize Clusters](#6-visualize-clusters)
  - [7. Correlation Heatmap](#7-correlation-heatmap)
  - [8. PCA Analysis](#8-pca-analysis)

---

## 🧰 Requirements

- R ≥ 4.0.0
- R packages:
  - `Seurat`
  - `ggpubr`
  - `factoextra`
  - `FactoMineR`

---

## 🛠️ Installation

1. Install R from [CRAN](https://cran.r-project.org/).

2. Install the required R packages:

```r
install.packages(c("Seurat", "ggpubr", "factoextra", "FactoMineR"))
```

3. Clone the repository:

```sh
git clone https://github.com/yourusername/scPseudobulk-analysis.git
cd scPseudobulk-analysis
```

---

## 🚀 Usage

1. Place your input `.rds` Seurat object into the project directory.
2. Edit `main.R` to specify:
   - Input file path
   - Cell types of interest
   - Assay type (e.g., "RNA")
   - Number of clusters
3. Run the pipeline:

```sh
Rscript main.R
```

---

## 🔄 Pipeline Overview

### 1. Load Data

```r
load_data <- function(file_path) {
  readRDS(file_path)
}
```

### 2. Prepare Data

Subsets the Seurat object to specific cell types and normalizes expression.

```r
prepare_data <- function(seurat_obj, cell_types, assay_type = "RNA") {
  Idents(seurat_obj) <- seurat_obj$CellType
  seurat_obj <- subset(seurat_obj, idents = cell_types)
  DefaultAssay(seurat_obj) <- assay_type
  NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
}
```

### 3. Create Pseudobulk

Aggregates expression per group (e.g., per clone or sample) into a pseudobulk object.

```r
create_pseudobulk <- function(seurat_obj, group_by = "Clone") {
  s.bulk <- AverageExpression(seurat_obj, return.seurat = TRUE, group.by = group_by)
  s.bulk <- FindVariableFeatures(s.bulk)
  s.bulk <- ScaleData(s.bulk)
  s.bulk <- RunPCA(s.bulk, features = VariableFeatures(s.bulk), npcs = 10)
  s.bulk
}
```

### 4. Save/Load Pseudobulk

```r
save_pseudobulk <- function(s.bulk, file_path) {
  saveRDS(s.bulk, file_path)
}

load_pseudobulk <- function(file_path) {
  readRDS(file_path)
}
```

### 5. K-means Clustering

```r
perform_kmeans <- function(s.bulk, assay = "RNA", num_clusters = 5) {
  exp <- GetAssayData(s.bulk, assay = assay)
  expression_matrix <- t(as.matrix(exp))
  kmeans(expression_matrix, centers = num_clusters)
}
```

### 6. Visualize Clusters

```r
plot_kmeans <- function(kmeans_result, expression_matrix) {
  fviz_cluster(kmeans_result, data = expression_matrix, stand = FALSE)
}
```

### 7. Correlation Heatmap

```r
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
```

### 8. PCA Analysis

```r
run_pca <- function(correlation_matrix) {
  pca_result <- prcomp(correlation_matrix, scale. = TRUE)
  pc_scores <- pca_result$x
  variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  cat("Variance explained by PC1:", variance_explained[1], "\n")
  cat("Variance explained by PC2:", variance_explained[2], "\n")
  list(pc_scores = pc_scores, variance_explained = variance_explained)
}
```

---

