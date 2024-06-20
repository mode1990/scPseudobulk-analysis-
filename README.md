# scPseudobulk-analysis-
Automated pipeline for psuedobulk analysis and downstream unsupervised analysis based on seurat v5
Sure, here is a README.md file for your automated pipeline.

```markdown
# Pseudobulk Analysis Pipeline with Seurat v5

This repository contains an automated pipeline for pseudobulk analysis and downstream unsupervised analysis using Seurat v5. The pipeline includes normalization, pseudobulk creation, clustering, heatmap generation, and principal component analysis (PCA).

## Table of Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Pipeline Steps](#pipeline-steps)
  - [1. Load Data](#1-load-data)
  - [2. Prepare Data](#2-prepare-data)
  - [3. Create Pseudobulk](#3-create-pseudobulk)
  - [4. Save and Load Pseudobulk](#4-save-and-load-pseudobulk)
  - [5. Perform K-means Clustering](#5-perform-k-means-clustering)
  - [6. Plot K-means Clustering Result](#6-plot-k-means-clustering-result)
  - [7. Calculate and Plot Heatmap](#7-calculate-and-plot-heatmap)
  - [8. Run PCA and Print Variance](#8-run-pca-and-print-variance)

## Requirements

- R (version >= 4.0.0)
- The following R packages:
  - Seurat
  - ggpubr
  - factoextra
  - FactoMineR

## Installation

1. Install R from [CRAN](https://cran.r-project.org/).
2. Install the required R packages:

    ```r
    install.packages(c("Seurat", "ggpubr", "factoextra", "FactoMineR"))
    ```

3. Clone this repository:

    ```sh
    git clone https://github.com/yourusername/pseudobulk-analysis-pipeline.git
    cd pseudobulk-analysis-pipeline
    ```

## Usage

1. Place your input Seurat object file (`.rds` format) in the repository directory.
2. Modify the `main.R` script to set the correct file paths and parameters (input file path, cell types, assay type, number of clusters).
3. Run the pipeline:

    ```sh
    Rscript main.R
    ```

## Pipeline Steps

### 1. Load Data

The `load_data` function reads the Seurat object from the specified file path.

```r
load_data <- function(file_path) {
  return(readRDS(file_path))
}
```

### 2. Prepare Data

The `prepare_data` function subsets the Seurat object based on specified cell types and normalizes the data.

```r
prepare_data <- function(seurat_obj, cell_types, assay_type = "RNA") {
  Idents(seurat_obj) <- seurat_obj$CellType
  seurat_obj <- subset(seurat_obj, idents = cell_types)
  DefaultAssay(seurat_obj) <- assay_type
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000, assay = assay_type)
  return(seurat_obj)
}
```

### 3. Create Pseudobulk

The `create_pseudobulk` function averages expression data to create a pseudobulk Seurat object.

```r
create_pseudobulk <- function(seurat_obj, group_by = 'Clone') {
  s.bulk <- AverageExpression(seurat_obj, return.seurat = TRUE, group.by = group_by)
  s.bulk <- FindVariableFeatures(s.bulk)
  s.bulk <- ScaleData(s.bulk)
  s.bulk <- RunPCA(s.bulk, features = VariableFeatures(object = s.bulk), npcs = 10)
  return(s.bulk)
}
```

### 4. Save and Load Pseudobulk

Functions to save and load the pseudobulk Seurat object.

```r
save_pseudobulk <- function(s.bulk, file_path) {
  saveRDS(s.bulk, file_path)
}

load_pseudobulk <- function(file_path) {
  return(readRDS(file_path))
}
```

### 5. Perform K-means Clustering

The `perform_kmeans` function performs K-means clustering on the pseudobulk data.

```r
perform_kmeans <- function(s.bulk, assay = "RNA", num_clusters = 5) {
  exp <- GetAssayData(object = s.bulk, assay = assay)
  expression_matrix <- t(as.matrix(exp))
  kmeans_result <- kmeans(expression_matrix, centers = num_clusters)
  return(kmeans_result)
}
```

### 6. Plot K-means Clustering Result

The `plot_kmeans` function visualizes the K-means clustering result.

```r
plot_kmeans <- function(kmeans_result, expression_matrix) {
  return(fviz_cluster(kmeans_result, data = expression_matrix, stand = FALSE))
}
```

### 7. Calculate and Plot Heatmap

The `plot_heatmap` function generates a heatmap of the correlation matrix.

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

### 8. Run PCA and Print Variance

The `run_pca` function performs PCA and prints the variance explained by the first two principal components.

```r
run_pca <- function(correlation_matrix) {
  pca_result <- prcomp(correlation_matrix, scale. = TRUE)
  pc_scores <- pca_result$x
  variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  cat("Variance explained by PC1:", variance_explained[1], "\n")
  cat("Variance explained by PC2:", variance_explained[2], "\n")
  return(list(pc_scores = pc_scores, variance_explained = variance_explained))
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
```

This README provides an overview of the pipeline, instructions for setup and usage, and a brief description of each function. Make sure to update any placeholders such as file paths or repository URLs with the correct information specific to your setup.
