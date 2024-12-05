suppressPackageStartupMessages({
  library(optparse)
  library(pheatmap)
  library(umap)
  library(ggplot2)
  library(Seurat)
  library(Signac)
})

##### command line options #####

# set options
option_list <- list(
  make_option(c("-n", "--motif_name"), type = "character", help = "ID and motif name matching table"),
  make_option(c("-m", "--sdmat"), type = "character", help = "Summit distance matrix path"),
  make_option(c("-o", "--fig_out_dir"), type = "character", 
              help = "Output directory for figures"),
  make_option(c("-c", "--thre_cell"), type = "numeric", default = 70, 
              help = "Threshold for cells with zero values (percentage) [default: 70]"),
  make_option(c("-t", "--thre_tf"), type = "numeric", default = 90, 
              help = "Threshold for TFs with zero values (percentage) [default: 90]"),
  make_option(c("-l", "--winsor_low"), type = "numeric", default = 0.05, 
              help = "Lower quantile for Winsorization [default: 0.05]"),
  make_option(c("-u", "--winsor_up"), type = "numeric", default = 0.95, 
              help = "Upper quantile for Winsorization [default: 0.95]"),
  make_option(c("-v", "--CV_thre"), type = "numeric", default = 0.5, 
              help = "Threshold for Coefficient of Variation (CV) [default: 0.5]"),
  make_option(c("-d", "--LSI_dims"), type = "numeric", default = 20, 
              help = "option 'dims' in function Signac::RunSVD [default: 20]"),
  make_option(c("-a", "--dim_low"), type = "numeric", default = 2, 
              help = "lower limit of option 'dims' in function Seurat::FindNeighbors [default: 2]"),
  make_option(c("-b", "--dim_up"), type = "numeric", default = 17, 
              help = "upper limit of option 'dims' in function Seurat::FindNeighbors [default: 17]"),
  make_option(c("-p", "--clustering_method"), type = "numeric", default = 1, 
              help = "Clustering method: 1 for graph-based, 2 for tree-based [default: 1]"),
  make_option(c("-k", "--tree_clusters"), type = "numeric", 
              help = "Number of clusters for tree-based clustering (must be provided if -p is 2)")
)


# Parse the command line options
opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

# required_args <- c("motif_name", "sdmat", "desired_number_of_clusters", "fig_out_dir")
required_args <- c("motif_name", "sdmat", "fig_out_dir")
missing_args <- required_args[!required_args %in% names(opts) | sapply(opts[required_args], is.null)]

if (length(missing_args) > 0) {
  stop(paste("Missing required arguments:", paste(missing_args, collapse = ", "), "\nUse -h for help."))
}

# Check if tree_clusters is provided when clustering_method is 2
if (opts$clustering_method == 2 && is.null(opts$tree_clusters)) {
  stop("Error: --tree_clusters must be specified when --clustering_method is set to 2.")
}

##### command line options #####


##### support functions #####
## (1) ## 
filter_zero_values <- function(df, n = 0, m = 0) {

  row_zero_percent <- apply(df, 1, function(x) sum(x == 0) / length(x) * 100)
  rows_to_remove <- which(row_zero_percent > (100 - n))
  filtered_df <- df[row_zero_percent <= (100 - n), ]
  
  col_zero_percent <- apply(filtered_df, 2, function(x) sum(x == 0) / length(x) * 100)
  cols_to_remove <- which(col_zero_percent > (100 - m))
  filtered_df <- filtered_df[, col_zero_percent <= (100 - m)]
  
  return(list(filtered_df = filtered_df, rows_removed = rows_to_remove, cols_removed = cols_to_remove))

}


filter_zero_values_old <- function(df, n = 0, m = 0) {

  row_zero_percent <- apply(df, 1, function(x) sum(x == 0) / length(x) * 100)
  rows_to_remove <- which(row_zero_percent > (100 - n))
  filtered_df <- df[row_zero_percent <= (100 - n), ]
  col_zero_percent <- apply(filtered_df, 2, function(x) sum(x == 0) / length(x) * 100)
  cols_to_remove <- which(col_zero_percent > (100 - m))
  filtered_df <- filtered_df[, col_zero_percent <= (100 - m)]
  
  return(list(filtered_df = filtered_df, rows_removed = rows_to_remove, cols_removed = cols_to_remove))
}


## (2) ## 
find_highVar_features_vst <- function(expression_df, top_percentage = 0.5) {

  feature_means <- colMeans(expression_df)
  feature_vars <- apply(expression_df, 2, var)

  fit <- loess(feature_vars ~ feature_means)
  expected_vars <- predict(fit, feature_means)
  standardized_vars <- feature_vars / expected_vars

  threshold <- quantile(standardized_vars, 1 - top_percentage)
  highly_variable_features <- names(standardized_vars[standardized_vars > threshold])

  highly_variable_features_indices <- match(highly_variable_features, colnames(expression_df))

  return(highly_variable_features_indices)
}


## (3) ##
reorderMatrix <- function(sourceMat) {
  
  rowNum <- nrow(sourceMat)
  colNum <- ncol(sourceMat)
  
  randomRowOrder <- sample(rowNum)

  randomMat <- sourceMat[randomRowOrder, , drop = FALSE]

  if (!is.null(rownames(sourceMat))) {
    rownames(randomMat) <- rownames(sourceMat)[randomRowOrder]
  }
  
  if (!is.null(colnames(sourceMat))) {
    colnames(randomMat) <- colnames(sourceMat)
  }
  
  return(randomMat)

}



## (4) ##
retain_common_columns <- function(df_list) {

  extract_prefix <- function(column_name) {
    parts <- unlist(strsplit(column_name, "_"))
    if (length(parts) >= 2) {
      return(paste(parts[1:2], collapse = "_"))
    } else {
      return(column_name)
    }
  }

  renamed_df_list <- lapply(df_list, function(df) {
    colnames(df) <- sapply(colnames(df), extract_prefix)
    return(df)
  })
  
  common_columns <- Reduce(intersect, lapply(renamed_df_list, colnames))
  
  df_list_common <- lapply(renamed_df_list, function(df) df[, common_columns, drop = FALSE])
  
  return(df_list_common)

}



## (5) Winsorization ##
winsorize_df <- function(df, lower_quantile = 0.05, upper_quantile = 0.95) {
  original_row_names <- rownames(df)
  
  df_winsorized <- as.data.frame(lapply(df, function(x) {
    lower_bound <- quantile(x, lower_quantile, na.rm = TRUE)
    upper_bound <- quantile(x, upper_quantile, na.rm = TRUE)
    x <- pmin(pmax(x, lower_bound, na.rm = TRUE), upper_bound, na.rm = TRUE)
    return(x)
  }))
  
  rownames(df_winsorized) <- original_row_names
  return(df_winsorized)
}



## (6) ##
process_column_name <- function(col_name, motifName) {
  parts <- strsplit(col_name, "_")[[1]]
  if (length(parts) > 2) {
    short_name <- paste(parts[1:2], collapse = ".")
  } else {
    short_name <- paste(parts, collapse = ".")
  }
  
  match_index <- match(short_name, motifName[, 1])
  
  if (!is.na(match_index)) {
    return(motifName[match_index, 2])
  } else {
    return(col_name)
  }
}

##### support functions #####


##### Main Script #####
# Read input files using the command line options
motifName <- read.table(opts$motif_name, header = F)
df <- read.table(opts$sdmat, header = T)

thre_cell <- opts$thre_cell
thre_tf <- opts$thre_tf
winsor_low <- opts$winsor_low
winsor_up <- opts$winsor_up
CV_thre <- opts$CV_thre
desired_number_of_clusters <- opts$desired_number_of_clusters
fig_out_dir <- opts$fig_out_dir
lsi_dims <- opts$LSI_dims
n_dim_low <- opts$dim_low
n_dim_up <- opts$dim_up


df <- filter_zero_values(df, n = thre_cell, m = thre_tf)$filtered_df
# df <- filter_zero_values_old(df, n = thre_cell)$filtered_df
# df <- filter_zero_values_old(df, m = thre_tf)$filtered_df

## Winsorization
df_norm <- winsorize_df(df, lower_quantile = winsor_low, upper_quantile = winsor_up)
rownames(df_norm) <- rownames(df)
dim(df_norm)

## name matching
name_motif <- unname(sapply(colnames(df_norm), process_column_name, motifName = motifName))
unique_motifs <- name_motif

counts <- table(name_motif)

seen <- list()

for (i in seq_along(name_motif)) {
  motif <- name_motif[i]
  
  if (counts[motif] > 1) {
    if (is.null(seen[[motif]])) {
      seen[[motif]] <- 1
    } else {
      seen[[motif]] <- seen[[motif]] + 1
    }
    unique_motifs[i] <- paste0(motif, "_", seen[[motif]])
  }
}

colnames(df_norm) <- unique_motifs

## CV 
sds <- apply(df_norm, 2, sd)
means <- apply(df_norm, 2, mean)

cv <- sds / means

sorted_cv <- sort(cv, decreasing = TRUE)

thre_type <- "cv"
filtered_motif <- names(sorted_cv[sorted_cv > CV_thre])
df_norm_filtered <- df_norm[, filtered_motif]

dim(df_norm_filtered)

df_norm_filtered[is.na(df_norm_filtered)] <- 0


## Based on clustering method
clustering_method <- opts$clustering_method

if (clustering_method == 1) {
  # Based on graph (SNN + Louvain)
  df_norm_filtered_clst <- df_norm_filtered
  df_norm_filtered_clst_TFIDF_topFeatures <- Signac::FindTopFeatures(df_norm_filtered_clst)
  df_norm_filtered_clst_TFIDF_SVD <- Signac::RunSVD(t(df_norm_filtered_clst), n = lsi_dims)

  svd_matrix <- df_norm_filtered_clst_TFIDF_SVD@cell.embeddings
  FindNeighbors_res <- Seurat::FindNeighbors(svd_matrix, reduction = 'lsi', dims = n_dim_low:n_dim_up)
  sample_clusters <- FindClusters(object = FindNeighbors_res$snn, verbose = FALSE, algorithm = 1)
  
  row_hclust <- pheatmap(df_norm_filtered,silent = TRUE)$tree_row
  sample_clusters <- cutree(row_hclust, k = 3)

} else if (clustering_method == 2) {
  # Based on tree
  row_hclust <- pheatmap(df_norm_filtered, silent = TRUE)$tree_row
  sample_clusters <- cutree(row_hclust, k = opts$tree_clusters)  # Use user-specified k value

} else {
  stop("Invalid clustering method. Please choose 1 (graph-based) or 2 (tree-based).")
}

sample_clusters <- as.data.frame(sample_clusters)


if (!dir.exists(fig_out_dir)) {
  dir.create(fig_out_dir, recursive = TRUE)
  print(paste("Directory created:", fig_out_dir))
}

colnames(sample_clusters) <- c("cluster")
write.table(sample_clusters, file = paste0(fig_out_dir, "/sample_clusters.txt"), sep = "\t", row.names = T, quote = F)

## UMAP
umap_res <- umap(df_norm_filtered, n_components = 2)

## cluster_label matching
cluster_label <- sample_clusters[rownames(umap_res$layout), 1]
UMAP_cluster <- cbind(umap_res$layout, cluster_label)
colnames(UMAP_cluster) <- c("umap1", "umap2", "cluster")

## p
p <- ggplot(UMAP_cluster, aes(x = umap1, y = umap2, color = as.factor(cluster))) +
  geom_point(size = 2.5) +
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme_classic() +
  labs(color = "clusters")

ggsave(paste0(fig_out_dir, "/umap_plot.png"), plot = p, dpi = 1200)
ggsave(paste0(fig_out_dir, "/umap_plot.pdf"), plot = p)
