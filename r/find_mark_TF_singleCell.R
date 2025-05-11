#!/usr/bin/env Rscript

# Load required packages
library(optparse)

options(warn = -1)

# Parse command-line arguments
option_list <- list(
  make_option(c("-c", "--center_count"), type = "character", help = "Path to center count file", metavar = "character"),
  make_option(c("-f", "--flank_count"), type = "character", help = "Path to flank count file", metavar = "character"),
  make_option(c("-b", "--center_count_bg"), type = "character", help = "Path to background center count file", metavar = "character"),
  make_option(c("-l", "--flank_count_bg"), type = "character", help = "Path to background flank count file", metavar = "character"),
  make_option(c("-s", "--sample_clusters"), type = "character", help = "Path to sample clusters file", metavar = "character"),
  make_option(c("-m", "--motif_name"), type = "character", help = "Path to motif name mapping file", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", help = "Output directory for results", metavar = "character"),
  make_option(c("-t", "--score_thre"), type = "numeric", default = 0.1, help = "Threshold for score (default = 0.1)", metavar = "numeric"),
  make_option(c("-p", "--chisq_p_thre"), type = "numeric", default = 0.01, help = "Threshold for chi-square p-value (default = 0.05)", metavar = "numeric")
)

# Set up option parser
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# Load input data
all_center_count <- read.table(args$center_count, header = TRUE, row.names = 1)
all_flank_count <- read.table(args$flank_count, header = TRUE, row.names = 1)
all_center_count_bg <- read.table(args$center_count_bg, header = TRUE, row.names = 1)
all_flank_count_bg <- read.table(args$flank_count_bg, header = TRUE, row.names = 1)
sample_clusters <- read.table(args$sample_clusters, header = TRUE, row.names = 1)
motifName <- read.table(args$motif_name, header = FALSE)
outDir <- args$output_dir
score_thre <- args$score_thre
chisq_p_thre <- args$chisq_p_thre


## Check if essential parameters are provided
if (is.null(all_center_count) || is.null(all_flank_count) || is.null(all_center_count_bg) || 
    is.null(all_flank_count_bg) || is.null(sample_clusters) || is.null(motifName) || is.null(outDir)) {
  stop("Error: 'all_center_count', 'all_flank_count', 'all_center_count_bg', 'all_flank_count_bg', 'sample_clusters', 'motifName', and 'outDir' are required parameters and must be specified.")
}



#### support functions ####

## (1)
contingency_tests_clst <- function(prefix) {
  
  center_count <- get(paste0(prefix, "_center_count_colSum"))
  flank_count <- get(paste0(prefix, "_flank_count_colSum"))
  center_count_bg <- get(paste0(prefix, "_center_count_bg_colSum"))
  flank_count_bg <- get(paste0(prefix, "_flank_count_bg_colSum"))
  
  p_values <- matrix(NA, nrow = nrow(center_count), ncol = ncol(center_count),
                     dimnames = list(rownames(center_count), colnames(center_count)))
  
  for (i in seq_len(nrow(center_count))) {
    for (j in seq_len(ncol(center_count))) {
      # 2x2 table
      center_motif <- center_count[i, j]
      center_bg <- center_count_bg[i, 1]
      flank_motif <- flank_count[i, j] - center_count[i, j]
      flank_bg <- flank_count_bg[i, 1] - center_count_bg[i, 1]
      
      contingency_table <- matrix(c(center_motif, flank_motif, center_bg, flank_bg), nrow = 2, byrow = TRUE)
      
      if (any(contingency_table == 0)) {
        fisher_test <- fisher.test(contingency_table)
        p_values[i, j] <- fisher_test$p.value
      }
      else{
        chi_test <- chisq.test(contingency_table, correct = F)
        if (all(chi_test$expected >= 5)) {
          # chisq
          p_values[i, j] <- chi_test$p.value
        } else {
          if (any(chi_test$expected < 5 & chi_test$expected > 1)) {
            # Yates
            yates_test <- chisq.test(contingency_table, correct = TRUE)
            p_values[i, j] <- yates_test$p.value
          } 
          else { # fisher
            fisher_test <- fisher.test(contingency_table)
            p_values[i, j] <- fisher_test$p.value
          }
        }
      }
    }
  }
  return(as.data.frame(p_values))
}



## (2) 
filter_significant_relations <- function(NS, contiTest_res, score_threshold = 0.1, pvalue_threshold = 0.05) {
  
  if (!all(rownames(NS) == rownames(contiTest_res)) || 
      !all(colnames(NS) == colnames(contiTest_res))) {
    stop("The row names and column names of the matrices must match.")
  }
  
  score_condition <- (NS > score_threshold)
  pvalue_condition <- (contiTest_res < pvalue_threshold)
  
  significant_indices <- which(score_condition & pvalue_condition, arr.ind = TRUE)
  
  if (nrow(significant_indices) == 0) {
    print("No significant cell-motif relationships found.")
    return("No significant cell-motif relationships found.")
  }
  
  significant_cells <- rownames(NS)[significant_indices[, 1]]
  significant_motifs <- colnames(NS)[significant_indices[, 2]]
  
  result <- data.frame(
    Cell = significant_cells,
    Motif = significant_motifs,
    Score = NS[significant_indices],
    P_value = contiTest_res[significant_indices]
  )
  
  return(result)

}


## (3)
split_by_cluster <- function(mydf, sample_clusters) {
  
  clusters <- unique(sample_clusters[, 1])
  split_dfs <- list()
  is_single_column <- ncol(mydf) == 1
  
  for (cluster in clusters) {
    cluster_cells <- rownames(sample_clusters[sample_clusters[, 1] == cluster, , drop = FALSE])
    
    sub_df <- mydf[rownames(mydf) %in% cluster_cells, , drop = FALSE]
    
    if (is_single_column) {
      split_dfs[[paste0("cluster_", cluster)]] <- sub_df
    } else {
      split_dfs[[paste0("cluster_", cluster)]] <- sub_df
    }
  }
  
  return(split_dfs)

}


## (4)
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


motif_name_match_column <- function(mydf, column_name, motifName, new_col_name = "name") {

  column_data <- mydf[[column_name]]
  name_motif <- unname(sapply(column_data, process_column_name, motifName = motifName))
  mydf[[new_col_name]] <- name_motif
  
  return(mydf)
}




## (5)
contingency_tests <- function(prefix) {
  
  center_count <- get(paste0(prefix, "_center_count"))
  flank_count <- get(paste0(prefix, "_flank_count"))
  center_count_bg <- get(paste0(prefix, "_center_count_bg"))
  flank_count_bg <- get(paste0(prefix, "_flank_count_bg"))
  
  p_values <- matrix(NA, nrow = nrow(center_count), ncol = ncol(center_count),
                     dimnames = list(rownames(center_count), colnames(center_count)))
  
  for (i in seq_len(nrow(center_count))) {
    for (j in seq_len(ncol(center_count))) {
      # 2x2 table
      center_motif <- center_count[i, j]
      center_bg <- center_count_bg[i, 1]
      flank_motif <- flank_count[i, j] - center_count[i, j]
      flank_bg <- flank_count_bg[i, 1] - center_count_bg[i, 1]
      
      contingency_table <- matrix(c(center_motif, flank_motif, center_bg, flank_bg), nrow = 2, byrow = TRUE)
      
      # chisq, Yates or Fisher's
      if (any(contingency_table == 0)) {
        fisher_test <- fisher.test(contingency_table)
        p_values[i, j] <- fisher_test$p.value
      }
      else{
        chi_test <- chisq.test(contingency_table, correct = F)
        if (all(chi_test$expected >= 5)) {
          p_values[i, j] <- chi_test$p.value
        } else {
          if (any(chi_test$expected < 5 & chi_test$expected > 1)) {
              yates_test <- chisq.test(contingency_table, correct = TRUE)
              p_values[i, j] <- yates_test$p.value
          } 
          else { # fisher
          fisher_test <- fisher.test(contingency_table)
          p_values[i, j] <- fisher_test$p.value
          }
        }
      }
    }
  }
  return(as.data.frame(p_values))
}


#### support functions ####


# Compute fractions and log ratios
all_frac <- all_center_count / all_flank_count
all_frac_bg <- all_center_count_bg / all_flank_count_bg
all_frac_ratio <- sweep(all_frac, 1, all_frac_bg[, 1], "/")
all_NS <- -log2(all_frac_ratio)
all_NS[] <- lapply(all_NS, function(x) { x[is.nan(x)] <- 0; return(x) })

all_contiTest_res <- contingency_tests("all")
rownames(all_NS) <- rownames(all_frac)

all_NS_by_cluster <- split_by_cluster(all_NS, sample_clusters)
all_contiTest_res_by_cluster <- split_by_cluster(all_contiTest_res, sample_clusters)


# Iterate through clusters and perform analysis
cluster_names <- names(all_NS_by_cluster)

motif_freq_list <- list()
for (cluster_name in cluster_names) {

  NS <- all_NS_by_cluster[[cluster_name]]
  contiTest_res <- all_contiTest_res_by_cluster[[cluster_name]]

  score_p_filRes <- filter_significant_relations(NS, contiTest_res, score_thre, chisq_p_thre)
  score_p_filRes <- score_p_filRes[!is.infinite(score_p_filRes$Score), ]
  
  tmp_motif_freq <- as.data.frame(table(score_p_filRes$Motif))
  motif_freq_list[[cluster_name]] <- tmp_motif_freq
  
}


## hypergeometric test
all_hyper_res_list <- list()
for (cluster_name in cluster_names) {
  ## merge other cluster's motif sequence
  df_list_left <- motif_freq_list[!names(motif_freq_list) %in% cluster_name]
  merged_df_left <- Reduce(function(x, y) merge(x, y, by = "Var1", all = TRUE), df_list_left) # checked
  
  n_cols <- length(merged_df_left) - 1
  colnames(merged_df_left) <- c("Var1", paste0("freq_", letters[1:n_cols]))
  
  merged_df_tmpClst_left <- merge(motif_freq_list[[cluster_name]], merged_df_left, by = "Var1", all = TRUE)
  tmp_freq_col_name <- paste0("freq_",cluster_name)
  colnames(merged_df_tmpClst_left)[1:2] <- c("motif", tmp_freq_col_name)

  merged_df_tmpClst_left[is.na(merged_df_tmpClst_left)] <- 0
  merged_df_tmpClst_left$freq_sum <- rowSums(merged_df_tmpClst_left[, c(tmp_freq_col_name, paste0("freq_", letters[1:n_cols]))], na.rm = TRUE)
  merged_df_tmpClst_left[is.na(merged_df_tmpClst_left)] <- 0

  tmp_hyper_res <- data.frame(target_motif = character(), p_value = numeric(), stringsAsFactors = FALSE)

  for (target_motif in merged_df_tmpClst_left$motif) {
    total_freq <- sum(merged_df_tmpClst_left$freq_sum)
    targetMotif_freq_in_all <- merged_df_tmpClst_left[which(merged_df_tmpClst_left$motif == target_motif),]$freq_sum
    total_freq_of_targetVec <- sum(merged_df_tmpClst_left[[tmp_freq_col_name]])
    targetMotif_freq_in_targetVec <- merged_df_tmpClst_left[which(merged_df_tmpClst_left$motif == target_motif), tmp_freq_col_name]

    N <- total_freq
    K <- targetMotif_freq_in_all
    n <- total_freq_of_targetVec
    k <- targetMotif_freq_in_targetVec

    p_value <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)

    tmp_hyper_res <- rbind(tmp_hyper_res, data.frame(target_motif = target_motif, p_value = p_value))

  }

  # tmp_hyper_res <- motif_name_match_column(tmp_hyper_res, "target_motif", motifName)
  tmp_hyper_res$p.adjust <- p.adjust(tmp_hyper_res$p_value)
  
  all_hyper_res_list[[cluster_name]] <- tmp_hyper_res
  
}


if (!dir.exists(outDir)) {
  dir.create(outDir)
}

for (cluster_name in names(all_hyper_res_list)) {
  all_hyper_res_list[[cluster_name]] <- motif_name_match_column(all_hyper_res_list[[cluster_name]], column_name = "target_motif", motifName = motifName)
  # print(colnames(all_hyper_res_list[[cluster_name]]))
  write.table(all_hyper_res_list[[cluster_name]], 
              file = paste0(outDir, "/", cluster_name, "_hypergeo.txt"), 
              quote = FALSE, row.names = FALSE)
}

