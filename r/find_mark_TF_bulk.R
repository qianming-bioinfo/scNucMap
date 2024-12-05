suppressPackageStartupMessages({
library(dplyr)
library(optparse)
})

# Define command line options
option_list <- list(
  make_option(c("-d", "--data_folder"), type = "character", help = "Data folder path (required)"),
  make_option(c("-s", "--sample_labels"), type = "character", help = "Comma-separated list of sample labels (required)"),
  make_option(c("-o", "--outDir"), type = "character", help = "Output directory (required)"),
  make_option(c("-t", "--score_thre"), type = "numeric", default = 0.1,
              help = "Score threshold [default %default]"),
  make_option(c("-p", "--p_thre"), type = "numeric", default = 0.01,
              help = "P-value threshold [default %default]")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

# Assign the parsed values to variables
outDir <- opts$outDir
sample_labels <- unlist(strsplit(opts$sample_labels, ","))
data_folder <- opts$data_folder
score_thre <- opts$score_thre
p_thre <- opts$p_thre

# Check if essential parameters are provided
if (is.null(outDir) || is.null(sample_labels) || is.null(data_folder)) {
  stop("Error: 'outDir', 'sample_labels', and 'data_folder' are required parameters and must be specified.")
}

#### support functions ####
## (1) ##
motif_name_match <- function(mydf, motifName) {
  
  name_motif <- unname(sapply(colnames(mydf), process_column_name, motifName = motifName))
  
  counts <- table(name_motif)
  unique_motifs <- name_motif
  
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
  
  colnames(mydf) <- unique_motifs
  
  return(mydf)
}



## (2) ##
motif_name_match_column <- function(mydf, column_name, motifName) {
  
  column_data <- mydf[[column_name]]
  name_motif <- unname(sapply(column_data, process_column_name, motifName = motifName))
  mydf[[column_name]] <- name_motif
  
  return(mydf)
}



## (3) ##
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



## (4) ##
contingency_tests_pooled <- function(center_count, flank_count, center_count_bg, flank_count_bg) {
  
  p_values <- numeric(length(center_count))
  names(p_values) <- names(center_count)
  
  for (i in seq_along(center_count)) {
    motif <- names(center_count)[i]
    
    center_motif <- center_count[i]
    flank_motif <- flank_count[i] - center_motif
    center_bg <- center_count_bg
    flank_bg <- flank_count_bg - center_count_bg
    
    contingency_table <- matrix(c(center_motif, flank_motif, center_bg, flank_bg), nrow = 2, byrow = TRUE)
    
    if (any(contingency_table == 0)) {
      fisher_test <- fisher.test(contingency_table)
      p_values[i] <- fisher_test$p.value
    } else {
      chi_test <- chisq.test(contingency_table, correct = FALSE)
      if (all(chi_test$expected >= 5)) {
        p_values[i] <- chi_test$p.value
      } else {
        if (any(chi_test$expected < 5 & chi_test$expected > 1)) {
          yates_test <- chisq.test(contingency_table, correct = TRUE)
          p_values[i] <- yates_test$p.value
        } else {
          fisher_test <- fisher.test(contingency_table)
          p_values[i] <- fisher_test$p.value
        }
      }
    }
  }
  return(as.data.frame(p_values))
}


#### support functions ####


center_count_paths <- list()
flank_count_paths <- list()
center_count_bg_paths <- list()
flank_count_bg_paths <- list()

for (sample in sample_labels) {
  center_count_paths[[sample]] <- file.path(data_folder, paste0(sample, "_center_count.txt"))
  flank_count_paths[[sample]] <- file.path(data_folder, paste0(sample, "_flank_count.txt"))
  center_count_bg_paths[[sample]] <- file.path(data_folder, paste0(sample, "_center_count_bg.txt"))
  flank_count_bg_paths[[sample]] <- file.path(data_folder, paste0(sample, "_flank_count_bg.txt"))
  
  if (!file.exists(center_count_paths[[sample]])) stop(paste0("Error: ", center_count_paths[[sample]], " not found"))
  if (!file.exists(flank_count_paths[[sample]])) stop(paste0("Error: ", flank_count_paths[[sample]], " not found"))
  if (!file.exists(center_count_bg_paths[[sample]])) stop(paste0("Error: ", center_count_bg_paths[[sample]], " not found"))
  if (!file.exists(flank_count_bg_paths[[sample]])) stop(paste0("Error: ", flank_count_bg_paths[[sample]], " not found"))
}


free_scores <- list()
conti_test_results <- list()

for (sample in sample_labels) {

  center_count <- read.table(center_count_paths[[sample]], header = T)
  flank_count <- read.table(flank_count_paths[[sample]], header = T)
  center_count_bg <- read.table(center_count_bg_paths[[sample]], header = T)
  flank_count_bg <- read.table(flank_count_bg_paths[[sample]], header = T)
  
  pooled_center_count <- colSums(center_count)
  pooled_flank_count <- colSums(flank_count)
  pooled_center_count_bg <- sum(center_count_bg)
  pooled_flank_count_bg <- sum(flank_count_bg)

  pooled_frac <- pooled_center_count / pooled_flank_count
  pooled_frac_bg <- pooled_center_count_bg / pooled_flank_count_bg
  frac_ratio <- pooled_frac_bg / pooled_frac
  free_score <- log2(frac_ratio)
  free_score[is.infinite(free_score)] <- 0
  free_scores[[sample]] <- as.data.frame(free_score)

  conti_test_res <- contingency_tests_pooled(pooled_center_count, pooled_flank_count, pooled_center_count_bg, pooled_flank_count_bg)
  conti_test_results[[sample]] <- conti_test_res
}


merged_results <- list()

for (sample in sample_labels) {
  free_score_data <- free_scores[[sample]]
  conti_test_res_data <- conti_test_results[[sample]]
  
  free_score_data$rownames_column <- rownames(free_score_data)
  conti_test_res_data$rownames_column <- rownames(conti_test_res_data)
  
  merged_res <- merge(free_score_data, conti_test_res_data, by = "rownames_column")
  merged_results[[sample]] <- merged_res
}

# # motif name matching
# for (sample in sample_labels) {
#   merged_results[[sample]] <- motif_name_match_column(merged_results[[sample]], column_name = "rownames_column", motifName = motifName)
# }


filtered_results <- list()

for (sample in sample_labels) {
  filtered_res <- merged_results[[sample]] %>%
    filter(free_score > score_thre & p_values < p_thre)
  filtered_results[[sample]] <- filtered_res
}


if (!dir.exists(outDir)) {
  dir.create(outDir)
}

for (sample in sample_labels) {
  output_path <- file.path(outDir, paste0(sample, "_filtered_results.txt"))
  write.table(filtered_results[[sample]], file = output_path, sep = "\t", row.names = FALSE)
}

