library(GOSemSim)

calculate_and_save_similarity <- function(type) {
  fn_hyper <- sprintf("Data/Processed/GoSemsimData/%s/mgeneSim_hyper.csv", type)
  
  hsGO <- godata('org.Hs.eg.db', keytype = "SYMBOL", ont=type, computeIC = FALSE)
  
  data <- read.csv("Data/Processed/Models/CrossFeature/cross_feature_selection.csv")
  genes <- data$gene
  genes_hyper <- data$gene[data$DNAm == "hyper"]
  
  gene_similarity_hyper <- mgeneSim(genes_hyper, semData = hsGO, measure = "Wang", combine = "BMA")
  
  result_hyper <- as.data.frame(gene_similarity_hyper)
  write.csv(result_hyper, file = fn_hyper, row.names = TRUE)
  
  return(gene_similarity_hyper)
}

combine_geometric_mean <- function(similarity_matrices) {
  combined_matrix <- similarity_matrices[[1]]

  for (i in 2:length(similarity_matrices)) {
    combined_matrix <- combined_matrix * similarity_matrices[[i]]
  }
  
  combined_matrix <- combined_matrix^(1/3)
  
  return(combined_matrix)
}

ontologies <- c("BP", "CC", "MF")
similarity_results <- list()

for (type in ontologies) {
  similarity_results[[type]] <- calculate_and_save_similarity(type)
}

fsim <- combine_geometric_mean(similarity_results)
write.csv(fsim, file = "Data/Processed/GoSemsimData/mean/mgeneSim_hyper.csv", row.names = TRUE)
