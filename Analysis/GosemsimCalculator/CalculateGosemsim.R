library(GOSemSim)

type <- "BP"
fn_o <- sprintf("Data/Processed/GoSemsimData/%s/mgeneSim_results.csv", type)
fn_hyper <- sprintf("Data/Processed/GoSemsimData/%s/mgeneSim_hyper.csv", type)
fn_hypo <- sprintf("Data/Processed/GoSemsimData/%s/mgeneSim_hypo.csv", type)

hsGO <- godata('org.Hs.eg.db', keytype = "SYMBOL", ont=type, computeIC = FALSE)

data <- read.csv("Data/Processed/Models/CrossFeature/cross_feature_selection.csv")
genes <- data$gene
genes_hyper <- data$gene[data$DNAm == "hyper"]
genes_hypo <- data$gene[data$DNAm == "hypo"]

gene_similarity <- mgeneSim(genes, semData = hsGO, measure = "Wang", combine = "BMA")
gene_similarity_hyper <- mgeneSim(genes_hyper, semData = hsGO, measure = "Wang", combine = "BMA")
genet_similarity_hypo <- mgeneSim(genes_hypo, semData = hsGO, measure = "Wang", combine = "BMA")

result_sim <- as.data.frame(gene_similarity)
write.csv(result_sim, file = fn_o, row.names = TRUE)

result_hyper <- as.data.frame(gene_similarity_hyper)
write.csv(result_hyper, file = fn_hyper, row.names = TRUE)

result_hypo <- as.data.frame(genet_similarity_hypo)
write.csv(result_hypo, file = fn_hypo, row.names = TRUE)
