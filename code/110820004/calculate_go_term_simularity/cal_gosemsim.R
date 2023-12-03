library(GOSemSim)

fn_o <- "C:/Users/user/Desktop/project/mgeneSim_results.csv"
fn_hyper <- "C:/Users/user/Desktop/project/mgeneSim_hyper.csv"
fn_hypo <- "C:/Users/user/Desktop/project/mgeneSim_hypo.csv"

hsGO <- godata('org.Hs.eg.db', keytype = "SYMBOL", ont="MF", computeIC = FALSE)

data <- read.csv("C:/Users/user/Desktop/project/comorbidity_group_auc.csv")
genes <- data$gene
genes_hyper <- data$gene[data$DNAm == "hyper"]
genes_hypo <- data$gene[data$DNAm == "hypo"]

gene_similarity <- mgeneSim(genes, semData = hsGO, measure = "Wang", combine = "BMA")
gene_similarity_hyper <- mgeneSim(genes_hyper, semData = hsGO, measure = "Wang", combine = "BMA")
genet_similarity_hypo <- mgeneSim(genes_hypo, semData = hsGO, measure = "Wang", combine = "BMA")

result_sim <- as.data.frame(gene_similarity)
write.csv(result_sim, file = fn_o, row.names = FALSE)

result_hyper <- as.data.frame(gene_similarity_hyper)
write.csv(result_hyper, file = fn_hyper, row.names = FALSE)

result_hypo <- as.data.frame(genet_similarity_hypo)
write.csv(result_hypo, file = fn_hypo, row.names = FALSE)
