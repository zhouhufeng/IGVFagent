# Reference: the authors' Open4Gene R package (hbliu/Open4Gene @ 6e7f36a) run
# unmodified on its bundled test data, as in the README.
.libPaths(c(file.path(Sys.getenv("IGVF_ROOT"), "Data/liu2025/Rlib"), .libPaths()))
suppressMessages({library(pscl); library(progress); library(GenomicRanges); library(Matrix); library(dplyr)})
root <- Sys.getenv("IGVF_ROOT"); dir.create(file.path(root, "Data/liu2025/ref_authors/open4gene"), recursive = TRUE, showWarnings = FALSE); setwd(file.path(root, "Data/liu2025/ref_authors/open4gene"))
source("../../repo_Open4Gene/R/Open4Gene.R")
load("../../repo_Open4Gene/inst/extdata/Open4Gene.Test.Data")
obj <- CreateOpen4GeneObj(RNA = RNA.Counts, ATAC = ATAC.Counts, Meta.data = Meta.Data,
                          Peak2Gene.Pairs = Peak.Gene, Covariates = c("lognCount_RNA","percent.mt"),
                          Celltypes = "Cell_Type")
for (ct in c("All", "Each")) {
  o <- Open4Gene(object = obj, Celltype = ct, Binary = FALSE, Method = "hurdle", MinNum.Cells = 5)
  write.table(o@Res, file = paste0("R_", ct, ".res.txt"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
}
# Export the same inputs for the Python port
writeMM(RNA.Counts, "rna.mtx"); writeLines(rownames(RNA.Counts), "rna_genes.txt")
writeMM(ATAC.Counts, "atac.mtx"); writeLines(rownames(ATAC.Counts), "atac_peaks.txt")
write.csv(Meta.Data, "meta.csv")
pg <- Peak.Gene; colnames(pg) <- c("Peak","Gene"); write.csv(pg, "pairs.csv", row.names=FALSE)
cat(dim(RNA.Counts), dim(ATAC.Counts), nrow(Peak.Gene), "\n")
