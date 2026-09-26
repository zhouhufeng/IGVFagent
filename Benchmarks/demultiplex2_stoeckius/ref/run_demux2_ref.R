#!/usr/bin/env Rscript
# Reference run of the authors' deMULTIplex2 (Gartner-Lab/deMULTIplex2 @ de48333)
# by sourcing its R/ files unmodified (em.R, classify.R). Avoids the heavy
# Bioconductor Imports (ShortRead, XVector) that only preprocessing needs.
# Usage: run_demux2_ref.R <pkg_repo> <tags.csv (cells x tags, first col = cell)> <out_prefix> [max.cell.fit] [max.iter]
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c("/n/holystore01/LABS/xlin/Lab/zhouhufeng/envs/Rlibs/R4.4_demux2", .libPaths()))
suppressMessages({library(MASS); library(Matrix); library(magrittr)})
pkg <- args[1]; inp <- args[2]; out <- args[3]
max.cell.fit <- if (length(args) >= 4) as.numeric(args[4]) else 1000   # benchmark_demultiplex2 default
max.iter     <- if (length(args) >= 5) as.numeric(args[5]) else 30     # benchmark_demultiplex2 default
for (f in c("em.R", "classify.R")) source(file.path(pkg, "R", f))
tag <- read.csv(inp, row.names = 1, check.names = FALSE)
tag <- as.matrix(tag)
res <- demultiplexTags(tag, init.cos.cut = 0.5, converge.threshold = 1e-3, max.iter = max.iter,
                       prob.cut = 0.5, min.cell.fit = 10, max.cell.fit = max.cell.fit,
                       min.quantile.fit = 0.05, max.quantile.fit = 0.95,
                       residual.type = "rqr", plot.umap = "none", plot.diagnostics = FALSE, seed = 1)
at <- res$assign_table
calls <- at$barcode_assign
calls[at$barcode_count == 0] <- "Negative"
calls[at$barcode_count > 1] <- "Multiplet"
write.csv(data.frame(cell = rownames(at), call = calls, barcode_count = at$barcode_count,
                     droplet_type = at$droplet_type), paste0(out, "_calls.csv"), row.names = FALSE, quote = FALSE)
write.csv(data.frame(cell = rownames(res$prob_mtx), res$prob_mtx, check.names = FALSE),
          paste0(out, "_prob.csv"), row.names = FALSE, quote = FALSE)
cat("cells", nrow(at), "\n"); print(table(at$droplet_type))
