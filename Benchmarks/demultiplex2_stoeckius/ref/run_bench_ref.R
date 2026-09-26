#!/usr/bin/env Rscript
# Authors' reference for a Table S3 dataset: deMULTIplex2 (Gartner-Lab/deMULTIplex2
# @ de48333, R/em.R + R/classify.R + R/benchmarking.R sourced unmodified) run per
# capture with benchmark_demultiplex2()'s defaults (max.cell.fit = 1000,
# max.iter = 30, seed = 1), calls combined across captures ("bc_cbn"), and scored
# with the authors' confusion_stats() against the genetic-donor truth.
# Usage: run_bench_ref.R <pkg_repo> <work_dir> <manifest.json> <dataset> <out_prefix> [max.cell.fit] [max.quantile.fit]
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c("/n/holystore01/LABS/xlin/Lab/zhouhufeng/envs/Rlibs/R4.4_demux2", .libPaths()))
suppressMessages({library(MASS); library(Matrix); library(magrittr); library(jsonlite)})
pkg <- args[1]; wd <- args[2]; man <- fromJSON(args[3]); ds <- args[4]; out <- args[5]
max.cell.fit <- if (length(args) >= 6) as.numeric(args[6]) else 1000
max.quantile.fit <- if (length(args) >= 7) as.numeric(args[7]) else 0.95
for (f in c("em.R", "classify.R", "benchmarking.R")) source(file.path(pkg, "R", f))
m <- man[[ds]]
calls_all <- c()
for (cap in m$captures) {
    tag <- as.matrix(read.csv(file.path(wd, cap), row.names = 1, check.names = FALSE))
    res <- demultiplexTags(tag, init.cos.cut = 0.5, converge.threshold = 1e-3, max.iter = 30,
                           prob.cut = 0.5, min.cell.fit = 10, max.cell.fit = max.cell.fit,
                           min.quantile.fit = 0.05, max.quantile.fit = max.quantile.fit,
                           residual.type = "rqr", plot.umap = "none", plot.diagnostics = FALSE, seed = 1)
    calls <- res$assign_table$barcode_assign
    calls[res$assign_table$barcode_count == 0] <- "Negative"
    calls[res$assign_table$barcode_count > 1] <- "Multiplet"
    names(calls) <- rownames(res$assign_table)
    calls_all <- c(calls_all, calls)
}
truth <- read.csv(file.path(wd, m$truth), stringsAsFactors = FALSE)
true_label <- setNames(truth$truth, truth$cell)
tag_mapping <- data.frame(tag = names(m$tag_mapping), true_label = unlist(m$tag_mapping))
cs <- confusion_stats(calls_all, true_label, tag_mapping, call.multiplet = "Multiplet", true.multiplet = "doublet")
write.csv(data.frame(cell = names(calls_all), call = calls_all), paste0(out, "_calls.csv"), row.names = FALSE, quote = FALSE)
res <- list(dataset = ds, n_cells = length(true_label),
            precision = unname(cs$singlet_avg_stats["precision"]),
            recall = unname(cs$singlet_avg_stats["recall"]),
            f_score = unname(cs$singlet_avg_stats["f_score"]),
            doublet = as.list(cs$doublet_avg_stats),
            per_tag = lapply(cs$tag_stats, as.list),
            max.cell.fit = max.cell.fit, max.quantile.fit = max.quantile.fit)
write(toJSON(res, auto_unbox = TRUE, pretty = TRUE, digits = NA), paste0(out, "_stats.json"))
cat(ds, "F =", res$f_score, "P =", res$precision, "R =", res$recall, "\n")
