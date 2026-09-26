#!/usr/bin/env Rscript
# Fig. 2B / Table S1 simulations 1-5, generated with the authors' simulateTags
# (Gartner-Lab/deMULTIplex2 @ de48333 R/simulation.R, sourced unmodified) using
# the parameter blocks of deMULTIplex2-benchmark @ ac0fa17 simulate_model_update.R
# (each block: set.seed(2023); ... simulateTags(seed = 2023, ...)). Truth labels
# follow that script, including homotypic doublets relabelled to their sample.
# Usage: make_simulations.R <pkg_repo> <out_dir>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c("/n/holystore01/LABS/xlin/Lab/zhouhufeng/envs/Rlibs/R4.4_demux2", .libPaths()))
suppressMessages({library(MASS); library(Matrix); library(jsonlite)})
source(file.path(args[1], "R", "simulation.R"))
out <- args[2]; dir.create(out, showWarnings = FALSE, recursive = TRUE)
P <- list(
  sim1 = list(n.bc = 5,  mean = 7, sd = .1, sdlog = .7,  min = 300, max = 500,  lm = function(a, b) log(mean(c(a, b))),     lsd = .1, theta = 2,  b0 = -5, amb = .1, dbl = .02, drop = 2),
  sim2 = list(n.bc = 10, mean = 5, sd = .1, sdlog = 1,   min = 50,  max = 2000, lm = function(a, b) log(mean(c(a, b))),     lsd = .8, theta = 7,  b0 = -5, amb = .1, dbl = .02, drop = 2),
  sim3 = list(n.bc = 10, mean = 5, sd = .1, sdlog = 1,   min = 50,  max = 2000, lm = function(a, b) log(mean(c(a, b))),     lsd = .8, theta = 7,  b0 = -4, amb = 1,  dbl = .1,  drop = 2),
  sim4 = list(n.bc = 30, mean = 5, sd = .1, sdlog = 1,   min = 50,  max = 3000, lm = function(a, b) log(mean(c(a, b)) / 2), lsd = 1,  theta = 7,  b0 = -4, amb = 1,  dbl = .1,  drop = 2),
  sim5 = list(n.bc = 30, mean = 5, sd = .5, sdlog = 1.5, min = 50,  max = 3000, lm = function(a, b) log(mean(c(a, b)) / 2), lsd = 1,  theta = 10, b0 = -4, amb = 3,  dbl = .2,  drop = .5))
man <- list()
for (s in names(P)) {
  p <- P[[s]]
  set.seed(2023)
  csm <- rnorm(p$n.bc, mean = p$mean, sd = p$sd)
  ncell <- round(rlnorm(p$n.bc, p$lm(p$min, p$max), p$lsd))
  ncell[ncell < p$min] <- p$min; ncell[ncell > p$max] <- p$max
  sim <- simulateTags(n.cell = ncell, n.bc = p$n.bc, seed = 2023, nb.theta = p$theta,
                      cell.staining.meanlog = csm, cell.staining.sdlog = p$sdlog,
                      ambient.meanlog = p$amb, b0 = p$b0, doublet.rate = p$dbl,
                      dropout.lambda = p$drop, cell.contam = TRUE, ambient.contam = TRUE,
                      return.all = TRUE)
  tag <- as.matrix(sim$final.umi.mtx)
  true_label <- sapply(strsplit(rownames(tag), "_"), function(x) x[1])
  cell_comp <- sapply(strsplit(rownames(tag), "[|]"), function(x) {
    if (length(x) >= 2) { x <- unlist(strsplit(x[[2]], "_")); if (x[1] == x[3]) x[1] else NA } else NA })
  true_label[!is.na(cell_comp)] <- cell_comp[!is.na(cell_comp)]
  write.csv(tag, file.path(out, paste0(s, "_c1_tags.csv")))
  write.csv(data.frame(cell = rownames(tag), truth = true_label), file.path(out, paste0(s, "_truth.csv")), row.names = FALSE)
  mp <- setNames(as.list(colnames(tag)), colnames(tag))
  write(toJSON(mp, auto_unbox = TRUE, pretty = TRUE), file.path(out, paste0(s, "_tag_mapping.json")))
  man[[s]] <- list(captures = list(paste0(s, "_c1_tags.csv")), truth = paste0(s, "_truth.csv"), tag_mapping = mp, n_cells = nrow(tag))
  cat(s, nrow(tag), "cells x", ncol(tag), "tags;", table(true_label == "doublet")[["TRUE"]], "doublets\n")
}
write(toJSON(man, auto_unbox = TRUE, pretty = TRUE), file.path(out, "sim_manifest.json"))
