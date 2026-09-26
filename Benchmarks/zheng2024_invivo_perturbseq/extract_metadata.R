#!/usr/bin/env Rscript
# Extract per-cell metadata from the paper's own deposited Seurat objects
# (GSE249416_Perturb_all.qs / GSE249416_Perturb_sg.qs) without needing the
# full Seurat package -- SeuratObject alone defines the S4 class needed for
# qs::qread() to deserialize it and for @meta.data / @assays access to work.
#
# Usage:
#   Rscript extract_metadata.R <raw_dir> <out_dir>
#
# <raw_dir> must contain the two *.qs files (gunzip'd from the GEO
# suppl *.qs.gz downloads -- qs's own serialization format is wrapped in an
# extra outer gzip by the depositor, so `gunzip -k GSE249416_Perturb_all.qs.gz`
# is required before qs::qread() will recognize the stream).
#
# See OPERATIONS.md for the one-time R package bootstrap this needs (qs and
# SeuratObject are not on the R module's default library and pull in a chain
# of compiled dependencies -- Rcpp, RcppEigen, RcppParallel, stringfish, BH,
# RApiSerialize -- plus a stringfish DOWNGRADE: CRAN's current stringfish
# (>=0.17) changed its internal sfstring API in a way that breaks qs 0.27.3,
# the last version of `qs` before it was superseded by `qs2` and archived).

args <- commandArgs(trailingOnly = TRUE)
raw_dir <- if (length(args) >= 1) args[1] else stop("usage: extract_metadata.R <raw_dir> <out_dir>")
out_dir <- if (length(args) >= 2) args[2] else stop("usage: extract_metadata.R <raw_dir> <out_dir>")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(qs)
  library(SeuratObject)
})

read_metadata <- function(qs_path, csv_path) {
  cat("Reading", qs_path, "...\n")
  obj <- qs::qread(qs_path)
  stopifnot(inherits(obj, "Seurat"))
  md <- obj@meta.data
  cat("  ", nrow(md), "cells x", ncol(md), "metadata columns\n")
  write.csv(md, csv_path, row.names = TRUE)
  invisible(md)
}

all_qs <- file.path(raw_dir, "GSE249416_Perturb_all.qs")
sg_qs  <- file.path(raw_dir, "GSE249416_Perturb_sg.qs")

if (file.exists(all_qs)) {
  read_metadata(all_qs, file.path(out_dir, "all_metadata.csv"))
} else {
  cat("Missing:", all_qs, "-- gunzip GSE249416_Perturb_all.qs.gz first\n")
}

if (file.exists(sg_qs)) {
  read_metadata(sg_qs, file.path(out_dir, "sg_metadata.csv"))
} else {
  cat("Missing:", sg_qs, "-- gunzip GSE249416_Perturb_sg.qs.gz first\n")
}

cat("DONE\n")
