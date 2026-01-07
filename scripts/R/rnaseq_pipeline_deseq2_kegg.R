#!/usr/bin/env Rscript

# ============================================================
# RNA-seq analysis pipeline (Salmon -> DESeq2 -> KEGG)
# Author: Askar Alshammari
#
# This script is designed for GitHub sharing.
# No metadata CSVs or raw data are included in the repository.
# ============================================================

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(apeglm)
  library(clusterProfiler)
  library(org.Cf.eg.db)
  library(readr)
  library(dplyr)
  library(tibble)
})

# -----------------------------
# User arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  args[i + 1]
}

quant_dir  <- get_arg("--quant_dir", "Quants")
gtf_file   <- get_arg("--gtf", "genomic.gtf")
meta_file  <- get_arg("--meta", "sample_sheet.csv")
outdir     <- get_arg("--out", "results")
contrast_s <- get_arg("--contrast", "condition,Virus,Mock")
kegg_org   <- get_arg("--kegg_org", "cfa")

padj_cutoff <- 0.05
lfc_cutoff  <- 1
top_var_genes <- 50

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Create metadata template if missing
# -----------------------------
if (!file.exists(meta_file)) {
  writeLines(
    c(
      "sample,condition",
      "SAMPLE_1,Mock",
      "SAMPLE_2,Virus",
      "SAMPLE_3,Virus"
    ),
    con = "sample_sheet_TEMPLATE.csv"
  )
  message("Metadata file not found.")
  message("Template 'sample_sheet_TEMPLATE.csv' created.")
  message("Fill it and rerun the script.")
  quit(status = 0)
}

# -----------------------------
# Load metadata
# -----------------------------
meta <- read_csv(meta_file, show_col_types = FALSE)

if (!all(c("sample", "condition") %in% colnames(meta))) {
  stop("Metadata must contain columns: sample, condition")
}

meta <- meta %>%
  mutate(
    sample = as.character(sample),
    condition = factor(condition)
  ) %>%
  as.data.frame()

rownames(meta) <- meta$sample

# -----------------------------
# Salmon quant files
# -----------------------------
files <- file.path(quant_dir, meta$sample, "quant.sf")
names(files) <- meta$sample

if (!all(file.exists(files))) {
  stop("Some quant.sf files are missing.")
}

# -----------------------------
# tx2gene mapping from GTF
# -----------------------------
gtf <- read_tsv(gtf_file, comment = "#", col_names = FALSE, show_col_types = FALSE)
colnames(gtf)[9] <- "attribute"

extract_attr <- function(x, key) {
  out <- sub(paste0('.*', key, ' "([^"]+)".*'), "\\1", x)
  out[grepl(";", out)] <- NA
  out
}

tx2gene <- tibble(
  TXNAME = extract_attr(gtf$attribute, "transcript_id"),
  GENEID = extract_attr(gtf$attribute, "gene_id")
) %>%
  filter(!is.na(TXNAME), !is.na(GENEID)) %>%
  distinct()

# -----------------------------
# tximport
# -----------------------------
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

# -----------------------------
# DESeq2
# -----------------------------
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ condition)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

contrast_vec <- strsplit(contrast_s, ",")[[1]]
res <- results(dds, contrast = contrast_vec)
res <- lfcShrink(dds, contrast = contrast_vec, res = res, type = "apeglm")

res_df <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  arrange(padj)

write_csv(res_df, file.path(outdir, "DESeq2_results.csv"))

# -----------------------------
# PCA
# -----------------------------
vsd <- vst(dds, blind = TRUE)
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))

p_pca <- ggplot(pca, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  theme_classic() +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%"))

ggsave(file.path(outdir, "PCA.png"), p_pca, width = 6, height = 5)

# -----------------------------
# Volcano plot
# -----------------------------
volcano <- res_df %>%
  mutate(
    sig = padj < padj_cutoff & abs(log2FoldChange) > lfc_cutoff,
    neglog10_padj = -log10(padj)
  )

p_volcano <- ggplot(volcano, aes(log2FoldChange, neglog10_padj, color = sig)) +
  geom_point(alpha = 0.7) +
  theme_classic() +
  labs(x = "log2 fold change", y = "-log10 adjusted p-value")

ggsave(file.path(outdir, "Volcano.png"), p_volcano, width = 6, height = 5)

# -----------------------------
# Heatmap (top variable genes)
# -----------------------------
mat <- assay(vsd)
vars <- apply(mat, 1, var)
top_genes <- names(sort(vars, decreasing = TRUE))[1:top_var_genes]

mat_scaled <- t(scale(t(mat[top_genes, ])))
anno <- data.frame(condition = meta$condition)
rownames(anno) <- rownames(meta)

png(file.path(outdir, "Heatmap_top_variable_genes.png"), 1600, 1400, res = 200)
pheatmap(mat_scaled, annotation_col = anno, show_rownames = FALSE)
dev.off()

# -----------------------------
# KEGG enrichment (GENE SYMBOLS)
# -----------------------------
sig_genes <- res_df %>%
  filter(padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff) %>%
  pull(gene_id)

mapped <- bitr(
  sig_genes,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Cf.eg.db
)

if (nrow(mapped) > 0) {
  ekegg <- enrichKEGG(
    gene = unique(mapped$ENTREZID),
    organism = kegg_org,
    pvalueCutoff = 0.05
  )

  write_csv(as.data.frame(ekegg),
            file.path(outdir, "KEGG_results.csv"))

  ggsave(
    file.path(outdir, "KEGG_dotplot.png"),
    dotplot(ekegg, showCategory = 20),
    width = 8, height = 6
  )
}

message("Pipeline finished successfully.")
