# Introductory public-data reanalysis, generated and executed with Codex assistance.
# Run from this directory: Rscript analyze.R
# This is a probe-level exploratory exercise, not a validated AD discovery study.
options(stringsAsFactors = FALSE)
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", ".UTF-8")
dir.create("rlib", showWarnings = FALSE)
.libPaths(c(normalizePath("rlib"), .libPaths()))
if (!requireNamespace("limma", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install("limma", ask = FALSE, update = FALSE)
}
library(limma)
dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
src <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE5nnn/GSE5281/matrix/GSE5281_series_matrix.txt.gz"
f <- "data/GSE5281_series_matrix.txt.gz"
if (!file.exists(f)) download.file(src, f, mode = "wb")
con <- gzfile(f, "rt", encoding = "UTF-8")
lines <- readLines(con, warn = FALSE)
close(con)
clean <- function(x) trimws(gsub("\u00a0", " ", gsub('^"|"$', "", x)))
fields <- function(line) clean(strsplit(line, "\t", fixed = TRUE)[[1]][-1])
get_one <- function(key) fields(lines[startsWith(lines, paste0(key, "\t"))][1])
titles <- get_one("!Sample_title")
ids <- get_one("!Sample_geo_accession")
chars <- lapply(lines[startsWith(lines, "!Sample_characteristics_ch1\t")], fields)
get_char <- function(key) {
  vapply(seq_along(ids), function(j) {
    vals <- vapply(chars, `[`, character(1), j)
    hit <- vals[startsWith(tolower(vals), paste0(tolower(key), ":"))]
    stopifnot(length(hit) == 1)
    trimws(sub("^[^:]+:", "", hit))
  }, character(1))
}
meta <- data.frame(sample = ids, title = titles,
  region = get_char("Organ Region"), disease = get_char("Disease State"),
  sex = tolower(get_char("Sex")), age_original = get_char("Age"))
meta <- meta[grepl("^HIP([ _])", meta$title), ]
meta$group <- factor(ifelse(grepl("control", meta$title), "control", "AD"),
                     levels = c("control", "AD"))
meta$sex <- factor(meta$sex)
meta$age <- ifelse(grepl("^[0-9.]+ years$", meta$age_original),
                  as.numeric(sub(" .*", "", meta$age_original)), NA_real_)
meta$primary_included <- is.finite(meta$age)
meta$exclusion_reason <- ifelse(meta$primary_included, "", "Age unit inconsistent with adult cohort; not silently corrected")
stopifnot(nrow(meta) == 23, !anyDuplicated(meta$sample),
          all(grepl("hippocampus", tolower(meta$region))),
          all(grepl("normal", tolower(meta$disease[meta$group == "control"]))),
          all(grepl("alzheimer", tolower(meta$disease[meta$group == "AD"]))))
start <- which(lines == "!series_matrix_table_begin")
end <- which(lines == "!series_matrix_table_end")
tab <- read.delim(text = paste(lines[(start + 1):(end - 1)], collapse = "\n"),
                  check.names = FALSE, row.names = 1)
x <- as.matrix(tab[, meta$sample, drop = FALSE])
stopifnot(all(is.finite(x)), all(x > 0), !anyDuplicated(rownames(x)))
raw_quantiles <- quantile(x, c(0, .25, .5, .75, .99, 1))
stopifnot(raw_quantiles[5] > 100) # This accession contains non-log signal values.
logx <- log2(x)
keep_probe <- !grepl("^AFFX", rownames(logx)) & apply(logx, 1, var) > 0
logx <- logx[keep_probe, , drop = FALSE]
primary <- meta$primary_included
primary_meta <- droplevels(meta[primary, ])
primary_x <- normalizeBetweenArrays(logx[, primary], method = "quantile")
run_fit <- function(expr, md, formula, label) {
  design <- model.matrix(formula, data = md)
  stopifnot(qr(design)$rank == ncol(design), ncol(expr) == nrow(md))
  fit <- eBayes(lmFit(expr, design), trend = TRUE)
  res <- topTable(fit, coef = "groupAD", number = Inf, adjust.method = "BH", sort.by = "P")
  res <- data.frame(probe_id = rownames(res), res, row.names = NULL)
  stopifnot(nrow(res) == nrow(expr), all(is.finite(res$P.Value)),
            all(res$adj.P.Val >= 0 & res$adj.P.Val <= 1))
  write.csv(res, paste0("results/", label, ".csv"), row.names = FALSE)
  res
}
res <- run_fit(primary_x, primary_meta, ~ group + age + sex, "primary_age_sex_adjusted")
no_norm <- run_fit(logx[, primary], primary_meta, ~ group + age + sex,
                   "sensitivity_no_additional_quantile_normalization")
same_sex <- run_fit(primary_x, primary_meta, ~ group + sex, "sensitivity_same_22_sex_adjusted")
all_sex <- run_fit(normalizeBetweenArrays(logx, method = "quantile"), meta,
                   ~ group + sex, "sensitivity_all_23_sex_adjusted")
write.csv(meta, "results/sample_metadata.csv", row.names = FALSE)
write.csv(head(res, 20), "results/top20_probes.csv", row.names = FALSE)
png("results/qc.png", width = 1800, height = 1300, res = 150)
par(mfrow = c(2, 2), mar = c(7, 4, 3, 1))
cols <- ifelse(primary_meta$group == "AD", "#d9534f", "#2878b5")
boxplot(logx[, primary], outline = FALSE, col = cols, las = 2, cex.axis = .5,
        main = "Submitted signal: log2", ylab = "log2 intensity")
boxplot(primary_x, outline = FALSE, col = cols, las = 2, cex.axis = .5,
        main = "After quantile normalization", ylab = "log2 intensity")
pca <- prcomp(t(primary_x), center = TRUE, scale. = FALSE)
vr <- 100 * pca$sdev^2 / sum(pca$sdev^2)
plot(pca$x[,1:2], col = cols, pch = 19,
     xlab = sprintf("PC1 (%.1f%%)", vr[1]), ylab = sprintf("PC2 (%.1f%%)", vr[2]),
     main = "PCA: all retained probes")
legend("topright", c("control", "AD"), col = c("#2878b5", "#d9534f"), pch = 19, bty = "n")
plot(res$logFC, -log10(res$P.Value), pch = 16, cex = .3,
     col = ifelse(res$adj.P.Val < .05, "#d9534f", "#999999"),
     xlab = "log2 fold change: AD - control", ylab = "-log10 unadjusted P",
     main = "Primary model; red = BH FDR < .05")
dev.off()
sig <- function(z) sum(z$adj.P.Val < .05)
ord <- match(res$probe_id, no_norm$probe_id)
summary <- c("GSE5281 hippocampal LCM neuronal expression: exploratory limma exercise",
  paste("Downloaded source:", src),
  paste("Input MD5:", unname(tools::md5sum(f))),
  paste("Raw signal quantiles:", paste(names(raw_quantiles), round(raw_quantiles, 3), collapse = "; ")),
  sprintf("Primary samples: %d control, %d AD", sum(primary_meta$group == "control"), sum(primary_meta$group == "AD")),
  paste("Excluded ambiguous age:", paste(meta$sample[!primary], meta$age_original[!primary], collapse = "; ")),
  paste("Tested probes:", nrow(res)),
  paste("Primary BH FDR < .05 probes:", sig(res)),
  paste("Without additional quantile normalization BH FDR < .05:", sig(no_norm)),
  paste("Same 22 samples, sex only BH FDR < .05:", sig(same_sex)),
  paste("All 23 samples, sex only BH FDR < .05:", sig(all_sex)),
  sprintf("LogFC Spearman correlation, with vs without additional normalization: %.4f", cor(res$logFC, no_norm$logFC[ord], method = "spearman")),
  "No gene-level aggregation, batch correction, external validation, or causal inference performed.",
  "Source processed intensities used; original normalization history is not documented in the matrix processing field.",
  "No claim of replicating Wan et al.; this is a separate, much smaller microarray exercise.")
writeLines(summary, "results/summary.txt")
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
cat(paste(summary, collapse = "\n"), "\n")
