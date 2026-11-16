# GSE5281 hippocampal expression: introductory limma exercise

This small public-data exercise was generated and executed to support preparation for an Alzheimer’s disease research internship. It is not prior independent research experience, a reproduction of Wan et al., a validated biomarker study, or evidence of proficiency in multi-omics. Review and understand the analysis before presenting it as a work sample.

## Research question and data

Which **probe-level expression associations** with Alzheimer’s disease appear in the hippocampal subset of GSE5281 after adjustment for age and sex, and how sensitive are estimates to additional quantile normalization?

- [NCBI GEO GSE5281](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5281)
- [Processed series matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE5nnn/GSE5281/matrix/GSE5281_series_matrix.txt.gz)
- Affymetrix GPL570 microarray, laser-capture microdissected neuronal samples; not RNA-seq counts or whole-brain bulk tissue.
- Only hippocampal samples are used. Pooling brain regions could reuse donors and confound region with disease. Donor identifiers are not independently validated in this exercise.
- The selected subset contains 13 control and 10 affected samples. The disease metadata agrees with these title-derived groups.
- GSM119628 has age `85 days` despite an adult cohort designation. The script preserves that value and excludes it from the age-adjusted model, without guessing a correction. A separate model includes all 23 samples without an age term.

## Run

From the directory containing `analyze.R`:

```r
source("analyze.R")
```

Or `Rscript analyze.R`. An internet connection is needed on first run. The script installs missing limma dependencies into `rlib/`, downloads the public matrix into `data/`, and writes results into `results/`. The delivered run used R 4.6.1 on Windows; see `results/sessionInfo.txt` for package versions. A fresh installation on other platforms has not been tested.

Downloaded input MD5: `15b785f0e1db5e223337fe8785313363`. The matrix contains positive, non-log intensity values. Its sample data-processing fields do not specify the original normalization history.

## Analysis choices

1. Parse sample annotations and select the hippocampus before fitting models.
2. Transform submitted intensities using log2; remove AFFX control probes and zero-variance probes.
3. Primary analysis: quantile normalization of the 22 retained arrays; linear model `~ group + age + sex`; limma empirical Bayes with `trend=TRUE`. Positive logFC means higher expression in AD relative to controls.
4. Benjamini-Hochberg adjustment across all 54,613 retained probes. Probe IDs are preserved; multiple probes can map to the same gene, and no gene-level conclusion is made.
5. Sensitivity analyses: omit additional quantile normalization with the same model and samples; omit age with the same 22 samples; and use all 23 samples with sex adjustment only. The last comparison changes sample inclusion as well as normalization inputs and is not a pure single-factor comparison.
6. Inspect distribution plots, PCA and a volcano plot; retain all samples except the prespecified ambiguous-age exclusion. PCA alone cannot establish whether separation is biological or technical.

## Actual run results

| Model | Samples | Probes with BH FDR < 0.05 |
|---|---:|---:|
| Primary: age + sex, added quantile normalization | 12 control / 10 AD | 6,171 |
| Age + sex, no additional quantile normalization | 12 control / 10 AD | 7,055 |
| Sex only, same primary sample set | 12 control / 10 AD | 6,115 |
| Sex only, all hippocampal samples | 13 control / 10 AD | 6,684 |

The Spearman correlation of logFC estimates with versus without additional quantile normalization was 0.9720. This is sensitivity evidence, not external replication. The number of significant probes is not a count of disease genes or a measure of analysis quality.

## Limitations and next steps

This is a small, observational, single-region reanalysis. Age and sex do not resolve unmeasured confounding. Batch, RNA quality, postmortem interval and donor provenance have not been sufficiently established here. Diagnostic-group separation in PCA does not rule out technical confounding. Quantile normalization assumes broadly comparable expression distributions and can remove real global differences; hence the sensitivity analysis. The source matrix may already reflect upstream processing, which was not reconstructed from raw CEL files.

Before making biological claims: resolve the anomalous age with original study documentation, obtain technical covariates and donor mapping, inspect raw-array quality, validate probe annotations, define a justified expression filter, and replicate results in an independent dataset. No batch correction, functional enrichment or external validation is claimed.

## Connection to Sage's work

[Wan et al. (2020), Cell Reports 32:107908](https://doi.org/10.1016/j.celrep.2020.107908) harmonizes AD definitions across ROSMAP, MSBB and Mayo before cross-cohort analysis. Its methods note that CERAD encodings differ across studies. Its discussion also cautions that bulk-tissue expression can reflect cell-proportion differences. Those lessons motivate careful annotation checks and constrained interpretation here; this small microarray exercise does not reproduce that multi-cohort study.

## Files

- `analyze.R`: reproducible analysis, download instructions and input checks.
- `results/summary.txt`: actual run counts and input checksum.
- `results/sample_metadata.csv`: all 23 annotations, inclusion flags and original age strings.
- `results/*adjusted.csv` and normalization sensitivity CSV: full probe-level statistics.
- `results/top20_probes.csv`, `results/qc.png`, `results/sessionInfo.txt`.

The distribution package intentionally excludes downloaded data and installed packages. GitHub publication has not been completed; no public repository URL is claimed.
