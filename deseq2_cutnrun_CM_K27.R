#!/usr/bin/env Rscript
# ============================================================================
# DESeq2 CUT&RUN Differential Peak Analysis
# ============================================================================
#
# Reads pre-made union consensus peaks, counts reads with featureCounts,
# runs DESeq2 with LFC shrinkage, and exports statistics at multiple
# thresholds for downstream GO analysis.
#
# ============================================================================
 
BASE_DIR   <- "/home/dpresasramos/gen_occ/mdcnr_final"
 
# Path to directory containing nf-core BAM files
BAM_DIR    <- "/home/dpresasramos/gen_occ/cutnrun_jan_md/results/02_alignment/bowtie2/target/markdup"
 
# Directory containing your pre-made union consensus BED files
UNION_PEAK_DIR <- file.path(BASE_DIR, "00_data/bed_peaks/consensus/union")
 
OUTPUT_DIR <- file.path(BASE_DIR, "01_DESeq2/CM_H3K27me3")
 
# ---- Experiment design ----
CELL_TYPE  <- "CM"
ANTIBODY   <- "K27"
 
# ---- Condition labels ----
CONDITION_LABELS <- list(
 cnd1 = "saline-sham",
 cnd2 = "malonate-sham",
 cnd3 = "saline-MI",
 cnd4 = "malonate-MI"
)
 
COMPARISONS <- list(
  c("cnd4", "cnd3"),
  c("cnd2", "cnd1"),
  c("cnd4", "cnd2"),
  c("cnd3", "cnd1")
)


# ---- Analysis parameters ----
N_THREADS  <- 8
PAIRED_END <- TRUE
 
# Minimum read count filter: keep peaks with >= MIN_COUNT reads in >= MIN_SAMPLES samples
MIN_COUNT   <- 5
MIN_SAMPLES <- 2
 
# ---- Multiple threshold sets ----
# Each list entry is one threshold combination that will be run and reported.
# raw_pval = TRUE uses pvalue instead of padj (exploratory)
# lfc_test = if > 0, tests H0:|LFC|<=lfc_test (built into DESeq2 test, more rigorous)
# lfc_filter = post-hoc absolute LFC filter applied to results table
THRESHOLD_SETS <- list(
 
  # --- Standard  ---
  list(name        = "padj05_lfc1",
       use_raw     = FALSE,
       p_threshold  = 0.05,
       correction  = "BH",         # Benjamini-Hochberg (DESeq2 default)
       lfc_test    = 1.0,           # test H0: |LFC| <= 1
       lfc_filter  = 1.0),
 
  # --- Relaxed adjusted p (low-n exploratory) ---
  list(name        = "padj10_lfc0",
       use_raw     = FALSE,
       p_threshold  = 0.10,
       correction  = "BH",
       lfc_test    = 0,
       lfc_filter  = 0),
 
  # --- Raw p-value ---
  list(name        = "rawp01_lfc1",
       use_raw     = TRUE,
       p_threshold  = 0.01,
       correction  = "none",
       lfc_test    = 1.0,
       lfc_filter  = 1.0),
 
  list(name        = "rawp05_lfc0",
       use_raw     = TRUE,
       p_threshold  = 0.05,
       correction  = "none",
       lfc_test    = 0,
       lfc_filter  = 0),
 
  # --- Bonferroni (most conservative; useful to show robustness) ---
  list(name        = "bonf05_lfc1",
       use_raw     = FALSE,
       p_threshold  = 0.05,
       correction  = "bonferroni",
       lfc_test    = 1.0,
       lfc_filter  = 1.0),
 

  list(name        = "padj05_lfc05_MLE",
       use_raw     = FALSE,
       p_threshold  = 0.05,
       correction  = "BH",
       lfc_test    = 0,       # no threshold in DESeq2 test itself
       lfc_filter  = 0.5,     # post-hoc: |MLE LFC| > 0.5
       use_mle     = TRUE),   # flag: filter on MLE not shrunken LFC
 
  list(name        = "padj05_lfc03_MLE",
       use_raw     = FALSE,
       p_threshold  = 0.05,
       correction  = "BH",
       lfc_test    = 0,
       lfc_filter  = 0.3,
       use_mle     = TRUE),
 
  # raw p < 0.1, MLE LFC >= 0.3, no correction
  list(name        = "rawp10_lfc03_MLE",
       use_raw     = TRUE,
       p_threshold  = 0.1,
       correction  = "none",
       lfc_test    = 0,
       lfc_filter  = 0.3,
       use_mle     = TRUE)
)
 
# ============================================================================
# SAMPLE METADATA + BAM LOOKUP
# ============================================================================
 
# Maps human-readable sample IDs -> original nf-core BAM filenames
# Hardcoded from sample_key_RID.xlsx - no BAM renaming needed
BAM_LOOKUP <- list(
  # H3K4me3 CM
  "CM2_K4_Rep1"  = "h3k4me3_R1.target.markdup.sorted.bam",
  "CM5_K4_Rep1"  = "h3k4me3_R2.target.markdup.sorted.bam",
  "CM8_K4_Rep1"  = "h3k4me3_R3.target.markdup.sorted.bam",
  "CM11_K4_Rep1" = "h3k4me3_R4.target.markdup.sorted.bam",
  "CM14_K4_Rep2" = "h3k4me3_R5.target.markdup.sorted.bam",
  "CM17_K4_Rep2" = "h3k4me3_R6.target.markdup.sorted.bam",
  "CM20_K4_Rep2" = "h3k4me3_R7.target.markdup.sorted.bam",
  "CM23_K4_Rep2" = "h3k4me3_R8.target.markdup.sorted.bam",
  # H3K4me3 CF
  "CF2_K4_Rep1"  = "h3k4me3_R9.target.markdup.sorted.bam",
  "CF5_K4_Rep1"  = "h3k4me3_R10.target.markdup.sorted.bam",
  "CF8_K4_Rep1"  = "h3k4me3_R11.target.markdup.sorted.bam",
  "CF11_K4_Rep1" = "h3k4me3_R12.target.markdup.sorted.bam",
  "CF14_K4_Rep2" = "h3k4me3_R13.target.markdup.sorted.bam",
  "CF17_K4_Rep2" = "h3k4me3_R14.target.markdup.sorted.bam",
  "CF20_K4_Rep2" = "h3k4me3_R15.target.markdup.sorted.bam",
  "CF23_K4_Rep2" = "h3k4me3_R16.target.markdup.sorted.bam",
  # H3K27me3 CM
  "CM3_K27_Rep1"  = "h3k27me3_R1.target.markdup.sorted.bam",
  "CM6_K27_Rep1"  = "h3k27me3_R2.target.markdup.sorted.bam",
  "CM9_K27_Rep1"  = "h3k27me3_R3.target.markdup.sorted.bam",
  "CM12_K27_Rep1" = "h3k27me3_R4.target.markdup.sorted.bam",
  "CM15_K27_Rep2" = "h3k27me3_R5.target.markdup.sorted.bam",
  "CM18_K27_Rep2" = "h3k27me3_R6.target.markdup.sorted.bam",
  "CM21_K27_Rep2" = "h3k27me3_R7.target.markdup.sorted.bam",
  "CM24_K27_Rep2" = "h3k27me3_R8.target.markdup.sorted.bam",
  # H3K27me3 CF
  "CF3_K27_Rep1"  = "h3k27me3_R9.target.markdup.sorted.bam",
  "CF6_K27_Rep1"  = "h3k27me3_R10.target.markdup.sorted.bam",
  "CF9_K27_Rep1"  = "h3k27me3_R11.target.markdup.sorted.bam",
  "CF12_K27_Rep1" = "h3k27me3_R12.target.markdup.sorted.bam",
  "CF15_K27_Rep2" = "h3k27me3_R13.target.markdup.sorted.bam",
  "CF18_K27_Rep2" = "h3k27me3_R14.target.markdup.sorted.bam",
  "CF21_K27_Rep2" = "h3k27me3_R15.target.markdup.sorted.bam",
  "CF24_K27_Rep2" = "h3k27me3_R16.target.markdup.sorted.bam"
)
 
build_sample_info <- function(ct, ab) {
  # NOTE: function arguments named 'ct' and 'ab' to avoid shadowing
  # the column names 'cell_type' and 'antibody' in subset()
  all_samples <- data.frame(
    sample_id = c(
      "CM2_K4_Rep1",  "CM5_K4_Rep1",  "CM8_K4_Rep1",  "CM11_K4_Rep1",
      "CM14_K4_Rep2", "CM17_K4_Rep2", "CM20_K4_Rep2", "CM23_K4_Rep2",
      "CM3_K27_Rep1", "CM6_K27_Rep1", "CM9_K27_Rep1", "CM12_K27_Rep1",
      "CM15_K27_Rep2","CM18_K27_Rep2","CM21_K27_Rep2","CM24_K27_Rep2",
      "CF2_K4_Rep1",  "CF5_K4_Rep1",  "CF8_K4_Rep1",  "CF11_K4_Rep1",
      "CF14_K4_Rep2", "CF17_K4_Rep2", "CF20_K4_Rep2", "CF23_K4_Rep2",
      "CF3_K27_Rep1", "CF6_K27_Rep1", "CF9_K27_Rep1", "CF12_K27_Rep1",
      "CF15_K27_Rep2","CF18_K27_Rep2","CF21_K27_Rep2","CF24_K27_Rep2"
    ),
    cell_type = c(rep("CM", 16), rep("CF", 16)),
    antibody  = c(rep("K4", 8), rep("K27", 8), rep("K4", 8), rep("K27", 8)),
    condition = rep(c("cnd1","cnd2","cnd3","cnd4",
                      "cnd1","cnd2","cnd3","cnd4"), 4),
    replicate = rep(c(rep("Rep1", 4), rep("Rep2", 4)), 4),
    stringsAsFactors = FALSE
  )
  result <- all_samples[all_samples$cell_type == ct & all_samples$antibody == ab, ]
  if (nrow(result) == 0) {
    stop(paste0("No samples found for cell_type='", ct, "' antibody='", ab,
                "'. Check CELL_TYPE and ANTIBODY settings."))
  }
  result
}
 
# ============================================================================
# LOAD LIBRARIES
# ============================================================================
 
message("Loading libraries...")
suppressPackageStartupMessages({
  library(DESeq2)
  library(Rsubread)
  library(GenomicRanges)
  library(rtracklayer)
  library(ashr)
  library(ggplot2)
  library(ggrepel)
})
 
# ============================================================================
# SETUP
# ============================================================================
 
analysis_name <- paste0(CELL_TYPE, "_", ANTIBODY)
analysis_dir  <- file.path(OUTPUT_DIR, analysis_name)
dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)
 
# Sub-directories for per-threshold outputs
for (ts in THRESHOLD_SETS) {
  dir.create(file.path(analysis_dir, ts$name), showWarnings = FALSE)
}
 
log_file <- file.path(analysis_dir, paste0(analysis_name, "_run.log"))
message(paste("Logging to:", log_file))
sink(log_file, append = FALSE, split = TRUE)   # tee to both console and file
 
message("============================================================")
message(paste("DESeq2 CUT&RUN Analysis:", analysis_name))
message(paste("Started:", Sys.time()))
message("============================================================\n")
 
# ============================================================================
# STEP 1: Build sample info
# ============================================================================
 
message("=== Step 1: Building sample metadata ===")
sample_info <- build_sample_info(ct = CELL_TYPE, ab = ANTIBODY)
sample_info$condition <- factor(sample_info$condition)
message(paste("Samples:", nrow(sample_info)))
print(sample_info[, c("sample_id", "condition", "replicate")])
 
# ============================================================================
# STEP 2: Load pre-made union consensus peaks
# ============================================================================
 
message("\n=== Step 2: Loading union consensus peaks ===")
 
# Build the consensus peak filename from the CONDITION_LABELS
# Since union peaks span all conditions, we load one union file per
# cell-type × antibody combination.
# Expected file: {CELL_TYPE}_{ANTIBODY}_union.bed  (covers all conditions)
# If you have per-condition union files, adjust the pattern below.
# Build pattern to match all union files for this cell type + antibody combination.
# e.g. with CELL_TYPE="CM", ANTIBODY="K4": matches
#   CM_saline-sham_H3K4me3_union.bed
#   CM_saline-MI_H3K4me3_union.bed
#   CM_treated-sham_H3K4me3_union.bed
#   CM_treated-MI_H3K4me3_union.bed
# (one file per condition; all four are loaded and merged into the counting universe)
ab_full      <- ifelse(ANTIBODY == "K4", "H3K4me3", "H3K27me3")
union_pattern <- paste0("^", CELL_TYPE, "_.*_", ab_full, "_union\\.bed$")
 
union_files <- list.files(UNION_PEAK_DIR,
                          pattern    = union_pattern,
                          full.names = TRUE)
 
if (length(union_files) == 0) {
  stop(paste0(
    "No union consensus BED files found in: ", UNION_PEAK_DIR,
    "\nExpected pattern: ", union_pattern,
    "\nFiles present:\n  ",
    paste(list.files(UNION_PEAK_DIR), collapse = "\n  ")
  ))
}
 
message(paste("Union BED files found (one per condition):", length(union_files)))
invisible(lapply(basename(union_files), function(f) message("  ", f)))
 
# Combine all four condition-specific union BEDs into one peak universe.
# This is the correct input for DESeq2: count reads at every peak called
# in any condition, then let DESeq2 determine which peaks change between them.
all_peaks <- GRanges()
for (f in union_files) {
  peaks <- tryCatch(
    import(f, format = "BED"),
    error = function(e) {
      message(paste("  WARNING: could not read", basename(f), "-", e$message))
      return(GRanges())
    }
  )
  if (length(peaks) > 0) {
    # Add chr prefix if chromosome names are bare numbers (e.g. "1" -> "chr1")
    if (any(grepl("^[0-9XYMxy]+$", as.character(seqnames(peaks))))) {
      seqlevels(peaks) <- paste0("chr", seqlevels(peaks))
      seqlevels(peaks) <- gsub("^chrMT$", "chrM", seqlevels(peaks))
    }
    all_peaks <- c(all_peaks, peaks)
    message(paste("  Loaded:", basename(f), "-", length(peaks), "peaks"))
  }
}
 
# reduce() merges any overlapping intervals across the four union files
# so each genomic region is counted exactly once
consensus <- reduce(all_peaks)
message(paste("Total peak universe after merging all conditions:", length(consensus)))
 
# Create SAF annotation for featureCounts
saf <- data.frame(
  GeneID = paste0("peak_", seq_along(consensus)),
  Chr    = as.character(seqnames(consensus)),
  Start  = start(consensus),
  End    = end(consensus),
  Strand = "+"
)
 
saf_file <- file.path(analysis_dir, paste0(analysis_name, "_consensus.saf"))
write.table(saf, saf_file, sep = "\t", quote = FALSE, row.names = FALSE)
 
consensus_bed_file <- file.path(analysis_dir, paste0(analysis_name, "_consensus.bed"))
export(consensus, consensus_bed_file, format = "BED")
message(paste("Saved consensus BED:", consensus_bed_file))
 
# ============================================================================
# STEP 3: Count reads in peaks
# ============================================================================
 
message("\n=== Step 3: Counting reads in peaks (featureCounts) ===")
 
# Resolve BAM filenames from lookup table (original nf-core naming)
bam_files <- sapply(sample_info$sample_id, function(sid) {
  bam <- BAM_LOOKUP[[sid]]
  if (is.null(bam)) stop(paste("No BAM mapping found for sample:", sid))
  file.path(BAM_DIR, bam)
})
names(bam_files) <- sample_info$sample_id
 
message("BAM files resolved:")
for (i in seq_along(bam_files)) {
  exists_str <- if (file.exists(bam_files[i])) "OK" else "MISSING"
  message(paste0("  [", exists_str, "] ", names(bam_files)[i],
                 " -> ", basename(bam_files[i])))
}
 
missing_bams <- bam_files[!file.exists(bam_files)]
if (length(missing_bams) > 0) {
  stop(paste("Missing BAM files:\n ",
             paste(paste0(names(missing_bams), " -> ", missing_bams),
                   collapse = "\n  ")))
}
 
fc <- featureCounts(
  files                   = bam_files,
  annot.ext               = saf_file,
  isGTFAnnotationFile     = FALSE,
  isPairedEnd             = PAIRED_END,
  nthreads                = N_THREADS,
  countMultiMappingReads  = FALSE,
  fraction                = FALSE,
  verbose                 = FALSE
)
 
counts <- fc$counts
colnames(counts) <- sample_info$sample_id
 
counts_file <- file.path(analysis_dir, paste0(analysis_name, "_raw_counts.csv"))
write.csv(counts, counts_file)
message(paste("Raw counts saved:", counts_file))
message(paste("Count matrix:", nrow(counts), "peaks x", ncol(counts), "samples"))
 
# ============================================================================
# STEP 4: DESeq2 object + QC
# ============================================================================
 
message("\n=== Step 4: DESeq2 setup ===")
 
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = sample_info,
  design    = ~ condition
)
 
# Filter: keep peaks with >= 10 reads in at least 2 samples
keep <- rowSums(counts(dds) >= MIN_COUNT) >= MIN_SAMPLES
message(paste("Low-count filter: >=", MIN_COUNT, "reads in >=", MIN_SAMPLES, "samples"))
dds  <- dds[keep, ]
message(paste("Peaks after low-count filter:", nrow(dds)))
 
dds <- DESeq(dds, quiet = FALSE)
 
# ---- PCA (VST-transformed) ----
message("Generating PCA plot...")
vsd      <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = c("condition", "replicate"), returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))
 
p_pca <- ggplot(pca_data,
                aes(PC1, PC2, color = condition, shape = replicate)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = name), size = 3) +
  xlab(paste0("PC1: ", pct_var[1], "% variance")) +
  ylab(paste0("PC2: ", pct_var[2], "% variance")) +
  ggtitle(paste("PCA -", CELL_TYPE, ANTIBODY)) +
  theme_bw()
 
ggsave(file.path(analysis_dir, paste0(analysis_name, "_PCA.pdf")),
       p_pca, width = 8, height = 6)
message("PCA saved.")
 
# ============================================================================
# STEP 5: Pairwise comparisons with LFC shrinkage + multiple thresholds
# ============================================================================
 
message("\n=== Step 5: Pairwise comparisons ===")
 
# Master summary table (one row per comparison × threshold)
master_summary <- data.frame()
 
apply_thresholds <- function(res_df, threshold_set) {
  ts <- threshold_set
 
  # use_mle flag: filter on raw MLE LFC instead of shrunken LFC
  # Used for the "old script equivalent" threshold sets
  use_mle <- isTRUE(ts$use_mle)
  lfc_col <- if (use_mle) "log2FoldChange_MLE" else "log2FoldChange_shrunk"
 
  p_col <- if (ts$use_raw) "pvalue" else "padj_corrected"
 
  if (ts$correction == "bonferroni") {
    res_df$padj_corrected <- p.adjust(res_df$pvalue, method = "bonferroni")
  } else if (ts$correction == "BH") {
    res_df$padj_corrected <- p.adjust(res_df$pvalue, method = "BH")
  } else {
    res_df$padj_corrected <- res_df$pvalue
  }
 
  sig_mask <- !is.na(res_df[[p_col]]) &
              res_df[[p_col]] < ts$p_threshold &
              abs(res_df[[lfc_col]]) >= ts$lfc_filter
 
  res_df$significant <- sig_mask
  res_df$direction   <- ifelse(!sig_mask, "NS",
                         ifelse(res_df[[lfc_col]] > 0, "Up", "Down"))
  return(res_df)
}
 
for (comp in COMPARISONS) {
  test_cnd <- comp[1]
  ref_cnd  <- comp[2]
  test_label <- CONDITION_LABELS[[test_cnd]]
  ref_label  <- CONDITION_LABELS[[ref_cnd]]
  comp_name  <- paste0(test_label, "_vs_", ref_label)
 
  message(paste("\n--- Comparison:", comp_name, "---"))
 
  # ---- Base results (lfc_test may vary per threshold set; run separately) ----
  # We run the core result once without an lfc_test for the shared base,
  # then re-run with lfcThreshold where needed.
 
  res_base <- tryCatch(
    results(dds,
            contrast          = c("condition", test_cnd, ref_cnd),
            independentFiltering = TRUE),
    error = function(e) {
      message(paste("  ERROR:", e$message)); return(NULL)
    }
  )
  if (is.null(res_base)) next
 
  # ---- LFC shrinkage with ashr (recommended for n=2) ----
  # ashr doesn't require a specific contrast format - works on any result
  res_shrunk <- tryCatch(
    lfcShrink(dds,
              contrast = c("condition", test_cnd, ref_cnd),
              type     = "ashr",
              res      = res_base,
              quiet    = TRUE),
    error = function(e) {
      message(paste("  WARNING: lfcShrink failed, using MLE LFC:", e$message))
      return(NULL)
    }
  )
 
  # Combine base stats with shrunk LFC
  res_df <- as.data.frame(res_base)
  res_df$peak_id <- rownames(res_df)
  res_df$log2FoldChange_MLE    <- res_df$log2FoldChange
  res_df$log2FoldChange_shrunk <- if (!is.null(res_shrunk)) {
    as.data.frame(res_shrunk)$log2FoldChange
  } else {
    res_df$log2FoldChange   # fall back to MLE
  }
  res_df$lfcSE_shrunk <- if (!is.null(res_shrunk)) {
    as.data.frame(res_shrunk)$lfcSE
  } else {
    res_df$lfcSE
  }
 
  # Add human-readable condition labels and peak coordinates
  res_df$test_condition <- test_label
  res_df$ref_condition  <- ref_label
  res_df$comparison     <- comp_name
  res_df$cell_type      <- CELL_TYPE
  res_df$antibody       <- ANTIBODY
  res_df <- merge(res_df, saf,
                  by.x = "peak_id", by.y = "GeneID", all.x = TRUE)
  res_df <- res_df[order(res_df$pvalue), ]
 
  # ---- Save full unfiltered results (always) ----
  full_file <- file.path(analysis_dir,
                         paste0(analysis_name, "_", comp_name, "_full_results.csv"))
  write.csv(res_df, full_file, row.names = FALSE)
  message(paste("  Full results saved:", basename(full_file)))
 
  # ---- Apply each threshold set ----
  for (ts in THRESHOLD_SETS) {
    ts_dir <- file.path(analysis_dir, ts$name)
 
    # Re-run results() with lfc_test built into the test if requested
    if (ts$lfc_test > 0) {
      res_lfc <- tryCatch(
        results(dds,
                contrast      = c("condition", test_cnd, ref_cnd),
                lfcThreshold  = ts$lfc_test,
                altHypothesis = "greaterAbs"),
        error = function(e) {
          message(paste("  WARNING lfc_test failed:", e$message))
          return(res_base)
        }
      )
      res_df$padj_lfc_test <- res_lfc$padj
      res_df$pvalue_lfc_test <- res_lfc$pvalue
    }
 
    res_ts <- apply_thresholds(res_df, ts)
 
    n_up   <- sum(res_ts$direction == "Up",   na.rm = TRUE)
    n_down <- sum(res_ts$direction == "Down",  na.rm = TRUE)
    n_sig  <- n_up + n_down
    message(paste0("  [", ts$name, "] UP=", n_up, "  DOWN=", n_down, "  TOTAL=", n_sig))
 
    # Save filtered significant-only results
    sig_df <- res_ts[res_ts$significant, ]
    sig_file <- file.path(ts_dir,
                          paste0(analysis_name, "_", comp_name, "_", ts$name, "_sig.csv"))
    write.csv(sig_df, sig_file, row.names = FALSE)
 
    # Save BED of significant peaks (for GO input / IGV)
    if (nrow(sig_df) > 0) {
      bed_up <- sig_df[sig_df$direction == "Up", ]
      bed_dn <- sig_df[sig_df$direction == "Down", ]
 
      save_bed <- function(df, suffix) {
        if (nrow(df) == 0) return()
        bed <- data.frame(
          chr    = df$Chr,
          start  = df$Start,
          end    = df$End,
          name   = df$peak_id,
          score  = round(-log10(pmax(df$pvalue, 1e-300)), 3),
          strand = ".",
          log2FC_shrunk = round(df$log2FoldChange_shrunk, 4)
        )
        out <- file.path(ts_dir,
                         paste0(analysis_name, "_", comp_name,
                                "_", ts$name, suffix))
        write.table(bed, out, sep = "\t", quote = FALSE,
                    row.names = FALSE, col.names = FALSE)
      }
      save_bed(bed_up, "_UP.bed")
      save_bed(bed_dn, "_DOWN.bed")
 
      # Gene list for GO: one gene symbol per line (requires ChIPseeker annotation
      # in next step - here we export peak IDs as placeholder)
      gene_file <- file.path(ts_dir,
                             paste0(analysis_name, "_", comp_name,
                                    "_", ts$name, "_peak_ids_for_GO.txt"))
      writeLines(sig_df$peak_id, gene_file)
    }
 
    # ---- Volcano plot (one per comparison × threshold) ----
    p_volcano <- ggplot(res_ts,
                        aes(log2FoldChange_shrunk,
                            -log10(pvalue),
                            color = direction)) +
      geom_point(alpha = 0.5, size = 0.8) +
      scale_color_manual(
        values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey70")) +
      geom_vline(xintercept = c(-ts$lfc_filter, ts$lfc_filter),
                 linetype = "dashed", color = "grey40", linewidth = 0.4) +
      labs(
        title    = paste0(CELL_TYPE, " ", ANTIBODY, " - ", comp_name),
        subtitle = paste0("Threshold: ", ts$name,
                          " | UP=", n_up, " DOWN=", n_down),
        x        = "Shrunken log2 Fold Change (ashr)",
        y        = "-log10(p-value)",
        color    = NULL
      ) +
      theme_bw(base_size = 11) +
      theme(plot.title    = element_text(face = "bold"),
            plot.subtitle = element_text(size = 9, color = "grey40"))
 
    ggsave(
      file.path(ts_dir, paste0(analysis_name, "_", comp_name,
                               "_", ts$name, "_volcano.pdf")),
      p_volcano, width = 7, height = 6)
 
    # Accumulate master summary row
    master_summary <- rbind(master_summary, data.frame(
      cell_type      = CELL_TYPE,
      antibody       = ANTIBODY,
      comparison     = comp_name,
      threshold_set  = ts$name,
      use_raw_pval   = ts$use_raw,
      p_threshold    = ts$p_threshold,
      correction     = ts$correction,
      lfc_test       = ts$lfc_test,
      lfc_filter     = ts$lfc_filter,
      total_peaks_tested = nrow(res_ts),
      sig_up         = n_up,
      sig_down       = n_down,
      sig_total      = n_sig,
      stringsAsFactors = FALSE
    ))
  }
}
 
# ============================================================================
# STEP 6: Master summary table
# ============================================================================
 
message("\n=== Step 6: Saving master summary ===")
 
summary_file <- file.path(analysis_dir,
                          paste0(analysis_name, "_threshold_comparison_summary.csv"))
write.csv(master_summary, summary_file, row.names = FALSE)
message(paste("Master summary saved:", summary_file))
 
message("\n--- Summary across all thresholds ---")
print(master_summary[, c("comparison", "threshold_set",
                          "sig_up", "sig_down", "sig_total")])
 
message("\n============================================================")
message("Analysis complete!")
message(paste("Finished:", Sys.time()))
message(paste("All outputs in:", analysis_dir))
message("============================================================")
sink()