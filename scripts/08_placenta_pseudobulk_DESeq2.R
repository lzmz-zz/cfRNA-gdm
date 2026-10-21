# ============================================================
# DESeq2 pseudobulk DEG analysis
# ============================================================
# 目的：
#   每个细胞类型内，把同一样本的细胞 counts 求和成 pseudobulk，
#   然后用 DESeq2 比较 GDM_Group vs healthy。
#
# 结果方向：
#   log2FoldChange > 0 代表 GDM_Group 高于 healthy
#   log2FoldChange < 0 代表 GDM_Group 低于 healthy
#
# 输出文件：
#   01_DESeq2_ALL_results.csv
#   02_DESeq2_summary_by_celltype.csv
#   03_DESeq2_skipped_celltypes.csv
#   04_DESeq2_celltype_sample_QC.csv
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(DESeq2)
})

# ============================================================
# 1. 参数区：以后主要只改这里
# ============================================================

rds_file <- file.path("data", "controlled", "placenta_scRNA_seurat.rds")

outdir <- file.path("results", "placenta", "pseudobulk_DESeq2")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

assay_name <- "RNA"

# metadata 列名
celltype_col <- "dalei_5"      # 细胞类型列
sample_col   <- "newsampleID"  # 样本 ID 列
group_col    <- "Group"        # 分组列：healthy / GDM_Group
batch_col    <- "dataname"     # 批次 / 数据来源列

# 分组名称
group_ref  <- "healthy"        # 对照组
group_case <- "GDM_Group"      # 疾病组

# 分组别名统一
ref_alias  <- c("healthy", "Healthy", "HC", "Control", "control")
case_alias <- c("GDM_Group", "GDM", "GDMM", "gdm", "GDM group")

# 是否在模型中加入 batch
use_batch <- TRUE

# pseudobulk 样本过滤
min_cells_per_sample_in_celltype <- 20
min_samples_per_group <- 3

# 模型安全阈值
min_residual_df <- 2

# DESeq2 低表达基因预过滤
# 含义：
#   某基因至少在 min_samples_expr 个 pseudobulk 样本中 counts >= min_gene_count 才保留
min_gene_count <- 10

# NULL 表示自动使用当前细胞类型中较小组的样本数
# 例如 healthy=5，GDM_Group=4，则 min_samples_expr 自动为 4
min_samples_expr <- NULL

# 四个独立阈值，不是组合阈值
p_cutoffs <- c(
  P005 = 0.05,
  P001 = 0.01
)

fdr_cutoffs <- c(
  FDR010 = 0.10,
  FDR005 = 0.05
)

# ============================================================
# 2. 工具函数
# ============================================================

standardize_group <- function(x) {
  x <- as.character(x)
  x[x %in% ref_alias]  <- group_ref
  x[x %in% case_alias] <- group_case
  return(x)
}

get_counts_safe <- function(obj, assay_name) {
  
  assay_name <- as.character(assay_name)
  
  # 优先使用 Seurat v5 的 LayerData
  mat <- tryCatch(
    {
      LayerData(
        object = obj,
        assay = assay_name,
        layer = "counts"
      )
    },
    error = function(e) {
      message("LayerData failed, trying GetAssayData(slot = 'counts')...")
      message("LayerData error: ", e$message)
      return(NULL)
    }
  )
  
  # 如果 LayerData 失败，再尝试 Seurat v4 写法
  if (is.null(mat)) {
    mat <- tryCatch(
      {
        GetAssayData(
          object = obj,
          assay = assay_name,
          slot = "counts"
        )
      },
      error = function(e) {
        message("GetAssayData also failed.")
        message("GetAssayData error: ", e$message)
        return(NULL)
      }
    )
  }
  
  if (is.null(mat)) {
    stop("Cannot extract raw counts from assay: ", assay_name)
  }
  
  return(mat)
}

make_skip <- function(celltype_name, reason, detail = NA_character_) {
  data.frame(
    celltype = celltype_name,
    reason = reason,
    detail = detail,
    stringsAsFactors = FALSE
  )
}

add_threshold_flags_DESeq2 <- function(res, p_cutoffs, fdr_cutoffs) {
  
  for (nm in names(p_cutoffs)) {
    cutoff <- p_cutoffs[[nm]]
    res[[paste0("sig_", nm)]] <- !is.na(res$pvalue) & res$pvalue < cutoff
  }
  
  for (nm in names(fdr_cutoffs)) {
    cutoff <- fdr_cutoffs[[nm]]
    res[[paste0("sig_", nm)]] <- !is.na(res$padj) & res$padj < cutoff
  }
  
  return(res)
}

make_summary_DESeq2 <- function(res, celltype_name, pb_meta, design_used, residual_df,
                                p_cutoffs, fdr_cutoffs,
                                min_gene_count, min_samples_expr_used) {
  
  summary_df <- data.frame(
    celltype = celltype_name,
    design_used = design_used,
    n_samples_total = nrow(pb_meta),
    n_ref = sum(pb_meta$group == group_ref),
    n_case = sum(pb_meta$group == group_case),
    n_batch = dplyr::n_distinct(pb_meta$batch),
    residual_df = residual_df,
    genes_tested = nrow(res),
    min_gene_count = min_gene_count,
    min_samples_expr_used = min_samples_expr_used,
    stringsAsFactors = FALSE
  )
  
  # pvalue 阈值统计
  for (nm in names(p_cutoffs)) {
    flag_col <- paste0("sig_", nm)
    
    summary_df[[paste0("DEG_", nm)]] <- sum(res[[flag_col]], na.rm = TRUE)
    summary_df[[paste0("up_", nm)]] <- sum(res[[flag_col]] & res$log2FoldChange > 0, na.rm = TRUE)
    summary_df[[paste0("down_", nm)]] <- sum(res[[flag_col]] & res$log2FoldChange < 0, na.rm = TRUE)
  }
  
  # padj / FDR 阈值统计
  for (nm in names(fdr_cutoffs)) {
    flag_col <- paste0("sig_", nm)
    
    summary_df[[paste0("DEG_", nm)]] <- sum(res[[flag_col]], na.rm = TRUE)
    summary_df[[paste0("up_", nm)]] <- sum(res[[flag_col]] & res$log2FoldChange > 0, na.rm = TRUE)
    summary_df[[paste0("down_", nm)]] <- sum(res[[flag_col]] & res$log2FoldChange < 0, na.rm = TRUE)
  }
  
  return(summary_df)
}

# ============================================================
# 3. 读取 Seurat 对象
# ============================================================

cat("========== Reading Seurat object ==========\n")

obj <- readRDS(rds_file)

# ------------------------------------------------------------
# 安全检查 assay 是否存在
# 不使用 Assays(obj)，避免某些环境下触发 match 报错
# ------------------------------------------------------------

assay_name <- as.character(assay_name)

assays_available <- as.character(names(obj@assays))

cat("Available assays:\n")
print(assays_available)

if (!assay_name %in% assays_available) {
  stop(
    "Assay not found: ", assay_name,
    "\nAvailable assays are: ",
    paste(assays_available, collapse = ", ")
  )
}

DefaultAssay(obj) <- assay_name

cat("Object loaded:\n")
print(obj)

cat("Default assay:\n")
print(DefaultAssay(obj))

# ============================================================
# 4. 整理 metadata
# ============================================================

cat("\n========== Preparing metadata ==========\n")

meta_raw <- obj@meta.data

need_cols <- c(celltype_col, sample_col, group_col, batch_col)
miss_cols <- setdiff(need_cols, colnames(meta_raw))

if (length(miss_cols) > 0) {
  stop("Missing metadata columns: ", paste(miss_cols, collapse = ", "))
}

meta <- meta_raw
meta$celltype <- as.character(meta[[celltype_col]])
meta$sample   <- as.character(meta[[sample_col]])
meta$group    <- standardize_group(meta[[group_col]])
meta$batch    <- as.character(meta[[batch_col]])

# 去掉关键信息缺失的细胞，只保留 healthy 和 GDM_Group
meta <- meta %>%
  dplyr::filter(
    !is.na(celltype),
    !is.na(sample),
    !is.na(group),
    !is.na(batch),
    group %in% c(group_ref, group_case)
  )

meta$group <- factor(meta$group, levels = c(group_ref, group_case))

cat("Cells kept:", nrow(meta), "\n")

cat("Group distribution:\n")
print(table(meta$group, useNA = "ifany"))

cat("Celltype distribution:\n")
print(table(meta$celltype))

cat("Batch x group:\n")
print(table(meta$batch, meta$group))

# 检查一个 sample 是否唯一对应一个 group 和一个 batch
sample_check <- meta %>%
  dplyr::distinct(sample, group, batch) %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(
    n_group = dplyr::n_distinct(group),
    n_batch = dplyr::n_distinct(batch),
    .groups = "drop"
  )

bad_samples <- sample_check %>%
  dplyr::filter(n_group > 1 | n_batch > 1)

if (nrow(bad_samples) > 0) {
  print(bad_samples)
  stop("Some samples map to more than one group or batch. Please fix metadata first.")
}

cat("Sample-to-group/batch uniqueness: PASS\n")

# ============================================================
# 5. 提取 raw counts 并检查
# ============================================================

cat("\n========== Extracting raw counts ==========\n")

counts_all <- get_counts_safe(obj, assay_name)

cat("Raw counts dimension genes x cells:\n")
print(dim(counts_all))

# 检查 metadata 中细胞是否都在 counts 里
missing_cells <- setdiff(rownames(meta), colnames(counts_all))

if (length(missing_cells) > 0) {
  cat("Example missing cells:\n")
  print(head(missing_cells))
  stop("Some metadata cells are not found in counts matrix.")
}

# 按 metadata 顺序重排 counts
counts <- counts_all[, rownames(meta), drop = FALSE]

if (!identical(colnames(counts), rownames(meta))) {
  stop("Counts columns and metadata rows are not aligned.")
}

cat("Counts dimension after metadata filtering genes x cells:\n")
print(dim(counts))

nonzero_vals <- if (inherits(counts, "sparseMatrix")) {
  counts@x
} else {
  counts[counts != 0]
}

if (length(nonzero_vals) == 0) {
  stop("Counts matrix has no non-zero values.")
}

is_integer_like <- all(abs(nonzero_vals - round(nonzero_vals)) < 1e-8)
has_negative <- any(nonzero_vals < 0)

cat("Counts integer-like:", is_integer_like, "\n")
cat("Counts contain negative values:", has_negative, "\n")

if (!is_integer_like) {
  stop("Counts are not integer-like. DESeq2 requires raw integer counts.")
}

if (has_negative) {
  stop("Counts contain negative values.")
}

lib_sizes <- Matrix::colSums(counts)

cat("Cell library size summary:\n")
print(summary(lib_sizes))

if (any(lib_sizes <= 0)) {
  cat("Warning: some cells have zero library size.\n")
}

# ============================================================
# 6. 单个细胞类型 DESeq2 分析函数
# ============================================================

run_one_celltype_DESeq2 <- function(celltype_name, meta, counts) {
  
  cat("\n--------------------------------------------------\n")
  cat("Running celltype:", celltype_name, "\n")
  
  # 6.1 提取当前细胞类型
  meta_ct <- meta %>%
    dplyr::filter(celltype == celltype_name)
  
  if (nrow(meta_ct) == 0) {
    return(list(ok = FALSE, skip = make_skip(celltype_name, "no_cells")))
  }
  
  # 6.2 统计每个样本在当前细胞类型中的细胞数
  sample_cell_counts <- meta_ct %>%
    dplyr::count(sample, group, batch, name = "n_cells") %>%
    dplyr::filter(n_cells >= min_cells_per_sample_in_celltype)
  
  qc_df <- sample_cell_counts %>%
    dplyr::mutate(celltype = celltype_name) %>%
    dplyr::select(celltype, sample, group, batch, n_cells)
  
  if (nrow(sample_cell_counts) == 0) {
    return(list(
      ok = FALSE,
      skip = make_skip(celltype_name, "no_sample_passed_min_cells"),
      qc = qc_df
    ))
  }
  
  # 6.3 检查两组样本数是否足够
  group_n <- sample_cell_counts %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(n_samples = dplyr::n_distinct(sample), .groups = "drop")
  
  if (!all(c(group_ref, group_case) %in% group_n$group)) {
    return(list(
      ok = FALSE,
      skip = make_skip(celltype_name, "not_both_groups_present"),
      qc = qc_df
    ))
  }
  
  if (any(group_n$n_samples[group_n$group %in% c(group_ref, group_case)] < min_samples_per_group)) {
    return(list(
      ok = FALSE,
      skip = make_skip(celltype_name, "too_few_samples_in_one_group"),
      qc = qc_df
    ))
  }
  
  # 6.4 保留合格样本中的细胞
  keep_samples <- sample_cell_counts$sample
  
  meta_ct2 <- meta_ct %>%
    dplyr::filter(sample %in% keep_samples)
  
  cells_use <- rownames(meta_ct2)
  counts_ct <- counts[, cells_use, drop = FALSE]
  
  # 6.5 pseudobulk：同一个 sample 内 counts 求和
  sample_ids <- sample_cell_counts$sample
  
  pb_list <- lapply(sample_ids, function(s) {
    cells_s <- rownames(meta_ct2)[meta_ct2$sample == s]
    Matrix::rowSums(counts_ct[, cells_s, drop = FALSE])
  })
  
  pb_mat <- do.call(cbind, pb_list)
  rownames(pb_mat) <- rownames(counts_ct)
  colnames(pb_mat) <- sample_ids
  
  # DESeq2 需要普通 matrix，且必须是整数
  pb_mat <- as.matrix(pb_mat)
  
  if (!all(abs(pb_mat - round(pb_mat)) < 1e-8)) {
    return(list(
      ok = FALSE,
      skip = make_skip(celltype_name, "pseudobulk_counts_not_integer"),
      qc = qc_df
    ))
  }
  
  pb_mat <- round(pb_mat)
  storage.mode(pb_mat) <- "integer"
  
  # 6.6 pseudobulk metadata
  pb_meta <- sample_cell_counts %>%
    dplyr::select(sample, group, batch, n_cells)
  
  pb_meta <- pb_meta[match(colnames(pb_mat), pb_meta$sample), , drop = FALSE]
  
  if (!identical(pb_meta$sample, colnames(pb_mat))) {
    stop("pb_meta and pb_mat are not aligned in celltype: ", celltype_name)
  }
  
  pb_meta$group <- factor(pb_meta$group, levels = c(group_ref, group_case))
  pb_meta$batch <- factor(pb_meta$batch)
  rownames(pb_meta) <- pb_meta$sample
  
  # 6.7 构建设计矩阵
  if (use_batch && nlevels(droplevels(pb_meta$batch)) >= 2) {
    design_formula <- ~ batch + group
    design_used <- "~ batch + group"
  } else {
    design_formula <- ~ group
    design_used <- "~ group"
  }
  
  design <- model.matrix(design_formula, data = pb_meta)
  rownames(design) <- pb_meta$sample
  
  # 6.8 检查 design 是否满秩
  design_rank <- qr(design)$rank
  
  if (design_rank < ncol(design)) {
    return(list(
      ok = FALSE,
      skip = make_skip(
        celltype_name,
        "design_not_full_rank",
        "batch and group may be confounded"
      ),
      qc = qc_df
    ))
  }
  
  residual_df <- nrow(design) - design_rank
  
  if (residual_df < min_residual_df) {
    return(list(
      ok = FALSE,
      skip = make_skip(
        celltype_name,
        "too_low_residual_df",
        paste0("residual_df=", residual_df)
      ),
      qc = qc_df
    ))
  }
  
  # 6.9 DESeq2 低表达基因预过滤
  group_table <- table(pb_meta$group)
  
  min_samples_expr_used <- if (is.null(min_samples_expr)) {
    min(as.numeric(group_table))
  } else {
    min_samples_expr
  }
  
  keep_gene <- rowSums(pb_mat >= min_gene_count) >= min_samples_expr_used
  pb_mat_filt <- pb_mat[keep_gene, , drop = FALSE]
  
  if (nrow(pb_mat_filt) == 0) {
    return(list(
      ok = FALSE,
      skip = make_skip(celltype_name, "no_genes_after_prefilter"),
      qc = qc_df
    ))
  }
  
  # 6.10 DESeq2 分析
  dds <- DESeqDataSetFromMatrix(
    countData = pb_mat_filt,
    colData = pb_meta,
    design = design_formula
  )
  
  dds <- DESeq(dds, quiet = TRUE)
  
  res <- results(
    dds,
    contrast = c("group", group_case, group_ref),
    alpha = max(fdr_cutoffs)
  )
  
  res <- as.data.frame(res)
  res$gene <- rownames(res)
  res$celltype <- celltype_name
  res$design_used <- design_used
  res$direction <- ifelse(
    res$log2FoldChange > 0,
    paste0("up_in_", group_case),
    paste0("down_in_", group_case)
  )
  
  # 添加四个独立阈值判断
  res <- add_threshold_flags_DESeq2(
    res = res,
    p_cutoffs = p_cutoffs,
    fdr_cutoffs = fdr_cutoffs
  )
  
  res <- res %>%
    dplyr::select(
      celltype,
      gene,
      baseMean,
      log2FoldChange,
      lfcSE,
      stat,
      pvalue,
      padj,
      direction,
      design_used,
      dplyr::starts_with("sig_"),
      dplyr::everything()
    )
  
  # 6.11 汇总
  summary_df <- make_summary_DESeq2(
    res = res,
    celltype_name = celltype_name,
    pb_meta = pb_meta,
    design_used = design_used,
    residual_df = residual_df,
    p_cutoffs = p_cutoffs,
    fdr_cutoffs = fdr_cutoffs,
    min_gene_count = min_gene_count,
    min_samples_expr_used = min_samples_expr_used
  )
  
  return(list(
    ok = TRUE,
    result = res,
    summary = summary_df,
    qc = qc_df
  ))
}

# ============================================================
# 7. 循环所有细胞类型
# ============================================================

cat("\n========== Running all celltypes ==========\n")

celltypes_to_run <- sort(unique(meta$celltype))

all_results <- list()
all_summary <- list()
all_skipped <- list()
all_qc <- list()

for (ct in celltypes_to_run) {
  
  out <- tryCatch(
    run_one_celltype_DESeq2(ct, meta, counts),
    error = function(e) {
      list(
        ok = FALSE,
        skip = make_skip(ct, "unexpected_error", e$message),
        qc = NULL
      )
    }
  )
  
  if (!is.null(out$qc)) {
    all_qc[[ct]] <- out$qc
  }
  
  if (isTRUE(out$ok)) {
    all_results[[ct]] <- out$result
    all_summary[[ct]] <- out$summary
  } else {
    all_skipped[[ct]] <- out$skip
  }
}

# ============================================================
# 8. 输出结果
# ============================================================

cat("\n========== Writing outputs ==========\n")

if (length(all_results) > 0) {
  res_all <- dplyr::bind_rows(all_results)
  write.csv(
    res_all,
    file.path(outdir, "01_DESeq2_ALL_results.csv"),
    row.names = FALSE
  )
} else {
  cat("No successful DESeq2 results.\n")
}

if (length(all_summary) > 0) {
  summary_all <- dplyr::bind_rows(all_summary)
  write.csv(
    summary_all,
    file.path(outdir, "02_DESeq2_summary_by_celltype.csv"),
    row.names = FALSE
  )
  print(summary_all)
} else {
  cat("No successful summary results.\n")
}

skipped_all <- if (length(all_skipped) > 0) {
  dplyr::bind_rows(all_skipped)
} else {
  data.frame(celltype = character(), reason = character(), detail = character())
}

write.csv(
  skipped_all,
  file.path(outdir, "03_DESeq2_skipped_celltypes.csv"),
  row.names = FALSE
)

qc_all <- if (length(all_qc) > 0) {
  dplyr::bind_rows(all_qc)
} else {
  data.frame(celltype = character(), sample = character(), group = character(), batch = character(), n_cells = numeric())
}

write.csv(
  qc_all,
  file.path(outdir, "04_DESeq2_celltype_sample_QC.csv"),
  row.names = FALSE
)

cat("\n========== Finished DESeq2 pseudobulk analysis ==========\n")
cat("Output directory:", outdir, "\n")