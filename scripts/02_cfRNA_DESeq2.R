# ============================================================
# DESeq2 最终简化版
# A/B/C 三期分别分析 GDM vs Healthy
#
# 低表达过滤：rowSums(count >= 10) >= 3
# 模型：
#   M3_basic          = ~ group
#   M4_age_collectGA = ~ age_z + collectGA_z + group
#
# 主阈值：
#   FDR < 0.10 & |log2FC| >= 0.5
#
# 输出文件：
#   01_DESeq2_ALL_all_results_for_GSEA.csv
#   02_DESeq2_DEG_FDR0p10_log2FC0p5.csv
# ============================================================


# ============================================================
# 0. 加载 R 包
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
})


# ============================================================
# 1. 参数区
# ============================================================

expr_file <- file.path("results", "cohort_selection", "19_FINAL_expression_matrix_matched.tsv")

meta_file <- file.path("results", "cohort_selection", "20_FINAL_sample_table_matched_to_matrix.csv")

outdir_base <- file.path("results", "cfRNA_DESeq2")

dir.create(outdir_base, recursive = TRUE, showWarnings = FALSE)


# ------------------------------
# 列名设置
# ------------------------------

gene_col      <- "sym_id"
sample_col    <- "sample"
group_col     <- "group2"
period_col    <- "period_new"

age_col       <- "age"
collectGA_col <- "collectGA"


# ------------------------------
# 分析时期
# ------------------------------

periods_to_run <- c(
  "A_pre_le16",
  "B_mid_gt16_le32",
  "C_late_gt32"
)


# ------------------------------
# 分组方向
# log2FoldChange > 0 表示 GDM 高于 Healthy
# ------------------------------

ref_group    <- "Healthy"
target_group <- "GDM"


# ------------------------------
# 低表达过滤：按你自己的标准
# ------------------------------

min_count <- 10
min_samples <- 3


# ------------------------------
# 主阈值
# ------------------------------

fdr_cutoff <- 0.10
logFC_cutoff <- 0.5


# ------------------------------
# RSEM expected counts 小数处理
# ------------------------------

count_noninteger_action <- "round"


# ------------------------------
# 每组最少样本数
# ------------------------------

min_group_n <- 3


# ============================================================
# 2. 工具函数
# ============================================================

read_table_auto <- function(file) {
  
  if (!file.exists(file)) {
    stop("文件不存在：", file)
  }
  
  if (grepl("\\.csv$", file, ignore.case = TRUE)) {
    x <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  } else if (grepl("\\.tsv$|\\.txt$", file, ignore.case = TRUE)) {
    x <- read.delim(file, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    stop("不支持的文件格式：", file)
  }
  
  return(x)
}


write_csv_safe <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, file, row.names = FALSE, fileEncoding = "UTF-8")
}


noninteger_prop_fun <- function(mat, tol = 1e-6) {
  x <- as.numeric(mat)
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  mean(abs(x - round(x)) > tol)
}


check_required_cols <- function(df, cols, df_name) {
  miss <- setdiff(cols, colnames(df))
  if (length(miss) > 0) {
    stop(df_name, " 缺少必要列：", paste(miss, collapse = ", "))
  }
}


make_direction <- function(log2fc) {
  out <- rep(NA_character_, length(log2fc))
  out[!is.na(log2fc) & log2fc > 0] <- "Up_in_GDM"
  out[!is.na(log2fc) & log2fc < 0] <- "Down_in_GDM"
  out[!is.na(log2fc) & log2fc == 0] <- "No_change"
  return(out)
}


check_design_full_rank <- function(design_formula, design_df) {
  mm <- model.matrix(design_formula, data = design_df)
  qr(mm)$rank == ncol(mm)
}


# ============================================================
# 3. 模型设置
# ============================================================

model_specs <- list(
  
  M3_basic = list(
    formula = as.formula("~ group"),
    description = "group only"
  ),
  
  M4_age_collectGA = list(
    formula = as.formula("~ age_z + collectGA_z + group"),
    description = "adjust age_z and collectGA_z"
  )
)


# ============================================================
# 4. 读取输入文件
# ============================================================

cat("\n============================================================\n")
cat("读取表达矩阵和样本表\n")
cat("============================================================\n")

expr_raw <- read_table_auto(expr_file)
meta_raw <- read_table_auto(meta_file)

cat("表达矩阵原始维度：", nrow(expr_raw), "行 ×", ncol(expr_raw), "列\n")
cat("样本表原始维度：", nrow(meta_raw), "行 ×", ncol(meta_raw), "列\n")


# ============================================================
# 5. 检查必要列
# ============================================================

required_meta_cols <- c(
  sample_col,
  group_col,
  period_col,
  age_col,
  collectGA_col
)

check_required_cols(meta_raw, required_meta_cols, "样本表")

if (!(gene_col %in% colnames(expr_raw))) {
  stop("表达矩阵中找不到基因列：", gene_col)
}


# ============================================================
# 6. 整理表达矩阵
# ============================================================

gene_vec <- as.character(expr_raw[[gene_col]])

expr_df <- expr_raw[, setdiff(colnames(expr_raw), gene_col), drop = FALSE]

expr_mat <- as.data.frame(expr_df, check.names = FALSE)

for (j in seq_len(ncol(expr_mat))) {
  expr_mat[[j]] <- suppressWarnings(as.numeric(expr_mat[[j]]))
}

expr_mat <- as.matrix(expr_mat)


# ------------------------------
# 去掉空基因名
# ------------------------------

valid_gene <- !is.na(gene_vec) & gene_vec != ""

expr_mat <- expr_mat[valid_gene, , drop = FALSE]
gene_vec <- gene_vec[valid_gene]

if (nrow(expr_mat) == 0) {
  stop("表达矩阵没有有效基因。")
}


# ------------------------------
# 检查 NA / 负值
# ------------------------------

if (any(is.na(expr_mat))) {
  stop("表达矩阵中存在 NA，请先检查输入矩阵。")
}

if (any(expr_mat < 0)) {
  stop("表达矩阵中存在负值，DESeq2 不接受负值。")
}


# ------------------------------
# 合并重复基因
# ------------------------------

if (any(duplicated(gene_vec))) {
  cat("检测到重复基因名，将按基因名求和合并。\n")
  expr_mat <- rowsum(expr_mat, group = gene_vec, reorder = FALSE)
} else {
  rownames(expr_mat) <- gene_vec
}


# ------------------------------
# RSEM expected count 小数 round
# ------------------------------

noninteger_prop <- noninteger_prop_fun(expr_mat)

cat("表达矩阵非整数比例：", noninteger_prop, "\n")

if (!is.na(noninteger_prop) && noninteger_prop > 0) {
  
  if (count_noninteger_action == "round") {
    cat("检测到非整数 count，按 RSEM expected count 处理：round 后运行 DESeq2。\n")
    expr_mat <- round(expr_mat)
  } else {
    stop("表达矩阵存在非整数。")
  }
}

storage.mode(expr_mat) <- "integer"

cat("整理后表达矩阵维度：", nrow(expr_mat), "基因 ×", ncol(expr_mat), "样本\n")


# ============================================================
# 7. 整理样本表并匹配样本
# ============================================================

meta <- meta_raw

meta[[sample_col]] <- as.character(meta[[sample_col]])
meta[[group_col]] <- as.character(meta[[group_col]])
meta[[period_col]] <- as.character(meta[[period_col]])

if (any(duplicated(meta[[sample_col]]))) {
  dup_samples <- unique(meta[[sample_col]][duplicated(meta[[sample_col]])])
  stop("样本表中存在重复 sample：", paste(head(dup_samples, 20), collapse = ", "))
}


# ------------------------------
# 只保留 Healthy / GDM
# ------------------------------

meta <- meta[meta[[group_col]] %in% c(ref_group, target_group), , drop = FALSE]

if (nrow(meta) == 0) {
  stop("样本表中没有 Healthy / GDM 样本，请检查 group_col。")
}


# ------------------------------
# 匹配共同样本
# ------------------------------

common_samples <- intersect(colnames(expr_mat), meta[[sample_col]])

cat("表达矩阵样本数：", ncol(expr_mat), "\n")
cat("样本表目标组样本数：", nrow(meta), "\n")
cat("共同样本数：", length(common_samples), "\n")

if (length(common_samples) < 10) {
  stop("表达矩阵和样本表共同样本数太少：", length(common_samples))
}

expr_mat <- expr_mat[, common_samples, drop = FALSE]
meta <- meta[match(common_samples, meta[[sample_col]]), , drop = FALSE]

rownames(meta) <- meta[[sample_col]]

if (!all(colnames(expr_mat) == meta[[sample_col]])) {
  stop("表达矩阵列名和样本表 sample 顺序没有对齐。")
}

cat("\n匹配后分组情况：\n")
print(table(meta[[group_col]], useNA = "ifany"))

cat("\n匹配后时期情况：\n")
print(table(meta[[period_col]], useNA = "ifany"))


# ============================================================
# 8. 循环分析：时期 × 模型
# ============================================================

all_results_list <- list()
deg_results_list <- list()

for (period_i in periods_to_run) {
  
  cat("\n############################################################\n")
  cat("开始分析时期：", period_i, "\n")
  cat("############################################################\n")
  
  meta_p <- meta[meta[[period_col]] == period_i, , drop = FALSE]
  
  if (nrow(meta_p) == 0) {
    warning("该时期没有样本，跳过：", period_i)
    next
  }
  
  meta_p$age_raw <- suppressWarnings(as.numeric(meta_p[[age_col]]))
  meta_p$collectGA_raw <- suppressWarnings(as.numeric(meta_p[[collectGA_col]]))
  
  meta_p$group <- factor(meta_p[[group_col]], levels = c(ref_group, target_group))
  
  samples_p <- meta_p[[sample_col]]
  count_p <- expr_mat[, samples_p, drop = FALSE]
  
  if (!all(colnames(count_p) == meta_p[[sample_col]])) {
    stop("时期 ", period_i, " 的 count 矩阵和样本表没有对齐。")
  }
  
  
  for (model_name in names(model_specs)) {
    
    cat("\n------------------------------------------------------------\n")
    cat("时期：", period_i, " | 模型：", model_name, "\n")
    cat("------------------------------------------------------------\n")
    
    model_formula <- model_specs[[model_name]]$formula
    
    # --------------------------------------------------------
    # 8.1 构建设计表
    # --------------------------------------------------------
    
    design_df <- data.frame(
      sample = meta_p[[sample_col]],
      group = meta_p$group,
      age_raw = meta_p$age_raw,
      collectGA_raw = meta_p$collectGA_raw,
      stringsAsFactors = FALSE
    )
    
    # M4 才需要 age 和 collectGA 完整
    if (model_name == "M4_age_collectGA") {
      
      design_df <- design_df[
        complete.cases(design_df[, c("group", "age_raw", "collectGA_raw")]),
        ,
        drop = FALSE
      ]
      
      if (nrow(design_df) == 0) {
        warning("M4 没有完整协变量样本，跳过：", period_i)
        next
      }
      
      age_sd <- sd(design_df$age_raw, na.rm = TRUE)
      collectGA_sd <- sd(design_df$collectGA_raw, na.rm = TRUE)
      
      if (is.na(age_sd) || age_sd == 0) {
        warning("age 在该时期无变异，跳过 M4：", period_i)
        next
      }
      
      if (is.na(collectGA_sd) || collectGA_sd == 0) {
        warning("collectGA 在该时期无变异，跳过 M4：", period_i)
        next
      }
      
      design_df$age_z <- as.numeric(scale(design_df$age_raw))
      design_df$collectGA_z <- as.numeric(scale(design_df$collectGA_raw))
    }
    
    # M3 不需要 age_z / collectGA_z，但为了列统一，补 NA
    if (!("age_z" %in% colnames(design_df))) {
      design_df$age_z <- NA_real_
    }
    
    if (!("collectGA_z" %in% colnames(design_df))) {
      design_df$collectGA_z <- NA_real_
    }
    
    design_df$group <- droplevels(factor(design_df$group, levels = c(ref_group, target_group)))
    
    rownames(design_df) <- design_df$sample
    
    
    # --------------------------------------------------------
    # 8.2 检查样本数
    # --------------------------------------------------------
    
    n_ref <- sum(design_df$group == ref_group, na.rm = TRUE)
    n_target <- sum(design_df$group == target_group, na.rm = TRUE)
    
    cat("模型使用样本数：", nrow(design_df), "\n")
    cat(ref_group, "=", n_ref, "；", target_group, "=", n_target, "\n")
    
    if (n_ref < min_group_n || n_target < min_group_n) {
      warning("某组样本数不足，跳过：", period_i, " | ", model_name)
      next
    }
    
    
    # --------------------------------------------------------
    # 8.3 匹配 count 矩阵
    # --------------------------------------------------------
    
    count_model <- count_p[, design_df$sample, drop = FALSE]
    
    if (!all(colnames(count_model) == rownames(design_df))) {
      stop("模型 ", period_i, " | ", model_name, " 的 count 和 design_df 没有对齐。")
    }
    
    
    # --------------------------------------------------------
    # 8.4 低表达过滤：按你自己的标准
    # 至少 3 个样本中 count >= 10
    # --------------------------------------------------------
    
    keep_gene <- rowSums(count_model >= min_count) >= min_samples
    
    count_filt <- count_model[keep_gene, , drop = FALSE]
    
    cat("过滤前基因数：", nrow(count_model), "\n")
    cat("过滤后基因数：", nrow(count_filt), "\n")
    
    if (nrow(count_filt) < 100) {
      warning("过滤后基因数太少，跳过：", period_i, " | ", model_name)
      next
    }
    
    
    # --------------------------------------------------------
    # 8.5 检查设计矩阵
    # --------------------------------------------------------
    
    if (!check_design_full_rank(model_formula, design_df)) {
      warning("design matrix 非满秩，跳过：", period_i, " | ", model_name)
      next
    }
    
    
    # --------------------------------------------------------
    # 8.6 运行 DESeq2
    # --------------------------------------------------------
    
    dds <- DESeqDataSetFromMatrix(
      countData = count_filt,
      colData = design_df,
      design = model_formula
    )
    
    dds$group <- relevel(dds$group, ref = ref_group)
    
    dds <- DESeq(
      dds,
      sfType = "poscounts",
      quiet = FALSE
    )
    
    
    # --------------------------------------------------------
    # 8.7 提取结果
    # 主阈值是 FDR < 0.10，所以 alpha 固定为 0.10
    # independentFiltering = TRUE
    # --------------------------------------------------------
    
    res_obj <- results(
      dds,
      contrast = c("group", target_group, ref_group),
      independentFiltering = TRUE,
      alpha = fdr_cutoff
    )
    
    res_df <- as.data.frame(res_obj)
    res_df$gene <- rownames(res_df)
    
    res_df <- res_df[, c(
      "gene",
      "baseMean",
      "log2FoldChange",
      "lfcSE",
      "stat",
      "pvalue",
      "padj"
    ), drop = FALSE]
    
    res_df$period <- period_i
    res_df$model <- model_name
    res_df$model_formula <- paste(deparse(model_formula), collapse = "")
    res_df$model_description <- model_specs[[model_name]]$description
    res_df$comparison <- paste0(target_group, "_vs_", ref_group)
    
    res_df$direction_by_log2FC <- make_direction(res_df$log2FoldChange)
    
    # GSEA 推荐使用 stat 排序
    res_df$rank_metric_for_GSEA <- res_df$stat
    
    # 主 DEG 标记
    res_df$DEG_FDR0p10_log2FC0p5 <- with(
      res_df,
      !is.na(padj) &
        padj < fdr_cutoff &
        !is.na(log2FoldChange) &
        abs(log2FoldChange) >= logFC_cutoff
    )
    
    # 排序
    res_df <- res_df[
      order(is.na(res_df$padj), res_df$padj, res_df$pvalue),
      ,
      drop = FALSE
    ]
    
    rownames(res_df) <- NULL
    
    
    # --------------------------------------------------------
    # 8.8 提取主 DEG
    # --------------------------------------------------------
    
    deg_df <- res_df[res_df$DEG_FDR0p10_log2FC0p5 %in% TRUE, , drop = FALSE]
    
    if (nrow(deg_df) > 0) {
      deg_df$threshold <- "FDR0p10_log2FC0p5"
    } else {
      deg_df$threshold <- character(0)
    }
    
    all_results_list[[paste(period_i, model_name, sep = "__")]] <- res_df
    deg_results_list[[paste(period_i, model_name, sep = "__")]] <- deg_df
    
    cat("完成：", period_i, " | ", model_name, "\n")
    cat("主 DEG 数量：", nrow(deg_df), "\n")
    cat("上调：", sum(deg_df$log2FoldChange > 0, na.rm = TRUE), "\n")
    cat("下调：", sum(deg_df$log2FoldChange < 0, na.rm = TRUE), "\n")
  }
}


# ============================================================
# 9. 合并并输出两个结果文件
# ============================================================

cat("\n============================================================\n")
cat("输出最终结果\n")
cat("============================================================\n")


# ------------------------------
# 9.1 全基因表：用于 GSEA
# ------------------------------

if (length(all_results_list) > 0) {
  all_results_final <- do.call(rbind, all_results_list)
  rownames(all_results_final) <- NULL
} else {
  stop("没有成功生成任何 DESeq2 全基因结果。")
}

outfile_all <- file.path(
  outdir_base,
  "01_DESeq2_ALL_all_results_for_GSEA.csv"
)

write_csv_safe(all_results_final, outfile_all)

cat("已输出全基因 GSEA 表：\n")
cat(outfile_all, "\n")


# ------------------------------
# 9.2 主 DEG 表
# ------------------------------

if (length(deg_results_list) > 0) {
  deg_final <- do.call(rbind, deg_results_list)
  rownames(deg_final) <- NULL
} else {
  deg_final <- all_results_final[0, , drop = FALSE]
  deg_final$threshold <- character(0)
}

outfile_deg <- file.path(
  outdir_base,
  "02_DESeq2_DEG_FDR0p10_log2FC0p5.csv"
)

write_csv_safe(deg_final, outfile_deg)

cat("已输出主 DEG 表：\n")
cat(outfile_deg, "\n")


# ============================================================
# 10. 输出结束信息
# ============================================================

cat("\n============================================================\n")
cat("DESeq2 最终简化版运行完成。\n")
cat("输出目录：", outdir_base, "\n\n")

cat("结果文件：\n")
cat("1. 01_DESeq2_ALL_all_results_for_GSEA.csv\n")
cat("   用于 GSEA，推荐使用 rank_metric_for_GSEA/stat 排序。\n\n")

cat("2. 02_DESeq2_DEG_FDR0p10_log2FC0p5.csv\n")
cat("   主差异基因：FDR < 0.10 & |log2FC| >= 0.5。\n\n")

cat("方向说明：\n")
cat("log2FoldChange > 0 表示 GDM 高于 Healthy。\n")
cat("log2FoldChange < 0 表示 GDM 低于 Healthy。\n")
cat("============================================================\n")