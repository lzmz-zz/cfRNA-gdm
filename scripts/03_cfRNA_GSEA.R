# cfRNA preranked GSEA (Figure 2 primary analysis)
# Run from repository root AFTER scripts/02_cfRNA_DESeq2.R.
# Source: GSEA (1)(7).ipynb, executed cells 0 / 1 / 2.
# M3_basic (~ group) produces the primary Figure 2 GSEA results.
# M4_age_collectGA is only used for covariate-adjusted sensitivity analysis
# (Supplementary Table 7 / Methods); it is optional and disabled by default.
# Only model dispatch is modified from the previously curated version.
# Gene ranking, gene sets, thresholds, and GSEA algorithm are unchanged.

# ============================================================
# Split DESeq2 all_results 总表
# 拆成：2 个模型 × 3 个时期 = 6 个 GSEA 输入文件
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ============================================================
# 1. 参数区
# ============================================================

# 你的 DESeq2 全基因总表
all_results_file <- file.path("results", "cfRNA_DESeq2", "01_DESeq2_ALL_all_results_for_GSEA.csv")

# 拆分后输出目录
split_outdir <- file.path("results", "cfRNA_DESeq2", "split_by_model_period")

dir.create(split_outdir, recursive = TRUE, showWarnings = FALSE)

# 主分析默认只跑 M3；如需复现校正敏感性分析，把 FALSE 改为 TRUE。
RUN_ADJUSTED_SENSITIVITY <- FALSE
models_to_run <- if (RUN_ADJUSTED_SENSITIVITY) {
  c("M3_basic", "M4_age_collectGA")
} else {
  "M3_basic"
}

# 只拆分本次实际运行的模型
models_to_split <- models_to_run

# 需要拆分的时期
periods_to_split <- c(
  "A_pre_le16",
  "B_mid_gt16_le32",
  "C_late_gt32"
)

# ============================================================
# 2. 读取总表
# ============================================================

all_df <- read.csv(
  all_results_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("读取总表：", all_results_file, "\n")
cat("总表维度：", nrow(all_df), "行 ×", ncol(all_df), "列\n")

# ============================================================
# 3. 检查必要列
# ============================================================

required_cols <- c(
  "gene",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj",
  "period",
  "model"
)

missing_cols <- setdiff(required_cols, colnames(all_df))

if (length(missing_cols) > 0) {
  stop("总表缺少必要列：", paste(missing_cols, collapse = ", "))
}

cat("\n总表中的 period：\n")
print(table(all_df$period, useNA = "ifany"))

cat("\n总表中的 model：\n")
print(table(all_df$model, useNA = "ifany"))

# ============================================================
# 4. 循环拆分并输出
# ============================================================

manifest_list <- list()

for (model_i in models_to_split) {
  
  for (period_i in periods_to_split) {
    
    sub_df <- all_df %>%
      dplyr::filter(
        model == model_i,
        period == period_i
      )
    
    if (nrow(sub_df) == 0) {
      warning("没有数据：", model_i, " | ", period_i)
      next
    }
    
    outdir_i <- file.path(split_outdir, model_i)
    dir.create(outdir_i, recursive = TRUE, showWarnings = FALSE)
    
    outfile_i <- file.path(
      outdir_i,
      paste0(
        "DESeq2_",
        period_i,
        "_",
        model_i,
        "_all_results_for_GSEA.csv"
      )
    )
    
    write.csv(
      sub_df,
      outfile_i,
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )
    
    cat("已输出：", outfile_i, " | 基因数：", nrow(sub_df), "\n")
    
    manifest_list[[paste(model_i, period_i, sep = "__")]] <- data.frame(
      model = model_i,
      period = period_i,
      n_genes = nrow(sub_df),
      file = outfile_i,
      stringsAsFactors = FALSE
    )
  }
}

manifest_df <- bind_rows(manifest_list)

manifest_file <- file.path(split_outdir, "00_split_manifest.csv")

write.csv(
  manifest_df,
  manifest_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("\n============================================================\n")
cat("拆分完成。\n")
cat("拆分目录：", split_outdir, "\n")
cat("文件清单：", manifest_file, "\n")
cat("============================================================\n")

run_original_gsea_model <- function(model_name) {
# ============================================================
# GSEA for DESeq2 all_results
# A/B/C 三期分别做 preranked GSEA
#
# 整体逻辑保持原代码不变：
#   1. 读取每个时期的 DESeq2 all_results
#   2. 用 stat 构建全基因排序
#   3. 使用 MSigDB gene sets
#   4. clusterProfiler::GSEA
#   5. 输出全部结果、筛选结果、数量汇总、NES 矩阵
#
# 本版修改：
#   1. 显著通路阈值：qvalue < 0.05
#   2. 基因集不变：HALLMARK / KEGG / REACTOME / GO_BP
#   3. 修复 GSEA 报错：object '.Random.seed' not found
#      方法：外部 set.seed(123)，GSEA 内部 seed = FALSE
# ============================================================


# ============================================================
# 0. 加载 R 包
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(enrichplot)
})


# ============================================================
# 1. 参数区
# ============================================================

# ------------------------------------------------------------
# 1.1 选择要跑 GSEA 的模型
#
# 你可以改成：
#   MODEL_TO_RUN <- "M3_basic"
# 或：
#   MODEL_TO_RUN <- "M4_age_collectGA"
# ------------------------------------------------------------

# The model name is supplied by run_original_gsea_model(model_name).
MODEL_TO_RUN <- model_name


# ------------------------------------------------------------
# 1.2 拆分后的 DESeq2 all_results 目录
# ------------------------------------------------------------

split_outdir <- file.path("results", "cfRNA_DESeq2", "split_by_model_period")


# ------------------------------------------------------------
# 1.3 输入文件：三个时期
# ------------------------------------------------------------

DEG_FILES <- c(
  A_pre_le16 = file.path(
    split_outdir,
    MODEL_TO_RUN,
    paste0("DESeq2_A_pre_le16_", MODEL_TO_RUN, "_all_results_for_GSEA.csv")
  ),
  
  B_mid_gt16_le32 = file.path(
    split_outdir,
    MODEL_TO_RUN,
    paste0("DESeq2_B_mid_gt16_le32_", MODEL_TO_RUN, "_all_results_for_GSEA.csv")
  ),
  
  C_late_gt32 = file.path(
    split_outdir,
    MODEL_TO_RUN,
    paste0("DESeq2_C_late_gt32_", MODEL_TO_RUN, "_all_results_for_GSEA.csv")
  )
)


# ------------------------------------------------------------
# 1.4 输出目录
# ------------------------------------------------------------

OUTDIR <- file.path(
  file.path("results", "cfRNA_GSEA"),
  MODEL_TO_RUN
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 1.5 排序指标
#
# 推荐 stat：
#   stat > 0：GDM 高于 Healthy
#   stat < 0：GDM 低于 Healthy
# ------------------------------------------------------------

RANK_METRIC <- "stat"

# 可选：
# RANK_METRIC <- "log2FoldChange"
# RANK_METRIC <- "signed_log10p"


# ------------------------------------------------------------
# 1.6 Gene sets
#
# 按你的要求：基因集不变
# ------------------------------------------------------------

GENESET_TYPES <- c(
  "HALLMARK",
  "KEGG",
  "REACTOME",
  "GO_BP"
)


# ------------------------------------------------------------
# 1.7 GSEA 参数
# ------------------------------------------------------------

MIN_GS_SIZE <- 10
MAX_GS_SIZE <- 500


# ------------------------------------------------------------
# 1.8 筛选阈值
#
# 按你的要求：
#   qvalue < 0.05
#
# 不额外限制 NES，因此 NES_CUTOFF = 0
# ------------------------------------------------------------

FILTER_Q_COL <- "qvalue"
Q_CUTOFF <- 0.05

NES_CUTOFF <- 0

USE_LEADING_EDGE_FILTER <- FALSE
LEADING_EDGE_N_CUTOFF <- 10


# ------------------------------------------------------------
# 1.9 是否自动把 ENSG 转为 SYMBOL
# ------------------------------------------------------------

AUTO_MAP_ENSEMBL_TO_SYMBOL <- TRUE


# ------------------------------------------------------------
# 1.10 是否保存基础 dotplot
# ------------------------------------------------------------

SAVE_BASIC_PLOT <- FALSE


# ------------------------------------------------------------
# 1.11 固定随机种子
# ------------------------------------------------------------

GSEA_SEED <- 123
set.seed(GSEA_SEED)


# ============================================================
# 2. 工具函数
# ============================================================

safe_read_csv <- function(file) {
  
  if (!file.exists(file)) {
    stop("文件不存在：", file)
  }
  
  read.csv(
    file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}


make_rank_score <- function(df, rank_metric = "stat") {
  
  if (rank_metric == "stat") {
    
    if (!("stat" %in% colnames(df))) {
      stop("输入表缺少 stat 列，不能使用 RANK_METRIC = 'stat'")
    }
    
    score <- df$stat
  }
  
  if (rank_metric == "log2FoldChange") {
    
    if (!("log2FoldChange" %in% colnames(df))) {
      stop("输入表缺少 log2FoldChange 列，不能使用 RANK_METRIC = 'log2FoldChange'")
    }
    
    score <- df$log2FoldChange
  }
  
  if (rank_metric == "signed_log10p") {
    
    if (!all(c("log2FoldChange", "pvalue") %in% colnames(df))) {
      stop("输入表缺少 log2FoldChange 或 pvalue 列，不能使用 signed_log10p")
    }
    
    p_use <- df$pvalue
    p_use[is.na(p_use)] <- 1
    
    if (any(p_use == 0, na.rm = TRUE)) {
      min_nonzero_p <- min(p_use[p_use > 0], na.rm = TRUE)
      p_use[p_use == 0] <- min_nonzero_p
    }
    
    score <- sign(df$log2FoldChange) * (-log10(p_use))
  }
  
  if (!exists("score")) {
    stop("RANK_METRIC 设置错误：", rank_metric)
  }
  
  return(score)
}


map_ensembl_to_symbol_if_needed <- function(gene_vec) {
  
  gene_vec <- as.character(gene_vec)
  
  is_ens <- grepl("^ENSG", gene_vec)
  
  if (!AUTO_MAP_ENSEMBL_TO_SYMBOL || sum(is_ens, na.rm = TRUE) == 0) {
    return(gene_vec)
  }
  
  ens_clean <- gene_vec
  ens_clean <- sub("\\..*$", "", ens_clean)
  
  mapped <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = ens_clean[is_ens],
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  
  out <- gene_vec
  out[is_ens] <- ifelse(
    is.na(mapped),
    gene_vec[is_ens],
    as.character(mapped)
  )
  
  return(out)
}


get_msigdb_term2gene <- function(geneset_type) {
  
  if (geneset_type == "HALLMARK") {
    
    msig <- msigdbr(
      species = "Homo sapiens",
      category = "H"
    )
  }
  
  if (geneset_type == "KEGG") {
    
    msig <- msigdbr(
      species = "Homo sapiens",
      category = "C2",
      subcategory = "CP:KEGG"
    )
  }
  
  if (geneset_type == "REACTOME") {
    
    msig <- msigdbr(
      species = "Homo sapiens",
      category = "C2",
      subcategory = "CP:REACTOME"
    )
  }
  
  if (geneset_type == "GO_BP") {
    
    msig <- msigdbr(
      species = "Homo sapiens",
      category = "C5",
      subcategory = "GO:BP"
    )
  }
  
  if (!exists("msig")) {
    stop("geneset_type 设置错误：", geneset_type)
  }
  
  term2gene <- msig %>%
    dplyr::select(gs_name, gene_symbol) %>%
    dplyr::distinct() %>%
    dplyr::filter(
      !is.na(gs_name),
      !is.na(gene_symbol),
      gs_name != "",
      gene_symbol != ""
    )
  
  if (nrow(term2gene) == 0) {
    stop("gene set 为空，请检查 msigdbr 版本或 geneset_type：", geneset_type)
  }
  
  return(term2gene)
}


parse_leading_edge_n <- function(core_enrichment) {
  
  if (is.na(core_enrichment) || core_enrichment == "") {
    return(NA_integer_)
  }
  
  length(unique(unlist(strsplit(core_enrichment, "/"))))
}


run_gsea_one <- function(gene_list, term2gene) {
  
  # 修复点：
  # 不再使用 seed = TRUE，避免 clusterProfiler/DOSE 内部调用 .Random.seed 报错。
  # 在函数外/函数内用 set.seed 固定随机性，GSEA 内部 seed = FALSE。
  
  set.seed(GSEA_SEED)
  
  gsea_res <- tryCatch({
    
    clusterProfiler::GSEA(
      geneList = gene_list,
      TERM2GENE = term2gene,
      pvalueCutoff = 1,
      minGSSize = MIN_GS_SIZE,
      maxGSSize = MAX_GS_SIZE,
      seed = FALSE,
      by = "fgsea",
      verbose = FALSE
    )
    
  }, error = function(e) {
    
    message("GSEA by='fgsea' 失败，尝试不使用 by 参数。错误信息：", e$message)
    
    set.seed(GSEA_SEED)
    
    clusterProfiler::GSEA(
      geneList = gene_list,
      TERM2GENE = term2gene,
      pvalueCutoff = 1,
      minGSSize = MIN_GS_SIZE,
      maxGSSize = MAX_GS_SIZE,
      seed = FALSE,
      verbose = FALSE
    )
  })
  
  return(gsea_res)
}


# ============================================================
# 3. 保存参数记录
# ============================================================

param_df <- data.frame(
  parameter = c(
    "MODEL_TO_RUN",
    "RANK_METRIC",
    "GENESET_TYPES",
    "MIN_GS_SIZE",
    "MAX_GS_SIZE",
    "FILTER_Q_COL",
    "Q_CUTOFF",
    "NES_CUTOFF",
    "USE_LEADING_EDGE_FILTER",
    "LEADING_EDGE_N_CUTOFF",
    "AUTO_MAP_ENSEMBL_TO_SYMBOL",
    "GSEA_SEED"
  ),
  value = c(
    MODEL_TO_RUN,
    RANK_METRIC,
    paste(GENESET_TYPES, collapse = ";"),
    MIN_GS_SIZE,
    MAX_GS_SIZE,
    FILTER_Q_COL,
    Q_CUTOFF,
    NES_CUTOFF,
    USE_LEADING_EDGE_FILTER,
    LEADING_EDGE_N_CUTOFF,
    AUTO_MAP_ENSEMBL_TO_SYMBOL,
    GSEA_SEED
  ),
  stringsAsFactors = FALSE
)

write.csv(
  param_df,
  file.path(OUTDIR, "00_GSEA_parameters_used.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ============================================================
# 4. 获取 gene sets
# ============================================================

cat("\n============================================================\n")
cat("读取 MSigDB gene sets\n")
cat("============================================================\n")

term2gene_list <- list()
geneset_summary_list <- list()

for (gs_type in GENESET_TYPES) {
  
  cat("读取 gene set：", gs_type, "\n")
  
  term2gene_i <- get_msigdb_term2gene(gs_type)
  
  term2gene_list[[gs_type]] <- term2gene_i
  
  geneset_summary_list[[gs_type]] <- data.frame(
    geneset_type = gs_type,
    n_terms = length(unique(term2gene_i$gs_name)),
    n_pairs = nrow(term2gene_i),
    n_unique_genes = length(unique(term2gene_i$gene_symbol)),
    stringsAsFactors = FALSE
  )
}

geneset_summary <- bind_rows(geneset_summary_list)

write.csv(
  geneset_summary,
  file.path(OUTDIR, "01_geneset_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

print(geneset_summary)


# ============================================================
# 5. 构建每个时期的 gene list 并运行 GSEA
# ============================================================

all_gsea_list <- list()
rank_check_list <- list()
overlap_report_list <- list()

for (period_i in names(DEG_FILES)) {
  
  cat("\n############################################################\n")
  cat("开始时期：", period_i, "\n")
  cat("输入文件：", DEG_FILES[[period_i]], "\n")
  cat("############################################################\n")
  
  df <- safe_read_csv(DEG_FILES[[period_i]])
  
  if (!("gene" %in% colnames(df))) {
    stop("输入文件缺少 gene 列：", DEG_FILES[[period_i]])
  }
  
  df$gene_raw <- as.character(df$gene)
  df$gene_for_gsea <- map_ensembl_to_symbol_if_needed(df$gene_raw)
  
  df$rank_score <- make_rank_score(df, RANK_METRIC)
  
  df_rank <- df %>%
    dplyr::filter(
      !is.na(gene_for_gsea),
      gene_for_gsea != "",
      !is.na(rank_score),
      is.finite(rank_score)
    ) %>%
    dplyr::group_by(gene_for_gsea) %>%
    dplyr::slice_max(
      order_by = abs(rank_score),
      n = 1,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup()
  
  gene_list <- df_rank$rank_score
  names(gene_list) <- df_rank$gene_for_gsea
  
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  if (length(gene_list) < 100) {
    stop("用于 GSEA 的基因数太少：", period_i, " | n = ", length(gene_list))
  }
  
  rank_check_i <- data.frame(
    period = period_i,
    model = MODEL_TO_RUN,
    input_rows = nrow(df),
    ranked_genes = length(gene_list),
    rank_metric = RANK_METRIC,
    max_rank = max(gene_list, na.rm = TRUE),
    min_rank = min(gene_list, na.rm = TRUE),
    top_gene = names(gene_list)[1],
    bottom_gene = names(gene_list)[length(gene_list)],
    stringsAsFactors = FALSE
  )
  
  rank_check_list[[period_i]] <- rank_check_i
  
  cat("用于 GSEA 的基因数：", length(gene_list), "\n")
  cat("Top gene：", names(gene_list)[1], " | score = ", gene_list[1], "\n")
  cat("Bottom gene：", names(gene_list)[length(gene_list)], " | score = ", gene_list[length(gene_list)], "\n")
  
  
  for (gs_type in GENESET_TYPES) {
    
    cat("运行 GSEA：", period_i, " | ", gs_type, "\n")
    
    term2gene_i <- term2gene_list[[gs_type]]
    
    overlap_genes <- intersect(names(gene_list), unique(term2gene_i$gene_symbol))
    
    overlap_report_list[[paste(period_i, gs_type, sep = "__")]] <- data.frame(
      period = period_i,
      model = MODEL_TO_RUN,
      geneset_type = gs_type,
      ranked_genes = length(gene_list),
      geneset_unique_genes = length(unique(term2gene_i$gene_symbol)),
      overlap_genes = length(overlap_genes),
      overlap_ratio_in_ranked = length(overlap_genes) / length(gene_list),
      stringsAsFactors = FALSE
    )
    
    if (length(overlap_genes) < 100) {
      warning("gene overlap 较少：", period_i, " | ", gs_type, " | overlap = ", length(overlap_genes))
    }
    
    gsea_res_i <- run_gsea_one(
      gene_list = gene_list,
      term2gene = term2gene_i
    )
    
    gsea_df_i <- as.data.frame(gsea_res_i)
    
    if (nrow(gsea_df_i) == 0) {
      warning("GSEA 无结果：", period_i, " | ", gs_type)
      next
    }
    
    gsea_df_i$period <- period_i
    gsea_df_i$model <- MODEL_TO_RUN
    gsea_df_i$geneset_type <- gs_type
    
    if ("core_enrichment" %in% colnames(gsea_df_i)) {
      gsea_df_i$leading_edge_n <- vapply(
        gsea_df_i$core_enrichment,
        parse_leading_edge_n,
        integer(1)
      )
    } else {
      gsea_df_i$leading_edge_n <- NA_integer_
    }
    
    all_gsea_list[[paste(period_i, gs_type, sep = "__")]] <- gsea_df_i
  }
}


# ============================================================
# 6. 输出输入检查表
# ============================================================

rank_check_df <- bind_rows(rank_check_list)

write.csv(
  rank_check_df,
  file.path(OUTDIR, "02_input_rank_check.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

overlap_report_df <- bind_rows(overlap_report_list)

write.csv(
  overlap_report_df,
  file.path(OUTDIR, "03_gene_overlap_report.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ============================================================
# 7. 合并 GSEA 结果
# ============================================================

if (length(all_gsea_list) == 0) {
  stop("没有任何 GSEA 结果。")
}

all_gsea_results <- bind_rows(all_gsea_list)

cat("\n全部 GSEA 结果数量：", nrow(all_gsea_results), "\n")


# ============================================================
# 8. 按 qvalue < 0.05 筛选
# ============================================================

if (!(FILTER_Q_COL %in% colnames(all_gsea_results))) {
  stop(
    "FILTER_Q_COL = ", FILTER_Q_COL, " 不在 GSEA 结果列中。\n",
    "可用列：", paste(colnames(all_gsea_results), collapse = ", ")
  )
}

all_gsea_results <- all_gsea_results %>%
  dplyr::mutate(
    q_for_filter = suppressWarnings(as.numeric(.data[[FILTER_Q_COL]])),
    NES = suppressWarnings(as.numeric(NES)),
    leading_edge_n = suppressWarnings(as.numeric(leading_edge_n)),
    
    pass_q = !is.na(q_for_filter) & q_for_filter < Q_CUTOFF,
    
    pass_NES = if (NES_CUTOFF > 0) {
      !is.na(NES) & abs(NES) > NES_CUTOFF
    } else {
      !is.na(NES)
    },
    
    pass_leading_edge = if (USE_LEADING_EDGE_FILTER) {
      !is.na(leading_edge_n) & leading_edge_n > LEADING_EDGE_N_CUTOFF
    } else {
      TRUE
    },
    
    pass_main_filter = pass_q & pass_NES & pass_leading_edge
  ) %>%
  dplyr::arrange(period, geneset_type, q_for_filter, dplyr::desc(abs(NES)))


# ------------------------------------------------------------
# 8.1 输出带筛选标记的全部 GSEA 结果
# ------------------------------------------------------------

write.csv(
  all_gsea_results,
  file.path(OUTDIR, "04_GSEA_all_results.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 8.2 输出 qvalue < 0.05 的筛选结果
# ------------------------------------------------------------

filtered_gsea_results <- all_gsea_results %>%
  dplyr::filter(pass_main_filter)

write.csv(
  filtered_gsea_results,
  file.path(OUTDIR, "05_GSEA_filtered_results_qvalue0p05.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("qvalue < 0.05 筛选后结果数量：", nrow(filtered_gsea_results), "\n")


# ============================================================
# 9. 数量汇总
# ============================================================

count_summary <- all_gsea_results %>%
  dplyr::group_by(period, model, geneset_type) %>%
  dplyr::summarise(
    n_all_terms = dplyr::n(),
    n_qvalue_lt_0p05 = sum(pass_main_filter, na.rm = TRUE),
    n_up_NES_positive = sum(pass_main_filter & NES > 0, na.rm = TRUE),
    n_down_NES_negative = sum(pass_main_filter & NES < 0, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  count_summary,
  file.path(OUTDIR, "06_GSEA_count_summary_qvalue0p05.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

print(count_summary)


# ============================================================
# 10. 生成 NES matrix
# ============================================================

if (nrow(filtered_gsea_results) > 0) {
  
  nes_matrix_df <- filtered_gsea_results %>%
    dplyr::select(period, model, geneset_type, ID, Description, NES, qvalue) %>%
    dplyr::mutate(
      pathway_label = paste(geneset_type, Description, sep = "__")
    ) %>%
    dplyr::select(pathway_label, period, NES) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(
      names_from = period,
      values_from = NES
    )
  
} else {
  
  nes_matrix_df <- data.frame()
}

write.csv(
  nes_matrix_df,
  file.path(OUTDIR, "07_GSEA_filtered_NES_matrix_qvalue0p05.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ============================================================
# 11. 可选：基础 dotplot
# ============================================================

if (SAVE_BASIC_PLOT && nrow(filtered_gsea_results) > 0) {
  
  pdf(
    file.path(OUTDIR, "08_GSEA_filtered_dotplot_qvalue0p05.pdf"),
    width = 12,
    height = 8
  )
  
  p <- filtered_gsea_results %>%
    dplyr::group_by(period, geneset_type) %>%
    dplyr::slice_min(order_by = q_for_filter, n = 10, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    ggplot(
      aes(
        x = NES,
        y = reorder(Description, NES),
        size = -log10(q_for_filter),
        color = NES
      )
    ) +
    geom_point() +
    facet_grid(geneset_type ~ period, scales = "free_y", space = "free_y") +
    theme_bw() +
    labs(
      x = "NES",
      y = NULL,
      size = "-log10(qvalue)",
      color = "NES"
    )
  
  print(p)
  
  dev.off()
}


# ============================================================
# 12. 结束信息
# ============================================================

cat("\n============================================================\n")
cat("GSEA 完成。\n")
cat("模型：", MODEL_TO_RUN, "\n")
cat("输出目录：", OUTDIR, "\n\n")

cat("核心输出文件：\n")
cat("00_GSEA_parameters_used.csv\n")
cat("01_geneset_summary.csv\n")
cat("02_input_rank_check.csv\n")
cat("03_gene_overlap_report.csv\n")
cat("04_GSEA_all_results.csv\n")
cat("05_GSEA_filtered_results_qvalue0p05.csv\n")
cat("06_GSEA_count_summary_qvalue0p05.csv\n")
cat("07_GSEA_filtered_NES_matrix_qvalue0p05.csv\n\n")

cat("筛选标准：qvalue < 0.05\n")
cat("NES 方向解释：\n")
cat("NES > 0：富集在 GDM 高表达端。\n")
cat("NES < 0：富集在 Healthy 高表达端，即 GDM 下调端。\n")
cat("============================================================\n")
}

# Run primary analysis; optionally also run adjusted sensitivity analysis.
for (model_name in models_to_run) {
  run_original_gsea_model(model_name)
}
