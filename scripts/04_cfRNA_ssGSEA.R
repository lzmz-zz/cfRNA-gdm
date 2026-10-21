# ============================================================
# GSEA-significant GO-BP full-gene-set ssGSEA
#
# 目的：
# 1. 输入 A/B/C 三期 GSEA 显著 GO-BP 通路表
# 2. 用 GSEA 结果决定哪些 GO-BP 通路进入 ssGSEA
# 3. 每条通路使用 MSigDB 完整 GO-BP gene set
# 4. 不再使用 core_enrichment genes
# 5. 表达矩阵使用 log2(DESeq2 normalized counts + 1)
# 6. 输出 sample × GO-BP ssGSEA 分数
# 7. A/B/C 每期分别比较 GDM vs Healthy
# 8. 输出 median difference、Cohen's d、Wilcoxon p、FDR、direction
# 9. 输出 GSEA NES 与 ssGSEA direction 是否一致
# 10. 输出 pathway × A/B/C effect matrix，用于后续动态模块聚类
#
# 重要：
# 这版已修复 GO_ID.x / GO_ID.y 导致 GO_ID 不存在的问题。
# ============================================================


# ============================================================
# 0. 加载 R 包
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(GSVA)
  library(msigdbr)
})

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("需要 readxl 包读取 xlsx 文件。请先安装：install.packages('readxl')")
}


# ============================================================
# 1. 参数区：主要改这里
# ============================================================

# ------------------------------------------------------------
# 1.1 表达矩阵
# 第一列为 gene symbol，后面为样本
# ------------------------------------------------------------

count_file <- file.path("results", "cohort_selection", "19_FINAL_expression_matrix_matched.tsv")


# ------------------------------------------------------------
# 1.2 样本表
# ------------------------------------------------------------

meta_file <- file.path("results", "cohort_selection", "20_FINAL_sample_table_matched_to_matrix.csv")


# ------------------------------------------------------------
# 1.3 GSEA 显著 GO-BP 通路表
#
# 你现在这个表应该是：
# GSEA-qvalue0.05-ABC-GOBP.xlsx
# 表中应该包含：
# ID / Description / NES / qvalue / period / geneset_type
# ------------------------------------------------------------

gsea_file <- file.path("data", "intermediate", "GSEA-qvalue0.05-ABC-GOBP.xlsx")

gsea_sheet <- 1


# ------------------------------------------------------------
# 1.4 输出目录
# ------------------------------------------------------------

outdir_base <- file.path("results", "cfRNA", "ssGSEA")

dir.create(outdir_base, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 1.5 列名设置
# ------------------------------------------------------------

gene_col <- "sym_id"

sample_col <- "sample"
group_col  <- "group2"
period_col <- "period_new"


# ------------------------------------------------------------
# 1.6 分组设置
# ------------------------------------------------------------

ref_group    <- "Healthy"
target_group <- "GDM"


# ------------------------------------------------------------
# 1.7 三个时期
# ------------------------------------------------------------

periods_to_use <- c(
  "A_pre_le16",
  "B_mid_gt16_le32",
  "C_late_gt32"
)

period_simple_map <- c(
  "A" = "A_pre_le16",
  "B" = "B_mid_gt16_le32",
  "C" = "C_late_gt32"
)


# ------------------------------------------------------------
# 1.8 GSEA 阈值
# ------------------------------------------------------------

gsea_qvalue_cutoff <- 0.05


# ------------------------------------------------------------
# 1.9 GO-BP gene set 大小过滤
# ------------------------------------------------------------

min_geneset_size_after_expr_filter <- 10
max_geneset_size_after_expr_filter <- 500


# ------------------------------------------------------------
# 1.10 表达矩阵类型
#
# count：raw count / expected count
# TPM：如果你的表达矩阵已经是 TPM，就改成 "TPM"
# ------------------------------------------------------------

expr_input_type <- "count"


# ------------------------------------------------------------
# 1.11 是否保存标准化表达矩阵
# ------------------------------------------------------------

save_expr_for_ssgsea <- TRUE


# ============================================================
# 2. 工具函数
# ============================================================

read_table_auto <- function(file, sheet = 1) {
  
  if (!file.exists(file)) {
    stop("文件不存在：", file)
  }
  
  if (grepl("\\.csv$", file, ignore.case = TRUE)) {
    
    x <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
    
  } else if (grepl("\\.tsv$|\\.txt$", file, ignore.case = TRUE)) {
    
    x <- read.delim(file, check.names = FALSE, stringsAsFactors = FALSE)
    
  } else if (grepl("\\.xlsx$|\\.xls$", file, ignore.case = TRUE)) {
    
    x <- readxl::read_excel(file, sheet = sheet)
    x <- as.data.frame(x, check.names = FALSE)
    
  } else {
    
    stop("不支持的文件格式：", file)
  }
  
  return(x)
}


write_csv_safe <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, file, row.names = FALSE, fileEncoding = "UTF-8")
}


sanitize_name <- function(x) {
  x <- as.character(x)
  x <- toupper(x)
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}


clean_pathway_key <- function(x) {
  
  x <- sanitize_name(x)
  
  x <- gsub("^GO_BP_FULL_", "", x)
  x <- gsub("^GO_BP_", "", x)
  x <- gsub("^GOBP_", "", x)
  x <- gsub("^MSIGDB_", "", x)
  
  return(x)
}


standardize_period_simple <- function(x) {
  
  x0 <- as.character(x)
  
  dplyr::case_when(
    grepl("A_pre|A_|^A$|le16|early|pre", x0, ignore.case = TRUE) ~ "A",
    grepl("B_mid|B_|^B$|gt16|le32|mid", x0, ignore.case = TRUE) ~ "B",
    grepl("C_late|C_|^C$|gt32|late", x0, ignore.case = TRUE) ~ "C",
    TRUE ~ NA_character_
  )
}


cohen_d_two_group <- function(x, g, target_group, ref_group) {
  
  x_target <- x[g == target_group]
  x_ref <- x[g == ref_group]
  
  x_target <- x_target[!is.na(x_target)]
  x_ref <- x_ref[!is.na(x_ref)]
  
  n1 <- length(x_target)
  n2 <- length(x_ref)
  
  if (n1 < 2 || n2 < 2) {
    return(NA_real_)
  }
  
  s1 <- sd(x_target)
  s2 <- sd(x_ref)
  
  pooled_sd <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  
  if (is.na(pooled_sd) || pooled_sd == 0) {
    return(NA_real_)
  }
  
  d <- (mean(x_target) - mean(x_ref)) / pooled_sd
  
  return(d)
}


safe_wilcox <- function(x, g, ref_group, target_group) {
  
  ok <- !is.na(x) & !is.na(g)
  x <- x[ok]
  g <- as.character(g[ok])
  
  n_ref <- sum(g == ref_group)
  n_target <- sum(g == target_group)
  
  if (n_ref < 3 || n_target < 3) {
    return(NA_real_)
  }
  
  x_ref <- x[g == ref_group]
  x_target <- x[g == target_group]
  
  p <- tryCatch({
    wilcox.test(x_target, x_ref, exact = FALSE)$p.value
  }, error = function(e) {
    NA_real_
  })
  
  return(p)
}


run_ssgsea_safe <- function(expr_mat, geneset_list) {
  
  if (length(geneset_list) == 0) {
    stop("geneset_list 为空，不能运行 ssGSEA。")
  }
  
  expr_mat <- as.matrix(expr_mat)
  storage.mode(expr_mat) <- "numeric"
  
  if (is.null(rownames(expr_mat))) {
    stop("表达矩阵没有 rownames。ssGSEA 要求 rownames 是基因名。")
  }
  
  if (is.null(colnames(expr_mat))) {
    stop("表达矩阵没有 colnames。ssGSEA 要求 colnames 是样本名。")
  }
  
  geneset_list <- lapply(geneset_list, unique)
  
  set.seed(123)
  
  ssgsea_mat <- tryCatch({
    
    if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
      
      param <- GSVA::ssgseaParam(
        exprData = expr_mat,
        geneSets = geneset_list,
        minSize = min_geneset_size_after_expr_filter,
        maxSize = max_geneset_size_after_expr_filter,
        alpha = 0.25,
        normalize = TRUE
      )
      
      GSVA::gsva(
        param,
        verbose = FALSE
      )
      
    } else {
      
      GSVA::gsva(
        expr = expr_mat,
        gset.idx.list = geneset_list,
        method = "ssgsea",
        kcdf = "Gaussian",
        min.sz = min_geneset_size_after_expr_filter,
        max.sz = max_geneset_size_after_expr_filter,
        ssgsea.norm = TRUE,
        verbose = FALSE
      )
    }
    
  }, error = function(e) {
    
    stop(
      "GSVA/ssGSEA 运行失败：", e$message, "\n\n",
      "请检查：\n",
      "1. 表达矩阵是否是 matrix，且行为基因、列为样本；\n",
      "2. rownames(expr_mat) 是否是 gene symbol；\n",
      "3. geneset_list 是否是 list；\n",
      "4. gene set 中的基因是否能匹配表达矩阵 rownames；\n",
      "5. GSVA 包版本是否与当前 R 环境兼容。"
    )
  })
  
  ssgsea_mat <- as.matrix(ssgsea_mat)
  
  return(ssgsea_mat)
}


# ============================================================
# 3. 读取 GSEA 显著 GO-BP 通路表，并整理成长表
# ============================================================

cat("\n============================================================\n")
cat("Step 1: 读取 GSEA 显著 GO-BP 通路表\n")
cat("============================================================\n")

gsea_raw <- read_table_auto(gsea_file, sheet = gsea_sheet)

cat("GSEA 输入表维度：", nrow(gsea_raw), "行 ×", ncol(gsea_raw), "列\n")
cat("GSEA 输入表列名：\n")
print(colnames(gsea_raw))


# ------------------------------------------------------------
# 3.1 如果表里有 gene set 类型列，筛 GO-BP
# ------------------------------------------------------------

type_col_candidates <- c(
  "geneset_type", "gene_set", "gene_set_type",
  "category", "gs_cat", "collection", "db", "source"
)

type_col <- intersect(type_col_candidates, colnames(gsea_raw))

if (length(type_col) > 0) {
  
  type_col <- type_col[1]
  
  before_n <- nrow(gsea_raw)
  
  gsea_raw <- gsea_raw %>%
    dplyr::filter(
      grepl("GO", .data[[type_col]], ignore.case = TRUE) |
        grepl("GOBP", .data[[type_col]], ignore.case = TRUE) |
        grepl("GO_BP", .data[[type_col]], ignore.case = TRUE) |
        grepl("GO:BP", .data[[type_col]], ignore.case = TRUE)
    )
  
  cat("根据", type_col, "筛选 GO-BP：", before_n, "→", nrow(gsea_raw), "\n")
}


# ------------------------------------------------------------
# 3.2 整理 GSEA 表为 long format
# ------------------------------------------------------------

make_gsea_long <- function(gsea_raw) {
  
  cols <- colnames(gsea_raw)
  
  pathway_col_candidates <- c(
    "ID", "GO_BP_name", "gs_name", "pathway", "term_name",
    "Description", "description", "NAME", "geneset_name"
  )
  
  pathway_col <- intersect(pathway_col_candidates, cols)
  
  if (length(pathway_col) == 0) {
    stop("GSEA 表中找不到通路名列。请检查是否有 ID / Description / gs_name / GO_BP_name。")
  }
  
  pathway_col <- pathway_col[1]
  
  period_col_candidates <- c("period", "Period", "stage", "Stage", "comparison", "contrast", "time")
  nes_col_candidates <- c("NES", "nes")
  q_col_candidates <- c("qvalue", "qvalues", "q_value", "FDR", "padj", "p.adjust")
  
  period_col <- intersect(period_col_candidates, cols)
  nes_col <- intersect(nes_col_candidates, cols)
  q_col <- intersect(q_col_candidates, cols)
  
  if (length(period_col) > 0 && length(nes_col) > 0 && length(q_col) > 0) {
    
    cat("识别为 long format GSEA 表。\n")
    
    period_col <- period_col[1]
    nes_col <- nes_col[1]
    q_col <- q_col[1]
    
    out <- gsea_raw %>%
      dplyr::transmute(
        pathway_raw = as.character(.data[[pathway_col]]),
        period_simple = standardize_period_simple(.data[[period_col]]),
        NES = suppressWarnings(as.numeric(.data[[nes_col]])),
        qvalue = suppressWarnings(as.numeric(.data[[q_col]]))
      ) %>%
      dplyr::filter(!is.na(period_simple))
    
    return(out)
  }
  
  cat("未识别为 long format，尝试识别为 wide dynamic format。\n")
  
  find_metric_col <- function(period_simple, metric_regex) {
    
    period_regex <- switch(
      period_simple,
      "A" = "A|A_pre|pre|le16|early",
      "B" = "B|B_mid|mid|gt16|le32",
      "C" = "C|C_late|late|gt32"
    )
    
    candidates <- cols[
      grepl(metric_regex, cols, ignore.case = TRUE) &
        grepl(period_regex, cols, ignore.case = TRUE)
    ]
    
    if (length(candidates) == 0) {
      return(NA_character_)
    }
    
    candidates[1]
  }
  
  out_list <- list()
  
  for (pp in c("A", "B", "C")) {
    
    nes_pp <- find_metric_col(pp, "NES")
    q_pp <- find_metric_col(pp, "qvalue|qvalues|q_value|FDR|padj|p.adjust")
    
    cat("Period", pp, ": NES列 =", nes_pp, "; qvalue列 =", q_pp, "\n")
    
    if (is.na(nes_pp) || is.na(q_pp)) {
      next
    }
    
    out_list[[pp]] <- gsea_raw %>%
      dplyr::transmute(
        pathway_raw = as.character(.data[[pathway_col]]),
        period_simple = pp,
        NES = suppressWarnings(as.numeric(.data[[nes_pp]])),
        qvalue = suppressWarnings(as.numeric(.data[[q_pp]]))
      )
  }
  
  if (length(out_list) == 0) {
    stop("无法从 GSEA 表中识别 A/B/C 的 NES 和 qvalue 列。")
  }
  
  out <- dplyr::bind_rows(out_list)
  
  return(out)
}


gsea_long <- make_gsea_long(gsea_raw) %>%
  dplyr::mutate(
    period_ssgsea = period_simple_map[period_simple],
    gsea_direction = dplyr::case_when(
      NES > 0 ~ "GDM_UP",
      NES < 0 ~ "GDM_DOWN",
      TRUE ~ "NO_DIRECTION"
    ),
    pathway_key = clean_pathway_key(pathway_raw)
  ) %>%
  dplyr::filter(
    !is.na(pathway_raw),
    pathway_raw != "",
    !is.na(period_ssgsea)
  )

gsea_sig <- gsea_long %>%
  dplyr::filter(
    !is.na(qvalue),
    qvalue < gsea_qvalue_cutoff
  )

cat("GSEA long 表记录数：", nrow(gsea_long), "\n")
cat("GSEA qvalue <", gsea_qvalue_cutoff, "显著记录数：", nrow(gsea_sig), "\n")
cat("GSEA 显著唯一通路数：", length(unique(gsea_sig$pathway_key)), "\n")

write_csv_safe(
  gsea_long,
  file.path(outdir_base, "00_input_GSEA_GO_BP_long_all_records.csv")
)

write_csv_safe(
  gsea_sig,
  file.path(outdir_base, "01_input_GSEA_GO_BP_long_qvalue_lt_0p05.csv")
)


# ============================================================
# 4. 读取 MSigDB 完整 GO-BP gene sets
# ============================================================

cat("\n============================================================\n")
cat("Step 2: 读取 MSigDB 完整 GO-BP gene sets\n")
cat("============================================================\n")

get_msigdb_gobp <- function() {
  
  msigdbr_args <- names(formals(msigdbr::msigdbr))
  
  if ("collection" %in% msigdbr_args) {
    
    gs <- msigdbr::msigdbr(
      species = "Homo sapiens",
      collection = "C5",
      subcollection = "GO:BP"
    )
    
  } else {
    
    gs <- msigdbr::msigdbr(
      species = "Homo sapiens",
      category = "C5",
      subcategory = "GO:BP"
    )
  }
  
  return(gs)
}

gobp_raw <- get_msigdb_gobp()

cat("MSigDB GO-BP 原始记录数：", nrow(gobp_raw), "\n")
cat("MSigDB GO-BP 原始通路数：", length(unique(gobp_raw$gs_name)), "\n")

required_msig_cols <- c("gs_name", "gene_symbol")
missing_msig_cols <- setdiff(required_msig_cols, colnames(gobp_raw))

if (length(missing_msig_cols) > 0) {
  stop("msigdbr 结果缺少必要列：", paste(missing_msig_cols, collapse = ", "))
}

desc_col <- intersect(c("gs_description", "description"), colnames(gobp_raw))

if (length(desc_col) == 0) {
  gobp_raw$gs_description_use <- gobp_raw$gs_name
} else {
  gobp_raw$gs_description_use <- gobp_raw[[desc_col[1]]]
}

goid_col <- intersect(c("gs_exact_source", "db_gene_set_id"), colnames(gobp_raw))

if (length(goid_col) == 0) {
  gobp_raw$GO_ID_use <- NA_character_
} else {
  gobp_raw$GO_ID_use <- as.character(gobp_raw[[goid_col[1]]])
}

gobp_def_all <- gobp_raw %>%
  dplyr::filter(!is.na(gs_name), !is.na(gene_symbol), gene_symbol != "") %>%
  dplyr::group_by(gs_name) %>%
  dplyr::summarise(
    geneset_type = "GO_BP_full",
    GO_ID = dplyr::first(GO_ID_use),
    Description = dplyr::first(gs_description_use),
    all_genes = paste(sort(unique(gene_symbol)), collapse = "/"),
    genes_before_expr_filter = length(unique(gene_symbol)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    geneset_name = paste0("GO_BP_FULL__", sanitize_name(gs_name)),
    pathway_key_gs_name = clean_pathway_key(gs_name),
    pathway_key_desc = clean_pathway_key(Description),
    pathway_key_goid = clean_pathway_key(GO_ID)
  )

write_csv_safe(
  gobp_def_all,
  file.path(outdir_base, "02_MSigDB_full_GO_BP_geneset_definition_all.csv")
)


# ============================================================
# 5. GSEA 显著通路与 MSigDB 完整 GO-BP 匹配
# ============================================================

cat("\n============================================================\n")
cat("Step 3: GSEA 显著 GO-BP 与 MSigDB 完整 GO-BP 匹配\n")
cat("============================================================\n")

gsea_keys <- unique(gsea_sig$pathway_key)

gobp_def <- gobp_def_all %>%
  dplyr::filter(
    pathway_key_gs_name %in% gsea_keys |
      pathway_key_desc %in% gsea_keys |
      pathway_key_goid %in% gsea_keys
  )

cat("GSEA 显著唯一通路数：", length(gsea_keys), "\n")
cat("匹配到 MSigDB 完整 GO-BP gene set 数：", nrow(gobp_def), "\n")

if (nrow(gobp_def) == 0) {
  
  cat("\n没有匹配到通路。请检查 GSEA 通路名和 MSigDB gs_name 是否一致。\n")
  cat("GSEA 通路名前 20 个：\n")
  print(head(unique(gsea_sig$pathway_raw), 20))
  cat("MSigDB GO-BP gs_name 前 20 个：\n")
  print(head(gobp_def_all$gs_name, 20))
  
  stop("GSEA 显著通路无法匹配 MSigDB 完整 GO-BP。")
}

write_csv_safe(
  gobp_def,
  file.path(outdir_base, "03_GSEA_sig_GO_BP_matched_full_MSigDB_genesets_before_expr_filter.csv")
)


# ------------------------------------------------------------
# 5.1 将 GSEA NES / qvalue 映射到 MSigDB geneset_name
# ------------------------------------------------------------

gsea_to_gobp_map <- gobp_def %>%
  dplyr::select(
    geneset_name,
    GO_BP_name = gs_name,
    GO_ID,
    Description,
    pathway_key_gs_name,
    pathway_key_desc,
    pathway_key_goid
  ) %>%
  tidyr::pivot_longer(
    cols = c(pathway_key_gs_name, pathway_key_desc, pathway_key_goid),
    names_to = "key_type",
    values_to = "pathway_key"
  ) %>%
  dplyr::filter(!is.na(pathway_key), pathway_key != "") %>%
  dplyr::distinct(geneset_name, GO_BP_name, GO_ID, Description, pathway_key)

gsea_sig_mapped <- gsea_sig %>%
  dplyr::left_join(
    gsea_to_gobp_map,
    by = "pathway_key"
  ) %>%
  dplyr::filter(!is.na(geneset_name)) %>%
  dplyr::group_by(geneset_name, period_ssgsea) %>%
  dplyr::arrange(qvalue, dplyr::desc(abs(NES)), .by_group = TRUE) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    geneset_name,
    GO_BP_name,
    GO_ID,
    Description,
    period_ssgsea,
    GSEA_NES = NES,
    GSEA_qvalue = qvalue,
    GSEA_direction = gsea_direction,
    GSEA_pathway_raw = pathway_raw
  )

write_csv_safe(
  gsea_sig_mapped,
  file.path(outdir_base, "04_GSEA_sig_GO_BP_mapped_to_MSigDB_full_genesets_by_period.csv")
)


# ============================================================
# 6. 读取表达矩阵和样本表
# ============================================================

cat("\n============================================================\n")
cat("Step 4: 读取表达矩阵和样本表\n")
cat("============================================================\n")

count_raw <- read_table_auto(count_file)
meta_raw <- read_table_auto(meta_file)

cat("表达矩阵原始维度：", nrow(count_raw), "行 ×", ncol(count_raw), "列\n")
cat("样本表原始维度：", nrow(meta_raw), "行 ×", ncol(meta_raw), "列\n")

required_meta_cols <- c(sample_col, group_col, period_col)
missing_meta_cols <- setdiff(required_meta_cols, colnames(meta_raw))

if (length(missing_meta_cols) > 0) {
  stop("样本表缺少必要列：", paste(missing_meta_cols, collapse = ", "))
}

meta <- meta_raw %>%
  dplyr::mutate(
    sample_ssgsea = as.character(.data[[sample_col]]),
    group_ssgsea = as.character(.data[[group_col]]),
    period_ssgsea = as.character(.data[[period_col]])
  ) %>%
  dplyr::filter(
    group_ssgsea %in% c(ref_group, target_group),
    period_ssgsea %in% periods_to_use
  )

if (any(duplicated(meta$sample_ssgsea))) {
  dup_samples <- unique(meta$sample_ssgsea[duplicated(meta$sample_ssgsea)])
  stop("样本表中存在重复 sample：", paste(head(dup_samples, 20), collapse = ", "))
}

if (!(gene_col %in% colnames(count_raw))) {
  warning("表达矩阵中找不到 gene_col = ", gene_col, "，将默认使用第一列作为基因列：", colnames(count_raw)[1])
  gene_col <- colnames(count_raw)[1]
}

gene_vec <- as.character(count_raw[[gene_col]])

sample_cols_in_count <- intersect(meta$sample_ssgsea, colnames(count_raw))

cat("表达矩阵中匹配到样本表的样本列数：", length(sample_cols_in_count), "\n")

if (length(sample_cols_in_count) < 10) {
  stop("表达矩阵中能匹配到样本表的样本列太少：", length(sample_cols_in_count))
}

count_df <- count_raw[, sample_cols_in_count, drop = FALSE]
count_num <- as.data.frame(count_df, check.names = FALSE)

for (j in seq_len(ncol(count_num))) {
  count_num[[j]] <- suppressWarnings(as.numeric(count_num[[j]]))
}

count_mat <- as.matrix(count_num)

if (any(is.na(count_mat))) {
  na_n <- sum(is.na(count_mat))
  stop("表达矩阵中存在 NA 或非数值转换失败，共 ", na_n, " 个。请检查是否有非表达列混入。")
}

if (any(count_mat < 0, na.rm = TRUE)) {
  stop("表达矩阵中存在负值，不符合 count/TPM 数据要求。")
}

valid_gene <- !is.na(gene_vec) & gene_vec != ""

count_mat <- count_mat[valid_gene, , drop = FALSE]
gene_vec <- gene_vec[valid_gene]

if (any(duplicated(gene_vec))) {
  
  cat("检测到重复基因名：", sum(duplicated(gene_vec)), "个重复行。将按 gene symbol 求和合并。\n")
  count_mat <- rowsum(count_mat, group = gene_vec, reorder = FALSE)
  
} else {
  
  rownames(count_mat) <- gene_vec
}

common_samples <- intersect(colnames(count_mat), meta$sample_ssgsea)

cat("表达矩阵样本数：", ncol(count_mat), "\n")
cat("样本表筛选后样本数：", nrow(meta), "\n")
cat("共同样本数：", length(common_samples), "\n")

if (length(common_samples) < 10) {
  stop("表达矩阵和样本表共同样本数太少：", length(common_samples))
}

count_mat <- count_mat[, common_samples, drop = FALSE]

meta <- meta %>%
  dplyr::filter(sample_ssgsea %in% common_samples) %>%
  dplyr::arrange(match(sample_ssgsea, common_samples))

if (!identical(meta$sample_ssgsea, colnames(count_mat))) {
  stop("样本表 sample 顺序和表达矩阵列名没有完全对齐。")
}

meta$group_ssgsea <- factor(meta$group_ssgsea, levels = c(ref_group, target_group))

cat("\n最终样本分布：\n")
print(table(meta$period_ssgsea, meta$group_ssgsea, useNA = "ifany"))

zero_sample <- colSums(count_mat, na.rm = TRUE) == 0

if (any(zero_sample)) {
  stop("存在总表达量为 0 的样本：", paste(colnames(count_mat)[zero_sample], collapse = ", "))
}

keep_gene_count <- rowSums(count_mat, na.rm = TRUE) > 0

cat("去除全 0 基因数：", sum(!keep_gene_count), "\n")

count_mat <- count_mat[keep_gene_count, , drop = FALSE]


# ============================================================
# 7. 表达矩阵标准化
# ============================================================

cat("\n============================================================\n")
cat("Step 5: 表达矩阵标准化\n")
cat("============================================================\n")

if (expr_input_type == "count") {
  
  count_mat_int <- round(count_mat)
  storage.mode(count_mat_int) <- "integer"
  
  coldata_tmp <- data.frame(
    sample = colnames(count_mat_int),
    row.names = colnames(count_mat_int)
  )
  
  dds_tmp <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_mat_int,
    colData = coldata_tmp,
    design = ~ 1
  )
  
  dds_tmp <- DESeq2::estimateSizeFactors(
    dds_tmp,
    type = "poscounts"
  )
  
  norm_counts <- DESeq2::counts(
    dds_tmp,
    normalized = TRUE
  )
  
  expr_for_ssgsea <- log2(norm_counts + 1)
  
  cat("ssGSEA 输入矩阵：log2(DESeq2 normalized counts + 1)\n")
  
} else if (expr_input_type == "TPM") {
  
  expr_for_ssgsea <- log2(count_mat + 1)
  
  cat("ssGSEA 输入矩阵：log2(TPM + 1)\n")
  
} else {
  
  stop("expr_input_type 只能是 'count' 或 'TPM'。")
}

expr_for_ssgsea <- as.matrix(expr_for_ssgsea)
storage.mode(expr_for_ssgsea) <- "numeric"

cat("ssGSEA 输入矩阵维度：", nrow(expr_for_ssgsea), "基因 ×", ncol(expr_for_ssgsea), "样本\n")

if (save_expr_for_ssgsea) {
  write_csv_safe(
    data.frame(
      gene = rownames(expr_for_ssgsea),
      expr_for_ssgsea,
      check.names = FALSE
    ),
    file.path(outdir_base, "05_expr_for_ssGSEA_log2_normalized_matrix.csv")
  )
}


# ============================================================
# 8. GO-BP gene sets 与表达矩阵匹配
# ============================================================

cat("\n============================================================\n")
cat("Step 6: GSEA 显著 GO-BP 的完整 gene set 与表达矩阵匹配\n")
cat("============================================================\n")

geneset_list_raw <- list()
match_report_list <- list()

for (ii in seq_len(nrow(gobp_def))) {
  
  gs_name <- gobp_def$geneset_name[ii]
  genes_raw <- unique(unlist(strsplit(gobp_def$all_genes[ii], "/")))
  genes_raw <- genes_raw[!is.na(genes_raw) & genes_raw != ""]
  
  genes_match <- intersect(genes_raw, rownames(expr_for_ssgsea))
  
  geneset_list_raw[[gs_name]] <- genes_match
  
  match_report_list[[gs_name]] <- data.frame(
    geneset_name = gs_name,
    n_genes_before_expr_filter = length(genes_raw),
    n_genes_after_expr_filter = length(genes_match),
    match_ratio = ifelse(length(genes_raw) > 0, length(genes_match) / length(genes_raw), NA_real_),
    matched_genes = paste(genes_match, collapse = "/"),
    stringsAsFactors = FALSE
  )
}

match_report <- dplyr::bind_rows(match_report_list)

write_csv_safe(
  match_report,
  file.path(outdir_base, "06_GSEA_sig_GO_BP_full_geneset_match_to_expr_matrix.csv")
)

keep_gs <- names(geneset_list_raw)[
  lengths(geneset_list_raw) >= min_geneset_size_after_expr_filter &
    lengths(geneset_list_raw) <= max_geneset_size_after_expr_filter
]

geneset_list <- geneset_list_raw[keep_gs]

cat("GSEA 显著 GO-BP 完整 gene set 数：", length(geneset_list_raw), "\n")
cat("表达矩阵过滤后保留 gene set 数：", length(geneset_list), "\n")

if (length(geneset_list) == 0) {
  stop("没有 GO-BP gene set 通过表达矩阵匹配过滤。")
}


# ------------------------------------------------------------
# 8.1 生成 geneset_def_keep
# 这里是修正版：match_report 不再带 GO_ID，避免 GO_ID.x / GO_ID.y
# ------------------------------------------------------------

geneset_def_keep <- gobp_def %>%
  dplyr::filter(geneset_name %in% keep_gs) %>%
  dplyr::left_join(
    match_report %>%
      dplyr::select(
        geneset_name,
        n_genes_before_expr_filter,
        n_genes_after_expr_filter,
        match_ratio,
        matched_genes
      ),
    by = "geneset_name"
  )


# ------------------------------------------------------------
# 8.2 保险修复：如果仍然出现 GO_ID.x / GO_ID.y，则统一回 GO_ID
# ------------------------------------------------------------

if (!("GO_ID" %in% colnames(geneset_def_keep))) {
  
  if ("GO_ID.x" %in% colnames(geneset_def_keep)) {
    geneset_def_keep$GO_ID <- geneset_def_keep$GO_ID.x
  } else if ("GO_ID.y" %in% colnames(geneset_def_keep)) {
    geneset_def_keep$GO_ID <- geneset_def_keep$GO_ID.y
  } else {
    geneset_def_keep$GO_ID <- NA_character_
  }
}

if (!("GO_BP_name" %in% colnames(geneset_def_keep))) {
  
  if ("gs_name" %in% colnames(geneset_def_keep)) {
    geneset_def_keep$GO_BP_name <- geneset_def_keep$gs_name
  } else {
    geneset_def_keep$GO_BP_name <- NA_character_
  }
}

if (!("Description" %in% colnames(geneset_def_keep))) {
  geneset_def_keep$Description <- geneset_def_keep$GO_BP_name
}

write_csv_safe(
  geneset_def_keep,
  file.path(outdir_base, "07_GSEA_sig_GO_BP_full_geneset_definition_after_expr_filter.csv")
)


# ============================================================
# 9. 运行 ssGSEA
# ============================================================

cat("\n============================================================\n")
cat("Step 7: 运行 GSEA 显著 GO-BP 的完整 gene set ssGSEA\n")
cat("============================================================\n")

ssgsea_mat <- run_ssgsea_safe(
  expr_mat = expr_for_ssgsea,
  geneset_list = geneset_list
)

cat("ssGSEA score 矩阵维度：", nrow(ssgsea_mat), "GO-BP gene sets ×", ncol(ssgsea_mat), "samples\n")

ssgsea_df <- as.data.frame(ssgsea_mat, check.names = FALSE)

ssgsea_df$geneset_name <- rownames(ssgsea_df)

ssgsea_df <- ssgsea_df %>%
  dplyr::select(geneset_name, dplyr::everything())

write_csv_safe(
  ssgsea_df,
  file.path(outdir_base, "08_GSEA_sig_GO_BP_full_ssgsea_scores_geneset_by_sample.csv")
)


# ------------------------------------------------------------
# 9.1 保存 sample × pathway 分数矩阵
# ------------------------------------------------------------

score_long <- ssgsea_df %>%
  tidyr::pivot_longer(
    cols = -geneset_name,
    names_to = "sample_ssgsea",
    values_to = "ssgsea_score"
  ) %>%
  dplyr::left_join(
    meta,
    by = "sample_ssgsea"
  ) %>%
  dplyr::left_join(
    geneset_def_keep %>%
      dplyr::select(
        geneset_name,
        GO_BP_name,
        GO_ID,
        Description,
        n_genes_after_expr_filter,
        match_ratio
      ),
    by = "geneset_name"
  )

score_wide_sample <- score_long %>%
  dplyr::select(
    sample_ssgsea,
    group_ssgsea,
    period_ssgsea,
    geneset_name,
    ssgsea_score
  ) %>%
  tidyr::pivot_wider(
    names_from = geneset_name,
    values_from = ssgsea_score
  )

write_csv_safe(
  score_wide_sample,
  file.path(outdir_base, "09_GSEA_sig_GO_BP_full_ssgsea_scores_sample_by_geneset_with_metadata.csv")
)


# ============================================================
# 10. A/B/C 三期内 GDM vs Healthy 比较
# ============================================================

cat("\n============================================================\n")
cat("Step 8: A/B/C 三期内 GDM vs Healthy 比较\n")
cat("============================================================\n")

wilcox_res <- score_long %>%
  dplyr::group_by(geneset_name, period_ssgsea) %>%
  dplyr::summarise(
    GO_BP_name = dplyr::first(GO_BP_name),
    GO_ID = dplyr::first(GO_ID),
    Description = dplyr::first(Description),
    n_genes_after_expr_filter = dplyr::first(n_genes_after_expr_filter),
    match_ratio = dplyr::first(match_ratio),
    
    n_Healthy = sum(group_ssgsea == ref_group & !is.na(ssgsea_score)),
    n_GDM = sum(group_ssgsea == target_group & !is.na(ssgsea_score)),
    
    median_Healthy = median(ssgsea_score[group_ssgsea == ref_group], na.rm = TRUE),
    median_GDM = median(ssgsea_score[group_ssgsea == target_group], na.rm = TRUE),
    mean_Healthy = mean(ssgsea_score[group_ssgsea == ref_group], na.rm = TRUE),
    mean_GDM = mean(ssgsea_score[group_ssgsea == target_group], na.rm = TRUE),
    
    median_diff_GDM_minus_Healthy = median_GDM - median_Healthy,
    mean_diff_GDM_minus_Healthy = mean_GDM - mean_Healthy,
    
    cohen_d_GDM_minus_Healthy = cohen_d_two_group(
      ssgsea_score,
      group_ssgsea,
      target_group,
      ref_group
    ),
    
    wilcoxon_pvalue = safe_wilcox(
      ssgsea_score,
      group_ssgsea,
      ref_group,
      target_group
    ),
    
    direction_ssGSEA = dplyr::case_when(
      median_diff_GDM_minus_Healthy > 0 ~ "GDM_UP",
      median_diff_GDM_minus_Healthy < 0 ~ "GDM_DOWN",
      TRUE ~ "NO_DIRECTION"
    ),
    
    .groups = "drop"
  ) %>%
  dplyr::group_by(period_ssgsea) %>%
  dplyr::mutate(
    wilcoxon_FDR_within_period = p.adjust(wilcoxon_pvalue, method = "BH")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    wilcoxon_FDR_all = p.adjust(wilcoxon_pvalue, method = "BH")
  )


# ------------------------------------------------------------
# 10.1 加入 GSEA NES / qvalue / direction，并判断方向一致性
# ------------------------------------------------------------

wilcox_res2 <- wilcox_res %>%
  dplyr::left_join(
    gsea_sig_mapped %>%
      dplyr::select(
        geneset_name,
        period_ssgsea,
        GSEA_NES,
        GSEA_qvalue,
        GSEA_direction,
        GSEA_pathway_raw
      ),
    by = c("geneset_name", "period_ssgsea")
  ) %>%
  dplyr::mutate(
    direction_consistency = dplyr::case_when(
      is.na(GSEA_direction) ~ "No_GSEA_sig_in_this_period",
      direction_ssGSEA == GSEA_direction ~ "Consistent",
      direction_ssGSEA != GSEA_direction ~ "Inconsistent",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::arrange(
    period_ssgsea,
    wilcoxon_FDR_within_period,
    dplyr::desc(abs(cohen_d_GDM_minus_Healthy))
  )

write_csv_safe(
  wilcox_res2,
  file.path(outdir_base, "10_ssGSEA_GDM_vs_Healthy_by_period_with_GSEA_direction_consistency.csv")
)

cat("每期 FDR < 0.10 的 ssGSEA GO-BP 数量：\n")
print(
  wilcox_res2 %>%
    dplyr::filter(!is.na(wilcoxon_FDR_within_period)) %>%
    dplyr::group_by(period_ssgsea, direction_ssGSEA) %>%
    dplyr::summarise(
      n_FDR_lt_0p10 = sum(wilcoxon_FDR_within_period < 0.10, na.rm = TRUE),
      .groups = "drop"
    )
)

cat("GSEA 与 ssGSEA 方向一致性统计：\n")
print(
  wilcox_res2 %>%
    dplyr::count(period_ssgsea, direction_consistency)
)


# ============================================================
# 11. 输出动态模块聚类用的三期效应矩阵
# ============================================================

cat("\n============================================================\n")
cat("Step 9: 输出动态模块聚类输入矩阵\n")
cat("============================================================\n")

make_effect_matrix <- function(df, value_col, out_file) {
  
  mat_df <- df %>%
    dplyr::select(
      geneset_name,
      GO_BP_name,
      GO_ID,
      Description,
      period_ssgsea,
      value = all_of(value_col)
    ) %>%
    tidyr::pivot_wider(
      names_from = period_ssgsea,
      values_from = value
    ) %>%
    dplyr::arrange(geneset_name)
  
  write_csv_safe(
    mat_df,
    file.path(outdir_base, out_file)
  )
  
  return(mat_df)
}


effect_cohenD <- make_effect_matrix(
  wilcox_res2,
  "cohen_d_GDM_minus_Healthy",
  "11_dynamic_effect_matrix_for_clustering_CohenD.csv"
)

effect_medianDiff <- make_effect_matrix(
  wilcox_res2,
  "median_diff_GDM_minus_Healthy",
  "12_dynamic_effect_matrix_for_clustering_medianDiff.csv"
)

wilcox_res_for_signed <- wilcox_res2 %>%
  dplyr::mutate(
    signed_log10FDR = sign(cohen_d_GDM_minus_Healthy) *
      (-log10(pmax(wilcoxon_FDR_within_period, 1e-300)))
  )

effect_signedFDR <- make_effect_matrix(
  wilcox_res_for_signed,
  "signed_log10FDR",
  "13_dynamic_effect_matrix_for_clustering_signed_log10FDR.csv"
)


# ============================================================
# 12. 输出最终宽表 summary
# ============================================================

cat("\n============================================================\n")
cat("Step 10: 输出最终 summary 表\n")
cat("============================================================\n")

summary_wide <- wilcox_res2 %>%
  dplyr::select(
    geneset_name,
    GO_BP_name,
    GO_ID,
    Description,
    n_genes_after_expr_filter,
    match_ratio,
    period_ssgsea,
    n_Healthy,
    n_GDM,
    median_Healthy,
    median_GDM,
    median_diff_GDM_minus_Healthy,
    mean_diff_GDM_minus_Healthy,
    cohen_d_GDM_minus_Healthy,
    wilcoxon_pvalue,
    wilcoxon_FDR_within_period,
    direction_ssGSEA,
    GSEA_NES,
    GSEA_qvalue,
    GSEA_direction,
    direction_consistency
  ) %>%
  tidyr::pivot_wider(
    names_from = period_ssgsea,
    values_from = c(
      n_Healthy,
      n_GDM,
      median_Healthy,
      median_GDM,
      median_diff_GDM_minus_Healthy,
      mean_diff_GDM_minus_Healthy,
      cohen_d_GDM_minus_Healthy,
      wilcoxon_pvalue,
      wilcoxon_FDR_within_period,
      direction_ssGSEA,
      GSEA_NES,
      GSEA_qvalue,
      GSEA_direction,
      direction_consistency
    )
  )

write_csv_safe(
  summary_wide,
  file.path(outdir_base, "14_ssGSEA_GSEA_sig_GO_BP_full_summary_wide_FINAL.csv")
)


# ============================================================
# 13. 输出方向一致性汇总表
# ============================================================

consistency_summary <- wilcox_res2 %>%
  dplyr::group_by(period_ssgsea, direction_consistency) %>%
  dplyr::summarise(
    n_pathways = dplyr::n(),
    .groups = "drop"
  )

write_csv_safe(
  consistency_summary,
  file.path(outdir_base, "15_GSEA_vs_ssGSEA_direction_consistency_summary.csv")
)


# ============================================================
# 14. 输出 session info
# ============================================================

sink(file.path(outdir_base, "99_sessionInfo.txt"))
print(sessionInfo())
sink()


# ============================================================
# 15. 完成提示
# ============================================================

cat("\n============================================================\n")
cat("GSEA 显著 GO-BP 的完整 gene set ssGSEA 分析完成。\n")
cat("输出目录：", outdir_base, "\n\n")

cat("重点结果文件：\n\n")

cat("08_GSEA_sig_GO_BP_full_ssgsea_scores_geneset_by_sample.csv\n")
cat("  GO-BP × sample 的 ssGSEA 分数矩阵。\n\n")

cat("09_GSEA_sig_GO_BP_full_ssgsea_scores_sample_by_geneset_with_metadata.csv\n")
cat("  sample × GO-BP 的 ssGSEA 分数矩阵，带分组和时期信息。\n\n")

cat("10_ssGSEA_GDM_vs_Healthy_by_period_with_GSEA_direction_consistency.csv\n")
cat("  A/B/C 三期 GDM vs Healthy 的 ssGSEA 差异结果，含 GSEA 方向一致性。\n\n")

cat("11_dynamic_effect_matrix_for_clustering_CohenD.csv\n")
cat("  后续动态模块聚类最推荐使用的输入矩阵：pathway × A/B/C Cohen's d。\n\n")

cat("12_dynamic_effect_matrix_for_clustering_medianDiff.csv\n")
cat("  后续动态模块聚类辅助输入矩阵：pathway × A/B/C median difference。\n\n")

cat("13_dynamic_effect_matrix_for_clustering_signed_log10FDR.csv\n")
cat("  后续动态模块聚类辅助输入矩阵：pathway × A/B/C signed -log10(FDR)。\n\n")

cat("14_ssGSEA_GSEA_sig_GO_BP_full_summary_wide_FINAL.csv\n")
cat("  最终宽表，每条 GO-BP 通路包含 A/B/C 三期的 ssGSEA 差异和 GSEA 方向信息。\n\n")

cat("15_GSEA_vs_ssGSEA_direction_consistency_summary.csv\n")
cat("  GSEA 与 ssGSEA 方向一致性统计。\n\n")

cat("文章解释口径：\n")
cat("本分析针对 A/B/C 任一时期 GSEA 显著的 GO-BP 通路，\n")
cat("基于 MSigDB 完整 GO-BP gene set 定义计算 ssGSEA 分数，\n")
cat("将 GSEA 提示的组间通路富集结果转化为单个样本层面的功能活性指标。\n")
cat("后续可基于 pathway × A/B/C effect matrix 聚类，归纳 GDM cfRNA 三期动态功能模块。\n")
cat("============================================================\n")