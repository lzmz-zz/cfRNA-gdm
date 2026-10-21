# ============================================================
# DESeq2 pseudobulk results -> GSEA
# ============================================================
# 输入：
#   01_DESeq2_ALL_results.csv
#
# 适用结果列：
#   celltype
#   gene
#   log2FoldChange
#   pvalue
#   padj
#   stat
#
# 排序方式：
#   默认优先使用 DESeq2 的 Wald statistic：stat
#   如果 stat 缺失或不可用，会自动退回：
#   sign(log2FoldChange) * -log10(pvalue)
#
# 结果方向：
#   NES > 0：该通路在 GDM_Group 上调方向富集
#   NES < 0：该通路在 healthy 上调方向富集
#
# 输出：
#   00_gene_set_summary.csv
#   01_GSEA_all_results.csv
#   02_GSEA_filtered_results.csv
#   03_GSEA_count_summary.csv
#   04_GSEA_rank_check.csv
#   05_GSEA_skipped.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(clusterProfiler)
  library(msigdbr)
  library(stringr)
})

# ============================================================
# 1. 参数区：以后主要只改这里
# ============================================================

deg_file <- file.path("results", "placenta", "pseudobulk_DESeq2", "01_DESeq2_ALL_results.csv")

outdir <- file.path("results", "placenta", "GSEA")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# 需要跑的基因集
GENESET_TYPES <- c("HALLMARK", "KEGG", "REACTOME", "GO_BP")

# 物种
SPECIES <- "Homo sapiens"

# 输入列名
CELLTYPE_COL <- "celltype"
GENE_COL     <- "gene"
LOGFC_COL    <- "log2FoldChange"
PVALUE_COL   <- "pvalue"
PADJ_COL     <- "padj"
STAT_COL     <- "stat"

# 排序模式
# 推荐：
#   stat：优先使用 DESeq2 Wald statistic
# 可选：
#   signed_log10p：sign(log2FoldChange) * -log10(pvalue)
#   logFC：只按 log2FoldChange 排序
RANK_MODE <- "stat"

# GSEA 参数
MIN_RANK_GENES <- 100
MIN_GS_SIZE    <- 10
MAX_GS_SIZE    <- 500
GSEA_P_CUTOFF  <- 1

# 结果筛选阈值
FILTER_PADJ <- 0.25
FILTER_ABS_NES <- 1
FILTER_LEADING_EDGE_MIN <- 10

set.seed(123)

# ============================================================
# 2. 工具函数
# ============================================================

write_csv_safe <- function(x, file) {
  write.csv(x, file = file, row.names = FALSE, fileEncoding = "UTF-8")
}

get_msig_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) {
    return(NA_character_)
  }
  return(hit[1])
}

prepare_msigdb_safe <- function(species, geneset_types) {
  
  cat("========== Preparing MSigDB gene sets ==========\n")
  
  msig_raw <- msigdbr(species = species)
  
  collection_col <- get_msig_col(msig_raw, c("gs_collection", "gs_cat"))
  subcollection_col <- get_msig_col(msig_raw, c("gs_subcollection", "gs_subcat"))
  name_col <- get_msig_col(msig_raw, c("gs_name"))
  gene_col <- get_msig_col(msig_raw, c("gene_symbol"))
  
  if (is.na(collection_col) || is.na(name_col) || is.na(gene_col)) {
    stop("Cannot identify required columns in msigdbr output.")
  }
  
  if (is.na(subcollection_col)) {
    msig_raw$.__subcollection__ <- NA_character_
    subcollection_col <- ".__subcollection__"
  }
  
  msig_raw <- msig_raw %>%
    dplyr::mutate(
      .collection = .data[[collection_col]],
      .subcollection = .data[[subcollection_col]],
      .gs_name = .data[[name_col]],
      .gene_symbol = .data[[gene_col]]
    )
  
  gene_sets_list <- list()
  
  for (gs_type in geneset_types) {
    
    cat("Preparing gene set:", gs_type, "\n")
    
    if (gs_type == "HALLMARK") {
      tmp <- msig_raw %>%
        dplyr::filter(.collection == "H")
      
    } else if (gs_type == "KEGG") {
      tmp <- msig_raw %>%
        dplyr::filter(
          .collection == "C2",
          grepl("KEGG", .subcollection, ignore.case = TRUE) |
            grepl("^KEGG_", .gs_name)
        )
      
    } else if (gs_type == "REACTOME") {
      tmp <- msig_raw %>%
        dplyr::filter(
          .collection == "C2",
          grepl("REACTOME", .subcollection, ignore.case = TRUE) |
            grepl("^REACTOME_", .gs_name)
        )
      
    } else if (gs_type == "GO_BP") {
      tmp <- msig_raw %>%
        dplyr::filter(
          .collection == "C5",
          grepl("GO:BP|BP", .subcollection, ignore.case = TRUE) |
            grepl("^GOBP_", .gs_name)
        )
      
    } else {
      stop("Unsupported gene set type: ", gs_type)
    }
    
    tmp2 <- tmp %>%
      dplyr::transmute(
        geneset_type = gs_type,
        term = .gs_name,
        gene = .gene_symbol
      ) %>%
      dplyr::filter(
        !is.na(term),
        !is.na(gene),
        term != "",
        gene != ""
      ) %>%
      dplyr::distinct()
    
    if (nrow(tmp2) == 0) {
      warning("No gene sets found for: ", gs_type)
    }
    
    gene_sets_list[[gs_type]] <- tmp2
  }
  
  gene_sets <- dplyr::bind_rows(gene_sets_list)
  
  gene_set_summary <- gene_sets %>%
    dplyr::group_by(geneset_type) %>%
    dplyr::summarise(
      n_terms = dplyr::n_distinct(term),
      n_unique_genes = dplyr::n_distinct(gene),
      .groups = "drop"
    )
  
  print(gene_set_summary)
  
  return(list(
    gene_sets = gene_sets,
    summary = gene_set_summary
  ))
}

make_rank_DESeq2 <- function(df,
                             gene_col,
                             logfc_col,
                             pvalue_col,
                             stat_col,
                             rank_mode = "stat") {
  
  need_cols <- c(gene_col, logfc_col, pvalue_col)
  miss_cols <- setdiff(need_cols, colnames(df))
  if (length(miss_cols) > 0) {
    stop("Missing columns for DESeq2 ranking: ", paste(miss_cols, collapse = ", "))
  }
  
  has_stat <- stat_col %in% colnames(df)
  
  tmp <- df %>%
    dplyr::transmute(
      gene = as.character(.data[[gene_col]]),
      logFC = suppressWarnings(as.numeric(.data[[logfc_col]])),
      pvalue = suppressWarnings(as.numeric(.data[[pvalue_col]])),
      stat = if (has_stat) suppressWarnings(as.numeric(.data[[stat_col]])) else NA_real_
    ) %>%
    dplyr::filter(
      !is.na(gene),
      gene != "",
      !is.na(logFC),
      is.finite(logFC)
    )
  
  if (nrow(tmp) == 0) {
    return(numeric(0))
  }
  
  # 处理 pvalue = 0
  positive_p <- tmp$pvalue[tmp$pvalue > 0 & is.finite(tmp$pvalue)]
  
  if (length(positive_p) == 0) {
    p_floor <- .Machine$double.xmin
  } else {
    p_floor <- min(positive_p, na.rm = TRUE) * 0.1
    p_floor <- max(p_floor, .Machine$double.xmin)
  }
  
  tmp$pvalue_safe <- ifelse(
    is.na(tmp$pvalue) | !is.finite(tmp$pvalue),
    NA_real_,
    ifelse(tmp$pvalue <= 0, p_floor, tmp$pvalue)
  )
  
  if (rank_mode == "stat") {
    
    # 如果 stat 可用，就用 stat
    if (has_stat && sum(!is.na(tmp$stat) & is.finite(tmp$stat)) >= MIN_RANK_GENES) {
      tmp$rank_score <- tmp$stat
    } else {
      message("stat is unavailable or insufficient; fallback to signed_log10p.")
      tmp$rank_score <- sign(tmp$logFC) * (-log10(tmp$pvalue_safe))
    }
    
  } else if (rank_mode == "signed_log10p") {
    tmp$rank_score <- sign(tmp$logFC) * (-log10(tmp$pvalue_safe))
    
  } else if (rank_mode == "logFC") {
    tmp$rank_score <- tmp$logFC
    
  } else {
    stop("Unsupported RANK_MODE for DESeq2: ", rank_mode)
  }
  
  tmp <- tmp %>%
    dplyr::filter(!is.na(rank_score), is.finite(rank_score))
  
  if (nrow(tmp) == 0) {
    return(numeric(0))
  }
  
  # 同一个 gene 多行时，保留 |rank_score| 最大的一行
  tmp2 <- tmp %>%
    dplyr::group_by(gene) %>%
    dplyr::slice_max(order_by = abs(rank_score), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  gene_list <- tmp2$rank_score
  names(gene_list) <- tmp2$gene
  
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  return(gene_list)
}

leading_edge_size <- function(x) {
  if (is.na(x) || x == "") {
    return(0)
  }
  length(unique(unlist(strsplit(x, "/"))))
}

run_one_gsea <- function(gene_list, term2gene, min_gs_size, max_gs_size) {
  
  gsea_obj <- tryCatch(
    {
      clusterProfiler::GSEA(
        geneList = gene_list,
        TERM2GENE = term2gene,
        pvalueCutoff = GSEA_P_CUTOFF,
        pAdjustMethod = "BH",
        minGSSize = min_gs_size,
        maxGSSize = max_gs_size,
        verbose = FALSE,
        seed = TRUE
      )
    },
    error = function(e) {
      message("GSEA failed: ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(gsea_obj)) {
    return(NULL)
  }
  
  res <- as.data.frame(gsea_obj@result)
  
  if (nrow(res) == 0) {
    return(NULL)
  }
  
  return(res)
}

# ============================================================
# 3. 读取 DESeq2 差异结果
# ============================================================

cat("========== Reading DESeq2 DEG results ==========\n")

deg <- read.csv(deg_file, check.names = FALSE, stringsAsFactors = FALSE)

need_cols <- c(CELLTYPE_COL, GENE_COL, LOGFC_COL, PVALUE_COL, PADJ_COL)
miss_cols <- setdiff(need_cols, colnames(deg))

if (length(miss_cols) > 0) {
  stop("Input DESeq2 result is missing columns: ", paste(miss_cols, collapse = ", "))
}

cat("Input rows:", nrow(deg), "\n")
cat("Celltypes:\n")
print(sort(unique(deg[[CELLTYPE_COL]])))

# ============================================================
# 4. 准备基因集
# ============================================================

msig <- prepare_msigdb_safe(SPECIES, GENESET_TYPES)
gene_sets <- msig$gene_sets
gene_set_summary <- msig$summary

write_csv_safe(
  gene_set_summary,
  file.path(outdir, "00_gene_set_summary.csv")
)

# ============================================================
# 5. 按 celltype × geneset_type 跑 GSEA
# ============================================================

cat("\n========== Running GSEA ==========\n")

all_gsea <- list()
rank_checks <- list()
skipped <- list()

celltypes <- sort(unique(deg[[CELLTYPE_COL]]))

for (ct in celltypes) {
  
  cat("\n--------------------------------------------------\n")
  cat("Celltype:", ct, "\n")
  
  deg_ct <- deg %>%
    dplyr::filter(.data[[CELLTYPE_COL]] == ct)
  
  gene_list <- make_rank_DESeq2(
    df = deg_ct,
    gene_col = GENE_COL,
    logfc_col = LOGFC_COL,
    pvalue_col = PVALUE_COL,
    stat_col = STAT_COL,
    rank_mode = RANK_MODE
  )
  
  rank_check <- data.frame(
    celltype = ct,
    n_input_rows = nrow(deg_ct),
    n_rank_genes = length(gene_list),
    n_positive_rank = sum(gene_list > 0),
    n_negative_rank = sum(gene_list < 0),
    n_zero_rank = sum(gene_list == 0),
    max_rank = ifelse(length(gene_list) > 0, max(gene_list), NA),
    min_rank = ifelse(length(gene_list) > 0, min(gene_list), NA),
    rank_mode = RANK_MODE,
    stringsAsFactors = FALSE
  )
  
  rank_checks[[ct]] <- rank_check
  
  if (length(gene_list) < MIN_RANK_GENES) {
    skipped[[paste0(ct, "__ALL")]] <- data.frame(
      celltype = ct,
      geneset_type = "ALL",
      reason = "too_few_ranked_genes",
      detail = paste0("n_rank_genes=", length(gene_list)),
      stringsAsFactors = FALSE
    )
    next
  }
  
  for (gs_type in GENESET_TYPES) {
    
    cat("Running gene set:", gs_type, "\n")
    
    term2gene <- gene_sets %>%
      dplyr::filter(geneset_type == gs_type) %>%
      dplyr::select(term, gene) %>%
      dplyr::distinct()
    
    overlap_genes <- intersect(names(gene_list), unique(term2gene$gene))
    
    if (length(overlap_genes) < MIN_RANK_GENES) {
      skipped[[paste0(ct, "__", gs_type)]] <- data.frame(
        celltype = ct,
        geneset_type = gs_type,
        reason = "too_few_overlap_genes",
        detail = paste0("n_overlap_genes=", length(overlap_genes)),
        stringsAsFactors = FALSE
      )
      next
    }
    
    gsea_res <- run_one_gsea(
      gene_list = gene_list,
      term2gene = term2gene,
      min_gs_size = MIN_GS_SIZE,
      max_gs_size = MAX_GS_SIZE
    )
    
    if (is.null(gsea_res) || nrow(gsea_res) == 0) {
      skipped[[paste0(ct, "__", gs_type)]] <- data.frame(
        celltype = ct,
        geneset_type = gs_type,
        reason = "no_gsea_result",
        detail = NA_character_,
        stringsAsFactors = FALSE
      )
      next
    }
    
    gsea_res$celltype <- ct
    gsea_res$geneset_type <- gs_type
    gsea_res$rank_mode <- RANK_MODE
    gsea_res$leading_edge_n <- vapply(gsea_res$core_enrichment, leading_edge_size, numeric(1))
    
    gsea_res <- gsea_res %>%
      dplyr::select(
        celltype,
        geneset_type,
        ID,
        Description,
        setSize,
        enrichmentScore,
        NES,
        pvalue,
        p.adjust,
        qvalue,
        rank,
        leading_edge_n,
        core_enrichment,
        rank_mode,
        dplyr::everything()
      )
    
    all_gsea[[paste0(ct, "__", gs_type)]] <- gsea_res
  }
}

# ============================================================
# 6. 汇总和输出
# ============================================================

cat("\n========== Writing outputs ==========\n")

gsea_all <- if (length(all_gsea) > 0) {
  dplyr::bind_rows(all_gsea)
} else {
  data.frame()
}

rank_check_all <- if (length(rank_checks) > 0) {
  dplyr::bind_rows(rank_checks)
} else {
  data.frame()
}

skipped_all <- if (length(skipped) > 0) {
  dplyr::bind_rows(skipped)
} else {
  data.frame(
    celltype = character(),
    geneset_type = character(),
    reason = character(),
    detail = character()
  )
}

if (nrow(gsea_all) > 0) {
  
  gsea_filtered <- gsea_all %>%
    dplyr::filter(
      !is.na(p.adjust),
      p.adjust < FILTER_PADJ,
      !is.na(NES),
      abs(NES) > FILTER_ABS_NES,
      leading_edge_n >= FILTER_LEADING_EDGE_MIN
    )
  
  gsea_summary <- gsea_all %>%
    dplyr::group_by(celltype, geneset_type) %>%
    dplyr::summarise(
      total_terms = dplyr::n(),
      filtered_terms = sum(
        !is.na(p.adjust) &
          p.adjust < FILTER_PADJ &
          !is.na(NES) &
          abs(NES) > FILTER_ABS_NES &
          leading_edge_n >= FILTER_LEADING_EDGE_MIN
      ),
      NES_positive_filtered = sum(
        !is.na(p.adjust) &
          p.adjust < FILTER_PADJ &
          !is.na(NES) &
          NES > FILTER_ABS_NES &
          leading_edge_n >= FILTER_LEADING_EDGE_MIN
      ),
      NES_negative_filtered = sum(
        !is.na(p.adjust) &
          p.adjust < FILTER_PADJ &
          !is.na(NES) &
          NES < -FILTER_ABS_NES &
          leading_edge_n >= FILTER_LEADING_EDGE_MIN
      ),
      .groups = "drop"
    )
  
} else {
  gsea_filtered <- data.frame()
  gsea_summary <- data.frame()
}

write_csv_safe(gsea_all, file.path(outdir, "01_GSEA_all_results.csv"))
write_csv_safe(gsea_filtered, file.path(outdir, "02_GSEA_filtered_results.csv"))
write_csv_safe(gsea_summary, file.path(outdir, "03_GSEA_count_summary.csv"))
write_csv_safe(rank_check_all, file.path(outdir, "04_GSEA_rank_check.csv"))
write_csv_safe(skipped_all, file.path(outdir, "05_GSEA_skipped.csv"))

cat("\n========== Finished DESeq2 GSEA ==========\n")
cat("Output directory:", outdir, "\n")