# ============================================================
# 06 | cfRNA cell type-associated signal analysis
# Source A: 细胞溯源-real (1) (2)(5).ipynb
# Source B: Untitled2(2).ipynb (participant permutation)
#
# Run at repository root: Rscript scripts/06_cell_origin_analysis.R
# Part I: 16-module x 63-cell Spearman + rank-residual adjustments
# Part II: 10,000 participant-stratified permutations
#
# Curation: path/metadata-reader changes; removed environment reset.
# No changes to statistics, model terms, permutation design or thresholds.
# ============================================================

# PART I | SPEARMAN + COVARIATE-ADJUSTED CORRELATIONS

# ============================================================
# 全部样本：16个功能模块 × 63种细胞来源 Spearman相关
#
# 不按A/B/C时期拆分
#
# 同时计算：
# 1. raw_rho：
#    全部样本直接计算Spearman相关
#
# 2. adjusted_rho：
#    不分时期，但校正Period和Group后计算相关
#
# 输入：
# 1. Step6_sample_module_scores_with_metadata.csv
# 2. 20_FINAL_sample_table_matched_to_matrix (11).csv
# 3. config/deconformer_celltypes_63.csv：固定的63种细胞名称（不含受试者信息）
#
# 输出：
# 1. 全部1008组原始相关结果
# 2. 全部1008组校正相关结果
# 3. 63细胞 × 16模块气泡热图
# 4. 每个模块最强相关细胞来源
# ============================================================

library(tidyverse)
library(data.table)


# ============================================================
# 1. 输入输出路径
# ============================================================

module_file <- file.path(
  "results", "cfRNA", "modules",
  "Step6_sample_module_scores_with_metadata.csv"
)

sample_file <- file.path(
  "results", "cohort_selection",
  "20_FINAL_sample_table_matched_to_matrix.csv"
)

celltype_file <- file.path(
  "config", "deconformer_celltypes_63.csv"
)

output_dir <- file.path(
  "results", "cfRNA", "cell_origin",
  "module_cellorigin_Spearman_all_samples"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. 读取模块分数
# ============================================================

module_dat <- fread(
  module_file,
  check.names = FALSE
)

sample_candidates <- c(
  "sample",
  "sample_id",
  "Sample",
  "SampleID"
)

module_sample_col <- intersect(
  sample_candidates,
  names(module_dat)
)[1]

if (is.na(module_sample_col)) {
  stop("模块分数文件中未找到样本列。")
}

if (module_sample_col != "sample") {
  setnames(
    module_dat,
    module_sample_col,
    "sample"
  )
}

module_dat[, sample := trimws(as.character(sample))]


# ============================================================
# 3. 识别16个模块
# ============================================================

module_cols_raw <- grep(
  "^Dynamic_Module_[0-9]{2}$",
  names(module_dat),
  value = TRUE
)

if (length(module_cols_raw) != 16) {
  stop(
    "识别到的模块数量不是16个，而是：",
    length(module_cols_raw)
  )
}

module_number <- as.integer(
  sub(
    "^Dynamic_Module_",
    "",
    module_cols_raw
  )
)

module_cols_raw <- module_cols_raw[
  order(module_number)
]

module_cols <- sprintf(
  "M%02d",
  sort(module_number)
)

setnames(
  module_dat,
  old = module_cols_raw,
  new = module_cols
)

module_dat[
  ,
  (module_cols) := lapply(
    .SD,
    function(x) suppressWarnings(as.numeric(x))
  ),
  .SDcols = module_cols
]


# ============================================================
# 4. 从公开配置表读取与原结果一致的63种细胞名称
#
# 注意：
# 此文件仅包含细胞名称，不包含rho、P值或受试者数据
# ============================================================

old_cor <- fread(
  celltype_file,
  check.names = FALSE
)

if (!"celltype" %in% names(old_cor)) {
  stop("细胞类型配置表中不存在celltype列。")
}

cell_cols <- unique(
  trimws(
    as.character(old_cor$celltype)
  )
)

cell_cols <- cell_cols[
  !is.na(cell_cols) &
    cell_cols != ""
]

message(
  "从细胞类型配置表中识别到细胞来源数量：",
  length(cell_cols)
)

if (length(cell_cols) != 63) {
  stop(
    "细胞来源数量不是63个，而是：",
    length(cell_cols),
    "。请检查细胞类型配置表。"
  )
}

write.csv(
  data.frame(celltype = cell_cols),
  file.path(
    output_dir,
    "00_celltype_list_63.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 5. 读取样本信息及63种细胞来源分数
# ============================================================

meta <- fread(
  sample_file,
  check.names = FALSE
)

meta_sample_col <- intersect(
  sample_candidates,
  names(meta)
)[1]

if (is.na(meta_sample_col)) {
  stop("样本信息文件中未找到样本列。")
}

if (meta_sample_col != "sample") {
  setnames(
    meta,
    meta_sample_col,
    "sample"
  )
}

meta[, sample := trimws(as.character(sample))]


# ============================================================
# 6. 检查必要列
# ============================================================

required_meta_cols <- c(
  "sample",
  "group2",
  "period_new",
  cell_cols
)

missing_meta_cols <- setdiff(
  required_meta_cols,
  names(meta)
)

if (length(missing_meta_cols) > 0) {
  stop(
    "样本信息文件缺少以下列：\n",
    paste(
      missing_meta_cols,
      collapse = "\n"
    )
  )
}

meta[
  ,
  (cell_cols) := lapply(
    .SD,
    function(x) suppressWarnings(as.numeric(x))
  ),
  .SDcols = cell_cols
]


# ============================================================
# 7. 合并模块分数和细胞来源分数
# ============================================================

analysis_dat <- merge(
  module_dat[
    ,
    c(
      "sample",
      module_cols
    ),
    with = FALSE
  ],
  meta[
    ,
    c(
      "sample",
      "group2",
      "period_new",
      cell_cols
    ),
    with = FALSE
  ],
  by = "sample",
  all = FALSE
)

message(
  "模块文件样本数：",
  nrow(module_dat)
)

message(
  "细胞来源文件样本数：",
  nrow(meta)
)

message(
  "最终共同样本数：",
  nrow(analysis_dat)
)

if (nrow(analysis_dat) == 0) {
  stop("模块文件和细胞来源文件没有共同样本。")
}


# ============================================================
# 8. 整理Group和Period
# ============================================================

analysis_dat <- as_tibble(analysis_dat) %>%
  mutate(
    Group = case_when(
      str_to_lower(trimws(group2)) %in%
        c("healthy", "control", "normal") ~ "Healthy",

      str_to_lower(trimws(group2)) %in%
        c("gdm", "case") ~ "GDM",

      TRUE ~ as.character(group2)
    ),

    Group = factor(
      Group,
      levels = c(
        "Healthy",
        "GDM"
      )
    ),

    Period = case_when(
      str_detect(period_new, "^A") ~ "A",
      str_detect(period_new, "^B") ~ "B",
      str_detect(period_new, "^C") ~ "C",
      TRUE ~ NA_character_
    ),

    Period = factor(
      Period,
      levels = c(
        "A",
        "B",
        "C"
      )
    )
  ) %>%
  filter(
    !is.na(Group),
    !is.na(Period)
  )

message(
  "完成Group和Period整理后的样本数：",
  nrow(analysis_dat)
)

print(
  table(
    analysis_dat$Period,
    analysis_dat$Group
  )
)


# ============================================================
# 9. 安全Spearman函数
# ============================================================

safe_spearman <- function(x, y) {

  keep <- is.finite(x) & is.finite(y)

  x <- x[keep]
  y <- y[keep]

  n_now <- length(x)

  if (
    n_now < 5 ||
      length(unique(x)) < 3 ||
      length(unique(y)) < 3
  ) {
    return(
      tibble(
        n = n_now,
        rho = NA_real_,
        pvalue = NA_real_
      )
    )
  }

  test_result <- suppressWarnings(
    cor.test(
      x,
      y,
      method = "spearman",
      exact = FALSE
    )
  )

  tibble(
    n = n_now,
    rho = unname(test_result$estimate),
    pvalue = test_result$p.value
  )
}


# ============================================================
# 10. 计算原始Spearman相关
#
# 全部样本直接计算，不分时期
# ============================================================

raw_result <- map_dfr(
  module_cols,
  function(module_now) {

    map_dfr(
      cell_cols,
      function(cell_now) {

        cor_now <- safe_spearman(
          analysis_dat[[module_now]],
          analysis_dat[[cell_now]]
        )

        tibble(
          module = module_now,
          celltype = cell_now,
          n = cor_now$n,
          raw_rho = cor_now$rho,
          raw_pvalue = cor_now$pvalue
        )
      }
    )
  }
) %>%
  mutate(
    raw_FDR_global = p.adjust(
      raw_pvalue,
      method = "BH"
    )
  ) %>%
  group_by(module) %>%
  mutate(
    raw_FDR_within_module = p.adjust(
      raw_pvalue,
      method = "BH"
    )
  ) %>%
  ungroup()


# ============================================================
# 11. 校正Period和Group
#
# 方法：
# 1. 对模块分数进行秩转换
# 2. 使用Period和Group解释模块秩
# 3. 保存模块残差
# 4. 对细胞分数采用相同处理
# 5. 对两组残差计算Spearman相关
#
# 这样仍然是不分时期分析，
# 但减少孕期和疾病分组共同驱动的伪相关
# ============================================================

get_adjusted_residual <- function(
    value,
    Period,
    Group
) {

  tmp <- tibble(
    value = suppressWarnings(
      as.numeric(value)
    ),
    Period = Period,
    Group = Group
  )

  valid <- is.finite(tmp$value) &
    !is.na(tmp$Period) &
    !is.na(tmp$Group)

  output <- rep(
    NA_real_,
    nrow(tmp)
  )

  if (sum(valid) < 10) {
    return(output)
  }

  tmp_valid <- tmp[valid, , drop = FALSE]

  if (length(unique(tmp_valid$value)) < 3) {
    return(output)
  }

  tmp_valid <- tmp_valid %>%
    mutate(
      value_rank = rank(
        value,
        ties.method = "average",
        na.last = "keep"
      )
    )

  fit <- tryCatch(
    lm(
      value_rank ~ Period + Group,
      data = tmp_valid
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(output)
  }

  output[valid] <- residuals(fit)

  output
}


# ============================================================
# 12. 为16模块和63种细胞计算校正残差
# ============================================================

adjusted_dat <- analysis_dat %>%
  select(
    sample,
    Period,
    Group,
    all_of(module_cols),
    all_of(cell_cols)
  )

all_score_cols <- c(
  module_cols,
  cell_cols
)

for (score_col in all_score_cols) {

  new_col <- paste0(
    score_col,
    "_adjusted"
  )

  adjusted_dat[[new_col]] <- get_adjusted_residual(
    value = adjusted_dat[[score_col]],
    Period = adjusted_dat$Period,
    Group = adjusted_dat$Group
  )
}


# ============================================================
# 13. 计算校正后的相关
# ============================================================

adjusted_result <- map_dfr(
  module_cols,
  function(module_now) {

    map_dfr(
      cell_cols,
      function(cell_now) {

        module_adjusted_col <- paste0(
          module_now,
          "_adjusted"
        )

        cell_adjusted_col <- paste0(
          cell_now,
          "_adjusted"
        )

        cor_now <- safe_spearman(
          adjusted_dat[[module_adjusted_col]],
          adjusted_dat[[cell_adjusted_col]]
        )

        tibble(
          module = module_now,
          celltype = cell_now,
          adjusted_n = cor_now$n,
          adjusted_rho = cor_now$rho,
          adjusted_pvalue = cor_now$pvalue
        )
      }
    )
  }
) %>%
  mutate(
    adjusted_FDR_global = p.adjust(
      adjusted_pvalue,
      method = "BH"
    )
  ) %>%
  group_by(module) %>%
  mutate(
    adjusted_FDR_within_module = p.adjust(
      adjusted_pvalue,
      method = "BH"
    )
  ) %>%
  ungroup()


# ============================================================
# 14. 合并原始和校正结果
# ============================================================

result_all <- raw_result %>%
  left_join(
    adjusted_result,
    by = c(
      "module",
      "celltype"
    )
  ) %>%
  mutate(
    same_direction = case_when(
      is.na(raw_rho) |
        is.na(adjusted_rho) ~ NA,

      sign(raw_rho) ==
        sign(adjusted_rho) ~ TRUE,

      TRUE ~ FALSE
    ),

    rho_change_after_adjustment =
      adjusted_rho - raw_rho,

    abs_rho_reduction =
      abs(raw_rho) -
      abs(adjusted_rho)
  ) %>%
  arrange(
    module,
    adjusted_FDR_global,
    desc(abs(adjusted_rho))
  )

write.csv(
  result_all,
  file.path(
    output_dir,
    "01_module_63cellorigin_Spearman_all_samples.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. 输出显著结果
# ============================================================

significant_result <- result_all %>%
  filter(
    !is.na(adjusted_FDR_global),
    adjusted_FDR_global < 0.10
  ) %>%
  arrange(
    adjusted_FDR_global,
    desc(abs(adjusted_rho))
  )

write.csv(
  significant_result,
  file.path(
    output_dir,
    "02_significant_adjusted_correlations_FDR0.10.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. 每个模块最强相关的10种细胞来源
# ============================================================

top10_each_module <- result_all %>%
  group_by(module) %>%
  arrange(
    desc(abs(adjusted_rho)),
    .by_group = TRUE
  ) %>%
  slice_head(n = 10) %>%
  ungroup()

write.csv(
  top10_each_module,
  file.path(
    output_dir,
    "03_top10_cellorigins_each_module.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. 每个模块显著相关数量汇总
# ============================================================

module_summary <- result_all %>%
  group_by(module) %>%
  summarise(
    n_celltypes = n(),

    n_raw_FDR005 = sum(
      raw_FDR_global < 0.05,
      na.rm = TRUE
    ),

    n_raw_FDR010 = sum(
      raw_FDR_global < 0.10,
      na.rm = TRUE
    ),

    n_adjusted_FDR005 = sum(
      adjusted_FDR_global < 0.05,
      na.rm = TRUE
    ),

    n_adjusted_FDR010 = sum(
      adjusted_FDR_global < 0.10,
      na.rm = TRUE
    ),

    n_adjusted_positive_FDR010 = sum(
      adjusted_FDR_global < 0.10 &
        adjusted_rho > 0,
      na.rm = TRUE
    ),

    n_adjusted_negative_FDR010 = sum(
      adjusted_FDR_global < 0.10 &
        adjusted_rho < 0,
      na.rm = TRUE
    ),

    strongest_positive_rho = max(
      adjusted_rho,
      na.rm = TRUE
    ),

    strongest_negative_rho = min(
      adjusted_rho,
      na.rm = TRUE
    ),

    .groups = "drop"
  )

write.csv(
  module_summary,
  file.path(
    output_dir,
    "04_module_correlation_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. 根据校正相关模式聚类63种细胞来源
# ============================================================

cluster_table <- result_all %>%
  select(
    celltype,
    module,
    adjusted_rho
  ) %>%
  distinct() %>%
  pivot_wider(
    names_from = module,
    values_from = adjusted_rho,
    values_fill = 0
  )

cluster_matrix <- cluster_table %>%
  column_to_rownames("celltype") %>%
  as.matrix()

cluster_matrix[
  !is.finite(cluster_matrix)
] <- 0

if (nrow(cluster_matrix) >= 2) {

  hc_cell <- hclust(
    dist(cluster_matrix),
    method = "ward.D2"
  )

  cell_order <- rownames(
    cluster_matrix
  )[hc_cell$order]

} else {

  cell_order <- rownames(
    cluster_matrix
  )
}


# ============================================================
# 19. 整理气泡热图数据
# ============================================================

if (length(cell_order) != length(unique(result_all$celltype))) {
  stop(
    "cell_order中的细胞数量与result_all中的细胞数量不一致：",
    "\ncell_order数量 = ", length(cell_order),
    "\nresult_all细胞数量 = ", length(unique(result_all$celltype))
  )
}

plot_data <- result_all %>%
  mutate(
    module = factor(
      module,
      levels = sprintf(
        "M%02d",
        1:16
      )
    ),

    Cell_source = factor(
      celltype,
      levels = rev(cell_order)
    ),

    Significance = case_when(
      !is.na(adjusted_FDR_global) &
        adjusted_FDR_global < 0.05 ~
        "Global FDR < 0.05",

      !is.na(adjusted_FDR_global) &
        adjusted_FDR_global < 0.10 ~
        "Global FDR < 0.10",

      !is.na(adjusted_pvalue) &
        adjusted_pvalue < 0.05 ~
        "P < 0.05",

      TRUE ~
        "Not significant"
    ),

    Significance = factor(
      Significance,
      levels = c(
        "Global FDR < 0.05",
        "Global FDR < 0.10",
        "P < 0.05",
        "Not significant"
      )
    ),

    FDR_label = case_when(
      !is.na(adjusted_FDR_global) &
        adjusted_FDR_global < 0.05 ~ "**",

      !is.na(adjusted_FDR_global) &
        adjusted_FDR_global < 0.10 ~ "*",

      TRUE ~ ""
    )
  )

write.csv(
  plot_data,
  file.path(
    output_dir,
    "05_heatmap_plot_data.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 20. 绘制校正相关气泡热图
# ============================================================

p_adjusted <- ggplot(
  plot_data,
  aes(
    x = module,
    y = Cell_source
  )
) +
  geom_point(
    aes(
      size = abs(adjusted_rho),
      fill = adjusted_rho,
      shape = Significance
    ),
    colour = "grey30",
    stroke = 0.35
  ) +
  geom_text(
    aes(
      label = FDR_label
    ),
    size = 2.2,
    fontface = "bold",
    colour = "black"
  ) +
  scale_y_discrete(
    labels = function(x) {
      x %>%
        str_replace_all("_", " ") %>%
        str_replace_all("-", " ") %>%
        str_squish()
    }
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Adjusted\nSpearman rho"
  ) +
  scale_size_continuous(
    range = c(0.4, 5),
    limits = c(0, 1),
    breaks = c(
      0.2,
      0.4,
      0.6,
      0.8
    ),
    name = "|rho|"
  ) +
  scale_shape_manual(
    values = c(
      "Global FDR < 0.05" = 21,
      "Global FDR < 0.10" = 22,
      "P < 0.05" = 24,
      "Not significant" = 1
    ),
    drop = FALSE,
    name = NULL
  ) +
  labs(
    x = "Functional module",
    y = "Cell-origin signal"
  ) +
  theme_bw(
    base_size = 10
  ) +
  theme(
    panel.grid.major = element_line(
      colour = "grey92",
      linewidth = 0.25
    ),

    panel.grid.minor = element_blank(),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 9
    ),

    axis.text.y = element_text(
      size = 6.5
    ),

    axis.title = element_text(
      size = 11
    ),

    legend.position = "bottom",
    legend.box = "horizontal"
  )

figure_height <- max(
  18,
  length(cell_cols) * 0.28 + 4
)

ggsave(
  file.path(
    output_dir,
    "Figure_module_63cellorigin_adjusted_bubble_heatmap.pdf"
  ),
  p_adjusted,
  width = 10,
  height = figure_height,
  limitsize = FALSE
)

ggsave(
  file.path(
    output_dir,
    "Figure_module_63cellorigin_adjusted_bubble_heatmap.png"
  ),
  p_adjusted,
  width = 10,
  height = figure_height,
  dpi = 600,
  limitsize = FALSE
)


# ============================================================
# 21. 绘制原始相关气泡热图
# ============================================================

plot_data_raw <- plot_data %>%
  mutate(
    Raw_significance = case_when(
      !is.na(raw_FDR_global) &
        raw_FDR_global < 0.05 ~
        "Global FDR < 0.05",

      !is.na(raw_FDR_global) &
        raw_FDR_global < 0.10 ~
        "Global FDR < 0.10",

      !is.na(raw_pvalue) &
        raw_pvalue < 0.05 ~
        "P < 0.05",

      TRUE ~
        "Not significant"
    ),

    Raw_significance = factor(
      Raw_significance,
      levels = c(
        "Global FDR < 0.05",
        "Global FDR < 0.10",
        "P < 0.05",
        "Not significant"
      )
    ),

    Raw_FDR_label = case_when(
      !is.na(raw_FDR_global) &
        raw_FDR_global < 0.05 ~ "**",

      !is.na(raw_FDR_global) &
        raw_FDR_global < 0.10 ~ "*",

      TRUE ~ ""
    )
  )

p_raw <- ggplot(
  plot_data_raw,
  aes(
    x = module,
    y = Cell_source
  )
) +
  geom_point(
    aes(
      size = abs(raw_rho),
      fill = raw_rho,
      shape = Raw_significance
    ),
    colour = "grey30",
    stroke = 0.35
  ) +
  geom_text(
    aes(
      label = Raw_FDR_label
    ),
    size = 2.2,
    fontface = "bold",
    colour = "black"
  ) +
  scale_y_discrete(
    labels = function(x) {
      x %>%
        str_replace_all("_", " ") %>%
        str_replace_all("-", " ") %>%
        str_squish()
    }
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Raw\nSpearman rho"
  ) +
  scale_size_continuous(
    range = c(0.4, 5),
    limits = c(0, 1),
    breaks = c(
      0.2,
      0.4,
      0.6,
      0.8
    ),
    name = "|rho|"
  ) +
  scale_shape_manual(
    values = c(
      "Global FDR < 0.05" = 21,
      "Global FDR < 0.10" = 22,
      "P < 0.05" = 24,
      "Not significant" = 1
    ),
    drop = FALSE,
    name = NULL
  ) +
  labs(
    x = "Functional module",
    y = "Cell-origin signal"
  ) +
  theme_bw(
    base_size = 10
  ) +
  theme(
    panel.grid.major = element_line(
      colour = "grey92",
      linewidth = 0.25
    ),

    panel.grid.minor = element_blank(),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 9
    ),

    axis.text.y = element_text(
      size = 6.5
    ),

    axis.title = element_text(
      size = 11
    ),

    legend.position = "bottom",
    legend.box = "horizontal"
  )

ggsave(
  file.path(
    output_dir,
    "Supplementary_raw_module_63cellorigin_bubble_heatmap.pdf"
  ),
  p_raw,
  width = 10,
  height = figure_height,
  limitsize = FALSE
)

ggsave(
  file.path(
    output_dir,
    "Supplementary_raw_module_63cellorigin_bubble_heatmap.png"
  ),
  p_raw,
  width = 10,
  height = figure_height,
  dpi = 600,
  limitsize = FALSE
)


# ============================================================
# 22. 完成提示
# ============================================================

message("分析完成。")
message("最终共同样本数：", nrow(analysis_dat))
message("模块数量：", length(module_cols))
message("细胞来源数量：", length(cell_cols))
message(
  "相关关系总数：",
  nrow(result_all)
)

message(
  "校正后Global FDR < 0.05的关系数：",
  sum(
    result_all$adjusted_FDR_global < 0.05,
    na.rm = TRUE
  )
)

message(
  "校正后Global FDR < 0.10的关系数：",
  sum(
    result_all$adjusted_FDR_global < 0.10,
    na.rm = TRUE
  )
)

p_adjusted

# ============================================================
# PART II | PARTICIPANT-STRATIFIED PERMUTATION (10,000 repeats)
# Uses the Spearman output generated immediately above.
# ============================================================

# ============================================================
# Subject-aware stratified permutation test
# for module × cell-origin correlations
#
# 数据：
# 1. results/cfRNA/modules/Step6_sample_module_scores_with_metadata.csv
#    211 samples × 16 dynamic modules
#
# 2. results/cohort_selection/20_FINAL_sample_table_matched_to_matrix.csv
#    sample metadata + Family + group2 + period_new
#    + Deconformer-p cell-origin scores
#
# 3. results/cfRNA/cell_origin/module_cellorigin_Spearman_all_samples/01_module_63cellorigin_Spearman_all_samples.csv
#    original 16 × 63 = 1008 Spearman results
#
# 统计设计：
# - 原始效应量：sample-level Spearman rho
# - 置换单位：Family / subject
# - 分层：group2 × sampling pattern
# - 同一受试者纵向cell-origin trajectory整体置换
# - 同一次置换中63种cell-origin同时移动
# - 10,000 permutations
# - 双侧empirical P
# - BH-FDR across all 1008 tests
#
# 目的：
# 检验module-cell association在控制
# group / stage structure / repeated measures后是否仍显著
# ============================================================


# ============================================================
# 0. 环境
# ============================================================

# NOTE: In the merged script, keep the objects created by Part I.
# Original standalone notebook had rm(list = ls()) and gc() here.
# Those environment-reset commands were removed; permutation logic is unchanged.

options(
  stringsAsFactors = FALSE,
  scipen = 999
)

# No readxl needed: cohort selection exports a CSV metadata table.


# ============================================================
# 1. 输入文件
# ============================================================

MODULE_FILE <- module_file
META_FILE <- sample_file
ORIGINAL_COR_FILE <- file.path(
  output_dir, "01_module_63cellorigin_Spearman_all_samples.csv"
)


# ============================================================
# 2. 输出目录
# ============================================================

OUT_DIR <- file.path(
  "results", "cfRNA", "cell_origin", "subject_permutation"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. 参数
# ============================================================

N_PERM <- 10000

SEED <- 20260914

GLOBAL_FDR_CUTOFF <- 0.05

set.seed(SEED)


# ============================================================
# 4. 检查输入文件
# ============================================================

input_files <- c(
  MODULE_FILE,
  META_FILE,
  ORIGINAL_COR_FILE
)

missing_files <-
  input_files[
    !file.exists(input_files)
  ]

if (length(missing_files) > 0) {

  stop(
    "以下文件不存在：\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}


# ============================================================
# 5. 读取数据
# ============================================================

module_df <-
  read.csv(
    MODULE_FILE,
    check.names = FALSE
  )

meta_df <-
  read.csv(
    META_FILE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

original_cor <-
  read.csv(
    ORIGINAL_COR_FILE,
    check.names = FALSE
  )


cat(
  "\n============================================================\n"
)

cat(
  "INPUT DATA\n"
)

cat(
  "============================================================\n"
)

cat(
  "Module score samples:",
  nrow(module_df),
  "\n"
)

cat(
  "Metadata samples:",
  nrow(meta_df),
  "\n"
)

cat(
  "Original correlations:",
  nrow(original_cor),
  "\n"
)


# ============================================================
# 6. 必需列检查
# ============================================================

required_module_cols <- "sample"

required_meta_cols <- c(
  "sample",
  "Family",
  "group2",
  "period_new"
)

required_cor_cols <- c(
  "module",
  "celltype",
  "n",
  "raw_rho",
  "raw_pvalue",
  "raw_FDR_global"
)


check_required <- function(
    x,
    required,
    object_name
) {

  missing <-
    setdiff(
      required,
      colnames(x)
    )

  if (length(missing) > 0) {

    stop(
      object_name,
      " 缺少列：",
      paste(
        missing,
        collapse = ", "
      )
    )
  }
}


check_required(
  module_df,
  required_module_cols,
  "module_df"
)

check_required(
  meta_df,
  required_meta_cols,
  "meta_df"
)

check_required(
  original_cor,
  required_cor_cols,
  "original_cor"
)


# ============================================================
# 7. 自动识别16个module
# ============================================================

module_cols <-
  grep(
    "^Dynamic_Module_[0-9]+$",
    colnames(module_df),
    value = TRUE
  )


if (length(module_cols) != 16) {

  stop(
    "检测到的module数量不是16，而是：",
    length(module_cols)
  )
}


module_numbers <-
  as.integer(
    sub(
      "^Dynamic_Module_",
      "",
      module_cols
    )
  )


module_order <-
  order(
    module_numbers
  )


module_cols <-
  module_cols[
    module_order
  ]


module_numbers <-
  module_numbers[
    module_order
  ]


module_names <-
  sprintf(
    "M%02d",
    module_numbers
  )


cat(
  "\nDetected modules:\n"
)

print(
  data.frame(
    column =
      module_cols,
    module =
      module_names
  )
)


# ============================================================
# 8. 自动识别原分析中的63种cell type
#
# 注意：
# 完全以原Spearman结果文件为准
# 不使用metadata中额外的3种cell type
# ============================================================

cell_types <-
  unique(
    as.character(
      original_cor$celltype
    )
  )


if (length(cell_types) != 63) {

  stop(
    "原相关性表中cell type不是63种，而是：",
    length(cell_types)
  )
}


missing_celltypes <-
  setdiff(
    cell_types,
    colnames(meta_df)
  )


if (length(missing_celltypes) > 0) {

  stop(
    "Metadata中缺少以下cell-origin列：\n",
    paste(
      missing_celltypes,
      collapse = "\n"
    )
  )
}


cat(
  "\nDetected cell types:",
  length(cell_types),
  "\n"
)


# ============================================================
# 9. 检查1008组合是否完整
# ============================================================

expected_tests <-
  length(module_names) *
  length(cell_types)


if (expected_tests != 1008) {

  stop(
    "预期检验数量不是1008。"
  )
}


if (nrow(original_cor) != 1008) {

  stop(
    "原Spearman表不是1008行。"
  )
}


pair_key <-
  paste(
    original_cor$module,
    original_cor$celltype,
    sep = "|||"
  )


if (anyDuplicated(pair_key) > 0) {

  stop(
    "原Spearman结果存在重复module-cell组合。"
  )
}


# ============================================================
# 10. sample ID完整性检查
# ============================================================

if (anyDuplicated(module_df$sample) > 0) {

  stop(
    "Module score文件sample存在重复。"
  )
}


if (anyDuplicated(meta_df$sample) > 0) {

  stop(
    "Metadata文件sample存在重复。"
  )
}


missing_in_meta <-
  setdiff(
    module_df$sample,
    meta_df$sample
  )


missing_in_module <-
  setdiff(
    meta_df$sample,
    module_df$sample
  )


if (
  length(missing_in_meta) > 0 |
  length(missing_in_module) > 0
) {

  stop(
    "两个样本文件sample ID不能完全匹配。\n",
    "Missing in metadata: ",
    length(missing_in_meta),
    "\n",
    "Missing in module score: ",
    length(missing_in_module)
  )
}


# ============================================================
# 11. 按module score文件顺序对齐metadata
# ============================================================

meta_match <-
  match(
    module_df$sample,
    meta_df$sample
  )


meta_use <-
  meta_df[
    meta_match,
    ,
    drop = FALSE
  ]


if (
  !all(
    module_df$sample ==
      meta_use$sample
  )
) {

  stop(
    "Sample alignment failed."
  )
}


# ============================================================
# 12. 检查211 samples / 132 subjects
# ============================================================

N_SAMPLE <-
  nrow(
    module_df
  )


N_SUBJECT <-
  length(
    unique(
      meta_use$Family
    )
  )


cat(
  "\n============================================================\n"
)

cat(
  "SAMPLE ALIGNMENT\n"
)

cat(
  "============================================================\n"
)

cat(
  "Samples:",
  N_SAMPLE,
  "\n"
)

cat(
  "Subjects:",
  N_SUBJECT,
  "\n"
)


if (N_SAMPLE != 211) {

  warning(
    "当前样本数不是预期的211。"
  )
}


if (N_SUBJECT != 132) {

  warning(
    "当前subject数量不是预期的132。"
  )
}


# ============================================================
# 13. group × period检查
# ============================================================

group_period_table <-
  table(
    meta_use$group2,
    meta_use$period_new
  )


cat(
  "\nGroup × period:\n"
)

print(
  group_period_table
)


write.csv(
  as.data.frame(
    group_period_table
  ),
  file.path(
    OUT_DIR,
    "01_group_period_sample_counts.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. 标准化period名称
# ============================================================

period_short <-
  ifelse(
    meta_use$period_new == "A_pre_le16",
    "A",
    ifelse(
      meta_use$period_new == "B_mid_gt16_le32",
      "B",
      ifelse(
        meta_use$period_new == "C_late_gt32",
        "C",
        NA_character_
      )
    )
  )


if (anyNA(period_short)) {

  stop(
    "period_new存在无法识别的时期。"
  )
}


Family <-
  as.character(
    meta_use$Family
  )


Group <-
  as.character(
    meta_use$group2
  )


Sample <-
  as.character(
    meta_use$sample
  )


# ============================================================
# 15. 检查每个Family在每个period最多一个样本
# ============================================================

family_period_key <-
  paste(
    Family,
    period_short,
    sep = "|||"
  )


if (anyDuplicated(family_period_key) > 0) {

  stop(
    "存在同一Family同一时期多个样本，",
    "无法按当前方案进行trajectory permutation。"
  )
}


# ============================================================
# 16. 为每个subject构建sampling pattern
#
# A
# B
# C
# AB
# AC
# BC
# ABC
# ============================================================

period_rank <-
  c(
    A = 1,
    B = 2,
    C = 3
  )


family_ids <-
  unique(
    Family
  )


family_pattern <-
  character(
    length(family_ids)
  )


family_group <-
  character(
    length(family_ids)
  )


for (
  i in seq_along(
    family_ids
  )
) {

  id_i <-
    family_ids[i]


  idx <-
    which(
      Family ==
        id_i
    )


  periods_i <-
    unique(
      period_short[idx]
    )


  periods_i <-
    periods_i[
      order(
        period_rank[
          periods_i
        ]
      )
    ]


  family_pattern[i] <-
    paste0(
      periods_i,
      collapse = ""
    )


  groups_i <-
    unique(
      Group[idx]
    )


  if (
    length(groups_i) != 1
  ) {

    stop(
      "Family ",
      id_i,
      " 存在多个group。"
    )
  }


  family_group[i] <-
    groups_i
}


family_info <-
  data.frame(

    Family =
      family_ids,

    group2 =
      family_group,

    sampling_pattern =
      family_pattern,

    stringsAsFactors =
      FALSE
  )


family_info$perm_stratum <-
  paste(
    family_info$group2,
    family_info$sampling_pattern,
    sep = "__"
  )


# ============================================================
# 17. sampling pattern summary
# ============================================================

sampling_summary <-
  as.data.frame(
    table(
      family_info$group2,
      family_info$sampling_pattern
    )
  )


colnames(
  sampling_summary
) <-
  c(
    "group2",
    "sampling_pattern",
    "N_subjects"
  )


sampling_summary <-
  sampling_summary[
    sampling_summary$N_subjects >
      0,
    ,
    drop = FALSE
  ]


cat(
  "\n============================================================\n"
)

cat(
  "SUBJECT SAMPLING PATTERNS\n"
)

cat(
  "============================================================\n"
)

print(
  sampling_summary
)


write.csv(
  sampling_summary,
  file.path(
    OUT_DIR,
    "02_sampling_pattern_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. 检查singleton strata
#
# 某个stratum只有1个subject时无法在该stratum内置换。
# 该subject保留原位置。
# ============================================================

stratum_sizes <-
  table(
    family_info$perm_stratum
  )


singleton_strata <-
  names(
    stratum_sizes[
      stratum_sizes ==
        1
    ]
  )


if (
  length(singleton_strata) > 0
) {

  cat(
    "\nNOTE: 以下置换层只有1个subject，",
    "这些subject在置换中将保持原位：\n"
  )

  print(
    singleton_strata
  )
}


# ============================================================
# 19. 给sample加入sampling pattern / stratum
# ============================================================

sample_pattern <-
  family_info$sampling_pattern[
    match(
      Family,
      family_info$Family
    )
  ]


sample_stratum <-
  family_info$perm_stratum[
    match(
      Family,
      family_info$Family
    )
  ]


alignment_check <-
  data.frame(

    sample =
      Sample,

    Family =
      Family,

    group2 =
      Group,

    period =
      period_short,

    sampling_pattern =
      sample_pattern,

    permutation_stratum =
      sample_stratum,

    stringsAsFactors =
      FALSE
  )


write.csv(
  alignment_check,
  file.path(
    OUT_DIR,
    "03_sample_alignment_and_permutation_strata.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 20. 提取module score矩阵
# ============================================================

module_matrix <-
  as.matrix(
    module_df[
      ,
      module_cols,
      drop = FALSE
    ]
  )


storage.mode(
  module_matrix
) <-
  "double"


colnames(
  module_matrix
) <-
  module_names


# ============================================================
# 21. 提取63种cell-origin矩阵
# ============================================================

cell_matrix <-
  as.matrix(
    meta_use[
      ,
      cell_types,
      drop = FALSE
    ]
  )


storage.mode(
  cell_matrix
) <-
  "double"


colnames(
  cell_matrix
) <-
  cell_types


# ============================================================
# 22. 检查缺失 / Inf
# ============================================================

if (
  any(
    !is.finite(
      module_matrix
    )
  )
) {

  stop(
    "Module matrix存在NA/Inf。"
  )
}


if (
  any(
    !is.finite(
      cell_matrix
    )
  )
) {

  stop(
    "Cell-origin matrix存在NA/Inf。"
  )
}


# ============================================================
# 23. Rank transform
#
# Spearman rho = Pearson correlation of ranks
#
# 因为当前所有1008组合均n=211，
# 可使用矩阵运算一次计算16×63。
# ============================================================

rank_matrix <- function(
    x
) {

  apply(
    x,
    2,
    rank,
    ties.method =
      "average"
  )
}


module_rank <-
  rank_matrix(
    module_matrix
  )


cell_rank <-
  rank_matrix(
    cell_matrix
  )


# ============================================================
# 24. Z-standardize ranked variables
# ============================================================

module_z <-
  scale(
    module_rank,
    center = TRUE,
    scale = TRUE
  )


cell_z <-
  scale(
    cell_rank,
    center = TRUE,
    scale = TRUE
  )


# ============================================================
# 25. 重新计算全部1008个原始Spearman rho
# ============================================================

rho_observed <-
  crossprod(
    module_z,
    cell_z
  ) /
  (
    N_SAMPLE -
      1
  )


rownames(
  rho_observed
) <-
  module_names


colnames(
  rho_observed
) <-
  cell_types


# ============================================================
# 26. 构建原始结果rho矩阵
# ============================================================

rho_original_file <-
  matrix(

    NA_real_,

    nrow =
      length(
        module_names
      ),

    ncol =
      length(
        cell_types
      ),

    dimnames =
      list(
        module_names,
        cell_types
      )
  )


for (
  i in seq_len(
    nrow(
      original_cor
    )
  )
) {

  m <-
    as.character(
      original_cor$module[i]
    )

  ctype <-
    as.character(
      original_cor$celltype[i]
    )


  rho_original_file[
    m,
    ctype
  ] <-
    original_cor$raw_rho[i]
}


# ============================================================
# 27. 强制验证：
# 新代码必须复现原Fig4的1008个rho
# ============================================================

rho_difference <-
  rho_observed -
  rho_original_file


MAX_RHO_DIFF <-
  max(
    abs(
      rho_difference
    ),
    na.rm = TRUE
  )


cat(
  "\n============================================================\n"
)

cat(
  "ORIGINAL SPEARMAN REPRODUCTION CHECK\n"
)

cat(
  "============================================================\n"
)

cat(
  "Maximum absolute rho difference:",
  format(
    MAX_RHO_DIFF,
    scientific = TRUE
  ),
  "\n"
)


if (
  MAX_RHO_DIFF >
    1e-10
) {

  stop(
    "重新计算的Spearman rho与原始结果不一致。\n",
    "停止置换分析，请检查输入数据。"
  )
}


cat(
  "PASS: all 1008 Spearman rho values reproduced.\n"
)


# ============================================================
# 28. 输出Spearman复现检查
# ============================================================

reproduction_check <-
  expand.grid(

    module =
      module_names,

    celltype =
      cell_types,

    KEEP.OUT.ATTRS =
      FALSE,

    stringsAsFactors =
      FALSE
  )


reproduction_check$rho_recalculated <-
  as.vector(
    rho_observed
  )


reproduction_check$rho_original <-
  as.vector(
    rho_original_file
  )


reproduction_check$absolute_difference <-
  abs(
    reproduction_check$rho_recalculated -
      reproduction_check$rho_original
  )


write.csv(
  reproduction_check,
  file.path(
    OUT_DIR,
    "04_original_Spearman_reproduction_check.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 29. 构建每个stratum中的subject列表
# ============================================================

families_by_stratum <-
  split(
    family_info$Family,
    family_info$perm_stratum
  )


# ============================================================
# 30. 创建 Family × period -> row index
# ============================================================

sample_row_key <-
  paste(
    Family,
    period_short,
    sep = "|||"
  )


row_lookup <-
  setNames(
    seq_len(
      N_SAMPLE
    ),
    sample_row_key
  )


# ============================================================
# 31. 置换函数
#
# 对每个stratum：
# target Family -> randomly selected source Family
#
# 因为sampling pattern完全相同，
# source Family拥有与target相同的period组合。
#
# 这样可将整条纵向cell-origin trajectory
# 作为一个整体移动。
# ============================================================

generate_source_row_index <-
  function() {

    # 默认自己映射自己
    source_family <-
      setNames(
        family_info$Family,
        family_info$Family
      )


    for (
      stratum_i in names(
        families_by_stratum
      )
    ) {

      fams <-
        families_by_stratum[[stratum_i]]


      if (
        length(fams) >
          1
      ) {

        source_family[
          fams
        ] <-
          sample(
            fams,
            size =
              length(
                fams
              ),
            replace =
              FALSE
          )
      }
    }


    # 对每一个target sample：
    # 找到source Family同一period对应的row
    source_keys <-
      paste(
        unname(
          source_family[
            Family
          ]
        ),
        period_short,
        sep =
          "|||"
      )


    source_idx <-
      unname(
        row_lookup[
          source_keys
        ]
      )


    if (
      anyNA(
        source_idx
      )
    ) {

      stop(
        "Permutation row mapping产生NA。"
      )
    }


    source_idx
  }


# ============================================================
# 32. 正式置换
#
# 重要：
# 每次置换对63种cell-origin使用相同subject mapping。
#
# 这样保留：
# - cell-origin之间的协方差结构
# - group
# - period
# - sampling pattern
# - longitudinal trajectory结构
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "STARTING SUBJECT-AWARE PERMUTATION\n"
)

cat(
  "============================================================\n"
)

cat(
  "N permutations:",
  N_PERM,
  "\n"
)

cat(
  "Tests:",
  expected_tests,
  "\n\n"
)


extreme_count <-
  matrix(

    0L,

    nrow =
      length(
        module_names
      ),

    ncol =
      length(
        cell_types
      ),

    dimnames =
      list(
        module_names,
        cell_types
      )
  )


set.seed(
  SEED
)


start_time <-
  Sys.time()


for (
  b in seq_len(
    N_PERM
  )
) {

  # ----------------------------------------------------------
  # subject-aware row permutation
  # ----------------------------------------------------------

  source_idx <-
    generate_source_row_index()


  # ----------------------------------------------------------
  # 63种cell-origin整体按照source_idx重排
  # ----------------------------------------------------------

  permuted_cell_z <-
    cell_z[
      source_idx,
      ,
      drop = FALSE
    ]


  # ----------------------------------------------------------
  # 一次计算16 × 63 correlations
  # ----------------------------------------------------------

  rho_perm <-
    crossprod(
      module_z,
      permuted_cell_z
    ) /
    (
      N_SAMPLE -
        1
    )


  # ----------------------------------------------------------
  # 双侧检验：
  # |rho_perm| >= |rho_observed|
  # ----------------------------------------------------------

  extreme_count <-
    extreme_count +
    (
      abs(
        rho_perm
      ) >=
      abs(
        rho_observed
      )
    )


  # ----------------------------------------------------------
  # progress
  # ----------------------------------------------------------

  if (
    b %% 500 ==
      0
  ) {

    elapsed <-
      difftime(
        Sys.time(),
        start_time,
        units =
          "mins"
      )


    cat(
      "[",
      b,
      "/",
      N_PERM,
      "] completed; elapsed ",
      round(
        as.numeric(
          elapsed
        ),
        2
      ),
      " min\n",
      sep = ""
    )
  }
}


# ============================================================
# 33. Empirical permutation P
#
# +1 correction避免P=0
# ============================================================

P_perm <-
  (
    extreme_count +
      1
  ) /
  (
    N_PERM +
      1
  )


rownames(
  P_perm
) <-
  module_names


colnames(
  P_perm
) <-
  cell_types


# ============================================================
# 34. Global BH correction across all 1008 tests
# ============================================================

FDR_perm_global <-
  matrix(

    p.adjust(
      as.vector(
        P_perm
      ),
      method =
        "BH"
    ),

    nrow =
      nrow(
        P_perm
      ),

    ncol =
      ncol(
        P_perm
      ),

    dimnames =
      dimnames(
        P_perm
      )
  )


# ============================================================
# 35. Within-module FDR
#
# 仅作为补充输出。
# 文章主要使用global FDR。
# ============================================================

FDR_perm_within_module <-
  t(
    apply(
      P_perm,
      1,
      function(x) {

        p.adjust(
          x,
          method =
            "BH"
        )
      }
    )
  )


rownames(
  FDR_perm_within_module
) <-
  module_names


colnames(
  FDR_perm_within_module
) <-
  cell_types


# ============================================================
# 36. 整理最终1008行结果
# ============================================================

perm_results <-
  expand.grid(

    module =
      module_names,

    celltype =
      cell_types,

    KEEP.OUT.ATTRS =
      FALSE,

    stringsAsFactors =
      FALSE
  )


perm_results$rho <-
  as.vector(
    rho_observed
  )


perm_results$permutation_p <-
  as.vector(
    P_perm
  )


perm_results$permutation_FDR_global <-
  as.vector(
    FDR_perm_global
  )


perm_results$permutation_FDR_within_module <-
  as.vector(
    FDR_perm_within_module
  )


perm_results$significant_global_FDR005 <-
  perm_results$permutation_FDR_global <
  GLOBAL_FDR_CUTOFF


perm_results$direction <-
  ifelse(
    perm_results$rho >
      0,
    "Positive",
    ifelse(
      perm_results$rho <
        0,
      "Negative",
      "Zero"
    )
  )


# ============================================================
# 37. 合并原始Spearman P/FDR
# ============================================================

orig_keep <-
  original_cor[
    ,
    c(
      "module",
      "celltype",
      "n",
      "raw_rho",
      "raw_pvalue",
      "raw_FDR_global",
      "raw_FDR_within_module"
    )
  ]


merge_key_result <-
  paste(
    perm_results$module,
    perm_results$celltype,
    sep =
      "|||"
  )


merge_key_orig <-
  paste(
    orig_keep$module,
    orig_keep$celltype,
    sep =
      "|||"
  )


orig_idx <-
  match(
    merge_key_result,
    merge_key_orig
  )


perm_results$n <-
  orig_keep$n[
    orig_idx
  ]


perm_results$raw_pvalue <-
  orig_keep$raw_pvalue[
    orig_idx
  ]


perm_results$raw_FDR_global <-
  orig_keep$raw_FDR_global[
    orig_idx
  ]


perm_results$raw_FDR_within_module <-
  orig_keep$raw_FDR_within_module[
    orig_idx
  ]


# ============================================================
# 38. 调整列顺序
# ============================================================

perm_results <-
  perm_results[
    ,
    c(
      "module",
      "celltype",
      "n",
      "rho",
      "raw_pvalue",
      "raw_FDR_global",
      "raw_FDR_within_module",
      "permutation_p",
      "permutation_FDR_global",
      "permutation_FDR_within_module",
      "significant_global_FDR005",
      "direction"
    )
  ]


# ============================================================
# 39. 排序
# ============================================================

perm_results <-
  perm_results[
    order(
      perm_results$permutation_FDR_global,
      -abs(
        perm_results$rho
      )
    ),
    ,
    drop = FALSE
  ]


# ============================================================
# 40. 输出全部1008组合
# ============================================================

write.csv(
  perm_results,
  file.path(
    OUT_DIR,
    "05_CellOrigin_16modules_63cells_subject_permutation_ALL.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 41. 输出global FDR < 0.05结果
# ============================================================

sig_results <-
  perm_results[
    perm_results$significant_global_FDR005,
    ,
    drop = FALSE
  ]


write.csv(
  sig_results,
  file.path(
    OUT_DIR,
    "06_CellOrigin_subject_permutation_globalFDR005.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 42. 原始显著 vs permutation显著比较
# ============================================================

perm_results$raw_significant_global_FDR005 <-
  perm_results$raw_FDR_global <
  0.05


perm_results$significance_change <-
  ifelse(

    perm_results$raw_significant_global_FDR005 &
      perm_results$significant_global_FDR005,

    "Significant_in_both",

    ifelse(

      perm_results$raw_significant_global_FDR005 &
        !perm_results$significant_global_FDR005,

      "Raw_only",

      ifelse(

        !perm_results$raw_significant_global_FDR005 &
          perm_results$significant_global_FDR005,

        "Permutation_only",

        "Neither"
      )
    )
  )


write.csv(
  perm_results,
  file.path(
    OUT_DIR,
    "07_raw_vs_permutation_significance_comparison.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 43. 重点结果：
# M9 / M13 肝胆来源
# M11 滋养层来源
# ============================================================

key_flag <-

  (
    perm_results$module %in%
      c(
        "M09",
        "M13"
      ) &
      perm_results$celltype %in%
      c(
        "hepatocyte",
        "intrahepatic_cholangiocyte"
      )
  ) |

  (
    perm_results$module ==
      "M11" &
      perm_results$celltype %in%
      c(
        "early_EVT",
        "late_EVT",
        "early_VCT",
        "late_VCT",
        "early_SCT",
        "late_SCT"
      )
  )


key_results <-
  perm_results[
    key_flag,
    ,
    drop = FALSE
  ]


key_results <-
  key_results[
    order(
      key_results$module,
      key_results$permutation_FDR_global
    ),
    ,
    drop = FALSE
  ]


write.csv(
  key_results,
  file.path(
    OUT_DIR,
    "08_key_M09_M13_M11_cellorigin_results.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 44. M09 / M13所有63种cell type
# ============================================================

M09_M13_all <-
  perm_results[
    perm_results$module %in%
      c(
        "M09",
        "M13"
      ),
    ,
    drop = FALSE
  ]


write.csv(
  M09_M13_all,
  file.path(
    OUT_DIR,
    "09_M09_M13_all_63_celltypes_permutation.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 45. M11所有63种cell type
# ============================================================

M11_all <-
  perm_results[
    perm_results$module ==
      "M11",
    ,
    drop = FALSE
  ]


write.csv(
  M11_all,
  file.path(
    OUT_DIR,
    "10_M11_all_63_celltypes_permutation.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 46. summary
# ============================================================

n_sig_raw <-
  sum(
    perm_results$raw_FDR_global <
      0.05,
    na.rm = TRUE
  )


n_sig_perm <-
  sum(
    perm_results$permutation_FDR_global <
      0.05,
    na.rm = TRUE
  )


comparison_table <-
  table(
    perm_results$significance_change
  )


end_time <-
  Sys.time()


runtime_min <-
  as.numeric(
    difftime(
      end_time,
      start_time,
      units =
        "mins"
    )
  )


summary_text <-
  c(

    "============================================================",
    "SUBJECT-AWARE CELL-ORIGIN PERMUTATION ANALYSIS",
    "============================================================",

    paste0(
      "Samples: ",
      N_SAMPLE
    ),

    paste0(
      "Subjects: ",
      N_SUBJECT
    ),

    paste0(
      "Modules: ",
      length(
        module_names
      )
    ),

    paste0(
      "Cell types: ",
      length(
        cell_types
      )
    ),

    paste0(
      "Total tests: ",
      nrow(
        perm_results
      )
    ),

    paste0(
      "Permutations: ",
      N_PERM
    ),

    paste0(
      "Maximum rho reproduction difference: ",
      format(
        MAX_RHO_DIFF,
        scientific =
          TRUE
      )
    ),

    "",

    paste0(
      "Raw Spearman global FDR < 0.05: ",
      n_sig_raw
    ),

    paste0(
      "Permutation global FDR < 0.05: ",
      n_sig_perm
    ),

    "",

    "Significance comparison:",

    capture.output(
      print(
        comparison_table
      )
    ),

    "",

    paste0(
      "Runtime: ",
      round(
        runtime_min,
        2
      ),
      " minutes"
    ),

    "",

    "Permutation design:",
    "Within group2 × sampling-pattern strata,",
    "the complete longitudinal cell-origin trajectory of each subject",
    "was permuted as a unit.",
    "All 63 cell-origin signatures were moved jointly.",
    "Empirical P values were two-sided.",
    "BH correction was performed across all 1008 module-cell tests.",

    "============================================================"
  )


writeLines(
  summary_text,
  file.path(
    OUT_DIR,
    "11_FINAL_PERMUTATION_SUMMARY.txt"
  )
)


# ============================================================
# 47. 控制台打印重点结果
# ============================================================

cat(
  "\n\n============================================================\n"
)

cat(
  "PERMUTATION ANALYSIS COMPLETED\n"
)

cat(
  "============================================================\n"
)

cat(
  paste(
    summary_text,
    collapse =
      "\n"
  )
)


cat(
  "\n\n============================================================\n"
)

cat(
  "KEY RESULTS\n"
)

cat(
  "============================================================\n"
)


print(
  key_results[
    ,
    c(
      "module",
      "celltype",
      "rho",
      "raw_FDR_global",
      "permutation_p",
      "permutation_FDR_global",
      "significant_global_FDR005"
    )
  ],
  row.names = FALSE
)


cat(
  "\n============================================================\n"
)

cat(
  "Output directory:\n",
  OUT_DIR,
  "\n"
)

cat(
  "============================================================\n"
)


# ============================================================
# 48. Session info
# ============================================================

sink(
  file.path(
    OUT_DIR,
    "12_sessionInfo.txt"
  )
)

print(
  sessionInfo()
)

sink()
