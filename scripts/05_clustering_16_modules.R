# ============================================================
# Step 1. 读取 A/B/C 三期 Cohen's d 动态效应矩阵
# ============================================================

library(tidyverse)

input_file <- file.path("results", "cfRNA", "ssGSEA", "11_dynamic_effect_matrix_for_clustering_CohenD.csv")

outdir <- file.path("results", "cfRNA", "modules")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

effect_df <- read.csv(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

period_cols <- c(
  "A_pre_le16",
  "B_mid_gt16_le32",
  "C_late_gt32"
)

stopifnot(
  "geneset_name" %in% colnames(effect_df),
  all(period_cols %in% colnames(effect_df))
)

effect_mat <- effect_df %>%
  select(geneset_name, all_of(period_cols)) %>%
  distinct(geneset_name, .keep_all = TRUE)

rownames(effect_mat) <- effect_mat$geneset_name

effect_mat <- effect_mat %>%
  select(-geneset_name) %>%
  as.matrix()

storage.mode(effect_mat) <- "numeric"

cat("通路数量：", nrow(effect_mat), "\n")
cat("时期数量：", ncol(effect_mat), "\n")
cat("缺失值数量：", sum(is.na(effect_mat)), "\n")

head(effect_mat)

# ============================================================
# Step 2. 基于 A/B/C 三期 Cohen's d 进行层次聚类
# ============================================================

# 计算通路之间的距离
dist_mat <- dist(effect_mat, method = "euclidean")

# Ward.D2 层次聚类
hc <- hclust(dist_mat, method = "ward.D2")

# 保存聚类结果
saveRDS(
  hc,
  file.path(outdir, "Step2_pathway_hclust.rds")
)

# 绘制聚类树
pdf(
  file.path(outdir, "Step2_pathway_dendrogram.pdf"),
  width = 12,
  height = 6
)

plot(
  hc,
  labels = FALSE,
  hang = -1,
  main = "Pathway clustering based on A/B/C Cohen's d",
  xlab = "Pathways",
  ylab = "Height"
)

dev.off()

cat("聚类通路数：", length(hc$order), "\n")
cat("聚类树已保存。\n")

# ============================================================
# Step 3. dynamicTreeCut 划分动态通路模块
# ============================================================

library(dynamicTreeCut)

module_id <- cutreeDynamic(
  dendro = hc,
  distM = as.matrix(dist_mat),
  method = "hybrid",
  deepSplit = 2,
  minClusterSize = 30,
  pamRespectsDendro = FALSE,
  verbose = 0
)

# 把数字模块改成统一名称
valid_id <- sort(unique(module_id[module_id != 0]))

module_name <- ifelse(
  module_id == 0,
  "Unassigned",
  sprintf(
    "Dynamic_Module_%02d",
    match(module_id, valid_id)
  )
)

# 保存每条通路属于哪个模块
module_assignment <- data.frame(
  geneset_name = rownames(effect_mat),
  module_id = module_id,
  dynamic_module = module_name,
  effect_mat,
  check.names = FALSE
)

write.csv(
  module_assignment,
  file.path(outdir, "Step3_pathway_module_assignment.csv"),
  row.names = FALSE
)

# 查看每个模块包含多少条通路
module_size <- as.data.frame(table(module_name))
colnames(module_size) <- c("dynamic_module", "n_pathways")

module_size <- module_size[order(-module_size$n_pathways), ]

print(module_size)

# ============================================================
# Step 4. 计算并绘制动态模块的三期效应
# ============================================================

module_profile <- module_assignment %>%
  filter(dynamic_module != "Unassigned") %>%
  group_by(dynamic_module) %>%
  summarise(
    n_pathways = n(),
    A_pre_le16 = mean(A_pre_le16),
    B_mid_gt16_le32 = mean(B_mid_gt16_le32),
    C_late_gt32 = mean(C_late_gt32),
    .groups = "drop"
  )

write.csv(
  module_profile,
  file.path(outdir, "Step4_dynamic_module_profile.csv"),
  row.names = FALSE
)

profile_long <- module_profile %>%
  pivot_longer(
    cols = all_of(period_cols),
    names_to = "period",
    values_to = "mean_cohen_d"
  ) %>%
  mutate(
    period = factor(period, levels = period_cols)
  )

p <- ggplot(
  profile_long,
  aes(
    x = period,
    y = dynamic_module,
    fill = mean_cohen_d
  )
) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = sprintf("%.2f", mean_cohen_d)),
    size = 3
  ) +
  scale_fill_gradient2(
    low = "#3B6FB6",
    mid = "white",
    high = "#C94C4C",
    midpoint = 0
  ) +
  labs(
    x = NULL,
    y = "Dynamic module",
    fill = "Mean Cohen's d"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(
  file.path(outdir, "Step4_dynamic_module_heatmap.pdf"),
  p,
  width = 7,
  height = max(5, nrow(module_profile) * 0.35)
)

print(module_profile)
print(p)

# ============================================================
# Step 5. 计算每个样本的动态模块分数
# ============================================================

ssgsea_file <- file.path("results", "cfRNA", "ssGSEA", "08_GSEA_sig_GO_BP_full_ssgsea_scores_geneset_by_sample.csv")

ssgsea_df <- read.csv(
  ssgsea_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# 只保留已分配模块的通路
pathway_module <- module_assignment %>%
  filter(dynamic_module != "Unassigned") %>%
  select(geneset_name, dynamic_module)

ssgsea_use <- ssgsea_df %>%
  inner_join(pathway_module, by = "geneset_name")

# 通路 × 样本矩阵
score_mat <- ssgsea_use %>%
  select(-geneset_name, -dynamic_module) %>%
  as.matrix()

rownames(score_mat) <- ssgsea_use$geneset_name
storage.mode(score_mat) <- "numeric"

# 每条通路在所有样本中做 Z-score
score_z <- t(scale(t(score_mat)))

# 删除标准差为0导致的缺失通路
keep <- apply(score_z, 1, function(x) all(is.finite(x)))
score_z <- score_z[keep, , drop = FALSE]

pathway_module <- pathway_module %>%
  filter(geneset_name %in% rownames(score_z))

# 模块内通路取平均
module_score <- sapply(
  unique(pathway_module$dynamic_module),
  function(m) {
    genes <- pathway_module$geneset_name[
      pathway_module$dynamic_module == m
    ]
    colMeans(score_z[genes, , drop = FALSE])
  }
)

module_score <- data.frame(
  sample = rownames(module_score),
  module_score,
  check.names = FALSE
)

write.csv(
  module_score,
  file.path(outdir, "Step5_sample_module_scores.csv"),
  row.names = FALSE
)

cat("样本数：", nrow(module_score), "\n")
cat("模块数：", ncol(module_score) - 1, "\n")

head(module_score)

# ============================================================
# Step 6. 比较各时期 GDM 与 Healthy 的模块分数
# ============================================================

meta_file <- file.path("results", "cohort_selection", "20_FINAL_sample_table_matched_to_matrix.csv")

meta <- read.csv(
  meta_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

score_meta <- module_score %>%
  inner_join(
    meta %>% select(sample, group2, period_new),
    by = "sample"
  ) %>%
  filter(
    group2 %in% c("GDM", "Healthy"),
    period_new %in% period_cols
  )

score_long <- score_meta %>%
  pivot_longer(
    cols = starts_with("Dynamic_Module_"),
    names_to = "dynamic_module",
    values_to = "module_score"
  )

module_test <- score_long %>%
  group_by(period_new, dynamic_module) %>%
  summarise(
    n_GDM = sum(group2 == "GDM"),
    n_Healthy = sum(group2 == "Healthy"),
    median_GDM = median(module_score[group2 == "GDM"]),
    median_Healthy = median(module_score[group2 == "Healthy"]),
    median_diff = median_GDM - median_Healthy,
    p_value = wilcox.test(
      module_score[group2 == "GDM"],
      module_score[group2 == "Healthy"],
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  group_by(period_new) %>%
  mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  ungroup()

write.csv(
  score_meta,
  file.path(outdir, "Step6_sample_module_scores_with_metadata.csv"),
  row.names = FALSE
)

write.csv(
  module_test,
  file.path(outdir, "Step6_module_GDM_vs_Healthy_statistics.csv"),
  row.names = FALSE
)

module_test %>%
  arrange(period_new, FDR) %>%
  print(n = Inf)

# ============================================================
# Step 7. 给动态模块补充 GO ID
# ============================================================

module_go <- module_assignment %>%
  select(
    geneset_name,
    dynamic_module,
    A_pre_le16,
    B_mid_gt16_le32,
    C_late_gt32
  ) %>%
  left_join(
    effect_df %>%
      select(
        geneset_name,
        GO_ID,
        GO_BP_name,
        Description
      ),
    by = "geneset_name"
  ) %>%
  filter(
    dynamic_module != "Unassigned",
    !is.na(GO_ID),
    grepl("^GO:", GO_ID)
  ) %>%
  distinct(dynamic_module, GO_ID, .keep_all = TRUE)

write.csv(
  module_go,
  file.path(outdir, "Step7_dynamic_modules_with_GOID.csv"),
  row.names = FALSE
)

cat("有 GO ID 的通路数：", nrow(module_go), "\n")
cat("动态模块数：", length(unique(module_go$dynamic_module)), "\n")
cat("缺失 GO ID 数：", sum(is.na(module_go$GO_ID)), "\n")

table(module_go$dynamic_module)

# ============================================================
# Step 8. 模块内 GO 语义相似性聚类
# ============================================================

library(GOSemSim)
library(org.Hs.eg.db)

sem_data <- godata(
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  computeIC = FALSE
)

assignment_list <- list()
representative_list <- list()

for (m in unique(module_go$dynamic_module)) {

  tmp <- module_go %>%
    filter(dynamic_module == m) %>%
    distinct(GO_ID, .keep_all = TRUE)

  go_ids <- tmp$GO_ID

  # GO term 太少时不再细分
  if (length(go_ids) < 3) {
    tmp$semantic_subcluster <- paste0(m, "_S01")
    tmp$is_representative <- FALSE
    tmp$is_representative[1] <- TRUE

    assignment_list[[m]] <- tmp
    representative_list[[m]] <- tmp[1, ]
    next
  }

  # Wang GO语义相似性
  sim_mat <- mgoSim(
    go_ids,
    go_ids,
    sem_data,
    measure = "Wang",
    combine = NULL
  )

  sim_mat[is.na(sim_mat)] <- 0
  diag(sim_mat) <- 1

  # 相似性转成距离后聚类
  hc_sem <- hclust(
    as.dist(1 - sim_mat),
    method = "average"
  )

  # 距离0.3，相当于大致按相似性0.7进行去冗余
  sub_id <- cutree(hc_sem, h = 0.3)

  tmp$semantic_subcluster <- paste0(
    m, "_S", sprintf("%02d", sub_id[tmp$GO_ID])
  )

  # 每个小簇选择平均语义相似性最高的代表term
  tmp$is_representative <- FALSE

  for (s in unique(tmp$semantic_subcluster)) {

    ids <- tmp$GO_ID[tmp$semantic_subcluster == s]

    if (length(ids) == 1) {
      representative_go <- ids
    } else {
      representative_go <- names(
        which.max(rowMeans(sim_mat[ids, ids, drop = FALSE]))
      )
    }

    tmp$is_representative[tmp$GO_ID == representative_go] <- TRUE
  }

  assignment_list[[m]] <- tmp
  representative_list[[m]] <- tmp %>%
    filter(is_representative)
}

semantic_assignment <- bind_rows(assignment_list)
representative_terms <- bind_rows(representative_list)

write.csv(
  semantic_assignment,
  file.path(outdir, "Step8_GO_semantic_subcluster_assignment.csv"),
  row.names = FALSE
)

write.csv(
  representative_terms,
  file.path(outdir, "Step8_GO_semantic_representative_terms.csv"),
  row.names = FALSE
)

semantic_summary <- semantic_assignment %>%
  group_by(dynamic_module) %>%
  summarise(
    n_GO_terms = n(),
    n_semantic_subclusters = n_distinct(semantic_subcluster),
    n_representative_terms = sum(is_representative),
    .groups = "drop"
  )

print(semantic_summary)

# ============================================================
# Step 8-1. 汇总每个模块的主要语义小组
# ============================================================

library(dplyr)

# 每个语义小组包含多少条原始 GO
semantic_size <- semantic_assignment %>%
  dplyr::count(
    dynamic_module,
    semantic_subcluster,
    name = "n_original_GO"
  )

# 给每条代表 GO 补上它代表了多少条原始 GO
representative_summary <- representative_terms %>%
  dplyr::left_join(
    semantic_size,
    by = c("dynamic_module", "semantic_subcluster")
  )

# 每个动态模块取最大的 10 个语义小组
top10_semantic_groups <- representative_summary %>%
  dplyr::group_by(dynamic_module) %>%
  dplyr::arrange(
    dplyr::desc(n_original_GO),
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    dynamic_module,
    semantic_subcluster,
    GO_ID,
    Description,
    n_original_GO,
    A_pre_le16,
    B_mid_gt16_le32,
    C_late_gt32
  )

write.csv(
  representative_summary,
  file.path(outdir, "Step8_1_all_representative_GO_with_size.csv"),
  row.names = FALSE
)

write.csv(
  top10_semantic_groups,
  file.path(outdir, "Step8_1_top10_semantic_groups_each_module.csv"),
  row.names = FALSE
)

print(top10_semantic_groups, n = Inf)

# ============================================================
# Step 8-2. 计算每个语义小组的平均 Cohen's d
# ============================================================

semantic_effect <- semantic_assignment %>%
  dplyr::group_by(dynamic_module, semantic_subcluster) %>%
  dplyr::summarise(
    n_original_GO = dplyr::n(),
    A_mean = mean(A_pre_le16, na.rm = TRUE),
    B_mean = mean(B_mid_gt16_le32, na.rm = TRUE),
    C_mean = mean(C_late_gt32, na.rm = TRUE),
    .groups = "drop"
  )

# 加入每个语义小组的代表 GO 名称
representative_info <- representative_terms %>%
  dplyr::select(
    dynamic_module,
    semantic_subcluster,
    GO_ID,
    Description
  ) %>%
  dplyr::rename(
    representative_GO_ID = GO_ID,
    representative_GO = Description
  )

semantic_effect <- semantic_effect %>%
  dplyr::left_join(
    representative_info,
    by = c("dynamic_module", "semantic_subcluster")
  ) %>%
  dplyr::select(
    dynamic_module,
    semantic_subcluster,
    representative_GO_ID,
    representative_GO,
    n_original_GO,
    A_mean,
    B_mean,
    C_mean
  )

write.csv(
  semantic_effect,
  file.path(outdir, "Step8_2_semantic_group_mean_effect.csv"),
  row.names = FALSE
)

head(semantic_effect)
