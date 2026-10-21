# ============================================================
# GDM cfRNA 样本筛选 + 协变量平衡 + 表达矩阵过滤【优化完整可复用版】
# 版本：v3_bugfix
#
# 目标：
# 1. 根据 collectGA 正确划分 A/B/C 孕期：A <= 16；B >16 & <=32；C >32
# 2. 尽量避免同一 Family 在同一时期重复采样
# 3. 默认保留全部 GDM；Healthy 重复样本和多余样本优先删除
# 4. 在样本表层面尽可能降低 age 和 collectGA 在 GDM/Healthy 之间的不平衡
# 5. 输出详细平衡指标：SMD、均值差、中位数差、t 检验、Wilcoxon、KS 检验、综合 balance score
# 6. 给每个样本计算协变量诊断指标：到对侧组中心距离、propensity score 等
# 7. 过滤表达矩阵，并严格检查最终样本表和矩阵列顺序完全一致
#
# 重要说明：
# - 样本层面平衡只能“降低” age 和 collectGA 的混杂风险，不能数学上保证完全消除。
# - 后续差异分析如果坚持 group-only 模型，需要在方法中报告本脚本输出的最终 SMD。
# ============================================================


# ============================================================
# 0. CONFIG 参数区：以后主要只改这里
# ============================================================

# ----------------------------
# 输入文件
# ----------------------------
sample_file <- file.path("data", "controlled", "最终样本表-信息最全-over.xlsx")
expr_file   <- file.path("data", "controlled", "最终表达矩阵-信息最全-over.tsv")

# ----------------------------
# 输出目录
# ----------------------------
outdir <- file.path("results", "cohort_selection")

# ----------------------------
# 样本表关键列名
# ----------------------------
sample_col <- "sample"
family_col <- "Family"
group_col  <- "group_main"
age_col    <- "age"
ga_col     <- "collectGA"

# ----------------------------
# 表达矩阵设置
# ----------------------------
gene_col <- "sym_id"     # 如果表达矩阵没有该列，会默认把第一列改成 gene_col
expr_sep <- "\t"         # tsv 用 "\t"，csv 用 ","

# ----------------------------
# 分组名称识别
# ----------------------------
healthy_patterns <- c("Healthy", "healthy", "Control", "control", "健康", "对照")
gdm_patterns     <- c("GDM", "GDM_zong", "gdm", "GDMM")

# ----------------------------
# 孕期划分：最终推荐方案
# A: collectGA <= 16
# B: 16 < collectGA <= 32
# C: collectGA > 32
# ----------------------------
A_upper <- 16
B_upper <- 32

period_labels <- c(
  A = "A_pre_le16",
  B = "B_mid_gt16_le32",
  C = "C_late_gt32"
)

# ----------------------------
# 协变量平衡设置
# ----------------------------
covariates_to_balance <- c("age", "collectGA")

# 权重：如果你认为 collectGA 比 age 更重要，可以设 collectGA = 1.5 或 2
covariate_weights <- c(age = 1, collectGA = 1)

# SMD 目标阈值：通常 |SMD| < 0.1 认为平衡较好
smd_cutoff <- 0.1

# Healthy 至少保留多少倍 GDM
# 1 = 至少 1:1；2 = Healthy 至少是 GDM 的 2 倍
min_control_ratio <- 1

# 如果 TRUE：一旦达到 SMD 阈值，就停止继续删除 Healthy，从而尽量保留更多样本
stop_when_smd_pass <- TRUE

# ----------------------------
# GDM 同一 Family 同一时期重复采样策略
# keep_report：默认，不删除 GDM，只输出重复报告和 flag
# stop：如果发现 GDM 重复，直接停止，让你手动决定
# dedup：也对 GDM 去重，保留最接近 Healthy 中心的一个样本
# ----------------------------
gdm_duplicate_strategy <- "keep_report"

# Family 重复去重时，age 和 collectGA 的距离权重
# 默认 age 和 collectGA 都参与；如果更重视孕周，可以 collectGA = 1.5
family_dedup_weights <- c(age = 1, collectGA = 1)

# ----------------------------
# 输出图设置
# ----------------------------
plot_width  <- 9
plot_height <- 5
plot_dpi    <- 300


# ============================================================
# 1. 加载 R 包
# ============================================================

required_pkgs <- c("dplyr", "tidyr", "ggplot2", "readr", "stringr", "purrr")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("缺少 R 包：", pkg, "。请先 install.packages('", pkg, "')")
  }
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)
library(purrr)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("开始运行：GDM cfRNA 样本筛选 + 协变量平衡 + 矩阵过滤【优化版】\n")
cat("样本表：", sample_file, "\n")
cat("表达矩阵：", expr_file, "\n")
cat("输出目录：", outdir, "\n")
cat("孕期划分：A <= ", A_upper, "; B > ", A_upper, " & <= ", B_upper, "; C > ", B_upper, "\n", sep = "")
cat("协变量平衡：", paste(covariates_to_balance, collapse = ", "), "\n")
cat("============================================================\n\n")


# ============================================================
# 2. 工具函数
# ============================================================

make_pattern <- function(patterns) {
  paste(patterns, collapse = "|")
}

safe_max_abs <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(Inf)
  max(abs(x))
}

safe_mean_abs <- function(x, w = NULL) {
  ok <- is.finite(x)
  x <- abs(x[ok])
  if (length(x) == 0) return(Inf)
  if (is.null(w)) return(mean(x))
  w <- w[ok]
  if (length(w) != length(x) || sum(w, na.rm = TRUE) == 0) return(mean(x))
  sum(x * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
}

# 解析孕周：支持 16.4、16、16+3、16周3天、16w3d 等格式
parse_ga_to_weeks <- function(x) {
  x_chr <- as.character(x)
  x_chr <- str_trim(x_chr)
  x_chr <- str_replace_all(x_chr, "＋", "+")
  x_chr <- str_replace_all(x_chr, "周", "w")
  x_chr <- str_replace_all(x_chr, "天", "d")
  x_chr <- str_replace_all(x_chr, "W", "w")
  x_chr <- str_replace_all(x_chr, "D", "d")
  
  out <- rep(NA_real_, length(x_chr))
  
  # 形式：16+3
  idx_plus <- str_detect(x_chr, "^\\s*\\d+(\\.\\d+)?\\s*\\+\\s*\\d+(\\.\\d+)?\\s*$")
  if (any(idx_plus, na.rm = TRUE)) {
    week <- as.numeric(str_extract(x_chr[idx_plus], "^\\s*\\d+(\\.\\d+)?"))
    day  <- as.numeric(str_extract(x_chr[idx_plus], "(?<=\\+)\\s*\\d+(\\.\\d+)?"))
    out[idx_plus] <- week + day / 7
  }
  
  # 形式：16w3d 或 16w 或 16w 3d
  idx_wd <- is.na(out) & str_detect(x_chr, "\\d+(\\.\\d+)?\\s*w")
  if (any(idx_wd, na.rm = TRUE)) {
    week <- as.numeric(str_extract(x_chr[idx_wd], "\\d+(\\.\\d+)?(?=\\s*w)"))
    day_txt <- str_extract(x_chr[idx_wd], "(?<=w)\\s*\\d+(\\.\\d+)?(?=\\s*d)")
    day <- suppressWarnings(as.numeric(day_txt))
    day[is.na(day)] <- 0
    out[idx_wd] <- week + day / 7
  }
  
  # 普通数字形式
  idx_num <- is.na(out)
  if (any(idx_num, na.rm = TRUE)) {
    out[idx_num] <- readr::parse_number(x_chr[idx_num])
  }
  
  return(out)
}

read_sample_any <- function(path) {
  if (str_detect(tolower(path), "\\.xlsx$|\\.xls$")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("样本表是 Excel 文件，但缺少 readxl 包。请先 install.packages('readxl')")
    }
    return(readxl::read_excel(path))
  }
  if (str_detect(tolower(path), "\\.csv$")) {
    return(readr::read_csv(path, show_col_types = FALSE))
  }
  return(readr::read_delim(path, delim = "\t", show_col_types = FALSE))
}

# 更安全的分组识别：先 trim，再小写，再精确匹配；避免 non-GDM 被误判成 GDM
normalize_group <- function(x) {
  x0 <- stringr::str_trim(as.character(x))
  x_low <- stringr::str_to_lower(x0)
  healthy_low <- stringr::str_to_lower(stringr::str_trim(healthy_patterns))
  gdm_low <- stringr::str_to_lower(stringr::str_trim(gdm_patterns))
  dplyr::case_when(
    x_low %in% gdm_low ~ "GDM",
    x_low %in% healthy_low ~ "Healthy",
    TRUE ~ NA_character_
  )
}

smd_cont <- function(x, g, case_group = "GDM", control_group = "Healthy") {
  x1 <- x[g == case_group]
  x0 <- x[g == control_group]
  x1 <- x1[!is.na(x1)]
  x0 <- x0[!is.na(x0)]
  if (length(x1) < 2 || length(x0) < 2) return(NA_real_)
  s_pool <- sqrt((var(x1) + var(x0)) / 2)
  if (!is.finite(s_pool) || s_pool == 0) return(NA_real_)
  (mean(x1) - mean(x0)) / s_pool
}

safe_t_p <- function(x, g) {
  x1 <- x[g == "GDM"]
  x0 <- x[g == "Healthy"]
  x1 <- x1[!is.na(x1)]
  x0 <- x0[!is.na(x0)]
  if (length(x1) < 2 || length(x0) < 2) return(NA_real_)
  tryCatch(t.test(x1, x0)$p.value, error = function(e) NA_real_)
}

safe_wilcox_p <- function(x, g) {
  x1 <- x[g == "GDM"]
  x0 <- x[g == "Healthy"]
  x1 <- x1[!is.na(x1)]
  x0 <- x0[!is.na(x0)]
  if (length(x1) < 2 || length(x0) < 2) return(NA_real_)
  tryCatch(suppressWarnings(wilcox.test(x1, x0, exact = FALSE)$p.value), error = function(e) NA_real_)
}

safe_ks_p <- function(x, g) {
  x1 <- x[g == "GDM"]
  x0 <- x[g == "Healthy"]
  x1 <- x1[!is.na(x1)]
  x0 <- x0[!is.na(x0)]
  if (length(x1) < 2 || length(x0) < 2) return(NA_real_)
  tryCatch(suppressWarnings(ks.test(x1, x0)$p.value), error = function(e) NA_real_)
}

calc_balance_metric_one <- function(df, variable) {
  x <- df[[variable]]
  g <- df$group2
  x_gdm <- x[g == "GDM"]
  x_ctl <- x[g == "Healthy"]
  
  tibble(
    variable = variable,
    n_GDM = sum(g == "GDM"),
    n_Healthy = sum(g == "Healthy"),
    GDM_mean = mean(x_gdm, na.rm = TRUE),
    Healthy_mean = mean(x_ctl, na.rm = TRUE),
    mean_diff_GDM_minus_Healthy = GDM_mean - Healthy_mean,
    GDM_sd = sd(x_gdm, na.rm = TRUE),
    Healthy_sd = sd(x_ctl, na.rm = TRUE),
    GDM_median = median(x_gdm, na.rm = TRUE),
    Healthy_median = median(x_ctl, na.rm = TRUE),
    median_diff_GDM_minus_Healthy = GDM_median - Healthy_median,
    GDM_IQR = IQR(x_gdm, na.rm = TRUE),
    Healthy_IQR = IQR(x_ctl, na.rm = TRUE),
    GDM_min = min(x_gdm, na.rm = TRUE),
    GDM_max = max(x_gdm, na.rm = TRUE),
    Healthy_min = min(x_ctl, na.rm = TRUE),
    Healthy_max = max(x_ctl, na.rm = TRUE),
    SMD = smd_cont(x, g),
    abs_SMD = abs(SMD),
    variance_ratio_GDM_over_Healthy = ifelse(
      is.finite(var(x_ctl, na.rm = TRUE)) && var(x_ctl, na.rm = TRUE) > 0,
      var(x_gdm, na.rm = TRUE) / var(x_ctl, na.rm = TRUE),
      NA_real_
    ),
    t_p = safe_t_p(x, g),
    wilcox_p = safe_wilcox_p(x, g),
    ks_p = safe_ks_p(x, g)
  )
}

calc_balance_long_one_period <- function(df_period, stage_name = "tmp") {
  if (nrow(df_period) == 0) return(tibble())
  map_dfr(covariates_to_balance, function(v) {
    calc_balance_metric_one(df_period, v)
  }) %>%
    mutate(
      stage = stage_name,
      period_new = as.character(unique(df_period$period_new)[1]),
      .before = 1
    )
}

calc_balance_long <- function(df, stage_name = "tmp") {
  if (nrow(df) == 0) return(tibble())
  periods <- unique(as.character(df$period_new))
  map_dfr(periods, function(p) {
    calc_balance_long_one_period(df %>% filter(as.character(period_new) == p), stage_name)
  })
}

calc_balance_wide <- function(df, stage_name = "tmp") {
  long <- calc_balance_long(df, stage_name)
  if (nrow(long) == 0) return(tibble())
  
  score_df <- long %>%
    mutate(weight = covariate_weights[variable]) %>%
    mutate(weight = ifelse(is.na(weight), 1, weight)) %>%
    group_by(stage, period_new) %>%
    summarise(
      n_GDM = max(n_GDM, na.rm = TRUE),
      n_Healthy = max(n_Healthy, na.rm = TRUE),
      max_abs_SMD = safe_max_abs(SMD),
      weighted_mean_abs_SMD = safe_mean_abs(SMD, weight),
      balance_score = max_abs_SMD + 0.25 * weighted_mean_abs_SMD,
      pass_all_SMD = all(is.finite(abs_SMD)) && all(abs_SMD <= smd_cutoff),
      .groups = "drop"
    )
  
  smd_wide <- long %>%
    select(stage, period_new, variable, SMD, abs_SMD, mean_diff_GDM_minus_Healthy, median_diff_GDM_minus_Healthy, t_p, wilcox_p, ks_p) %>%
    pivot_wider(
      names_from = variable,
      values_from = c(SMD, abs_SMD, mean_diff_GDM_minus_Healthy, median_diff_GDM_minus_Healthy, t_p, wilcox_p, ks_p),
      names_sep = "_"
    )
  
  score_df %>% left_join(smd_wide, by = c("stage", "period_new"))
}

balance_objective <- function(df_period) {
  b <- calc_balance_long_one_period(df_period, "tmp")
  if (nrow(b) == 0) return(Inf)
  w <- covariate_weights[b$variable]
  w[is.na(w)] <- 1
  max_abs <- safe_max_abs(b$SMD)
  w_mean <- safe_mean_abs(b$SMD, w)
  max_abs + 0.25 * w_mean
}

calc_distance_to_reference_center <- function(candidates, reference_df) {
  if (nrow(candidates) == 0) return(numeric(0))
  out <- rep(0, nrow(candidates))
  for (v in names(family_dedup_weights)) {
    if (!v %in% colnames(candidates) || !v %in% colnames(reference_df)) next
    ref_med <- median(reference_df[[v]], na.rm = TRUE)
    pooled_sd <- sd(c(candidates[[v]], reference_df[[v]]), na.rm = TRUE)
    if (!is.finite(pooled_sd) || pooled_sd == 0) pooled_sd <- 1
    out <- out + family_dedup_weights[[v]] * abs(candidates[[v]] - ref_med) / pooled_sd
  }
  out
}

dedup_group_by_family <- function(target_df, reference_df, target_group_name) {
  if (nrow(target_df) == 0) return(target_df)
  if (nrow(reference_df) == 0) {
    target_df <- target_df %>%
      group_by(Family_id) %>%
      arrange(collectGA, age, .by_group = TRUE) %>%
      mutate(family_dedup_rank = row_number()) %>%
      ungroup()
  } else {
    target_df <- target_df %>%
      mutate(
        family_dedup_distance_to_opposite = calc_distance_to_reference_center(., reference_df)
      ) %>%
      group_by(Family_id) %>%
      arrange(family_dedup_distance_to_opposite, collectGA, age, .by_group = TRUE) %>%
      mutate(family_dedup_rank = row_number()) %>%
      ungroup()
  }
  
  target_df %>%
    mutate(
      family_dedup_status = ifelse(family_dedup_rank == 1, "kept_by_family_dedup", paste0("removed_", target_group_name, "_family_duplicate"))
    )
}

dedup_within_period <- function(df_period) {
  p <- as.character(unique(df_period$period_new)[1])
  gdm <- df_period %>% filter(group2 == "GDM")
  healthy <- df_period %>% filter(group2 == "Healthy")
  
  gdm_dup <- gdm %>% count(Family_id, name = "n_gdm_in_period") %>% filter(n_gdm_in_period > 1)
  
  if (nrow(gdm_dup) > 0 && gdm_duplicate_strategy == "stop") {
    stop("发现 GDM 同一 Family 在同一时期重复采样，时期：", p,
         "。请查看输出的 GDM 重复表后手动决定，或将 gdm_duplicate_strategy 改为 keep_report / dedup。")
  }
  
  healthy_tagged <- dedup_group_by_family(healthy, gdm, "Healthy")
  healthy_kept <- healthy_tagged %>% filter(family_dedup_rank == 1)
  
  if (gdm_duplicate_strategy == "dedup") {
    gdm_tagged <- dedup_group_by_family(gdm, healthy, "GDM")
    gdm_kept <- gdm_tagged %>% filter(family_dedup_rank == 1)
  } else {
    # 重要修复：距离必须先在未分组数据上计算，得到长度 = nrow(gdm) 的向量。
    # 不能在 group_by(Family_id) 后用 calc_distance_to_reference_center(., healthy)，
    # 否则 grouped mutate 会在每个 Family 小组里收到整列距离，导致 size 不匹配报错。
    gdm_kept <- gdm %>%
      mutate(
        family_dedup_distance_to_opposite = calc_distance_to_reference_center(gdm, healthy)
      ) %>%
      group_by(Family_id) %>%
      mutate(
        family_dedup_rank = row_number(),
        family_dedup_status = ifelse(n() > 1, "kept_GDM_duplicate_report_only", "kept_by_family_dedup")
      ) %>%
      ungroup()
  }
  
  bind_rows(gdm_kept, healthy_kept)
}

choose_balanced_healthy <- function(df_period, smd_cutoff = 0.1, min_control_ratio = 1) {
  p <- as.character(unique(df_period$period_new)[1])
  gdm <- df_period %>% filter(group2 == "GDM")
  healthy <- df_period %>% filter(group2 == "Healthy")
  n_gdm <- nrow(gdm)
  min_ctrl_n <- ceiling(n_gdm * min_control_ratio)
  
  start_df <- bind_rows(gdm, healthy)
  start_score <- balance_objective(start_df)
  
  if (nrow(gdm) < 2 || nrow(healthy) < 2) {
    return(list(
      data = start_df %>% mutate(selection_status = "too_few_samples", balance_selection_rank = NA_integer_),
      status = "too_few_samples",
      removal_log = tibble(),
      final_score = start_score
    ))
  }
  
  current_ctrl <- healthy
  best_ctrl <- current_ctrl
  best_df <- bind_rows(gdm, best_ctrl)
  best_score <- balance_objective(best_df)
  best_status <- "initial"
  removal_log <- tibble()
  iter <- 0
  
  repeat {
    iter <- iter + 1
    current_df <- bind_rows(gdm, current_ctrl)
    current_long <- calc_balance_long_one_period(current_df, "current")
    current_score <- balance_objective(current_df)
    current_max_smd <- safe_max_abs(current_long$SMD)
    
    if (is.finite(current_score) && current_score < best_score) {
      best_score <- current_score
      best_ctrl <- current_ctrl
      best_df <- current_df
      best_status <- "best_available"
    }
    
    if (stop_when_smd_pass && is.finite(current_max_smd) && current_max_smd <= smd_cutoff) {
      final_df <- current_df %>%
        mutate(selection_status = "PASS_smd_less_than_cutoff", balance_selection_rank = NA_integer_)
      return(list(
        data = final_df,
        status = "PASS_smd_less_than_cutoff",
        removal_log = removal_log,
        final_score = current_score
      ))
    }
    
    if (nrow(current_ctrl) <= min_ctrl_n) {
      final_df <- best_df %>%
        mutate(selection_status = "NOT_PASS_keep_best_available", balance_selection_rank = NA_integer_)
      return(list(
        data = final_df,
        status = "NOT_PASS_keep_best_available",
        removal_log = removal_log,
        final_score = best_score
      ))
    }
    
    remove_score <- map_dbl(seq_len(nrow(current_ctrl)), function(i) {
      candidate_ctrl <- current_ctrl[-i, , drop = FALSE]
      candidate_df <- bind_rows(gdm, candidate_ctrl)
      balance_objective(candidate_df)
    })
    
    remove_score[!is.finite(remove_score)] <- Inf
    if (all(is.infinite(remove_score))) {
      final_df <- best_df %>%
        mutate(selection_status = "NOT_PASS_no_valid_removal", balance_selection_rank = NA_integer_)
      return(list(
        data = final_df,
        status = "NOT_PASS_no_valid_removal",
        removal_log = removal_log,
        final_score = best_score
      ))
    }
    
    remove_i <- which.min(remove_score)
    removed_one <- current_ctrl[remove_i, , drop = FALSE] %>%
      mutate(
        period_new = p,
        balance_selection_iteration = iter,
        score_before_removal = current_score,
        predicted_score_after_removal = remove_score[remove_i],
        removal_reason = "removed_Healthy_to_improve_age_collectGA_balance"
      )
    removal_log <- bind_rows(removal_log, removed_one)
    current_ctrl <- current_ctrl[-remove_i, , drop = FALSE]
  }
}

add_covariate_diagnostics_one_period <- function(df_period) {
  if (nrow(df_period) == 0) return(df_period)
  gdm <- df_period %>% filter(group2 == "GDM")
  healthy <- df_period %>% filter(group2 == "Healthy")
  
  out <- df_period
  
  for (v in covariates_to_balance) {
    s <- sd(out[[v]], na.rm = TRUE)
    if (!is.finite(s) || s == 0) s <- 1
    m <- mean(out[[v]], na.rm = TRUE)
    out[[paste0(v, "_z_within_period")]] <- (out[[v]] - m) / s
  }
  
  out <- out %>%
    mutate(
      distance_to_opposite_center = case_when(
        group2 == "Healthy" ~ calc_distance_to_reference_center(., gdm),
        group2 == "GDM" ~ calc_distance_to_reference_center(., healthy),
        TRUE ~ NA_real_
      )
    )
  
  if (nrow(gdm) >= 2 && nrow(healthy) >= 2) {
    tmp <- out %>% mutate(gdm01 = ifelse(group2 == "GDM", 1, 0))
    for (v in covariates_to_balance) {
      s <- sd(tmp[[v]], na.rm = TRUE)
      if (!is.finite(s) || s == 0) s <- 1
      tmp[[paste0(v, "_model_z")]] <- (tmp[[v]] - mean(tmp[[v]], na.rm = TRUE)) / s
    }
    form <- as.formula(paste("gdm01 ~", paste(paste0(covariates_to_balance, "_model_z"), collapse = " + ")))
    fit <- tryCatch(
      suppressWarnings(glm(form, data = tmp, family = binomial())),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      out$propensity_score_age_collectGA <- tryCatch(
        suppressWarnings(as.numeric(predict(fit, newdata = tmp, type = "response"))),
        error = function(e) NA_real_
      )
      out$propensity_score_extreme_flag <- ifelse(
        is.na(out$propensity_score_age_collectGA),
        NA,
        out$propensity_score_age_collectGA < 0.05 | out$propensity_score_age_collectGA > 0.95
      )
    } else {
      out$propensity_score_age_collectGA <- NA_real_
      out$propensity_score_extreme_flag <- NA
    }
  } else {
    out$propensity_score_age_collectGA <- NA_real_
    out$propensity_score_extreme_flag <- NA
  }
  
  out
}

add_covariate_diagnostics <- function(df) {
  if (nrow(df) == 0) return(df)
  # 重要修复：不用 group_modify()，因为 group_modify 默认不把 period_new 放进 .x，
  # 容易出现 Unknown or uninitialised column: period_new。
  split(df, as.character(df$period_new), drop = TRUE) %>%
    purrr::map_dfr(add_covariate_diagnostics_one_period) %>%
    mutate(period_new = factor(period_new, levels = unname(period_labels)))
}

write_csv_utf8 <- function(x, file) {
  write.csv(x, file = file, row.names = FALSE, fileEncoding = "UTF-8")
}


# ============================================================
# 3. 读取并清洗样本表
# ============================================================

sample_raw <- read_sample_any(sample_file)

required_cols <- c(sample_col, family_col, group_col, age_col, ga_col)
missing_cols <- setdiff(required_cols, colnames(sample_raw))
if (length(missing_cols) > 0) {
  stop("样本表缺少以下列，请检查 CONFIG 参数区列名：\n", paste(missing_cols, collapse = ", "))
}

write_csv_utf8(data.frame(colname = colnames(sample_raw)), file.path(outdir, "00_sample_column_names.csv"))

healthy_pattern <- make_pattern(healthy_patterns)
gdm_pattern <- make_pattern(gdm_patterns)

sample_clean0 <- sample_raw %>%
  mutate(
    sample_id = stringr::str_trim(as.character(.data[[sample_col]])),
    Family_id = stringr::str_trim(as.character(.data[[family_col]])),
    group_original = stringr::str_trim(as.character(.data[[group_col]])),
    age = readr::parse_number(as.character(.data[[age_col]])),
    collectGA_raw = as.character(.data[[ga_col]]),
    collectGA = parse_ga_to_weeks(.data[[ga_col]]),
    group2 = normalize_group(.data[[group_col]])
  )

sample_excluded_basic <- sample_clean0 %>%
  filter(
    is.na(sample_id) | is.na(Family_id) | is.na(group2) |
      is.na(age) | is.na(collectGA)
  ) %>%
  rowwise() %>%
  mutate(
    exclude_reason = paste(
      na.omit(c(
        ifelse(is.na(sample_id), "missing_sample_id", NA),
        ifelse(is.na(Family_id), "missing_Family_id", NA),
        ifelse(is.na(group2), "unrecognized_group", NA),
        ifelse(is.na(age), "missing_or_unparsed_age", NA),
        ifelse(is.na(collectGA), "missing_or_unparsed_collectGA", NA)
      )),
      collapse = ";"
    )
  ) %>%
  ungroup()

sample_clean <- sample_clean0 %>%
  filter(!is.na(sample_id), !is.na(Family_id), !is.na(group2), !is.na(age), !is.na(collectGA)) %>%
  filter(group2 %in% c("Healthy", "GDM"))

if (anyDuplicated(sample_clean$sample_id) > 0) {
  dup_sample <- sample_clean %>% count(sample_id, name = "n") %>% filter(n > 1)
  write_csv_utf8(dup_sample, file.path(outdir, "00_duplicate_sample_ids_in_sample_table.csv"))
  stop("样本表中存在重复 sample_id，请查看 00_duplicate_sample_ids_in_sample_table.csv")
}

range_problem <- sample_clean %>%
  filter(age < 15 | age > 55 | collectGA < 5 | collectGA > 45)
if (nrow(range_problem) > 0) {
  write_csv_utf8(range_problem, file.path(outdir, "01_WARNING_age_collectGA_out_of_reasonable_range.csv"))
  warning("发现 age 或 collectGA 超出常见范围，请查看 01_WARNING_age_collectGA_out_of_reasonable_range.csv")
}

write_csv_utf8(sample_clean, file.path(outdir, "01_sample_clean.csv"))
write_csv_utf8(sample_excluded_basic, file.path(outdir, "01_sample_excluded_missing_or_unrecognized.csv"))

cat("清洗后样本数量：\n")
print(table(sample_clean$group2))


# ============================================================
# 4. 重新划分孕期
# ============================================================

sample_period <- sample_clean %>%
  mutate(
    period_new = case_when(
      collectGA <= A_upper ~ period_labels[["A"]],
      collectGA > A_upper & collectGA <= B_upper ~ period_labels[["B"]],
      collectGA > B_upper ~ period_labels[["C"]],
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(period_new)) %>%
  mutate(period_new = factor(period_new, levels = unname(period_labels)))

sample_period <- add_covariate_diagnostics(sample_period)

period_summary_before <- sample_period %>%
  group_by(period_new, group2) %>%
  summarise(
    n_sample = n(),
    n_family = n_distinct(Family_id),
    age_mean = mean(age, na.rm = TRUE),
    age_sd = sd(age, na.rm = TRUE),
    age_median = median(age, na.rm = TRUE),
    age_IQR = IQR(age, na.rm = TRUE),
    collectGA_mean = mean(collectGA, na.rm = TRUE),
    collectGA_sd = sd(collectGA, na.rm = TRUE),
    collectGA_median = median(collectGA, na.rm = TRUE),
    collectGA_IQR = IQR(collectGA, na.rm = TRUE),
    collectGA_min = min(collectGA, na.rm = TRUE),
    collectGA_max = max(collectGA, na.rm = TRUE),
    .groups = "drop"
  )

balance_before_long <- calc_balance_long(sample_period, "before_selection")
balance_before_wide <- calc_balance_wide(sample_period, "before_selection")

write_csv_utf8(sample_period, file.path(outdir, "02_sample_with_period_and_covariate_diagnostics.csv"))
write_csv_utf8(period_summary_before, file.path(outdir, "03_period_summary_before_selection.csv"))
write_csv_utf8(balance_before_long, file.path(outdir, "04_balance_before_selection_long.csv"))
write_csv_utf8(balance_before_wide, file.path(outdir, "04_balance_before_selection_wide.csv"))

cat("\n原始分期后样本数量：\n")
print(table(sample_period$period_new, sample_period$group2))
cat("\n原始分期后平衡情况：\n")
print(balance_before_wide)


# ============================================================
# 5. 检查 Family-period 重复采样
# ============================================================

family_dup_before <- sample_period %>%
  count(period_new, group2, Family_id, name = "n_in_this_period") %>%
  filter(n_in_this_period > 1) %>%
  arrange(period_new, group2, desc(n_in_this_period), Family_id)

gdm_dup_before <- sample_period %>%
  filter(group2 == "GDM") %>%
  count(period_new, Family_id, name = "n_gdm_in_this_period") %>%
  filter(n_gdm_in_this_period > 1) %>%
  arrange(period_new, desc(n_gdm_in_this_period), Family_id)

write_csv_utf8(family_dup_before, file.path(outdir, "05_family_duplicates_before_dedup.csv"))
write_csv_utf8(gdm_dup_before, file.path(outdir, "05_GDM_family_duplicates_in_period.csv"))

if (nrow(gdm_dup_before) > 0 && gdm_duplicate_strategy == "keep_report") {
  warning(
    "发现 GDM 同一 Family 在同一时期重复采样。当前设置 gdm_duplicate_strategy = 'keep_report'，不会删除 GDM，只在最终表中保留 flag。\n",
    "请查看 05_GDM_family_duplicates_in_period.csv。"
  )
}


# ============================================================
# 6. Family-period 去重
# ============================================================

# 重要修复：不用 group_modify()，避免 period_new 在 .x 中被默认移除。
sample_dedup <- split(sample_period, as.character(sample_period$period_new), drop = TRUE) %>%
  purrr::map_dfr(dedup_within_period) %>%
  mutate(period_new = factor(period_new, levels = unname(period_labels)))

sample_dedup <- add_covariate_diagnostics(sample_dedup)

removed_by_family_dedup <- sample_period %>%
  anti_join(sample_dedup %>% select(sample_id), by = "sample_id") %>%
  mutate(removal_reason = "removed_by_family_period_dedup")

family_dup_after_dedup <- sample_dedup %>%
  count(period_new, group2, Family_id, name = "n_in_this_period") %>%
  filter(n_in_this_period > 1) %>%
  arrange(period_new, group2, desc(n_in_this_period), Family_id)

balance_after_dedup_long <- calc_balance_long(sample_dedup, "after_family_dedup")
balance_after_dedup_wide <- calc_balance_wide(sample_dedup, "after_family_dedup")

write_csv_utf8(sample_dedup, file.path(outdir, "06_sample_after_family_dedup.csv"))
write_csv_utf8(removed_by_family_dedup, file.path(outdir, "06_removed_by_family_dedup.csv"))
write_csv_utf8(family_dup_after_dedup, file.path(outdir, "06_family_duplicates_after_dedup.csv"))
write_csv_utf8(balance_after_dedup_long, file.path(outdir, "07_balance_after_family_dedup_long.csv"))
write_csv_utf8(balance_after_dedup_wide, file.path(outdir, "07_balance_after_family_dedup_wide.csv"))

cat("\nFamily 去重后样本数量：\n")
print(table(sample_dedup$period_new, sample_dedup$group2))


# ============================================================
# 7. Healthy 平衡筛选：尽可能降低 age 和 collectGA 的 SMD
# ============================================================

periods <- levels(droplevels(sample_dedup$period_new))
balanced_list <- setNames(vector("list", length(periods)), periods)

for (p in periods) {
  cat("\n开始平衡筛选时期：", p, "\n")
  df_p <- sample_dedup %>% filter(as.character(period_new) == p)
  balanced_list[[p]] <- choose_balanced_healthy(
    df_p,
    smd_cutoff = smd_cutoff,
    min_control_ratio = min_control_ratio
  )
  cat("完成：", p, "；状态：", balanced_list[[p]]$status, "；final_score=", balanced_list[[p]]$final_score, "\n")
}

sample_final <- map_dfr(names(balanced_list), function(p) {
  balanced_list[[p]]$data %>%
    mutate(period_new = factor(period_new, levels = unname(period_labels)))
})

sample_final <- add_covariate_diagnostics(sample_final)

selection_iteration_log <- map_dfr(names(balanced_list), function(p) {
  balanced_list[[p]]$removal_log
})

removed_by_balance <- sample_dedup %>%
  anti_join(sample_final %>% select(sample_id), by = "sample_id") %>%
  mutate(removal_reason = "removed_by_balance_selection")

removed_all <- bind_rows(
  removed_by_family_dedup,
  removed_by_balance
) %>%
  distinct(sample_id, .keep_all = TRUE)

balance_final_long <- calc_balance_long(sample_final, "final_balanced")
balance_final_wide <- calc_balance_wide(sample_final, "final_balanced") %>%
  left_join(
    sample_final %>% distinct(period_new, selection_status) %>% mutate(period_new = as.character(period_new)),
    by = "period_new"
  )

final_summary <- sample_final %>%
  group_by(period_new, group2) %>%
  summarise(
    n_sample = n(),
    n_family = n_distinct(Family_id),
    age_mean = mean(age, na.rm = TRUE),
    age_sd = sd(age, na.rm = TRUE),
    age_median = median(age, na.rm = TRUE),
    age_IQR = IQR(age, na.rm = TRUE),
    collectGA_mean = mean(collectGA, na.rm = TRUE),
    collectGA_sd = sd(collectGA, na.rm = TRUE),
    collectGA_median = median(collectGA, na.rm = TRUE),
    collectGA_IQR = IQR(collectGA, na.rm = TRUE),
    collectGA_min = min(collectGA, na.rm = TRUE),
    collectGA_max = max(collectGA, na.rm = TRUE),
    median_distance_to_opposite_center = median(distance_to_opposite_center, na.rm = TRUE),
    .groups = "drop"
  )

write_csv_utf8(sample_final, file.path(outdir, "08_FINAL_sample_table_balanced_with_covariate_metrics.csv"))
write_csv_utf8(balance_final_long, file.path(outdir, "09_FINAL_balance_after_selection_long.csv"))
write_csv_utf8(balance_final_wide, file.path(outdir, "09_FINAL_balance_after_selection_wide.csv"))
write_csv_utf8(final_summary, file.path(outdir, "10_FINAL_sample_summary.csv"))
write_csv_utf8(removed_all, file.path(outdir, "11_removed_samples_all_reasons.csv"))
write_csv_utf8(selection_iteration_log, file.path(outdir, "12_balance_selection_iteration_log.csv"))

cat("\n最终筛选后样本数量：\n")
print(table(sample_final$period_new, sample_final$group2))
cat("\n最终平衡情况：\n")
print(balance_final_wide)


# ============================================================
# 8. 汇总所有平衡指标并画图
# ============================================================

balance_all_long <- bind_rows(
  balance_before_long,
  balance_after_dedup_long,
  balance_final_long
)

balance_all_wide <- bind_rows(
  balance_before_wide,
  balance_after_dedup_wide,
  balance_final_wide %>% select(-any_of("selection_status"))
)

write_csv_utf8(balance_all_long, file.path(outdir, "13_balance_before_after_all_long.csv"))
write_csv_utf8(balance_all_wide, file.path(outdir, "13_balance_before_after_all_wide.csv"))

balance_plot_df <- balance_all_long %>%
  mutate(
    abs_SMD = abs(SMD),
    variable = recode(variable, age = "Age", collectGA = "collectGA"),
    stage = factor(stage, levels = c("before_selection", "after_family_dedup", "final_balanced"))
  )

p_smd <- ggplot(balance_plot_df, aes(x = period_new, y = abs_SMD, fill = stage)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = smd_cutoff, linetype = "dashed") +
  facet_wrap(~ variable, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    title = "Age and collectGA balance before/after selection",
    subtitle = paste0("Dashed line: |SMD| = ", smd_cutoff),
    x = "Period",
    y = "|SMD|",
    fill = "Stage"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(outdir, "14_SMD_before_after.png"), p_smd, width = plot_width, height = plot_height, dpi = plot_dpi)
ggsave(file.path(outdir, "14_SMD_before_after.pdf"), p_smd, width = plot_width, height = plot_height)

p_count <- sample_final %>%
  count(period_new, group2) %>%
  ggplot(aes(x = period_new, y = n, fill = group2)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  theme_bw(base_size = 12) +
  labs(title = "Final sample size by period and group", x = "Period", y = "Sample number", fill = "Group") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(outdir, "15_FINAL_sample_count_by_period_group.png"), p_count, width = plot_width, height = plot_height, dpi = plot_dpi)
ggsave(file.path(outdir, "15_FINAL_sample_count_by_period_group.pdf"), p_count, width = plot_width, height = plot_height)

p_age <- ggplot(sample_final, aes(x = group2, y = age, fill = group2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.12, alpha = 0.45, size = 1) +
  facet_wrap(~ period_new, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(title = "Final age distribution after selection", x = "Group", y = "Age") +
  theme(legend.position = "none")

ggsave(file.path(outdir, "16_FINAL_age_boxplot_by_period.png"), p_age, width = plot_width, height = plot_height, dpi = plot_dpi)
ggsave(file.path(outdir, "16_FINAL_age_boxplot_by_period.pdf"), p_age, width = plot_width, height = plot_height)

p_ga <- ggplot(sample_final, aes(x = group2, y = collectGA, fill = group2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.12, alpha = 0.45, size = 1) +
  facet_wrap(~ period_new, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(title = "Final collectGA distribution after selection", x = "Group", y = "collectGA, weeks") +
  theme(legend.position = "none")

ggsave(file.path(outdir, "17_FINAL_collectGA_boxplot_by_period.png"), p_ga, width = plot_width, height = plot_height, dpi = plot_dpi)
ggsave(file.path(outdir, "17_FINAL_collectGA_boxplot_by_period.pdf"), p_ga, width = plot_width, height = plot_height)

p_distance <- ggplot(sample_final, aes(x = group2, y = distance_to_opposite_center, fill = group2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.12, alpha = 0.45, size = 1) +
  facet_wrap(~ period_new, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    title = "Covariate distance to opposite-group center",
    subtitle = "Lower value means closer to the opposite group in age and collectGA",
    x = "Group",
    y = "Distance to opposite-group center"
  ) +
  theme(legend.position = "none")

ggsave(file.path(outdir, "18_FINAL_covariate_distance_to_opposite_center.png"), p_distance, width = plot_width, height = plot_height, dpi = plot_dpi)
ggsave(file.path(outdir, "18_FINAL_covariate_distance_to_opposite_center.pdf"), p_distance, width = plot_width, height = plot_height)


# ============================================================
# 9. 读取表达矩阵并按最终样本表过滤
# ============================================================

cat("\n开始读取表达矩阵...\n")

expr_raw <- readr::read_delim(
  expr_file,
  delim = expr_sep,
  col_types = cols(),
  progress = FALSE
)

if (!gene_col %in% colnames(expr_raw)) {
  warning("表达矩阵中没有找到 gene_col = ", gene_col, "，将默认第一列作为 gene 列。")
  colnames(expr_raw)[1] <- gene_col
}

expr_sample_cols <- setdiff(colnames(expr_raw), gene_col)

if (anyDuplicated(expr_sample_cols) > 0) {
  dup_expr_sample <- data.frame(sample_id = expr_sample_cols[duplicated(expr_sample_cols)])
  write_csv_utf8(dup_expr_sample, file.path(outdir, "18_duplicate_sample_ids_in_expression_matrix.csv"))
  stop("表达矩阵中存在重复样本列名，请查看 18_duplicate_sample_ids_in_expression_matrix.csv")
}

gene_dup_n <- sum(duplicated(expr_raw[[gene_col]]))
gene_check <- tibble(
  check_item = c("gene_col", "n_gene", "n_duplicate_gene_id"),
  value = c(gene_col, as.character(nrow(expr_raw)), as.character(gene_dup_n))
)
write_csv_utf8(gene_check, file.path(outdir, "18_expression_gene_id_check.csv"))

final_sample_ids <- sample_final$sample_id
samples_in_both <- intersect(final_sample_ids, expr_sample_cols)
samples_missing_in_expr <- setdiff(final_sample_ids, expr_sample_cols)
extra_expr_samples <- setdiff(expr_sample_cols, final_sample_ids)

match_report <- tibble(
  item = c(
    "final_sample_table_n",
    "expr_matrix_sample_n",
    "samples_in_both_n",
    "samples_missing_in_expr_n",
    "extra_expr_samples_n"
  ),
  n = c(
    length(final_sample_ids),
    length(expr_sample_cols),
    length(samples_in_both),
    length(samples_missing_in_expr),
    length(extra_expr_samples)
  )
)

write_csv_utf8(match_report, file.path(outdir, "18_matrix_sample_match_report.csv"))
write_csv_utf8(data.frame(sample_id = samples_missing_in_expr), file.path(outdir, "18_samples_in_final_table_but_missing_in_expr.csv"))
write_csv_utf8(data.frame(sample_id = extra_expr_samples), file.path(outdir, "18_samples_in_expr_but_not_in_final_table.csv"))

if (length(samples_missing_in_expr) > 0) {
  warning(
    "有最终样本表中的样本不在表达矩阵中。请查看 18_samples_in_final_table_but_missing_in_expr.csv。\n",
    "后续表达矩阵只能保留二者交集样本。"
  )
}

# 按最终样本表顺序排列表达矩阵列
final_sample_ids_in_expr <- final_sample_ids[final_sample_ids %in% expr_sample_cols]

expr_final <- expr_raw %>%
  select(all_of(gene_col), all_of(final_sample_ids_in_expr))

# 同步过滤样本表，只保留表达矩阵中存在的样本，并按表达矩阵列顺序排列
sample_final_matched <- sample_final %>%
  filter(sample_id %in% final_sample_ids_in_expr) %>%
  mutate(sample_id = factor(sample_id, levels = final_sample_ids_in_expr)) %>%
  arrange(sample_id) %>%
  mutate(sample_id = as.character(sample_id))

check_order <- identical(sample_final_matched$sample_id, colnames(expr_final)[-1])

check_df <- tibble(
  check_item = c(
    "sample_order_identical",
    "n_sample_table_matched",
    "n_expr_matrix_matched",
    "n_gene",
    "n_duplicate_gene_id"
  ),
  value = c(
    as.character(check_order),
    as.character(nrow(sample_final_matched)),
    as.character(ncol(expr_final) - 1),
    as.character(nrow(expr_final)),
    as.character(gene_dup_n)
  )
)

write_csv_utf8(check_df, file.path(outdir, "18_FINAL_matrix_order_check.csv"))

if (!check_order) {
  stop("最终样本表和表达矩阵列顺序不一致，请检查。")
}

# 矩阵匹配后重新计算最终 balance，这是后续 DEG 真正使用的样本平衡情况
balance_matched_long <- calc_balance_long(sample_final_matched, "final_matched_to_matrix")
balance_matched_wide <- calc_balance_wide(sample_final_matched, "final_matched_to_matrix")

write_csv_utf8(balance_matched_long, file.path(outdir, "21_FINAL_balance_after_matrix_match_long.csv"))
write_csv_utf8(balance_matched_wide, file.path(outdir, "21_FINAL_balance_after_matrix_match_wide.csv"))

write.table(
  expr_final,
  file = file.path(outdir, "19_FINAL_expression_matrix_matched.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write_csv_utf8(
  sample_final_matched,
  file.path(outdir, "20_FINAL_sample_table_matched_to_matrix.csv")
)

cat("\n表达矩阵过滤完成。\n")
cat("最终表达矩阵：", file.path(outdir, "19_FINAL_expression_matrix_matched.tsv"), "\n")
cat("最终匹配样本表：", file.path(outdir, "20_FINAL_sample_table_matched_to_matrix.csv"), "\n")
cat("顺序是否一致：", check_order, "\n")


# ============================================================
# 10. 输出 Excel 汇总，可选
# ============================================================

if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(
    list(
      parameters = tibble(
        parameter = c(
          "A_definition", "B_definition", "C_definition", "smd_cutoff", "min_control_ratio",
          "gdm_duplicate_strategy", "covariates_to_balance", "covariate_weights"
        ),
        value = c(
          paste0("collectGA <= ", A_upper),
          paste0(A_upper, " < collectGA <= ", B_upper),
          paste0("collectGA > ", B_upper),
          as.character(smd_cutoff),
          as.character(min_control_ratio),
          gdm_duplicate_strategy,
          paste(covariates_to_balance, collapse = ", "),
          paste(names(covariate_weights), covariate_weights, sep = "=", collapse = ", ")
        )
      ),
      sample_clean = sample_clean,
      sample_with_period = sample_period,
      balance_before_long = balance_before_long,
      balance_before_wide = balance_before_wide,
      family_duplicates_before = family_dup_before,
      GDM_duplicates_before = gdm_dup_before,
      sample_after_family_dedup = sample_dedup,
      balance_after_dedup_wide = balance_after_dedup_wide,
      FINAL_sample_balanced = sample_final,
      FINAL_balance_long = balance_final_long,
      FINAL_balance_wide = balance_final_wide,
      FINAL_summary = final_summary,
      removed_samples = removed_all,
      selection_iteration_log = selection_iteration_log,
      matrix_match_report = match_report,
      matrix_order_check = check_df,
      FINAL_sample_matched_to_matrix = sample_final_matched,
      FINAL_balance_after_matrix_match = balance_matched_wide
    ),
    path = file.path(outdir, "FINAL_all_results_summary_optimized.xlsx")
  )
  cat("Excel 汇总文件：", file.path(outdir, "FINAL_all_results_summary_optimized.xlsx"), "\n")
} else {
  cat("未安装 writexl，跳过 Excel 输出。CSV/TSV 结果已完整输出。\n")
}


# ============================================================
# 11. 最终提示
# ============================================================

cat("\n============================================================\n")
cat("全部完成！重点查看以下文件：\n\n")
cat("【后续 DEG 真正使用】表达矩阵：\n")
cat(file.path(outdir, "19_FINAL_expression_matrix_matched.tsv"), "\n\n")
cat("【后续 DEG 真正使用】样本表：\n")
cat(file.path(outdir, "20_FINAL_sample_table_matched_to_matrix.csv"), "\n\n")
cat("【最重要】矩阵匹配后的最终平衡指标：\n")
cat(file.path(outdir, "21_FINAL_balance_after_matrix_match_wide.csv"), "\n")
cat(file.path(outdir, "21_FINAL_balance_after_matrix_match_long.csv"), "\n\n")
cat("【PPT 可用图】SMD 前后变化：\n")
cat(file.path(outdir, "14_SMD_before_after.png"), "\n\n")
cat("【PPT 可用图】最终 age / collectGA 分布：\n")
cat(file.path(outdir, "16_FINAL_age_boxplot_by_period.png"), "\n")
cat(file.path(outdir, "17_FINAL_collectGA_boxplot_by_period.png"), "\n\n")
cat("【删除样本原因】\n")
cat(file.path(outdir, "11_removed_samples_all_reasons.csv"), "\n")
cat(file.path(outdir, "12_balance_selection_iteration_log.csv"), "\n")
cat("============================================================\n")
