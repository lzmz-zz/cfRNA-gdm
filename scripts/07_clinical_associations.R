
# ============================================================
# 不分时期：16个模块与非OGTT临床指标整体关联分析
#
# 1. 普通 Spearman
# 2. 未校正秩模型
# 3. 校正秩模型：age + collectGA + group2
#
# 有重复受试者且数据允许时使用线性混合效应模型；
# 无重复受试者时使用普通线性模型；
# 单个模型失败不会中断全部分析。
# ============================================================

library(readxl)
library(tidyverse)
library(nlme)

# ============================================================
# 1. 输入与输出
# ============================================================

module_file <- file.path("results", "cfRNA", "modules", "Step6_sample_module_scores_with_metadata.csv")

clinical_file <- file.path("data", "controlled", "临床数据表.xlsx")

output_dir <- file.path("results", "clinical", "overall_three_methods")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

expected_n <- 211
min_n <- 20
min_family <- 10
min_repeated_family_for_lmm <- 5

clinical_vars <- c(
  "WBC", "NLR", "SII",
  "Hb", "RDW", "PLT", "MPV", "PDW",
  "HBA1C",
  "ALT", "AST", "ALB", "TBA", "UA",
  "HPL", "UE3", "TSH2", "FT4",
  "FIB-C", "APTT", "INR"
)

# ============================================================
# 2. 辅助函数
# ============================================================

repair_column_names <- function(df) {
  nm <- names(df)

  bad <- is.na(nm) | trimws(nm) == ""
  if (any(bad)) {
    nm[bad] <- paste0("unnamed_", which(bad))
  }

  names(df) <- make.unique(nm, sep = "_")
  df
}

drop_empty_columns <- function(df) {
  keep <- vapply(
    df,
    function(x) {
      any(!is.na(x) & trimws(as.character(x)) != "")
    },
    logical(1)
  )

  df[, keep, drop = FALSE]
}


# 合并Excel中因重复列名产生的 sample...1、sample...17 等列
collapse_sample_columns <- function(df) {
  sample_cols <- names(df)[
    stringr::str_detect(
      names(df),
      stringr::regex(
        "^sample(?:\\.\\.\\.[0-9]+|\\.[0-9]+|_[0-9]+)?$",
        ignore_case = TRUE
      )
    )
  ]

  if (length(sample_cols) == 0) {
    stop(
      "临床表中没有识别到样本编号列。当前列名：\n",
      paste(names(df), collapse = ", ")
    )
  }

  sample_value_df <- as.data.frame(
    lapply(
      df[sample_cols],
      function(x) {
        x <- trimws(as.character(x))
        x[is.na(x) | x == ""] <- NA_character_
        x
      }
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  conflict_rows <- which(
    apply(
      sample_value_df,
      1,
      function(z) {
        length(unique(z[!is.na(z)])) > 1
      }
    )
  )

  if (length(conflict_rows) > 0) {
    stop(
      "多个sample列在以下行不一致：",
      paste(conflict_rows, collapse = ", "),
      "\n识别到的sample列：",
      paste(sample_cols, collapse = ", ")
    )
  }

  sample_value <- apply(
    sample_value_df,
    1,
    function(z) {
      z <- z[!is.na(z)]
      if (length(z) == 0) NA_character_ else z[1]
    }
  )

  other_cols <- setdiff(names(df), sample_cols)

  dplyr::bind_cols(
    tibble::tibble(sample = as.character(sample_value)),
    df[, other_cols, drop = FALSE]
  )
}

to_numeric <- function(x) {
  if (is.numeric(x)) {
    return(as.numeric(x))
  }

  x <- trimws(as.character(x))
  x[x %in% c(
    "", "NA", "N/A", "NaN",
    "NULL", "-", "--", ".", "None"
  )] <- NA_character_

  suppressWarnings(
    readr::parse_number(
      x,
      locale = readr::locale(grouping_mark = ",")
    )
  )
}

normalize_group <- function(x) {
  x <- toupper(trimws(as.character(x)))

  case_when(
    x %in% c("GDM", "CASE", "DISEASE", "1") ~ "GDM",
    x %in% c(
      "HEALTHY", "CTRL", "CONTROL",
      "NORMAL", "0"
    ) ~ "Healthy",
    TRUE ~ NA_character_
  )
}

rank_z <- function(x) {
  r <- rank(
    x,
    ties.method = "average",
    na.last = "keep"
  )

  if (
    sum(!is.na(r)) < 3 ||
    is.na(sd(r, na.rm = TRUE)) ||
    sd(r, na.rm = TRUE) == 0
  ) {
    return(rep(NA_real_, length(x)))
  }

  as.numeric(scale(r))
}

empty_result <- function(
    n = 0,
    n_family = 0,
    n_repeated_family = 0,
    model_type = NA_character_,
    status = "FAILED",
    message = NA_character_) {

  tibble(
    n = n,
    n_family = n_family,
    n_repeated_family = n_repeated_family,
    model_type = model_type,
    effect = NA_real_,
    se = NA_real_,
    df = NA_real_,
    statistic = NA_real_,
    ci_low = NA_real_,
    ci_high = NA_real_,
    p_value = NA_real_,
    status = status,
    message = message
  )
}

safe_fdr <- function(p) {
  out <- rep(NA_real_, length(p))
  ok <- is.finite(p)

  if (any(ok)) {
    out[ok] <- p.adjust(p[ok], method = "BH")
  }

  out
}

# ============================================================
# 3. 读取模块分数
# ============================================================

if (!file.exists(module_file)) {
  stop("模块分数文件不存在：", module_file)
}

module_df <- read.csv(
  module_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

module_df <- module_df %>%
  repair_column_names() %>%
  drop_empty_columns()

if (!("sample" %in% names(module_df))) {
  stop(
    "模块表中没有sample列。当前列名：\n",
    paste(names(module_df), collapse = ", ")
  )
}

module_df$sample <- trimws(as.character(module_df$sample))

module_cols <- grep(
  "^Dynamic_Module_[0-9]+$",
  names(module_df),
  value = TRUE
)

module_cols <- module_cols[
  order(
    as.numeric(
      str_extract(module_cols, "[0-9]+$")
    )
  )
]

if (length(module_cols) != 16) {
  stop(
    "识别到的模块列不是16个，当前为：",
    length(module_cols),
    "\n识别到：",
    paste(module_cols, collapse = ", ")
  )
}

if (any(is.na(module_df$sample) | module_df$sample == "")) {
  stop("模块表中存在空sample")
}

if (anyDuplicated(module_df$sample) > 0) {
  dup <- unique(
    module_df$sample[duplicated(module_df$sample)]
  )

  stop(
    "模块表中存在重复sample：\n",
    paste(dup, collapse = "\n")
  )
}

if (nrow(module_df) != expected_n) {
  stop(
    "模块表样本数不是", expected_n,
    "，当前为：", nrow(module_df)
  )
}

module_df <- module_df %>%
  mutate(
    across(all_of(module_cols), to_numeric)
  )

if (any(
  vapply(
    module_df[module_cols],
    function(x) all(is.na(x)),
    logical(1)
  )
)) {
  stop("至少一个模块列全部为NA，请检查模块表")
}

# ============================================================
# 4. 读取临床数据表
# ============================================================

if (!file.exists(clinical_file)) {
  stop("临床数据表不存在：", clinical_file)
}

clinical_df <- read_excel(
  clinical_file,
  sheet = 1,
  na = c(
    "", "NA", "N/A", "NaN",
    "NULL", "-", "--"
  ),
  .name_repair = "unique"
) %>%
  as.data.frame(check.names = FALSE) %>%
  repair_column_names() %>%
  drop_empty_columns() %>%
  collapse_sample_columns()

cat("临床表维度：", nrow(clinical_df), "×", ncol(clinical_df), "\n")
cat("样本编号列已统一为：sample\n")

required_cols <- c(
  "sample",
  "group2",
  "period_new",
  "age",
  "collectGA",
  "WBC",
  "NEUT%",
  "LY%",
  "Hb",
  "RDW",
  "PLT",
  "MPV",
  "PDW"
)

missing_required <- setdiff(
  required_cols,
  names(clinical_df)
)

if (length(missing_required) > 0) {
  stop(
    "临床表缺少必要列：",
    paste(missing_required, collapse = ", "),
    "\n当前列名：\n",
    paste(names(clinical_df), collapse = ", ")
  )
}

clinical_df <- clinical_df %>%
  mutate(
    sample = trimws(as.character(sample)),
    Family = str_remove(sample, "^PW-PR-[ABC]-"),
    group2 = normalize_group(group2),
    period_new = trimws(as.character(period_new)),
    age = to_numeric(age),
    collectGA = to_numeric(collectGA)
  )

# 原始Notebook包含一条经核实的年龄更正。为避免公开样本标识符，
# 请将相同的更正记录放在受控目录的CSV文件中；不要跳过该更正。
# CSV 必须包含 sample,corrected_age 两列，且更正记录须与原始研究一致。
age_correction_file <- file.path("data", "controlled", "clinical_age_corrections.csv")

if (!file.exists(age_correction_file)) {
  stop(
    "缺少受控年龄更正表：", age_correction_file,
    "\n请从原Notebook已核实的更正记录生成 sample,corrected_age 两列表格，不要上传GitHub。"
  )
}

age_corrections <- read.csv(
  age_correction_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!all(c("sample", "corrected_age") %in% names(age_corrections))) {
  stop("年龄更正表必须包含 sample,corrected_age 两列。")
}

age_corrections <- age_corrections %>%
  transmute(
    sample = trimws(as.character(sample)),
    corrected_age = suppressWarnings(as.numeric(corrected_age))
  )

if (
  nrow(age_corrections) < 1 ||
  any(is.na(age_corrections$sample) | age_corrections$sample == "") ||
  any(!is.finite(age_corrections$corrected_age)) ||
  anyDuplicated(age_corrections$sample) > 0 ||
  !all(age_corrections$sample %in% clinical_df$sample)
) {
  stop("年龄更正表存在无效、重复或未匹配的样本，请核实原始记录。")
}

clinical_df <- clinical_df %>%
  left_join(age_corrections, by = "sample") %>%
  mutate(
    age = if_else(
      !is.na(corrected_age),
      corrected_age,
      age
    )
  ) %>%
  select(-corrected_age)

if (any(is.na(clinical_df$sample) | clinical_df$sample == "")) {
  stop("临床表中存在空sample")
}

if (anyDuplicated(clinical_df$sample) > 0) {
  dup <- unique(
    clinical_df$sample[
      duplicated(clinical_df$sample)
    ]
  )

  stop(
    "临床表中存在重复sample：\n",
    paste(dup, collapse = "\n")
  )
}

if (nrow(clinical_df) != expected_n) {
  stop(
    "临床表样本数不是", expected_n,
    "，当前为：", nrow(clinical_df)
  )
}

if (any(is.na(clinical_df$Family) | clinical_df$Family == "")) {
  stop("部分样本无法从sample中提取Family")
}

if (any(is.na(clinical_df$group2))) {
  bad <- clinical_df %>%
    filter(is.na(group2)) %>%
    pull(sample)

  stop(
    "以下样本的group2无法识别：\n",
    paste(bad, collapse = "\n")
  )
}

if (any(is.na(clinical_df$age))) {
  bad <- clinical_df %>%
    filter(is.na(age)) %>%
    pull(sample)

  stop(
    "以下样本age缺失或无法转为数值：\n",
    paste(bad, collapse = "\n")
  )
}

if (any(is.na(clinical_df$collectGA))) {
  bad <- clinical_df %>%
    filter(is.na(collectGA)) %>%
    pull(sample)

  stop(
    "以下样本collectGA缺失或无法转为数值：\n",
    paste(bad, collapse = "\n")
  )
}

# ============================================================
# 5. 转换临床指标并计算NLR/SII
# ============================================================

direct_numeric_vars <- intersect(
  c(
    "WBC",
    "NEUT%",
    "LY%",
    "Hb",
    "RDW",
    "PLT",
    "MPV",
    "PDW",
    "HBA1C",
    "ALT",
    "AST",
    "ALB",
    "TBA",
    "UA",
    "HPL",
    "UE3",
    "TSH2",
    "FT4",
    "FIB-C",
    "APTT",
    "INR"
  ),
  names(clinical_df)
)

clinical_df <- clinical_df %>%
  mutate(
    across(all_of(direct_numeric_vars), to_numeric)
  ) %>%
  mutate(
    NLR = if_else(
      !is.na(`NEUT%`) &
        !is.na(`LY%`) &
        `LY%` > 0,
      `NEUT%` / `LY%`,
      NA_real_
    ),
    SII = if_else(
      !is.na(PLT) &
        !is.na(NLR),
      PLT * NLR,
      NA_real_
    )
  )

missing_clinical <- setdiff(
  clinical_vars,
  names(clinical_df)
)

if (length(missing_clinical) > 0) {
  stop(
    "以下临床指标不存在：",
    paste(missing_clinical, collapse = ", ")
  )
}

clinical_df <- clinical_df %>%
  mutate(
    across(all_of(clinical_vars), to_numeric)
  )

# ============================================================
# 6. 样本匹配与分组检查
# ============================================================

missing_in_clinical <- setdiff(
  module_df$sample,
  clinical_df$sample
)

extra_in_clinical <- setdiff(
  clinical_df$sample,
  module_df$sample
)

if (length(missing_in_clinical) > 0) {
  stop(
    "以下模块样本不在临床表中：\n",
    paste(missing_in_clinical, collapse = "\n")
  )
}

if (length(extra_in_clinical) > 0) {
  stop(
    "以下临床样本不在模块表中：\n",
    paste(extra_in_clinical, collapse = "\n")
  )
}

clinical_df <- clinical_df[
  match(module_df$sample, clinical_df$sample),
  ,
  drop = FALSE
]

if (!identical(
  module_df$sample,
  clinical_df$sample
)) {
  stop("临床表按模块表排序失败")
}

if ("group2" %in% names(module_df)) {
  module_group <- normalize_group(
    module_df$group2
  )

  bad <- which(
    is.na(module_group) |
      module_group != clinical_df$group2
  )

  if (length(bad) > 0) {
    stop(
      "模块表与临床表的group2不一致，样本：\n",
      paste(module_df$sample[bad], collapse = "\n")
    )
  }
}

if ("period_new" %in% names(module_df)) {
  module_period <- trimws(
    as.character(module_df$period_new)
  )

  bad <- which(
    is.na(module_period) |
      module_period != clinical_df$period_new
  )

  if (length(bad) > 0) {
    stop(
      "模块表与临床表的period_new不一致，样本：\n",
      paste(module_df$sample[bad], collapse = "\n")
    )
  }
}

# ============================================================
# 7. 合并数据
# ============================================================

dat <- module_df %>%
  select(
    sample,
    all_of(module_cols)
  ) %>%
  left_join(
    clinical_df %>%
      select(
        sample,
        Family,
        group2,
        period_new,
        age,
        collectGA,
        all_of(clinical_vars)
      ),
    by = "sample"
  ) %>%
  mutate(
    Family = factor(Family),
    group2 = factor(
      group2,
      levels = c("Healthy", "GDM")
    )
  )

if (nrow(dat) != expected_n) {
  stop(
    "合并后样本数不是", expected_n,
    "，当前为：", nrow(dat)
  )
}

if (anyDuplicated(dat$sample) > 0) {
  stop("合并后出现重复sample")
}

if (any(is.na(dat$Family))) {
  stop("合并后Family存在NA")
}

cat("\n========================================\n")
cat("样本匹配成功\n")
cat("样本数：", nrow(dat), "\n")
cat("受试者数：", n_distinct(dat$Family), "\n")
cat("GDM样本数：", sum(dat$group2 == "GDM"), "\n")
cat("Healthy样本数：", sum(dat$group2 == "Healthy"), "\n")
cat("\n分期与分组：\n")
print(table(dat$period_new, dat$group2))
cat("========================================\n")

# ============================================================
# 8. 临床数据完整性
# ============================================================

availability <- map_dfr(
  clinical_vars,
  function(v) {
    d <- dat %>%
      filter(!is.na(.data[[v]]))

    family_count <- table(d$Family)

    tibble(
      clinical = v,
      n_samples = nrow(d),
      n_family = n_distinct(d$Family),
      n_repeated_family = sum(
        family_count >= 2
      ),
      missing_n = sum(is.na(dat[[v]])),
      missing_percent = round(
        mean(is.na(dat[[v]])) * 100,
        2
      ),
      unique_values = n_distinct(
        d[[v]]
      )
    )
  }
)

write_csv(
  availability,
  file.path(
    output_dir,
    "00_clinical_availability.csv"
  )
)

# ============================================================
# 9. 普通Spearman
# ============================================================

run_spearman <- function(
    data,
    module_name,
    clinical_name) {

  d <- data %>%
    transmute(
      x = .data[[module_name]],
      y = .data[[clinical_name]],
      Family = Family
    ) %>%
    drop_na(x, y)

  n <- nrow(d)
  n_family <- n_distinct(d$Family)
  n_repeated <- sum(
    table(d$Family) >= 2
  )

  if (n < min_n) {
    return(
      empty_result(
        n,
        n_family,
        n_repeated,
        "Spearman",
        "FAILED",
        "有效样本数不足"
      )
    )
  }

  if (
    n_distinct(d$x) < 3 ||
    n_distinct(d$y) < 3
  ) {
    return(
      empty_result(
        n,
        n_family,
        n_repeated,
        "Spearman",
        "FAILED",
        "变量有效取值不足"
      )
    )
  }

  test <- tryCatch(
    suppressWarnings(
      cor.test(
        d$x,
        d$y,
        method = "spearman",
        exact = FALSE
      )
    ),
    error = function(e) e
  )

  if (inherits(test, "error")) {
    return(
      empty_result(
        n,
        n_family,
        n_repeated,
        "Spearman",
        "FAILED",
        conditionMessage(test)
      )
    )
  }

  tibble(
    n = n,
    n_family = n_family,
    n_repeated_family = n_repeated,
    model_type = "Spearman",
    effect = unname(test$estimate),
    se = NA_real_,
    df = NA_real_,
    statistic = unname(test$statistic),
    ci_low = NA_real_,
    ci_high = NA_real_,
    p_value = test$p.value,
    status = "OK",
    message = ""
  )
}

# ============================================================
# 10. 秩模型
#
# 有足够重复受试者：LMM
# 无重复或重复不足：LM
# ============================================================

run_rank_model <- function(
    data,
    module_name,
    clinical_name,
    adjusted = FALSE) {

  if (!adjusted) {
    d <- data %>%
      transmute(
        x = .data[[module_name]],
        y = .data[[clinical_name]],
        Family = Family
      ) %>%
      drop_na()
  } else {
    d <- data %>%
      transmute(
        x = .data[[module_name]],
        y = .data[[clinical_name]],
        Family = Family,
        age = age,
        collectGA = collectGA,
        group2 = group2
      ) %>%
      drop_na() %>%
      droplevels()
  }

  n <- nrow(d)
  n_family <- n_distinct(d$Family)
  n_repeated <- sum(
    table(d$Family) >= 2
  )

  model_label <- ifelse(
    adjusted,
    "Adjusted",
    "Crude"
  )

  if (n < min_n) {
    return(
      empty_result(
        n,
        n_family,
        n_repeated,
        model_label,
        "FAILED",
        "有效样本数不足"
      )
    )
  }

  if (n_family < min_family) {
    return(
      empty_result(
        n,
        n_family,
        n_repeated,
        model_label,
        "FAILED",
        "有效受试者数不足"
      )
    )
  }

  if (
    n_distinct(d$x) < 3 ||
    n_distinct(d$y) < 3
  ) {
    return(
      empty_result(
        n,
        n_family,
        n_repeated,
        model_label,
        "FAILED",
        "变量有效取值不足"
      )
    )
  }

  if (
    adjusted &&
    n_distinct(d$group2) < 2
  ) {
    return(
      empty_result(
        n,
        n_family,
        n_repeated,
        model_label,
        "FAILED",
        "有效数据中只有一个疾病分组"
      )
    )
  }

  d <- d %>%
    mutate(
      x_rank = rank_z(x),
      y_rank = rank_z(y)
    )

  if (
    any(is.na(d$x_rank)) ||
    any(is.na(d$y_rank))
  ) {
    return(
      empty_result(
        n,
        n_family,
        n_repeated,
        model_label,
        "FAILED",
        "秩转换失败"
      )
    )
  }

  if (adjusted) {
    if (
      is.na(sd(d$age)) ||
      sd(d$age) == 0 ||
      is.na(sd(d$collectGA)) ||
      sd(d$collectGA) == 0
    ) {
      return(
        empty_result(
          n,
          n_family,
          n_repeated,
          model_label,
          "FAILED",
          "年龄或采样孕周没有变异"
        )
      )
    }

    d <- d %>%
      mutate(
        age_z = as.numeric(scale(age)),
        collectGA_z = as.numeric(
          scale(collectGA)
        )
      )
  }

  use_lmm <- n_repeated >=
    min_repeated_family_for_lmm

  if (use_lmm) {
    model_type <- ifelse(
      adjusted,
      "Adjusted_LMM",
      "Crude_LMM"
    )

    fit <- tryCatch(
      {
        if (!adjusted) {
          nlme::lme(
            fixed = y_rank ~ x_rank,
            random = ~ 1 | Family,
            data = d,
            method = "REML",
            na.action = na.omit,
            control = nlme::lmeControl(
              opt = "optim",
              msMaxIter = 200,
              msMaxEval = 400,
              niterEM = 50,
              tolerance = 1e-6,
              returnObject = FALSE
            )
          )
        } else {
          nlme::lme(
            fixed = y_rank ~
              x_rank +
              age_z +
              collectGA_z +
              group2,
            random = ~ 1 | Family,
            data = d,
            method = "REML",
            na.action = na.omit,
            control = nlme::lmeControl(
              opt = "optim",
              msMaxIter = 200,
              msMaxEval = 400,
              niterEM = 50,
              tolerance = 1e-6,
              returnObject = FALSE
            )
          )
        }
      },
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      return(
        empty_result(
          n,
          n_family,
          n_repeated,
          model_type,
          "FAILED",
          conditionMessage(fit)
        )
      )
    }

    tt <- summary(fit)$tTable

    if (!("x_rank" %in% rownames(tt))) {
      return(
        empty_result(
          n,
          n_family,
          n_repeated,
          model_type,
          "FAILED",
          "模型结果中未找到x_rank"
        )
      )
    }

    estimate <- unname(
      tt["x_rank", "Value"]
    )
    standard_error <- unname(
      tt["x_rank", "Std.Error"]
    )
    df_value <- unname(
      tt["x_rank", "DF"]
    )
    statistic <- unname(
      tt["x_rank", "t-value"]
    )
    p_value <- unname(
      tt["x_rank", "p-value"]
    )
  } else {
    model_type <- ifelse(
      adjusted,
      "Adjusted_LM",
      "Crude_LM"
    )

    fit <- tryCatch(
      {
        if (!adjusted) {
          lm(
            y_rank ~ x_rank,
            data = d
          )
        } else {
          lm(
            y_rank ~
              x_rank +
              age_z +
              collectGA_z +
              group2,
            data = d
          )
        }
      },
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      return(
        empty_result(
          n,
          n_family,
          n_repeated,
          model_type,
          "FAILED",
          conditionMessage(fit)
        )
      )
    }

    tt <- summary(fit)$coefficients

    if (!("x_rank" %in% rownames(tt))) {
      return(
        empty_result(
          n,
          n_family,
          n_repeated,
          model_type,
          "FAILED",
          "模型结果中未找到x_rank"
        )
      )
    }

    estimate <- unname(
      tt["x_rank", "Estimate"]
    )
    standard_error <- unname(
      tt["x_rank", "Std. Error"]
    )
    df_value <- df.residual(fit)
    statistic <- unname(
      tt["x_rank", "t value"]
    )
    p_value <- unname(
      tt["x_rank", "Pr(>|t|)"]
    )
  }

  if (
    !is.finite(df_value) ||
    df_value <= 0
  ) {
    critical_value <- qnorm(0.975)
  } else {
    critical_value <- qt(
      0.975,
      df = df_value
    )
  }

  tibble(
    n = n,
    n_family = n_family,
    n_repeated_family = n_repeated,
    model_type = model_type,
    effect = estimate,
    se = standard_error,
    df = df_value,
    statistic = statistic,
    ci_low = estimate -
      critical_value * standard_error,
    ci_high = estimate +
      critical_value * standard_error,
    p_value = p_value,
    status = "OK",
    message = ""
  )
}

# ============================================================
# 11. 批量运行
# ============================================================

analysis_grid <- expand_grid(
  module_col = module_cols,
  clinical = clinical_vars
)

result_spearman <- analysis_grid %>%
  mutate(
    result = map2(
      module_col,
      clinical,
      ~ run_spearman(
        dat,
        .x,
        .y
      )
    )
  ) %>%
  unnest(result) %>%
  mutate(Method = "Spearman")

result_crude <- analysis_grid %>%
  mutate(
    result = map2(
      module_col,
      clinical,
      ~ run_rank_model(
        dat,
        .x,
        .y,
        adjusted = FALSE
      )
    )
  ) %>%
  unnest(result) %>%
  mutate(Method = "Crude")

result_adjusted <- analysis_grid %>%
  mutate(
    result = map2(
      module_col,
      clinical,
      ~ run_rank_model(
        dat,
        .x,
        .y,
        adjusted = TRUE
      )
    )
  ) %>%
  unnest(result) %>%
  mutate(Method = "Adjusted")

# ============================================================
# 12. 合并结果与FDR
# ============================================================

all_results <- bind_rows(
  result_spearman,
  result_crude,
  result_adjusted
) %>%
  mutate(
    module_number = as.numeric(
      str_extract(
        module_col,
        "[0-9]+$"
      )
    ),
    Module = paste0(
      "M",
      sprintf(
        "%02d",
        module_number
      )
    )
  ) %>%
  group_by(Method) %>%
  mutate(
    FDR = safe_fdr(p_value)
  ) %>%
  ungroup() %>%
  mutate(
    direction = case_when(
      effect > 0 ~ "Positive",
      effect < 0 ~ "Negative",
      TRUE ~ NA_character_
    ),
    significance = case_when(
      FDR < 0.001 ~ "FDR<0.001",
      FDR < 0.01 ~ "FDR<0.01",
      FDR < 0.05 ~ "FDR<0.05",
      FDR < 0.10 ~ "FDR<0.10",
      TRUE ~ "NS"
    )
  ) %>%
  select(
    Method,
    model_type,
    Module,
    module_col,
    clinical,
    n,
    n_family,
    n_repeated_family,
    effect,
    se,
    df,
    statistic,
    ci_low,
    ci_high,
    p_value,
    FDR,
    direction,
    significance,
    status,
    message
  ) %>%
  arrange(
    Method,
    clinical,
    Module
  )

write_csv(
  all_results,
  file.path(
    output_dir,
    "01_three_methods_all_results.csv"
  )
)

# ============================================================
# 13. 横向比较结果
# ============================================================

comparison <- all_results %>%
  select(
    Module,
    clinical,
    Method,
    model_type,
    n,
    n_family,
    n_repeated_family,
    effect,
    p_value,
    FDR,
    status
  ) %>%
  pivot_wider(
    names_from = Method,
    values_from = c(
      model_type,
      n,
      n_family,
      n_repeated_family,
      effect,
      p_value,
      FDR,
      status
    ),
    names_glue = "{Method}_{.value}"
  ) %>%
  mutate(
    direction_consistent = case_when(
      is.na(Spearman_effect) |
        is.na(Crude_effect) |
        is.na(Adjusted_effect) ~ NA,
      sign(Spearman_effect) ==
        sign(Crude_effect) &
        sign(Crude_effect) ==
        sign(Adjusted_effect) ~ TRUE,
      TRUE ~ FALSE
    ),
    robust_association = case_when(
      direction_consistent %in% TRUE &
        Adjusted_FDR < 0.05 ~
        "High confidence",
      direction_consistent %in% TRUE &
        Adjusted_FDR < 0.10 ~
        "Trend",
      TRUE ~
        "Not robust"
    )
  )

write_csv(
  comparison,
  file.path(
    output_dir,
    "02_three_methods_comparison.csv"
  )
)

write_csv(
  all_results %>%
    filter(!is.na(FDR), FDR < 0.05) %>%
    arrange(Method, FDR),
  file.path(
    output_dir,
    "03_significant_FDR0.05.csv"
  )
)

write_csv(
  all_results %>%
    filter(!is.na(FDR), FDR < 0.10) %>%
    arrange(Method, FDR),
  file.path(
    output_dir,
    "04_trend_FDR0.10.csv"
  )
)

write_csv(
  comparison %>%
    filter(
      robust_association !=
        "Not robust"
    ) %>%
    arrange(Adjusted_FDR),
  file.path(
    output_dir,
    "05_robust_associations.csv"
  )
)

# ============================================================
# 14. 热图
# ============================================================

module_order <- sprintf(
  "M%02d",
  1:16
)

method_order <- c(
  "Spearman",
  "Crude",
  "Adjusted"
)

method_labels <- c(
  Spearman = "Ordinary Spearman",
  Crude = "Crude rank model",
  Adjusted = "Adjusted rank model"
)

plot_data <- all_results %>%
  mutate(
    Method = factor(
      Method,
      levels = method_order
    ),
    Module = factor(
      Module,
      levels = rev(module_order)
    ),
    clinical = factor(
      clinical,
      levels = clinical_vars
    ),
    star = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      FDR < 0.10 ~ "·",
      TRUE ~ ""
    ),
    label = if_else(
      is.na(effect),
      "",
      paste0(
        sprintf("%.2f", effect),
        star
      )
    )
  )

finite_effects <- abs(
  plot_data$effect[
    is.finite(plot_data$effect)
  ]
)

if (length(finite_effects) == 0) {
  effect_limit <- 1
} else {
  effect_limit <- max(
    0.30,
    as.numeric(
      quantile(
        finite_effects,
        probs = 0.98,
        na.rm = TRUE
      )
    )
  )
}

for (method_i in method_order) {
  p <- plot_data %>%
    filter(Method == method_i) %>%
    ggplot(
      aes(
        x = clinical,
        y = Module,
        fill = effect
      )
    ) +
    geom_tile(
      color = "white",
      linewidth = 0.4
    ) +
    geom_text(
      aes(label = label),
      size = 2.5
    ) +
    scale_fill_gradient2(
      low = "#3B4CC0",
      mid = "white",
      high = "#B40426",
      midpoint = 0,
      limits = c(
        -effect_limit,
        effect_limit
      ),
      oob = scales::squish,
      na.value = "grey90"
    ) +
    labs(
      title = method_labels[[method_i]],
      subtitle = paste0(
        "*** FDR<0.001; ** FDR<0.01; ",
        "* FDR<0.05; · FDR<0.10"
      ),
      x = NULL,
      y = NULL,
      fill = ifelse(
        method_i == "Spearman",
        "Spearman rho",
        "Rank beta"
      )
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 9
      )
    )

  ggsave(
    filename = file.path(
      output_dir,
      paste0(
        "06_",
        method_i,
        "_heatmap.pdf"
      )
    ),
    plot = p,
    width = 14,
    height = 9
  )
}

p_combined <- ggplot(
  plot_data,
  aes(
    x = clinical,
    y = Module,
    fill = effect
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.35
  ) +
  geom_text(
    aes(label = label),
    size = 2.1
  ) +
  facet_grid(
    rows = vars(Method),
    labeller = as_labeller(
      method_labels
    )
  ) +
  scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0,
    limits = c(
      -effect_limit,
      effect_limit
    ),
    oob = scales::squish,
    na.value = "grey90"
  ) +
  labs(
    title = "Overall module-clinical associations",
    subtitle = paste0(
      "*** FDR<0.001; ** FDR<0.01; ",
      "* FDR<0.05; · FDR<0.10"
    ),
    x = NULL,
    y = NULL,
    fill = "Effect"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    strip.text.y = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    )
  )

ggsave(
  filename = file.path(
    output_dir,
    "07_three_methods_combined_heatmap.pdf"
  ),
  plot = p_combined,
  width = 14,
  height = 23
)

# ============================================================
# 15. 模型运行状态
# ============================================================

status_summary <- all_results %>%
  count(
    Method,
    model_type,
    status,
    name = "n_tests"
  )

write_csv(
  status_summary,
  file.path(
    output_dir,
    "08_model_status_summary.csv"
  )
)

print(status_summary)

failed_models <- all_results %>%
  filter(status != "OK")

write_csv(
  failed_models,
  file.path(
    output_dir,
    "09_failed_models.csv"
  )
)

cat("\n分析完成\n")
cat("样本数：", nrow(dat), "\n")
cat("受试者数：", n_distinct(dat$Family), "\n")
cat("模块数：", length(module_cols), "\n")
cat("临床指标数：", length(clinical_vars), "\n")
cat("输出目录：", output_dir, "\n")
