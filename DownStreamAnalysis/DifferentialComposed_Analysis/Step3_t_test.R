suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

# 解析 "(...,...)" 字符串为数值向量
parse_enrich <- function(s) {
  s <- gsub("[()]", "", s)
  x <- strsplit(s, ",\\s*")[[1]]
  as.numeric(x)
}

# 对某个 CellType 的两组向量做 t 检验，返回一行数据框
do_one_ttest <- function(cell_type, a, b) {
  if (length(a) >= 2 && length(b) >= 2 && all(is.finite(c(a, b)))) {
    t_left  <- tryCatch(t.test(a, b, alternative = "less",    var.equal = TRUE)$p.value, error = function(e) NA_real_)
    t_right <- tryCatch(t.test(a, b, alternative = "greater", var.equal = TRUE)$p.value, error = function(e) NA_real_)
    data.frame(CellType = cell_type,
               p_value_left = t_left,
               p_value_right = t_right,
               stringsAsFactors = FALSE)
  } else {
    NULL
  }
}

# 三个对比与映射
comparisons <- list(
  `_t0-t1` = c("t0", "t1"),
  `_t0-t2` = c("t0", "t2"),
  `_t1-t2` = c("t1", "t2")
)

for (i in 1:10) {
  root <- "./data/EnrichScoreMatrix/t_test/"
  input_file <- paste0(root, "TCN", i, ".csv")
  
  if (!file.exists(input_file)) {
    warning("Missing: ", input_file)
    next
  }
  
  # 读入，并把 0/1/2 标准化为 "t0/t1/t2"
  data_raw <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
  data <- data_raw %>%
    mutate(
      Condition = ifelse(Condition %in% c(0, 1, 2),
                         paste0("t", Condition),
                         as.character(Condition))
    ) %>%
    filter(Condition %in% c("t0", "t1", "t2"))
  
  # 逐个对比分别输出到各自文件夹
  for (obj in names(comparisons)) {
    condA <- comparisons[[obj]][1]
    condB <- comparisons[[obj]][2]
    
    output_dir  <- file.path(root, "bias_result", obj)
    output_file <- file.path(output_dir, paste0("TCN", i, obj, ".csv"))
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # 分 CellType 处理
    results_list <- lapply(split(data, data$CellType), function(df_ct) {
      vecs <- lapply(df_ct$EnrichSco, parse_enrich)
      cond <- df_ct$Condition
      a <- unlist(vecs[cond == condA], use.names = FALSE)
      b <- unlist(vecs[cond == condB], use.names = FALSE)
      do_one_ttest(unique(df_ct$CellType), a, b)
    })
    
    # 丢掉 NULL，若为空则不写文件
    results_list <- results_list[!vapply(results_list, is.null, logical(1))]
    if (length(results_list) == 0) {
      message("No results for TCN", i, " ", obj, " — skipped writing.")
      next
    }
    
    results <- do.call(rbind, results_list)
    # 如果出于稳妥再检查一下行数
    if (is.data.frame(results) && nrow(results) > 0) {
      write.csv(results, output_file, row.names = FALSE)
      message("Wrote: ", output_file)
    } else {
      message("No results for TCN", i, " ", obj, " — skipped writing.")
    }
  }
}

