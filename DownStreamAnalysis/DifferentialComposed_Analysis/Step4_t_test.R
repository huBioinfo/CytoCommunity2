library(dplyr)
library(stringr)

# 解析每个字符串里的所有数字为一个数向量
parse_nums <- function(x) {
  v <- str_extract_all(as.character(x), "-?\\d*\\.?\\d+(?:[eE][+-]?\\d+)?")
  as.numeric(unlist(v))
}

safe_t <- function(x, y, alt = c("less","greater")) {
  alt <- match.arg(alt)
  x <- x[is.finite(x) & !is.na(x)]
  y <- y[is.finite(y) & !is.na(y)]
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  if (sd(x) == 0 || sd(y) == 0)       return(NA_real_)
  t.test(x, y, alternative = alt, var.equal = TRUE)$p.value
}

for (i in 1:10){
  input_file  <- paste0('../../data/EnrichScoreMatrix/t_test/TCN',as.character(i),'.csv')
  output_file <- paste0('../../data/EnrichScoreMatrix/t_test/bias_result/TCN',as.character(i),'.csv')
  
  data <- read.csv(input_file, stringsAsFactors = FALSE)
  data$Condition <- suppressWarnings(as.numeric(data$Condition))
  
  results <- data.frame(CellType = character(),
                        p_value_left = numeric(),
                        p_value_right = numeric(),
                        stringsAsFactors = FALSE)
  
  grouped_data <- data %>% group_by(CellType)
  
  for (cell_type in unique(data$CellType)) {
    current_group <- grouped_data %>% filter(CellType == cell_type)
    
    # 提取两个 Condition 下的样本，逐元素正则抽数值并合并
    t0_vector <- current_group$EnrichSco[current_group$Condition == 0]
    t1_vector <- current_group$EnrichSco[current_group$Condition == 1]
    
    t0 <- parse_nums(t0_vector)
    t1 <- parse_nums(t1_vector)
    
    # 计算 p 值（不足/零方差则返回 NA，不报错）
    p_left  <- safe_t(t0, t1, "less")
    p_right <- safe_t(t0, t1, "greater")
    
    # 仅当两侧 p 值都不是 NA 时保留
    if (!is.na(p_left) && !is.na(p_right)) {
      results <- dplyr::bind_rows(results, data.frame(
        CellType      = cell_type,
        p_value_left  = p_left,
        p_value_right = p_right,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  write.csv(results, output_file, row.names = FALSE)
  message("[SAVE] ", output_file)
}

