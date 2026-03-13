suppressPackageStartupMessages({
  library(reshape2)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
})

# ---------------------- 路径配置（按你的目录） ----------------------
score_folder <- "./data/EnrichScoreMatrix/Score"                               # 原始长表目录
final_folder <- "./data/EnrichScoreMatrix/t_test/bias_result/final"      # Python筛选后的final目录

# ---------------------- 小工具：清洗与映射 ----------------------
# 去空格并转字符
norm_str <- function(x) trimws(as.character(x))

# 统一 Condition -> t0/t1/t2
to_t012 <- function(v) {
  v <- norm_str(v)
  v[v %in% c("0","1","2")] <- paste0("t", v[v %in% c("0","1","2")])
  v
}

# 从 final 读取某个CN的显著 CellType（并按最小p排序）
final_celltypes_for_cn <- function(i, final_dir) {
  cand_files <- sprintf("TCN%d_%s.csv", i, c("t0-t1","t0-t2","t1-t2"))
  cand_paths <- file.path(final_dir, cand_files)
  cand_paths <- cand_paths[file.exists(cand_paths)]
  if (length(cand_paths) == 0) return(character(0))

  dfs <- lapply(cand_paths, function(p){
    df <- tryCatch(read.csv(p, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    if (is.null(df)) return(NULL)
    names(df) <- norm_str(names(df))
    if (!all(c("CellType","p_value") %in% names(df))) return(NULL)
    df$CellType <- norm_str(df$CellType)
    df$p_value  <- as.numeric(df$p_value)
    df
  })
  dfs <- Filter(Negate(is.null), dfs)
  if (length(dfs) == 0) return(character(0))

  all_df <- do.call(rbind, dfs)
  all_df <- all_df[is.finite(all_df$p_value), , drop = FALSE]
  if (nrow(all_df) == 0) return(character(0))

  agg <- aggregate(p_value ~ CellType, data = all_df, FUN = min, na.rm = TRUE)
  agg <- agg[order(agg$p_value, decreasing = FALSE), ]
  agg$CellType
}

# ---------------------- 生成图形 ----------------------
plots <- list()
plot_idx <- 0

# 检查目录
if (!dir.exists(score_folder)) stop("找不到 score_folder: ", score_folder)
if (!dir.exists(final_folder)) stop("找不到 final_folder: ", final_folder)

for (i in 1:10) {
  # 1) 取该 CN 在 final 中的显著细胞类型（按p从小到大）
  keep_ct <- final_celltypes_for_cn(i, final_folder)
  message(sprintf("[CN%d] final中显著细胞类型数: %d", i, length(keep_ct)))
  if (length(keep_ct) == 0) {
    message(sprintf("  -> 跳过 CN%d（final 中无显著 CellType 或文件未找到）。", i))
    next
  }

  # 2) 读取该 CN 的原始长表
  fn <- sprintf("TCN%d_EnrichmentScore.csv", i)
  fp <- file.path(score_folder, fn)
  if (!file.exists(fp)) {
    message(sprintf("  -> 缺少原始文件：%s，跳过。", fp))
    next
  }
  df <- read.csv(fp, stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- norm_str(names(df))
  # 关键列校验
  need_cols <- c("CellType","Condition","EnrichSco")
  if (!all(need_cols %in% names(df))) {
    message(sprintf("  -> 原始文件缺列(%s)中的一部分，跳过。", paste(need_cols, collapse=",")))
    next
  }

  # 清洗
  df$CellType  <- norm_str(df$CellType)
  df$Condition <- to_t012(df$Condition)
  # 如果原始就是 t0/t1/t2，则保持；如果有别的值，先过滤
  df <- df %>% filter(Condition %in% c("t0","t1","t2"))

  # 只保留出现在 final 的细胞类型
  df <- df %>% filter(CellType %in% keep_ct)

  # 如果清洗后一个都匹配不上，提示并跳过
  if (nrow(df) == 0) {
    message(sprintf("  -> CN%d: final中的CellType与原始文件不匹配，可能是大小写/空格差异；请检查。", i))
    next
  }

  # x轴顺序：按显著性（final中的顺序）
  df$CellType <- factor(df$CellType, levels = keep_ct)

  # 确保 EnrichSco 数值
  df$EnrichSco <- suppressWarnings(as.numeric(df$EnrichSco))

  p <- ggplot(df, aes(
    x      = CellType,
    y      = EnrichSco,
    fill   = Condition,
    group  = Condition
  )) +
    labs(
      x    = paste0("CN", i),
      y    = "EnrichSco",
      fill = NULL
    ) +
    geom_jitter(
      position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8),
      shape    = 21,
      colour   = "black",
      stroke   = 0.3,
      aes(fill = Condition),
      size     = 2,
      alpha    = 0.8,
      show.legend = FALSE
    ) +
    stat_summary(
      aes(group = Condition),
      fun         = mean,
      geom        = "crossbar",
      position    = position_dodge(width = 0.8),
      width       = 0.6,
      size        = 0.3,
      color       = "black",
      show.legend = FALSE
    ) +
    scale_fill_manual(
      values = c("t0" = "#E41A1C", "t1" = "#4DAF4A", "t2" = "#377EB8")
    ) +
    scale_x_discrete(
      limits = keep_ct,
      breaks = keep_ct,
      labels = keep_ct,
      drop   = FALSE,
      expand = c(0, 0)
    ) +
    # y 轴根据数据自适应一个合理范围（更稳）
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
    theme_bw() +
    theme(
      axis.line        = element_line(colour = "black"),
      panel.grid       = element_blank(),
      panel.border     = element_blank(),
      panel.background = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.text.x      = element_text(size = 10, angle = 90, hjust = 1),
      axis.title.x     = element_text(size = 12),
      axis.text.y      = element_text(size = 10),
      legend.position  = "none"
    )

  plot_idx <- plot_idx + 1
  plots[[plot_idx]] <- p
  names(plots)[plot_idx] <- paste0("CN", i)
}


aligned <- align_plots(plotlist = plots, align = "h", axis = "l")

outdir <- "./plot/EnrichScoreMatrix/ScorePlots"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

for (i in seq_along(aligned)) {
  img_name <- paste0(names(plots)[i], ".pdf") 
  ggsave(
    filename = file.path(outdir, img_name),
    plot     = ggdraw(aligned[[i]]),
    width    = 6,
    height   = 6,
    dpi      = 1200
  )
}

message("完成：图像已输出到 ", outdir)

