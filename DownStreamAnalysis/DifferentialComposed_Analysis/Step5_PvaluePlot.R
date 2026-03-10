suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  library(rlang)
  library(showtext)
})

font_add("Arial", "C:/Windows/Fonts/arial.ttf")
showtext_auto()

# ---------------------- 路径配置 ----------------------
score_folder <- "../../data/EnrichScoreMatrix/Score"                         
final_folder <- "../../data/EnrichScoreMatrix/t_test/bias_result/final"      
outdir       <- "../../plot/EnrichScoreMatrix/ScorePlots"

norm_str <- function(x) trimws(as.character(x))

# 从 final 读取某个 TCN 的显著 CellType（按 p_value 从小到大）
final_celltypes_for_tcn <- function(i, final_dir) {
  fp <- file.path(final_dir, sprintf("TCN%d.csv", i))
  if (!file.exists(fp)) return(character(0))
  
  df <- read.csv(fp, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(df) == 0) return(character(0))
  
  df$CellType <- norm_str(df$CellType)
  df$p_value  <- suppressWarnings(as.numeric(df$p_value))
  
  df <- df %>%
    filter(.data$CellType != "", is.finite(.data$p_value)) %>%
    arrange(.data$p_value)
  
  unique(df$CellType)
}

# ---------------------- 生成图形 ----------------------
plots <- list()

for (i in 1:10) {
  keep_ct <- final_celltypes_for_tcn(i, final_folder)

  if (length(keep_ct) == 0) {
    message(sprintf("  -> 跳过 TCN%d（final 文件缺失或无有效 CellType/p_value）。", i))
    next
  }
  
  # 读取原始 EnrichmentScore
  score_fp <- file.path(score_folder, sprintf("TCN%d_EnrichmentScore.csv", i))
  if (!file.exists(score_fp)) {
    message(sprintf("  -> 缺少原始文件：%s，跳过。", score_fp))
    next
  }
  
  df <- read.csv(score_fp, stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- norm_str(names(df))
  
  need_cols <- c("CellType","Condition","EnrichSco")
  
  # 清洗 + 只保留 0/1
  df$CellType  <- norm_str(df$CellType)
  df$Condition <- norm_str(df$Condition)
  df <- df %>% filter(Condition %in% c("0","1"))
  
  # 只保留 final 中的 cell types
  df <- df %>% filter(CellType %in% keep_ct)
  df$Condition <- factor(df$Condition, levels = c("0","1"), labels = c("0","1"))
  df$CellType <- factor(df$CellType, levels = keep_ct)
  df$EnrichSco <- suppressWarnings(as.numeric(df$EnrichSco))
  df <- df %>% filter(is.finite(EnrichSco))
  
  if (nrow(df) == 0) {
    message(sprintf("  -> TCN%d: EnrichSco 全部非数值/异常，跳过。", i))
    next
  }
  
  p <- ggplot(df, aes(x = CellType, y = EnrichSco, fill = Condition, group = Condition)) +
    labs(
      x    = paste0("TCN", i),
      y    = "EnrichSco",
      fill = NULL
    ) +
    geom_jitter(
      position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8),
      shape    = 21,
      colour   = "black",
      stroke   = 0.3,
      size     = 2,
      alpha    = 0.8,
      show.legend = FALSE
    ) +
    stat_summary(
      aes(group = Condition),
      fun      = mean,
      geom     = "crossbar",
      position = position_dodge(width = 0.8),
      width    = 0.6,
      size     = 0.3,
      color    = "black",
      show.legend = FALSE
    ) +
    scale_fill_manual(values = c("0" = "#E41A1C", "1" = "#4DAF4A")) +
    scale_x_discrete(drop = TRUE, expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
    coord_cartesian(ylim = c(0, 20)) +
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
  
  plots[[paste0("TCN", i)]] <- p
}

# ---------------------- 对齐并保存 ----------------------
aligned <- align_plots(plotlist = plots, align = "h", axis = "l")

idx <- 0
for (nm in names(plots)) {
  idx <- idx + 1
  ggsave(
    filename = file.path(outdir, paste0(nm, ".pdf")),
    plot     = ggdraw(aligned[[idx]]),
    width    = 6,
    height   = 7,
    dpi      = 1200
  )
}

message("完成：图像已输出到 ", outdir)
