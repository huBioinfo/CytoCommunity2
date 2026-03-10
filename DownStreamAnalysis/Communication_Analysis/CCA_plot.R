suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(stringr)
  library(showtext)
})

font_add("Arial", "C:/Windows/Fonts/arial.ttf")
showtext_auto()

LOAD_DIR <- "../../data/Communication/config"
OUT_DIR  <- "../../data/Communication/SpearmanBetweenCNs_plot/CCA_plots"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 只读取 coefficients 输出
csvs <- list.files(
  LOAD_DIR,
  pattern = "^CCA_coefficients_CN[0-9]+_vs_CN[0-9]+\\.csv$",
  full.names = TRUE
)

for (f in csvs) {
  # 提取 CN 编号
  m <- stringr::str_match(basename(f), "^CCA_coefficients_CN([0-9]+)_vs_CN([0-9]+)\\.csv$")
  if (any(is.na(m))) next
  cnA <- m[2]; cnB <- m[3]
  
  # 读取数据
  df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  
  # 相关系数（可选）
  rho1 <- if ("rho1" %in% names(df)) round(df$rho1[1], 2) else NA_real_
  
  # 绘图（坐标轴固定在 [-1, 1]）
  p <- ggplot(df, aes(x = CN_A_coef, y = CN_B_coef)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey80", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey80", linewidth = 0.5) +
    geom_point(color = "#2E86AB", alpha = 0.8, size = 3) +
    geom_text_repel(
      aes(label = CellType),
      size = 3.2,
      max.overlaps = Inf,
      min.segment.length = 0,
      segment.size = 0.25,
      segment.color = "grey50",
      box.padding = 0.35,
      point.padding = 0.3
    ) +
    coord_fixed(ratio = 1, xlim = c(-1, 1), ylim = c(-1, 1)) +  # 固定范围
    theme_classic(base_size = 12, base_family = "Arial") +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      axis.title  = element_text(face = "bold", size = 10),
      axis.text   = element_text(size = 10),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white", color = NA),
      plot.margin = margin(t = 15, r = 10, b = 10, l = 10)
    ) +
    labs(
      title = sprintf("First canonical variate pair\n(correlation coefficient = %.2f)", rho1),
      x = sprintf("CN%s canonical coefficient", cnA),
      y = sprintf("CN%s canonical coefficient", cnB)
    )
  
  out_pdf <- file.path(
    OUT_DIR,
    sub("\\.csv$", "_scatter.pdf", basename(f))
  )
  ggsave(out_pdf, p, width = 6, height = 6, dpi = 1200, bg = "white")
  message("✅ [COEFFICIENTS] CN", cnA, " vs CN", cnB, " 已保存：", out_png)
}

message("\n✅ 典型系数散点图已全部生成")

