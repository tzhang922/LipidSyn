# ============================================================
# PCA analysis and visualization
# 数据结构：
# 第一列为脂质名称 Lipids
# 后续各列为不同样本
# ============================================================

# 1. 安装并加载需要的R包
packages <- c("ggplot2", "dplyr")

new_packages <- packages[
  !(packages %in% rownames(installed.packages()))
]

if (length(new_packages) > 0) {
  install.packages(new_packages)
}

library(ggplot2)
library(dplyr)


# ============================================================
# 2. 设置工作目录和输入文件
# ============================================================

# 将CSV文件放在当前工作目录，并修改为实际文件名
input_file <- "PCA_data.csv"

# 读取数据
raw_data <- read.csv(
  input_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# 查看数据结构
print(dim(raw_data))
print(colnames(raw_data))
head(raw_data)


# ============================================================
# 3. 数据检查与清理
# ============================================================

# 第一列作为脂质名称
lipid_names <- raw_data[[1]]

# 其余列为样本数据
numeric_data <- raw_data[, -1, drop = FALSE]

# 转换为数值型，避免字符格式影响PCA
numeric_data[] <- lapply(
  numeric_data,
  function(x) as.numeric(as.character(x))
)

# 检查转换过程中是否产生缺失值
if (anyNA(numeric_data)) {
  warning("数据中存在NA，将使用每个脂质的中位数填补。")
  
  for (i in seq_len(nrow(numeric_data))) {
    
    row_values <- as.numeric(numeric_data[i, ])
    
    if (anyNA(row_values)) {
      
      row_median <- median(
        row_values,
        na.rm = TRUE
      )
      
      # 如果整行均为NA，则填充为0
      if (!is.finite(row_median)) {
        row_median <- 0
      }
      
      row_values[is.na(row_values)] <- row_median
      numeric_data[i, ] <- row_values
    }
  }
}


# ============================================================
# 4. 删除零方差脂质
# ============================================================

# PCA要求变量具有一定变异性
lipid_sd <- apply(
  numeric_data,
  1,
  sd,
  na.rm = TRUE
)

keep_lipids <- is.finite(lipid_sd) & lipid_sd > 0

numeric_data <- numeric_data[keep_lipids, , drop = FALSE]
lipid_names <- lipid_names[keep_lipids]

cat(
  "保留用于PCA的脂质数量：",
  nrow(numeric_data),
  "\n"
)


# ============================================================
# 5. 数据转置
# ============================================================

# 原始数据：
# 行 = 脂质
# 列 = 样本
#
# PCA需要：
# 行 = 样本
# 列 = 脂质

pca_matrix <- t(as.matrix(numeric_data))

# 保留原始样本名
sample_original <- rownames(pca_matrix)

# 检查矩阵
print(dim(pca_matrix))
print(sample_original)


# ============================================================
# 6. PCA分析
# ============================================================

# center = TRUE：中心化
# scale. = TRUE：标准化
#
# 该设置得到：
# PC1约为45.66%
# PC2约为18.41%

pca_result <- prcomp(
  pca_matrix,
  center = TRUE,
  scale. = TRUE
)

# 计算各主成分解释率
variance_percent <- (
  pca_result$sdev^2 /
    sum(pca_result$sdev^2)
) * 100

pc1_percent <- variance_percent[1]
pc2_percent <- variance_percent[2]

cat(
  sprintf(
    "PC1解释率：%.2f%%\nPC2解释率：%.2f%%\n",
    pc1_percent,
    pc2_percent
  )
)


# ============================================================
# 7. 提取样本得分
# ============================================================

pca_scores <- as.data.frame(
  pca_result$x[, 1:2]
)

colnames(pca_scores) <- c("PC1", "PC2")

pca_scores$Original_sample <- rownames(pca_scores)

# 根据样本名前缀区分两组
pca_scores$Group <- ifelse(
  grepl("^OW450", pca_scores$Original_sample),
  "OW450",
  "OW40"
)


# ============================================================
# 8. 按原始数据中的排列顺序重新编号
# ============================================================

# 原始CSV样本编号并不连续：
# OW40：1、2、3、5、6、7
# OW450：2、3、7、1、5、6
#
# 为匹配示例图，分别按出现顺序重新编号为1–6

pca_scores <- pca_scores %>%
  group_by(Group) %>%
  mutate(
    Sample_order = row_number(),
    Sample = paste0(Group, "_", Sample_order)
  ) %>%
  ungroup()


# ============================================================
# 9. 固定PCA坐标方向
# ============================================================

# PCA坐标正负号没有固定生物学含义。
# 为使OW40位于左侧、OW450位于右侧，统一PC1方向。

mean_pc1_ow450 <- mean(
  pca_scores$PC1[pca_scores$Group == "OW450"]
)

mean_pc1_ow40 <- mean(
  pca_scores$PC1[pca_scores$Group == "OW40"]
)

if (mean_pc1_ow450 < mean_pc1_ow40) {
  pca_scores$PC1 <- -pca_scores$PC1
}

# 为使OW450_1位于图形上方，统一PC2方向
ow450_1_pc2 <- pca_scores$PC2[
  pca_scores$Sample == "OW450_1"
]

ow450_2_pc2 <- pca_scores$PC2[
  pca_scores$Sample == "OW450_2"
]

if (
  length(ow450_1_pc2) == 1 &&
  length(ow450_2_pc2) == 1 &&
  ow450_1_pc2 < ow450_2_pc2
) {
  pca_scores$PC2 <- -pca_scores$PC2
}

print(pca_scores)


# ============================================================
# 10. 计算各组凸包
# ============================================================

get_convex_hull <- function(data) {
  
  if (nrow(data) < 3) {
    return(data)
  }
  
  hull_index <- chull(
    data$PC1,
    data$PC2
  )
  
  hull_data <- data[hull_index, , drop = FALSE]
  
  # 闭合多边形
  rbind(
    hull_data,
    hull_data[1, , drop = FALSE]
  )
}

hull_data <- pca_scores %>%
  group_by(Group) %>%
  group_modify(
    ~ get_convex_hull(.x)
  ) %>%
  ungroup()


# ============================================================
# 11. 手动设置标签位置
# ============================================================

# 以下位置按照示例图布局设置。
# 标签位置采用相对于样本点的偏移量。

label_offset <- data.frame(
  Sample = c(
    "OW40_1",
    "OW40_2",
    "OW40_3",
    "OW40_4",
    "OW40_5",
    "OW40_6",
    "OW450_1",
    "OW450_2",
    "OW450_3",
    "OW450_4",
    "OW450_5",
    "OW450_6"
  ),
  
  offset_x = c(
    2.8,
    0.7,
    0.0,
    3.0,
    3.2,
    3.4,
    -0.4,
    0.0,
    -1.8,
    -0.1,
    0.0,
    -0.2
  ),
  
  offset_y = c(
    -2.9,
    4.2,
    -4.8,
    1.1,
    4.2,
    5.4,
    3.0,
    2.0,
    -0.2,
    2.2,
    -2.7,
    -3.5
  )
)

plot_data <- pca_scores %>%
  left_join(
    label_offset,
    by = "Sample"
  ) %>%
  mutate(
    label_x = PC1 + offset_x,
    label_y = PC2 + offset_y
  )


# ============================================================
# 12. 设置配色
# ============================================================

# 颜色接近示例图
point_colors <- c(
  "OW40"  = "#46A9DF",
  "OW450" = "#E69F00"
)

fill_colors <- c(
  "OW40"  = "#A9DAF2",
  "OW450" = "#F3D58B"
)


# ============================================================
# 13. 绘制PCA图
# ============================================================

pca_plot <- ggplot() +
  
  # 分组凸包
  geom_polygon(
    data = hull_data,
    aes(
      x = PC1,
      y = PC2,
      group = Group,
      fill = Group
    ),
    alpha = 0.52,
    color = NA
  ) +
  
  # 标签与样本点之间的虚线
  geom_segment(
    data = plot_data,
    aes(
      x = PC1,
      y = PC2,
      xend = label_x,
      yend = label_y
    ),
    color = "#777777",
    linewidth = 0.45,
    linetype = "dashed"
  ) +
  
  # PCA样本点
  geom_point(
    data = plot_data,
    aes(
      x = PC1,
      y = PC2,
      color = Group
    ),
    size = 3.4,
    alpha = 1
  ) +
  
  # 样本标签
  geom_text(
    data = plot_data,
    aes(
      x = label_x,
      y = label_y,
      label = Sample
    ),
    family = "sans",
    size = 4.8,
    color = "#202020",
    hjust = 0.5,
    vjust = 0.5
  ) +
  
  scale_color_manual(
    values = point_colors
  ) +
  
  scale_fill_manual(
    values = fill_colors
  ) +
  
  # 与示例图相近的坐标刻度
  scale_x_continuous(
    breaks = c(-10, 0, 10),
    expand = expansion(mult = c(0.06, 0.06))
  ) +
  
  scale_y_continuous(
    breaks = c(-10, 0, 10),
    expand = expansion(mult = c(0.06, 0.06))
  ) +
  
  labs(
    x = sprintf(
      "PCA1: %.2f %%",
      pc1_percent
    ),
    y = sprintf(
      "PCA2: %.2f %%",
      pc2_percent
    ),
    tag = "B"
  ) +
  
  coord_cartesian(
    xlim = c(-15, 19),
    ylim = c(-16.5, 16),
    clip = "off"
  ) +
  
  theme_classic(
    base_size = 16,
    base_family = "sans"
  ) +
  
  theme(
    # 不显示图例
    legend.position = "none",
    
    # 坐标轴标题
    axis.title.x = element_text(
      size = 17,
      color = "black",
      margin = margin(t = 8)
    ),
    
    axis.title.y = element_text(
      size = 17,
      color = "black",
      margin = margin(r = 8)
    ),
    
    # 坐标刻度文字
    axis.text = element_text(
      size = 14,
      color = "#555555"
    ),
    
    # 坐标刻度线
    axis.ticks = element_line(
      color = "#666666",
      linewidth = 0.6
    ),
    
    axis.ticks.length = unit(
      0.16,
      "cm"
    ),
    
    # 四周边框
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    
    # B面板编号
    plot.tag = element_text(
      size = 25,
      face = "plain",
      color = "black"
    ),
    
    plot.tag.position = c(
      -0.16,
      1.04
    ),
    
    # 为外部面板编号留出空间
    plot.margin = margin(
      t = 18,
      r = 18,
      b = 12,
      l = 42
    )
  )


# ============================================================
# 14. 显示图形
# ============================================================

print(pca_plot)


# ============================================================
# 15. 保存高分辨率图片
# ============================================================

ggsave(
  filename = "PCA_OW40_OW450.png",
  plot = pca_plot,
  width = 6.3,
  height = 5.6,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "PCA_OW40_OW450.pdf",
  plot = pca_plot,
  width = 6.3,
  height = 5.6,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)


# ============================================================
# 16. 导出PCA坐标和解释率
# ============================================================

write.csv(
  plot_data,
  file = "PCA_sample_scores.csv",
  row.names = FALSE
)

variance_table <- data.frame(
  Principal_component = paste0(
    "PC",
    seq_along(variance_percent)
  ),
  Explained_variance_percent = variance_percent,
  Cumulative_variance_percent = cumsum(variance_percent)
)

write.csv(
  variance_table,
  file = "PCA_explained_variance.csv",
  row.names = FALSE
)

cat(
  "\n分析完成，已输出：\n",
  "1. PCA_OW40_OW450.png\n",
  "2. PCA_OW40_OW450.pdf\n",
  "3. PCA_sample_scores.csv\n",
  "4. PCA_explained_variance.csv\n"
)

