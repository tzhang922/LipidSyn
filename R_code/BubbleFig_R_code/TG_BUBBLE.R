library(tidyverse)
library(ggplot2)
library(ggrepel)
library(grid)

# 不要在 dplyr 后加载 plyr
# library(plyr)

#------------------------------------------------------------
# 1. 读取数据
#------------------------------------------------------------

df_bind2 <- read.csv(
  "df_bind2.csv",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!"Lipids" %in% names(df_bind2)) {
  stop("数据中不存在 Lipids 列，请检查列名。")
}

#------------------------------------------------------------
# 2. 提取 TG 数据
#------------------------------------------------------------

subsetS_TG <- df_bind2 %>%
  dplyr::filter(
    stringr::str_detect(Lipids, "^TG\\(")
  )

if (nrow(subsetS_TG) == 0) {
  stop("没有找到以 TG( 开头的脂质。")
}

#------------------------------------------------------------
# 3. 识别样本列
# 假定 Lipids 后面的前12列为样本列
#------------------------------------------------------------

lipid_col_index <- match("Lipids", names(subsetS_TG))

sample_cols <- names(subsetS_TG)[
  seq(
    from = lipid_col_index + 1,
    length.out = min(12, ncol(subsetS_TG) - lipid_col_index)
  )
]

if (length(sample_cols) < 12) {
  warning("检测到的样本列少于12列，请检查 sample_cols。")
}

# 转换为数值
subsetS_TG <- subsetS_TG %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(sample_cols),
      ~ suppressWarnings(
        as.numeric(
          gsub(",", "", trimws(as.character(.x)))
        )
      )
    )
  )

#------------------------------------------------------------
# 4. 从脂质名称提取总碳数和总双键数
#------------------------------------------------------------

extract_tg_structure <- function(x) {
  
  fatty_acids <- stringr::str_extract_all(
    x,
    "\\d+:\\d+"
  )[[1]]
  
  if (length(fatty_acids) < 3) {
    return(
      tibble::tibble(
        CabonNumber = NA_real_,
        DoubleBonds = NA_real_
      )
    )
  }
  
  fatty_acids <- fatty_acids[1:3]
  
  carbon_numbers <- as.numeric(
    stringr::str_extract(
      fatty_acids,
      "^\\d+"
    )
  )
  
  double_bonds <- as.numeric(
    stringr::str_extract(
      fatty_acids,
      "(?<=:)\\d+"
    )
  )
  
  tibble::tibble(
    CabonNumber = sum(carbon_numbers, na.rm = TRUE),
    DoubleBonds = sum(double_bonds, na.rm = TRUE)
  )
}

tg_structure <- purrr::map_dfr(
  subsetS_TG$Lipids,
  extract_tg_structure
)

subsetS_TG <- dplyr::bind_cols(
  subsetS_TG,
  tg_structure
)

subsetS_TG <- subsetS_TG %>%
  dplyr::filter(
    !is.na(CabonNumber),
    !is.na(DoubleBonds)
  )

#------------------------------------------------------------
# 5. 相同碳数和双键数的 TG 合并
#------------------------------------------------------------

TG_CB <- subsetS_TG %>%
  dplyr::select(
    dplyr::all_of(sample_cols),
    CabonNumber,
    DoubleBonds
  )

TG_CB_2 <- TG_CB %>%
  dplyr::group_by(
    CabonNumber,
    DoubleBonds
  ) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(sample_cols),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

#------------------------------------------------------------
# 6. 每个样本内换算为相对百分比
#------------------------------------------------------------

TG_CB_3 <- TG_CB_2 %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(sample_cols),
      ~ {
        total_value <- sum(.x, na.rm = TRUE)
        
        if (total_value == 0) {
          rep(0, length(.x))
        } else {
          .x / total_value * 100
        }
      }
    )
  )

# 前6列视为OW40，后6列视为OW450
OW40_cols <- sample_cols[1:min(6, length(sample_cols))]
OW450_cols <- sample_cols[
  7:min(12, length(sample_cols))
]

TG_CB_3$OW40 <- rowMeans(
  TG_CB_3[, OW40_cols, drop = FALSE],
  na.rm = TRUE
)

TG_CB_3$OW450 <- rowMeans(
  TG_CB_3[, OW450_cols, drop = FALSE],
  na.rm = TRUE
)

# 防止全部为NA时得到NaN
TG_CB_3$OW40[is.nan(TG_CB_3$OW40)] <- NA
TG_CB_3$OW450[is.nan(TG_CB_3$OW450)] <- NA

TG_CB_4 <- TG_CB_3 %>%
  dplyr::select(
    CabonNumber,
    DoubleBonds,
    OW40,
    OW450
  ) %>%
  dplyr::mutate(
    Structure = paste0(
      "TG(",
      CabonNumber,
      ":",
      DoubleBonds,
      ")"
    )
  )

#------------------------------------------------------------
# 7. 气泡图
#------------------------------------------------------------

p_bubble <- ggplot(
  TG_CB_4,
  aes(
    x = CabonNumber,
    y = DoubleBonds
  )
) +
  # OW450 / YFP
  geom_point(
    aes(
      size = OW450,
      colour = "YFP"
    ),
    shape = 16,
    alpha = 0.22
  ) +
  geom_point(
    aes(
      size = OW450,
      colour = "YFP"
    ),
    shape = 1,
    alpha = 0.9,
    stroke = 0.7
  ) +
  
  # OW40 / αSyn-YFP
  geom_point(
    aes(
      size = OW40,
      colour = "αSyn-YFP"
    ),
    shape = 16,
    alpha = 0.22
  ) +
  geom_point(
    aes(
      size = OW40,
      colour = "αSyn-YFP"
    ),
    shape = 1,
    alpha = 0.9,
    stroke = 0.7
  ) +
  
  scale_x_continuous(
    limits = c(46, 60),
    breaks = seq(46, 60, 2),
    expand = expansion(mult = c(0.03, 0.08))
  ) +
  
  scale_y_continuous(
    limits = c(0, 15),
    breaks = seq(0, 15, 5),
    expand = expansion(mult = c(0.03, 0.08))
  ) +
  
  scale_size_continuous(
    range = c(1, 15),
    guide = "none"
  ) +
  
  scale_colour_manual(
    values = c(
      "YFP" = "#E69F00",
      "αSyn-YFP" = "#56B4E9"
    )
  ) +
  
  labs(
    x = "Carbon Number",
    y = "Double Bonds",
    colour = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.line = element_line(
      colour = "black",
      linewidth = 0.5
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    ),
    legend.text = element_text(size = 14),
    legend.background = element_blank(),
    legend.direction = "horizontal",
    legend.position = c(0.29, 0.92),
    legend.key.size = grid::unit(1, "lines")
  )

#------------------------------------------------------------
# 8. 添加标签
#------------------------------------------------------------

p_final <- p_bubble +
  
  geom_text_repel(
    data = TG_CB_4 %>%
      dplyr::filter(
        Structure %in% c(
          "TG(58:10)",
          "TG(58:9)",
          "TG(56:8)",
          "TG(56:7)",
          "TG(56:6)"
        )
      ),
    aes(
      label = Structure,
      colour = "αSyn-YFP"
    ),
    fontface = "bold",
    size = 3.5,
    box.padding = grid::unit(0.35, "lines"),
    point.padding = grid::unit(0.3, "lines"),
    segment.color = "#56B4E9",
    nudge_x = 1.5,
    nudge_y = 1,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  
  geom_text_repel(
    data = TG_CB_4 %>%
      dplyr::filter(
        Structure %in% c(
          "TG(49:2)",
          "TG(50:2)",
          "TG(48:2)",
          "TG(49:3)",
          "TG(50:3)"
        )
      ),
    aes(
      label = Structure,
      colour = "YFP"
    ),
    fontface = "bold",
    size = 3.5,
    box.padding = grid::unit(0.35, "lines"),
    point.padding = grid::unit(0.3, "lines"),
    segment.color = "#E69F00",
    nudge_x = -1.5,
    nudge_y = -1.2,
    show.legend = FALSE,
    max.overlaps = Inf
  )

print(p_final)

ggsave(
  filename = "TG_carbon_doublebond_bubble_plot.pdf",
  plot = p_final,
  width = 7,
  height = 6,
  device = cairo_pdf
)

ggsave(
  filename = "TG_carbon_doublebond_bubble_plot.png",
  plot = p_final,
  width = 7,
  height = 6,
  dpi = 600,
  bg = "white"
)




#-----------------------------------------------------------------------------------------------以下为All PE饱和度计算并做气泡图
#================================================================================
# All PE：总碳数和总双键数计算及气泡图
#================================================================================

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(grid)
library(openxlsx)

# 重要：
# 不建议加载 plyr，因为 plyr 会覆盖 dplyr::summarise() 和 dplyr::mutate()
# 如果前面已经运行 library(plyr)，下面代码通过 dplyr:: 明确指定函数来源

#--------------------------------------------------------------------------------
# 1. 检查数据
#--------------------------------------------------------------------------------

if (!exists("df_bind2")) {
  stop("当前环境中不存在 df_bind2，请先读取 df_bind2.csv。")
}

if (!"Lipids" %in% names(df_bind2)) {
  stop("df_bind2 中不存在 Lipids 列，请检查列名。")
}

#--------------------------------------------------------------------------------
# 2. 提取包含 PE、但不包含 LPE 的脂质
#--------------------------------------------------------------------------------

subseth_PE <- df_bind2[
  !grepl("LPE", df_bind2$Lipids, ignore.case = TRUE) &
    grepl("PE", df_bind2$Lipids, ignore.case = TRUE),
]

subsetS_PE <- subseth_PE

if (nrow(subsetS_PE) == 0) {
  stop("未提取到 PE 数据，请检查 Lipids 列中的脂质命名。")
}

cat("提取到的 PE 数量：", nrow(subsetS_PE), "\n")

# 查看提取到的脂质名称
print(head(subsetS_PE$Lipids, 20))

#--------------------------------------------------------------------------------
# 3. 设置样本列
#
# 根据原始代码，默认第2至第13列为12个样本：
# 第2至第7列：OW40，共6个重复
# 第8至第13列：OW450，共6个重复
#--------------------------------------------------------------------------------

if (ncol(subsetS_PE) < 13) {
  stop("数据列数少于13列，无法按照第2至第13列提取12个样本。")
}

sample_cols <- names(subsetS_PE)[2:13]

OW40_cols <- sample_cols[1:6]
OW450_cols <- sample_cols[7:12]

cat("\n样本列：\n")
print(sample_cols)

cat("\nOW40样本列：\n")
print(OW40_cols)

cat("\nOW450样本列：\n")
print(OW450_cols)

#--------------------------------------------------------------------------------
# 4. 将样本列转换为数值型
#--------------------------------------------------------------------------------

subsetS_PE <- subsetS_PE %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(sample_cols),
      ~ {
        x <- trimws(as.character(.x))
        
        # 去掉千位分隔符
        x <- gsub(",", "", x, fixed = TRUE)
        
        # 将常见非数值标记转换为NA
        x[x %in% c(
          "", "-", "--", "NA", "N/A",
          "ND", "N.D.", "nd", "n.d.",
          "NaN", "NULL"
        )] <- NA_character_
        
        suppressWarnings(as.numeric(x))
      }
    )
  )

# 查看各样本列中的缺失值数量
PE_NA_check <- data.frame(
  Sample = sample_cols,
  NA_number = colSums(
    is.na(subsetS_PE[, sample_cols, drop = FALSE])
  ),
  stringsAsFactors = FALSE
)

print(PE_NA_check)

#--------------------------------------------------------------------------------
# 5. 定义 PE 总碳数和总双键数提取函数
#
# 可识别：
# PE(16:0/18:1)       -> PE(34:1)
# PE(16:0_18:1)       -> PE(34:1)
# PE(16:0/18:1(9Z))   -> PE(34:1)
# PE(34:1)            -> PE(34:1)
#--------------------------------------------------------------------------------

extract_pe_structure <- function(x) {
  
  if (is.na(x) || trimws(x) == "") {
    return(
      data.frame(
        CabonNumber = NA_real_,
        DoubleBonds = NA_real_
      )
    )
  }
  
  # 提取全部“数字:数字”形式
  lipid_parts <- stringr::str_extract_all(
    x,
    "\\d+\\s*:\\s*\\d+"
  )[[1]]
  
  if (length(lipid_parts) == 0) {
    return(
      data.frame(
        CabonNumber = NA_real_,
        DoubleBonds = NA_real_
      )
    )
  }
  
  # 提取碳数
  carbon_values <- as.numeric(
    stringr::str_extract(
      lipid_parts,
      "^\\d+"
    )
  )
  
  # 提取双键数
  double_bond_values <- as.numeric(
    stringr::str_extract(
      lipid_parts,
      "(?<=:)\\s*\\d+"
    )
  )
  
  # 判断是总组成形式还是脂肪酸链形式
  # 如PE(34:1)只使用第一个值
  # 如PE(16:0/18:1)将前两个值相加
  lipid_inside <- stringr::str_match(
    x,
    "\\((.*)\\)"
  )[, 2]
  
  has_chain_separator <- !is.na(lipid_inside) &&
    stringr::str_detect(
      lipid_inside,
      "[/_]"
    )
  
  if (has_chain_separator && length(carbon_values) >= 2) {
    
    carbon_number <- sum(
      carbon_values[1:2],
      na.rm = TRUE
    )
    
    double_bonds <- sum(
      double_bond_values[1:2],
      na.rm = TRUE
    )
    
  } else {
    
    carbon_number <- carbon_values[1]
    double_bonds <- double_bond_values[1]
  }
  
  data.frame(
    CabonNumber = carbon_number,
    DoubleBonds = double_bonds
  )
}

#--------------------------------------------------------------------------------
# 6. 批量提取总碳数和总双键数
#--------------------------------------------------------------------------------

PE_structure <- purrr::map_dfr(
  subsetS_PE$Lipids,
  extract_pe_structure
)

subsetS_PE <- dplyr::bind_cols(
  subsetS_PE,
  PE_structure
)

# 检查无法提取的脂质名称
PE_unrecognized <- subsetS_PE %>%
  dplyr::filter(
    is.na(CabonNumber) |
      is.na(DoubleBonds)
  ) %>%
  dplyr::select(
    Lipids,
    CabonNumber,
    DoubleBonds
  )

if (nrow(PE_unrecognized) > 0) {
  
  warning(
    paste0(
      "有 ",
      nrow(PE_unrecognized),
      " 个PE名称未能识别，请检查对象 PE_unrecognized。"
    )
  )
  
  print(PE_unrecognized)
}

# 保留成功识别的PE
subsetS_PE_valid <- subsetS_PE %>%
  dplyr::filter(
    !is.na(CabonNumber),
    !is.na(DoubleBonds)
  )

if (nrow(subsetS_PE_valid) == 0) {
  stop("没有任何PE名称成功提取总碳数和双键数。")
}

# 保存结构提取检查结果
PE_structure_check <- subsetS_PE_valid %>%
  dplyr::select(
    Lipids,
    CabonNumber,
    DoubleBonds
  )

openxlsx::write.xlsx(
  PE_structure_check,
  file = "PE_structure_check.xlsx",
  rowNames = FALSE,
  overwrite = TRUE
)

#--------------------------------------------------------------------------------
# 7. 整理样本丰度数据
#--------------------------------------------------------------------------------

PE_CB <- subsetS_PE_valid %>%
  dplyr::select(
    CabonNumber,
    DoubleBonds,
    dplyr::all_of(sample_cols)
  )

#--------------------------------------------------------------------------------
# 8. 合并具有相同总碳数和总双键数的PE
#
# 必须明确使用dplyr::summarise，避免plyr冲突
#--------------------------------------------------------------------------------

PE_CB_2 <- PE_CB %>%
  dplyr::group_by(
    CabonNumber,
    DoubleBonds
  ) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(sample_cols),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

# 查看合并结果
print(PE_CB_2)

#--------------------------------------------------------------------------------
# 9. 每个样本列内部标准化为百分比
#
# 每个样本中所有PE总和为100%
#--------------------------------------------------------------------------------

PE_CB_3 <- PE_CB_2 %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(sample_cols),
      ~ {
        column_sum <- sum(.x, na.rm = TRUE)
        
        if (is.na(column_sum) || column_sum == 0) {
          rep(0, length(.x))
        } else {
          .x / column_sum * 100
        }
      }
    )
  )

# 检查各列百分比之和，正常情况下应接近100
PE_percentage_check <- colSums(
  PE_CB_3[, sample_cols, drop = FALSE],
  na.rm = TRUE
)

cat("\n各样本PE百分比之和：\n")
print(PE_percentage_check)

#--------------------------------------------------------------------------------
# 10. 计算OW40和OW450组平均值
#--------------------------------------------------------------------------------

PE_CB_3$OW40 <- rowMeans(
  PE_CB_3[, OW40_cols, drop = FALSE],
  na.rm = TRUE
)

PE_CB_3$OW450 <- rowMeans(
  PE_CB_3[, OW450_cols, drop = FALSE],
  na.rm = TRUE
)

# 全部为NA时rowMeans可能得到NaN
PE_CB_3$OW40[is.nan(PE_CB_3$OW40)] <- NA_real_
PE_CB_3$OW450[is.nan(PE_CB_3$OW450)] <- NA_real_

#--------------------------------------------------------------------------------
# 11. 生成最终气泡图数据
#--------------------------------------------------------------------------------

PE_CB_4 <- PE_CB_3 %>%
  dplyr::select(
    CabonNumber,
    DoubleBonds,
    OW40,
    OW450
  ) %>%
  dplyr::mutate(
    CabonNumber = as.numeric(CabonNumber),
    DoubleBonds = as.numeric(DoubleBonds),
    OW40 = as.numeric(OW40),
    OW450 = as.numeric(OW450),
    Structure = paste0(
      "PE(",
      CabonNumber,
      ":",
      DoubleBonds,
      ")"
    )
  ) %>%
  dplyr::arrange(
    CabonNumber,
    DoubleBonds
  )

print(PE_CB_4)

# 保存绘图源数据
openxlsx::write.xlsx(
  PE_CB_4,
  file = "PE_bubble_plot_source_data.xlsx",
  rowNames = FALSE,
  overwrite = TRUE
)

#--------------------------------------------------------------------------------
# 12. 设置需要标注的PE
#
# 修正了原代码中的：
# "PE38:6" -> "PE(38:6)"
#--------------------------------------------------------------------------------

alphaSyn_labels <- c(
  "PE(40:1)",
  "PE(38:5)",
  "PE(38:4)",
  "PE(38:3)",
  "PE(38:6)"
)

YFP_labels <- c(
  "PE(33:1)",
  "PE(32:0)",
  "PE(34:0)",
  "PE(33:0)"
)

alphaSyn_label_data <- PE_CB_4 %>%
  dplyr::filter(
    Structure %in% alphaSyn_labels
  )

YFP_label_data <- PE_CB_4 %>%
  dplyr::filter(
    Structure %in% YFP_labels
  )

# 检查哪些标签实际存在
cat("\nαSyn-YFP组实际存在的标签：\n")
print(alphaSyn_label_data$Structure)

cat("\nYFP组实际存在的标签：\n")
print(YFP_label_data$Structure)

#--------------------------------------------------------------------------------
# 13. 绘制PE气泡图
#
# OW450：YFP，橙色
# OW40：αSyn-YFP，蓝色
#--------------------------------------------------------------------------------

p_PE <- ggplot(
  PE_CB_4,
  aes(
    x = CabonNumber,
    y = DoubleBonds
  )
) +
  
  # OW450 / YFP：半透明实心晕影
  geom_point(
    aes(
      size = OW450,
      colour = "YFP"
    ),
    shape = 16,
    alpha = 0.20,
    na.rm = TRUE
  ) +
  
  # OW450 / YFP：外部圆环
  geom_point(
    aes(
      size = OW450,
      colour = "YFP"
    ),
    shape = 1,
    alpha = 0.82,
    stroke = 0.6,
    na.rm = TRUE
  ) +
  
  # OW40 / αSyn-YFP：半透明实心晕影
  geom_point(
    aes(
      size = OW40,
      colour = "αSyn-YFP"
    ),
    shape = 16,
    alpha = 0.20,
    na.rm = TRUE
  ) +
  
  # OW40 / αSyn-YFP：外部圆环
  geom_point(
    aes(
      size = OW40,
      colour = "αSyn-YFP"
    ),
    shape = 1,
    alpha = 0.82,
    stroke = 0.6,
    na.rm = TRUE
  ) +
  
  scale_x_continuous(
    limits = c(26, 42),
    breaks = seq(26, 42, 2),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  
  scale_y_continuous(
    limits = c(0, 7),
    breaks = seq(0, 7, 2),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  
  scale_size_continuous(
    range = c(0, 20),
    breaks = c(0, 5, 10, 15, 20),
    guide = "none"
  ) +
  
  scale_colour_manual(
    values = c(
      "YFP" = "#E69F00",
      "αSyn-YFP" = "#56B4E9"
    )
  ) +
  
  labs(
    x = "Carbon Number",
    y = "Double Bonds",
    colour = NULL
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    axis.text.x = element_text(
      size = 14,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 14,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 14,
      colour = "black"
    ),
    axis.title.y = element_text(
      size = 14,
      colour = "black"
    ),
    axis.line = element_line(
      colour = "black",
      linewidth = 0.5
    ),
    panel.border = element_rect(
      fill = NA,
      colour = "black",
      linewidth = 0.5
    ),
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    legend.text = element_text(
      size = 14
    ),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.direction = "horizontal",
    legend.position = c(0.29, 0.92),
    legend.key.size = grid::unit(
      1,
      "lines"
    ),
    plot.margin = margin(
      t = 10,
      r = 15,
      b = 10,
      l = 10
    )
  )

#--------------------------------------------------------------------------------
# 14. 添加αSyn-YFP组标签
#--------------------------------------------------------------------------------

p_PE_label1 <- p_PE +
  ggrepel::geom_text_repel(
    data = alphaSyn_label_data,
    aes(
      label = Structure,
      colour = "αSyn-YFP"
    ),
    fontface = "bold",
    size = 3.5,
    box.padding = grid::unit(
      0.35,
      "lines"
    ),
    point.padding = grid::unit(
      0.30,
      "lines"
    ),
    segment.color = "#56B4E9",
    segment.size = 0.4,
    min.segment.length = 0,
    nudge_x = 1,
    nudge_y = 0.5,
    max.overlaps = Inf,
    show.legend = FALSE,
    seed = 123
  )

#--------------------------------------------------------------------------------
# 15. 添加YFP组标签
#--------------------------------------------------------------------------------

p_PE_final <- p_PE_label1 +
  ggrepel::geom_text_repel(
    data = YFP_label_data,
    aes(
      label = Structure,
      colour = "YFP"
    ),
    fontface = "bold",
    size = 3.5,
    box.padding = grid::unit(
      0.35,
      "lines"
    ),
    point.padding = grid::unit(
      0.30,
      "lines"
    ),
    segment.color = "#E69F00",
    segment.size = 0.4,
    min.segment.length = 0,
    nudge_x = -0.5,
    nudge_y = 0.5,
    max.overlaps = Inf,
    show.legend = FALSE,
    seed = 456
  )

# 显示图形
print(p_PE_final)

#--------------------------------------------------------------------------------
# 16. 保存图片
#--------------------------------------------------------------------------------

ggsave(
  filename = "All_PE_carbon_doublebond_bubble_plot.pdf",
  plot = p_PE_final,
  width = 7,
  height = 6,
  device = grDevices::cairo_pdf,
  bg = "white"
)

ggsave(
  filename = "All_PE_carbon_doublebond_bubble_plot.png",
  plot = p_PE_final,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "All_PE_carbon_doublebond_bubble_plot.tiff",
  plot = p_PE_final,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

cat(
  "\nPE分析完成，生成以下文件：\n",
  "1. All_PE_carbon_doublebond_bubble_plot.pdf\n",
  "2. All_PE_carbon_doublebond_bubble_plot.png\n",
  "3. All_PE_carbon_doublebond_bubble_plot.tiff\n",
  "4. PE_bubble_plot_source_data.xlsx\n",
  "5. PE_structure_check.xlsx\n"
)