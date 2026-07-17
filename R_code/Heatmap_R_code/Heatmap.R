##############################################################
# All Lipids Heatmap
# Blue gradient Level 0-1
##############################################################


##############################################################
# 1. 加载R包
##############################################################

library(readxl)
library(magrittr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(viridis)
library(MASS)
library(cowplot)
library(grid)
library(gridExtra)
library(ggrepel)
library(Cairo)
library(openxlsx)
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(tidyselect)
library(FactoMineR)
library(circlize)
library(plyr)


##############################################################
# 2. 数据读取与标准化
##############################################################

setwd("C:/R_Project/Revision/Heatmap")

# 提取所有Lipids数据

subseth2_all <- read.csv(
  "Heatmap_data.csv",
  header = TRUE,
  check.names = FALSE
)

subseth_all <- subseth2_all

# 转换为数值矩阵

subseth2_all <- apply(
  
  as.matrix(
    subseth2_all[,2:13]
  ),
  
  2,
  
  as.numeric
  
)


#-------------------------------------------------------------
# Log标准化函数
#-------------------------------------------------------------

normalizeh_all <- function(x){
  
  min_x <- min(x)
  
  max_x <- max(x)
  
  
  scaled <- 
    log(x - min_x + 1) /
    log(max_x - min_x + 1)
  
  
  return(scaled)
  
}



# 行标准化

h_all <- t(
  
  apply(
    
    subseth2_all,
    
    1,
    
    normalizeh_all
    
  )
)



# 设置行名

rownames(h_all) <- subseth_all$Lipids




##############################################################
# 3. 主热图颜色设置
##############################################################

color_gradient <- colorRampPalette(
  
  colors = c(
    
    "#FFFFFF",   # Level 0
    
    "#DCEEFF",   # Level 0.25
    
    "#4F9BD5",   # Level 0.5
    
    "#2871A9",   # Level 0.75
    
    "#00467D"    # Level 1
    
  )
  
)




##############################################################
# 4. 热图单元格绘制函数
##############################################################

my_cell_fun <- function(
    j,
    i,
    x,
    y,
    width,
    height,
    fill
){
  
  
  data_value <- h_all[i,j]
  
  
  grid.rect(
    
    x=x,
    y=y,
    
    width=width,
    height=height,
    
    gp=gpar(
      
      fill=fill,
      
      col=NA
      
    )
    
  )
  
}




##############################################################
# 5. 顶部实验组注释
##############################################################

class <- anno_block(
  
  gp=gpar(
    
    fill=c(
      
      "#FFFFFF",
      "#FFFFFF"
      
    ),
    
    col="white"
    
  ),
  
  
  height=unit(
    
    4,
    
    "mm"
    
  ),
  
  
  labels=c(
    
    "αSyn-YFP",
    "YFP"
    
  ),
  
  
  labels_gp=gpar(
    
    col="#4B4B4B",
    
    fontsize=10
    
  )
  
)



strains <- HeatmapAnnotation(
  
  group=class
  
)




##############################################################
# 6. 主Heatmap绘制
##############################################################

heatmap1 <- Heatmap(
  
  
  h_all,
  
  
  name="Level",
  
  
  top_annotation=strains,
  
  
  bottom_annotation=NULL,
  
  
  column_title=NULL,
  
  
  col=color_gradient(200),
  
  
  cell_fun=my_cell_fun,
  
  
  # 聚类设置
  
  column_km=2,
  
  
  column_dend_height=
    unit(
      0.5,
      "cm"
    ),
  
  
  cluster_rows=FALSE,
  
  
  cluster_columns=TRUE,
  
  
  # 标签设置
  
  show_column_names=FALSE,
  
  
  show_row_names=FALSE,
  
  
  # 图例设置
  
  heatmap_legend_param=list(
    
    title="Level",
    
    at=c(
      
      0,
      0.5,
      1
      
    ),
    
    labels=c(
      
      "0",
      "0.5",
      "1"
      
    )
    
  )
  
)



heatmap1





##############################################################
# 7. 脂质类别侧边注释
##############################################################

letter <- data.frame(
  
  Lipids=rownames(h_all)
  
)



# 提取脂质类别

letter$Lipids <- substr(
  
  letter$Lipids,
  
  1,
  
  3
  
)



letter$Lipids <- gsub(
  
  "\\(|\\)",
  
  "",
  
  letter$Lipids
  
)



letter$Lipids <- gsub(
  
  "\\d+",
  
  "",
  
  letter$Lipids
  
)



letter <- letter$Lipids




##############################################################
# 8. 脂质类别标签
##############################################################

row_anno <- rowAnnotation(
  
  
  mark_gene =
    
    anno_mark(
      
      at=c(
        
        19,
        29,
        38,
        47,
        86,
        126,
        145,
        161,
        182,
        194,
        199,
        204,
        210,
        219,
        247,
        268,
        272
        
      ),
      
      
      labels=c(
        
        "Cer",
        "DAG",
        "FA",
        "Glc",
        "GPC",
        "GPE",
        "GPI",
        "GPS",
        "LPC",
        "LPE",
        "MAG",
        "NAE",
        "PA",
        "PC",
        "PE",
        "GPG",
        "SM",
        "TAG"
        
      )
      
    )
  
)





##############################################################
# 9. 脂质类别颜色
##############################################################

colors <- c(
  
  "#0072B2",
  "#E69F00",
  "#56B4E9",
  "#009E73"
  
)



color_vector <- rep(
  
  colors,
  
  length.out=30
  
)



color_vector <- setNames(
  
  color_vector,
  
  unique(as.vector(letter))
  
)




##############################################################
# 10. 构建脂质类别Heatmap
##############################################################

heatmap2 <- Heatmap(
  
  
  letter,
  
  
  show_column_names=FALSE,
  
  
  col=color_vector,
  
  
  show_heatmap_legend=FALSE,
  
  
  right_annotation=row_anno
  
)



heatmap2





##############################################################
# 11. 合并热图
##############################################################

heatmap1 + heatmap2
