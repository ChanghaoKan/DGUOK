####设置多线程运行 <- 每次使用前运行####
#加载包
  library(future)
#查看可用核
  availableCores()
#当前可用核
  nbrOfWorkers() 
#八核处理
  plan(multisession, workers=40)
#扩大存储内存
  options(future.globals.maxSize= 200000*1024^2)    #当前限制为200GB 1000*1024^2为1GB
  
  
  
  
 ####创建目的Seurat对象####
#设置项目工作目录
  setwd("~/DGUOK")
#加载包
  library(dplyr)
  library(Seurat)
  library(tidyverse)
#读取DGUOK/DGUOK_Mutant
  Mutant <- Read10X(data.dir = "~/DGUOK/DGUOK_Mutant/raw_cell_gene_matrix/")
#读取DGUOK_Control
  Control <- Read10X(data.dir = "~/DGUOK/DGUOK_Control/raw_cell_gene_matrix/")
#创建Seurat对象
  seurat_Mutant <- CreateSeuratObject(counts = Mutant,        # counts：原始UMI计数矩阵（行是基因，列是细胞）   
                                      project = "Mutant",     # project：项目名称             
                                      min.cells = 3,          # min.cells：基因至少在多少个细胞中表达才保留                               
                                      min.features = 200)     # min.features：细胞至少检测到多少个基因才保留
  seurat_Control <- CreateSeuratObject(counts = Control,      # counts：原始UMI计数矩阵（行是基因，列是细胞）   
                                       project = "Control",   # project：项目名称             
                                       min.cells = 3,         # min.cells：基因至少在多少个细胞中表达才保留                               
                                       min.features = 200)    # min.features：细胞至少检测到多少个基因才保留
#保存原始Seurat对象
  saveRDS(seurat_Mutant,"seurat_Mutant.rds")
  saveRDS(seurat_Control,"seurat_Control.rds")
  
  
  
  
####对数据进行QC####
#读取数据
  seurat_Mutant <- readRDS("seurat_Mutant.rds") 
  seurat_Control <- readRDS("seurat_Control.rds") 
#对数据分别进行质控
#计算线粒体、核糖体和血红细胞基因比例
  #线粒体基因以"mt-"开头（小鼠基因命名规则）
  seurat_Mutant[["percent.mt"]] <- PercentageFeatureSet(seurat_Mutant,
                                                        pattern = "^mt-")
  seurat_Control[["percent.mt"]] <- PercentageFeatureSet(seurat_Control, 
                                                         pattern = "^mt-")
  #核糖体基因以"Rp[sl]"开头（大写Rp表示小鼠基因）
  
  seurat_Mutant[["percent.ribo"]] <- PercentageFeatureSet(seurat_Mutant,
                                                          pattern = "^Rp[sl]")
  seurat_Control[["percent.ribo"]] <- PercentageFeatureSet(seurat_Control, 
                                                           pattern = "^Rp[sl]")
  #血红细胞基因以"Hb[ab]"开头（小鼠中常见血红蛋白基因）
  seurat_Mutant[["percent.hb"]] <- PercentageFeatureSet(seurat_Mutant, 
                                                        pattern = "^Hb[ab]")
  seurat_Control[["percent.hb"]] <- PercentageFeatureSet(seurat_Control, 
                                                         pattern = "^Hb[ab]")
#绘制小提琴图查看nFeature_RNA，nCount_RNA，percent.mt
  VlnPlot(seurat_Mutant, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
          ncol = 3) 
  VlnPlot(seurat_Control, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
          ncol = 3) 
#定义最低标准
  min_nCount <- 1000      # UMI counts最低值
  min_nFeature <- 500     # 基因数最低值
  max_percent_mt <- 20    # 线粒体比例最大值
#应用最低标准
  seurat_Mutant <- subset(seurat_Mutant, 
                                   subset = nCount_RNA > min_nCount & 
                                     nFeature_RNA > min_nFeature & 
                                     percent.mt < max_percent_mt)
  seurat_Control <- subset(seurat_Control, 
                                    subset = nCount_RNA > min_nCount & 
                                      nFeature_RNA > min_nFeature & 
                                      percent.mt < max_percent_mt)
#动态质控函数
#祖传代码别动（^_^）
  dynamic_filter <- function(seurat_obj) {
    meta_data <- seurat_obj@meta.data %>%
      rownames_to_column(var = "cell_id") %>%
      group_by(orig.ident) %>%  # 使用orig.ident作为样本标识
      mutate(
        nCount_median = median(nCount_RNA),
        nCount_mad = mad(nCount_RNA),
        nFeature_median = median(nFeature_RNA),
        nFeature_mad = mad(nFeature_RNA),
        percent_mt_median = median(percent.mt),
        percent_mt_mad = mad(percent.mt),
        percent_ribo_median = median(percent.ribo),
        percent_ribo_mad = mad(percent.ribo),
        percent_hb_median = median(percent.hb),
        percent_hb_mad = mad(percent.hb)
      ) %>%
      filter(
        nCount_RNA >= nCount_median - 3 * nCount_mad,
        nCount_RNA <= nCount_median + 3 * nCount_mad,
        nFeature_RNA >= nFeature_median - 3 * nFeature_mad,
        nFeature_RNA <= nFeature_median + 3 * nFeature_mad,
        percent.mt <= percent_mt_median + 3 * percent_mt_mad,        # 只设上限，去除高线粒体含量细胞
        percent.ribo <= percent_ribo_median + 3 * percent_ribo_mad,  # 只设上限, 去除高核糖体含量细胞
        percent.hb <= percent_hb_median + 3 * percent_hb_mad         # 只设上限, 去除高血红细胞基因含量细胞
      ) %>%
      pull(cell_id)
    return(meta_data)
  }
#对Mutant和Control分别进行动态质控
  cells_filtered_Mutant <- dynamic_filter(seurat_Mutant)
  seurat_Mutant_filtered <- subset(seurat_Mutant, 
                                   cells = cells_filtered_Mutant)
  cells_filtered_Control <- dynamic_filter(seurat_Control)
  seurat_Control_filtered <- subset(seurat_Control, 
                                    cells = cells_filtered_Control)
#去除线粒体、核糖体和血红细胞基因
  mt_genes <- rownames(seurat_Mutant_filtered)[grepl("^mt-", 
                                                     rownames(seurat_Mutant_filtered))]
  ribo_genes <- rownames(seurat_Mutant_filtered)[grepl("^Rp[sl]", 
                                                       rownames(seurat_Mutant_filtered))]
  hb_genes <- rownames(seurat_Mutant_filtered)[grepl("^Hb[ab]",
                                                     rownames(seurat_Mutant_filtered))]
  genes_to_remove <- unique(c(mt_genes,
                              ribo_genes, 
                              hb_genes))
#获得处理后Seurat对象
  seurat_Mutant_filtered <- subset(seurat_Mutant_filtered, 
                                   features = setdiff(rownames(seurat_Mutant_filtered), 
                                                      genes_to_remove))
  seurat_Control_filtered <- subset(seurat_Control_filtered, 
                                    features = setdiff(rownames(seurat_Control_filtered), 
                                                       genes_to_remove))  
#保存处理后Seurat对象
  saveRDS(seurat_Mutant_filtered,"seurat_Mutant_filtered.rds")
  saveRDS(seurat_Control_filtered,"seurat_Control_filtered.rds")
  
  
  
  
####合并Seurat对象####
#读取数据
  seurat_Mutant_filtered <- readRDS("seurat_Mutant_filtered.rds") 
  seurat_Control_filtered <- readRDS("seurat_Control_filtered.rds") 
#保持细胞名单一
#为每个对象添加前缀以避免重复
  seurat_Mutant_filtered <- RenameCells(seurat_Mutant_filtered, 
                                        add.cell.id = "Mutant")
  seurat_Control_filtered <- RenameCells(seurat_Control_filtered, 
                                         add.cell.id = "Control")
#保持基因一致性
#提取Mutant和Control的共有基因
  common_genes <- intersect(rownames(seurat_Mutant_filtered), 
                            rownames(seurat_Control_filtered))
  seurat_Mutant_filtered <- subset(seurat_Mutant_filtered,
                                   features = common_genes)
  seurat_Control_filtered <- subset(seurat_Control_filtered, 
                                    features = common_genes)
#合并数据
  seurat_combined <- merge(seurat_Mutant_filtered, seurat_Control_filtered)
#合并Counts矩阵 
  seurat_combined <- JoinLayers(seurat_combined)
#归一化前细胞总表达量
  hist(colSums(seurat_combined@assays[["RNA"]]@layers[["counts"]]),
       breaks = 100,
       main = "Total expression before normalisation",
       xlab = "Sum of expression")
#进行归一化 <- 细胞标准化
#使用全局缩放归一化方法“LogNormalize”
#每个细胞每个基因的特征计数除以该细胞（一列）的特征总计数，再乘以scale.factor(默认10,000)，然后对结果进行log1p对数转换（log1p=log(n+1)）
  seurat_combined <- NormalizeData(seurat_combined, 
                                  normalization.method = "LogNormalize", 
                                  scale.factor = 10000)
#不调整，为默认值则使用这个
  #seurat <- NormalizeData(seurat)
#归一化后细胞总表达量
  hist(colSums(seurat_combined@assays[["RNA"]]@layers[["data"]]),
       breaks = 100,
       main = "Total expression after normalisation",
       xlab = "Sum of expression")  
#查看高变基因数量
  seurat_combined <- FindVariableFeatures(seurat_combined, 
                                         selection.method = "vst", 
                                         nfeatures = 2000) #nfeatures <- 高变基因数量
#高变基因标准化
#做PCA时，最好使用高变基因，否则会引入噪声。低丰度，变化低的基因
#使用高变基因标准化，
  seurat_combined <- ScaleData(seurat_combined)
#使用高变基因进行PCA
  seurat_combined <- RunPCA(seurat_combined, 
                           features = VariableFeatures(object = seurat_combined))  
#保存Seurat对象
  saveRDS(seurat_combined,"seurat_combined.rds")  





####使用harmony去批次####
#加载seurat对象
  seurat_combined <- readRDS("seurat_combined.rds")
  library(Seurat)
  library(harmony)
  set.seed(999)
#根据90%方差原理选择PCs
  xx <- cumsum(seurat_combined[["pca"]]@stdev^2)
  xx <- xx / max(xx)
  which(xx > 0.9) # 可以查看多少PCs解释了90%的方差，假设10%的方差来自于噪声，然后就可以选择相应的PCs
#dims出现的第一个参数记录
  ndim <- which(xx > 0.9)[1]
#使用harmony进行去批次
  seurat_combined <- RunHarmony(seurat_combined,
                               reduction = "pca",
                               group.by.vars = "orig.ident",
                               reduction.save = "harmony",
                               max_iter = 10)
#进行UMAP降维
  seurat_combined <- RunUMAP(seurat_combined, 
                            reduction = "harmony", 
                            dims = 1:ndim,
                            reduction.name = "umap")
#保存Seurat对象
  saveRDS(seurat_combined,"seurat_combined_harmony.rds")  
  
  
  
  
####绘制UMAP聚类图-origin.ident####
#加载包
  library(tidydr)
#加载seurat对象
  seurat_combined_harmony <- readRDS("seurat_combined_harmony.rds")
#提取原始身份信息
  orig_ident <- seurat_combined_harmony@meta.data$orig.ident
#合并UMAP数据和原始身份信息
  UMAP_with_orig_ident <- as.data.frame(seurat_combined_harmony@reductions$umap@cell.embeddings) %>%
    cbind(orig_ident = orig_ident)  # 确保这里使用的是正确的列名
#建立自定义主题
  mytheme <- theme_minimal() + 
    theme(plot.margin = margin(5.5, 15, 5.5, 5.5)) # 画布空白页缘调整
  colours <- c("#008695","#EE7C79")
#创建UMAP散点图
  p_orig_ident <- ggplot(UMAP_with_orig_ident, aes(x = umap_1, 
                                                   y = umap_2, 
                                                   color = orig_ident)) +
    geom_point(alpha = 0.5, size = 0.7) +  # 根据需要调整点的大小
    mytheme +
    theme(axis.text = element_text(size = 12), # 调整坐标轴文字大小
          axis.title = element_text(size = 14)) + # 调整坐标轴标题大小
    labs(x = "UMAP 1", y = "UMAP 2") +  # 命名坐标轴名称
    theme(legend.position = "right") +  # 根据需要调整图例位置
    scale_color_manual(values = colours)
#添加UMAP坐标轴箭头
  p_orig_ident <- p_orig_ident +
    theme_dr(xlength = 0.2, ylength = 0.2,
             arrow = grid::arrow(length = unit(0.1, "inches"),
                                 ends = 'last', type = "closed"))+
    labs(x = "UMAP 1",
         y = "UMAP 2", 
         color = "orig ident") +#命名坐标轴名称
    theme(panel.grid = element_blank())#去除背景条带
#增大注释大小
  p_orig_ident <- p_orig_ident + 
    guides(color = guide_legend(override.aes = list(size = 3)))
#绘图
  print(p_orig_ident)

  
  
  
  
####使用CHOIR进行细胞分群####
#下载专用
# remotes::install_github("corceslab/CHOIR", 
#                         ref="main",
#                         repos = BiocManager::repositories(),
#                         upgrade = "never")
#加载CHOIR包
  library(CHOIR)
#加载seurat对象
  seurat_combined_harmony <- readRDS("seurat_combined_harmony.rds")
#进行CHOIR
#我劝你别跑，主播的高性能服务器跑了一周（^_^#）
  object <- CHOIR(seurat_combined_harmony, 
                  n_cores = 40)
  saveRDS(object,"seurat_combined_CHOIR.rds")  
#显示聚类结果的UMAP图
  DimPlot(seurat_combined_harmony,
          reduction = "umap",
          group.by = "seurat_clusters",
          label=T ) 
#寻找亚型Marker
  seurat.markers <- FindAllMarkers(seurat_combined_harmony, 
                                   only.pos = T,     #是否选择只上调的基因
                                   min.pct = 0.25,   #某基因在细胞中表达的细胞数占相应cluster细胞数最低25%
                                   logfc.threshold = 0.25)
#保存上调的marker基因
  write.csv(seurat.markers,"seurat.markers_resolution_0.01.csv")
#保存之前Seurat对象
  saveRDS(seurat_combined_harmony,"seurat_combined_harmony_before.rds")

  
  
  
####根据ACT进行手动注释####
#绘图所需包
  library(ggsci)
  library(ggtext)
  library(ggrepel)
  library(tidyverse)
  library(paletteer)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(circlize)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(Seurat)
#读取Seurat对象
  seurat_combined_harmony <- readRDS("seurat_combined_harmony_before.rds")
#读取上调的marker基因
  seurat.markers <- read.csv("seurat.markers_resolution_0.01.csv", 
                             header = T,
                             row.names = 1)
#根据数据库进行手动命名
  new.cluster.ids <- c("Keratinocyte",                                #0
                       "Fibroblast",                                  #1
                       "Marcophage",                                  #2
                       "Endothelial cell of vascular tree",           #3
                       "Pericyte",                                    #4
                       "Melanocyte",                                  #5
                       "Endothelial cell"                             #6
                       )
#将新聚类ID名称映射到Seurat对象上
  names(new.cluster.ids) <- levels(seurat_combined_harmony)
#使用RenameIdents函数更新Seurat对象中聚类ID
  seurat_combined_harmony <- RenameIdents(seurat_combined_harmony, new.cluster.ids)
#将当前身份（重命名后的群集标识符）保存到新的元数据列 "new_cluster_idents"
  seurat_combined_harmony$new_cluster_idents <- Idents(seurat_combined_harmony)
#将new_cluster_idents命名为marker对应的cluster，便于画图
  seurat.markers$cluster <- new.cluster.ids[seurat.markers$cluster + 1]
#显示注释结果的UMAP图
  DimPlot(seurat_combined_harmony,reduction = "umap",label=T ) 
#保存之前Seurat对象
  saveRDS(seurat_combined_harmony,"seurat_combined_comment.rds")
#seurat绘图之前运载结束 
#读取之前保存的Seurat对象
  seurat_combined_harmony <- readRDS("seurat_combined_comment.rds") 
  
  
  
  
####绘制细胞比例图####
#加载包
  library(dplyr)
  library(ggplot2)
#生成细胞类型和样本的频数表
  tb <- data.frame(table(
    seurat_combined_harmony@meta.data$new_cluster_idents,
    seurat_combined_harmony@meta.data$orig.ident
  ))
  colnames(tb) <- c("Celltype", "Sample", "Freq")
#计算每个样本的总细胞数并合并
  sample_totals <- tb %>%
    group_by(Sample) %>%
    summarise(Total = sum(Freq)) %>%
    ungroup()
#更直观的百分比计算
  tb <- tb %>%
    left_join(sample_totals, by = "Sample") %>%
    mutate(Percentage = round(Freq / Total * 100, 1))  
#计算总体比例
  total_freq <- tb %>%
    group_by(Celltype) %>%
    summarise(Freq = sum(Freq))
  total_cells <- sum(total_freq$Freq)
  tb_overall <- total_freq %>%
    mutate(
      Sample = "Overall",
      Total = total_cells,
      Percentage = round(Freq / Total * 100, 1)
    )
#合并数据
  tb_combined <- bind_rows(tb, tb_overall)
#保存上调的marker基因
  write.csv(tb_combined,"tb_combined.csv")
  library(RColorBrewer)
  colours <- brewer.pal(n = n_celltypes, name = "Set3")  # 使用 Set3 调色板
#甜甜圈图
  ggplot(tb_combined, aes(x = 2, y = Percentage, fill = Celltype)) +  # 调整x值
    geom_col(width = 0.5, color = "white") +  # 调整宽度
    facet_wrap(~Sample, nrow = 1) +
    coord_polar(theta = "y") +
    xlim(c(0.5, 2.5)) +  # 调整xlim以控制环形大小
    scale_fill_manual(values = colours) +
    theme_void() +
    theme(
      strip.text.x = element_text(size = 14),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 15)
    )

  
  

####绘制各细胞类型marker基因展示####
#提取注释后筛选所得marker
  top5 <- read.csv("Dguok.marker.Primary.top5.csv",row.names = 1) # <- 自己手动获得
#提取基因名称，并计算平均表达值
  genes <- unique(top5$gene)
  aver_dt <- AverageExpression(seurat_combined_harmony, features = genes, 
                               group.by = 'new_cluster_idents', layer = 'data')
  aver_dt <- as.data.frame(aver_dt$RNA)
#提取具体的分组名称
  cluster_names <- levels(seurat_combined_harmony$new_cluster_idents)
#确认列数与分组名称数量匹配
  if (length(cluster_names) == ncol(aver_dt)) {
    colnames(aver_dt) <- cluster_names
  } else {
    stop("The number of clusters does not match the number of columns in the averaged data")
  }
#创建基因和细胞注释数据框
  gene_anno <- data.frame(gene_anno = top5$cluster, row.names = top5$gene)
  cell_anno <- data.frame(cell_anno = colnames(aver_dt), row.names = colnames(aver_dt))
#确保 gene_anno 的行名与 aver_dt 的行名一致并按顺序排列
  gene_anno <- gene_anno[rownames(aver_dt), , drop = FALSE]
#按 gene_anno 的值对 aver_dt 进行排序
  order_index <- order(gene_anno$gene_anno)
  aver_dt <- aver_dt[order_index, ]
  gene_anno <- gene_anno[order_index, , drop = FALSE]
#确保 gene_anno 的 gene_anno 列是因子类型
  gene_anno$gene_anno <- factor(gene_anno$gene_anno,
                                levels = unique(gene_anno$gene_anno))
#确保 cell_anno 的 cell_anno 列是因子类型
  cell_anno$cell_anno <- factor(cell_anno$cell_anno, 
                                levels = unique(cell_anno$cell_anno))
#确认 aver_dt 的列数
  ncol_aver_dt <- ncol(aver_dt)
  print(ncol_aver_dt)  # 输出 aver_dt 的列数
#生成正确数量的列名
  new_colnames <- paste('cluster ', 0:(ncol_aver_dt - 1), sep = '')
#将新列名赋给 aver_dt
  colnames(aver_dt) <- new_colnames
#Z-score
  htdf <- t(scale(t(aver_dt), scale = TRUE, center = TRUE))
#使用正确的变量来生成足够的颜色
  num_clusters <- length(cluster_names)
  colors <- brewer.pal(min(num_clusters, 12), "Set3")
  if (num_clusters > 12) {
    colors <- colorRampPalette(brewer.pal(12, "Set3"))(num_clusters)
  }
#创建 Z-Score 颜色映射
  col_fun <- colorRamp2(c(-2, 0, 2), c("#0099CC", "white", "#CC0033"))
#确保列名与具体分组名称匹配
  colnames(htdf) <- cluster_names
#创建 HeatmapAnnotation 并使用正确的颜色映射
  column_ha <- HeatmapAnnotation(cluster = colnames(htdf), 
                                 col = list(cluster = setNames(colors, cluster_names)))
#创建行注释
  row_cols <- setNames(rep(colors, each = 4), rownames(htdf))
  row_anno <- rowAnnotation(foo = anno_text(rownames(htdf), 
                                            location = 0, 
                                            just = "left", 
                                            gp = gpar(fill = row_cols, 
                                                      col = "black", 
                                                      fontface = 'italic'),
                                            width = max_text_width(rownames(htdf)) * 1.5))
#绘制热图
  Heatmap(htdf, 
          name = "Z-score", 
          cluster_columns = FALSE, #列聚类
          cluster_rows = FALSE,    #行聚类
          #左侧标题调整
          row_title = "Cell Type Marker genes",#标题名称
          row_title_gp = gpar(fontsize = 12, fontface = "bold"),#调整行标题的字体大小为 12，并设置字体为粗体。
          row_title_side = "left", #将行标题放置在热图的左侧
          #上方标题调整
          column_title = "Cell Type", #标题名称
          row_names_gp = gpar(col = NA),#隐藏左侧基因名称
          #row_names_gp = gpar(fontface = 'italic', fontsize = 10), #如需调整
          row_names_side = 'left',#将行名发放在左侧
          #设置单元格
          rect_gp = gpar(col = "white", lwd = 1), #设置单个外框
          column_names_side = 'top', column_names_rot = 45,#设置列名称
          top_annotation = column_ha, #设置上方注释
          col = col_fun, #设置颜色
          #分组
          row_split = rep(LETTERS[1:7], each = 4),#按照每组五行分组
          row_gap = unit(1.5, "mm"),#设置总外框
          border = T, 
          border_gp = gpar(col = "#151515", lwd = 1.2)) + row_anno
  
  
  
  
#绘制UMAP聚类图-celltype
#加载所需要的包
  library(tidydr)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(cols4all)
  library(mascarade)
  library(SeuratData)
  library(data.table)
#提取细胞类型和UMAP数据并合并
  celltype <- Idents(seurat_combined_harmony)
  UMAP <- as.data.frame(seurat_combined_harmony@reductions$umap@cell.embeddings) %>%
    cbind(celltype = celltype)
#设置配色
  colours <- c("Keratinocyte" = "#6E8FB2",
               "Fibroblast" = "#EAB67A",
               "Marcophage" = "#7DA494",
               "Endothelial cell of vascular tree" = "#C16E71",
               "Pericyte" = "#ABC8E5",
               "Melanocyte" = "#90C67C",
               "Endothelial cell" = "#D8ABC1")
#建立自定义主题：
  mytheme <- theme_void() + #空白主题，便于我们后期添加UMAP箭头
    theme(plot.margin = margin(5.5,15,5.5,5.5)) #画布空白页缘调整
#创建UMAP散点图并添加置信区间
  p <- ggplot(data = UMAP, aes(x = umap_1, y = umap_2)) +
    geom_point(aes(color = celltype), size = 0.4, alpha = 0.8) +
    #旧版添加范围圈
    #    stat_ellipse(aes(color = celltype, fill = celltype),
    #                 level = 0.95, linetype = 2, show.legend = FALSE,
    #                 geom = 'polygon', alpha = 0.1) +
    mytheme+
    scale_color_manual(values = colours)
#添加UMAP坐标轴箭头
  p <- p +
    theme_dr(xlength = 0.2, ylength = 0.2,
             arrow = grid::arrow(length = unit(0.1, "inches"),
                                 ends = 'last', type = "closed"))+
    labs(x = "UMAP 1",
         y = "UMAP 2",
         color = "cell type") +#命名坐标轴名称
    theme(panel.grid = element_blank())#去除背景条带
#在图中增加亚群标签
  label <- UMAP %>%
    group_by(celltype) %>%
    summarise(umap_1 = median(umap_1), umap_2 = median(umap_2))
  p <- p +
    geom_text(data = label,
              aes(x = umap_1, y = umap_2, label = celltype),
              fontface = "bold", color = 'black', size = 4)+
    guides(color = guide_legend(override.aes = list(size = 3)))
#为UMAP图添加范围
  maskTable <- generateMask(dims=UMAP, 
                            cluster=celltype, 
                            minDensity = 100 , 
                            smoothSigma = 0.05 )
#绘图
  print(p)
  
  
  
  
#### 提取黑色素细胞子集 ####
#加载已注释的Seurat对象
  seurat_combined_harmony <- readRDS("seurat_combined_comment.rds")
#提取黑色素细胞
  melanocyte_subset <- subset(seurat_combined_harmony, 
                              subset = new_cluster_idents == "Melanocyte")
#展示提取后的细胞
  table(melanocyte_subset@meta.data$new_cluster_idents)
#归一化前细胞总表达量
  hist(colSums(melanocyte_subset@assays[["RNA"]]@layers[["counts"]]),
       breaks = 100,
       main = "Total expression before normalisation",
       xlab = "Sum of expression")
#进行归一化 <- 细胞标准化
#使用全局缩放归一化方法“LogNormalize”
#每个细胞每个基因的特征计数除以该细胞（一列）的特征总计数，再乘以scale.factor(默认10,000)，然后对结果进行log1p对数转换（log1p=log(n+1)）
  melanocyte_subset <- NormalizeData(melanocyte_subset, 
                                   normalization.method = "LogNormalize", 
                                   scale.factor = 10000)
  #不调整，为默认值则使用这个
  #seurat <- NormalizeData(seurat)
#归一化后细胞总表达量
  hist(colSums(melanocyte_subset@assays[["RNA"]]@layers[["data"]]),
       breaks = 100,
       main = "Total expression after normalisation",
       xlab = "Sum of expression")  
#查看高变基因数量
  melanocyte_subset <- FindVariableFeatures(melanocyte_subset, 
                                          selection.method = "vst", 
                                          nfeatures = 2000) #nfeatures <- 高变基因数量
#查看最高变的10个基因
  top10 <- head(VariableFeatures(melanocyte_subset), 10)
#高变基因标准化
#做PCA时，最好使用高变基因，否则会引入噪声。低丰度，变化低的基因
#使用全局基因进行标准化，
  melanocyte_subset <- ScaleData(melanocyte_subset,
                               features = rownames(melanocyte_subset))
#使用高变基因进行PCA
  melanocyte_subset <- RunPCA(melanocyte_subset, features = VariableFeatures(object = melanocyte_subset))  
  
  
  
  
#保存之前Seurat对象
  saveRDS(melanocyte_subset,"melanocyte_subset_filter.rds")
#seurat去批次之前运载结束 
#读取之前保存的Seurat对象
  melanocyte_subset <- readRDS("melanocyte_subset_filter.rds") 
  
  
  
  
####对非免疫细胞亚群进行harmony去批次####
#加载seurat
  library(Seurat)
  library(harmony)
  set.seed(999)
  #UMAP可视化（未进行批次校正）
  #  seurat <- RunUMAP(seurat,reduction = "pca", dims = 1:10, reduction.name = "umap_naive")
  #  p1 <- DimPlot(seurat, reduction = "umap_naive",group.by = "orig.ident")
  #  p1
#选择不同范围的dims判断分为两类的显著特点
  DimHeatmap(melanocyte_subset, dims = 25:30, cells = 200, balanced = TRUE)
#绘制肘部图以确定重要的维度数量
  ElbowPlot(melanocyte_subset, reduction = "pca", ndims = 25)
#根据90%方差原理选择PCs
  xx <- cumsum(melanocyte_subset[["pca"]]@stdev^2)
  xx <- xx / max(xx)
  which(xx > 0.9) # 可以查看多少PCs解释了90%的方差，假设10%的方差来自于噪声，然后就可以选择相应的PCs
#dims出现的第一个参数记录
  ndim <- which(xx > 0.9)[1]
#使用harmony进行去批次
  melanocyte_subset <- RunHarmony(melanocyte_subset,
                                reduction = "pca",
                                group.by.vars = "orig.ident",
                                reduction.save = "harmony",
                                max.iter.harmony = 10)
#为排除harmony聚类可能掩盖误差，再次使用pca进行二次运算，结果类似，故选择harmony
  melanocyte_subset <- RunUMAP(melanocyte_subset, reduction = "harmony", dims = 1:ndim,reduction.name = "umap") 
#保存Seurat对象
  saveRDS(melanocyte_subset,"melanocyte_subset_harmony.rds")  
  
  
  
  
####绘制UMAP聚类图-origin.ident####
#加载包
  library(tidydr)
#加载seurat对象
  melanocyte_subset <- readRDS("melanocyte_subset_harmony.rds")
#提取原始身份信息
  orig_ident <- melanocyte_subset@meta.data$orig.ident
#合并UMAP数据和原始身份信息
  UMAP_with_orig_ident <- as.data.frame(melanocyte_subset@reductions$umap@cell.embeddings) %>%
    cbind(orig_ident = orig_ident)  # 确保这里使用的是正确的列名
#建立自定义主题
  mytheme <- theme_minimal() + 
    theme(plot.margin = margin(5.5, 15, 5.5, 5.5)) # 画布空白页缘调整
  colours <- c("#008695","#EE7C79")
#创建UMAP散点图
  p_orig_ident <- ggplot(UMAP_with_orig_ident, aes(x = umap_1, 
                                                   y = umap_2, 
                                                   color = orig_ident)) +
    geom_point(alpha = 1, size = 1) +  # 根据需要调整点的大小
    mytheme +
    theme(axis.text = element_text(size = 12), # 调整坐标轴文字大小
          axis.title = element_text(size = 14)) + # 调整坐标轴标题大小
    labs(x = "UMAP 1", y = "UMAP 2") +  # 命名坐标轴名称
    theme(legend.position = "right") +  # 根据需要调整图例位置
    scale_color_manual(values = colours)
#添加UMAP坐标轴箭头
  p_orig_ident <- p_orig_ident +
    theme_dr(xlength = 0.2, ylength = 0.2,
             arrow = grid::arrow(length = unit(0.1, "inches"),
                                 ends = 'last', type = "closed"))+
    labs(x = "UMAP 1",
         y = "UMAP 2", 
         color = "orig ident") +#命名坐标轴名称
    theme(panel.grid = element_blank())#去除背景条带
#增大注释大小
  p_orig_ident <- p_orig_ident + 
    guides(color = guide_legend(override.aes = list(size = 3)))
#绘图
  print(p_orig_ident)

  
  
  
#绘制UMAP聚类图-new.cluster.ids
#提取原始身份信息
  new_cluster_idents <- melanocyte_subset@meta.data$new_cluster_idents
#合并UMAP数据和原始身份信息
  UMAP_with_new_cluster_idents <- as.data.frame(melanocyte_subset@reductions$umap@cell.embeddings) %>%
    cbind(new_cluster_idents = new_cluster_idents)  # 确保这里使用的是正确的列名
#建立自定义主题
  mytheme <- theme_minimal() + 
    theme(plot.margin = margin(5.5, 15, 5.5, 5.5)) # 画布空白页缘调整
#设置配色
  colours <- c("Keratinocyte" = "#6E8FB2",
               "Fibroblast" = "#EAB67A",
               "Marcophage" = "#7DA494",
               "Endothelial cell of vascular tree" = "#C16E71",
               "Pericyte" = "#ABC8E5",
               "Melanocyte" = "#90C67C",
               "Endothelial cell" = "#D8ABC1")
#创建UMAP散点图
  p_orig_ident <- ggplot(UMAP_with_new_cluster_idents, aes(x = umap_1, y = umap_2, color = new_cluster_idents)) +
    geom_point(size = 0.5) +  # 根据需要调整点的大小
    mytheme +
    theme(axis.text = element_text(size = 12), # 调整坐标轴文字大小
          axis.title = element_text(size = 14)) + # 调整坐标轴标题大小
    labs(x = "UMAP 1", y = "UMAP 2") +  # 命名坐标轴名称
    theme(legend.position = "right")+  # 根据需要调整图例位置
    scale_color_manual(values = colours) 
#添加UMAP坐标轴箭头
  p_orig_ident <- p_orig_ident +
    theme_dr(xlength = 0.2, ylength = 0.2,
             arrow = grid::arrow(length = unit(0.1, "inches"),
                                 ends = 'last', type = "closed"))+
    labs(x = "UMAP 1",
         y = "UMAP 2", 
         color = "cell type") +#命名坐标轴名称
    theme(panel.grid = element_blank())#去除背景条带
#在图中增加亚群标签
  label <- UMAP_with_new_cluster_idents %>%
    group_by(new_cluster_idents) %>%
    summarise(umap_1 = median(umap_1), umap_2 = median(umap_2))
  p_orig_ident <- p_orig_ident +
    geom_text(data = label,
              aes(x = umap_1, y = umap_2, label = new_cluster_idents),
              fontface = "bold", color = 'black', size = 4)+
    guides(color = guide_legend(override.aes = list(size = 3)))
#增大注释大小
  p_orig_ident <- p_orig_ident + guides(color = guide_legend(override.aes = list(size = 3)))
#绘图
  print(p_orig_ident)  
  
 
  
  
####进行基因集评分####
#加载必要的包
  library(Seurat)
  library(homologene)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  
#加载黑色素细胞子集
  melanocyte_subset <- readRDS("melanocyte_subset_harmony.rds")
#提取log归一化后的表达矩阵
  expr_matrix <- GetAssayData(melanocyte_subset, assay = "RNA", layer = "data")
  expr_matrix <- as.matrix(expr_matrix)
  dim(expr_matrix)
#设置基因集
  geneSet_human <- read.csv("ROS Top50.csv")
  human_genes <- geneSet_human$Gene.Symbol
  mapped_genes <- homologene(human_genes, inTax = 9606, outTax = 10090)
  geneSet_mouse <- unique(mapped_genes$`10090`)
  present_genes <- geneSet_mouse[geneSet_mouse %in% rownames(expr_matrix)]
  missing_genes <- setdiff(geneSet_mouse, rownames(expr_matrix))
  if (length(missing_genes) > 0) {
    warning(paste("以下小鼠基因在表达矩阵中未找到：", paste(missing_genes, collapse = ", ")))
  }
#计算基因集评分
  gene_sets <- list(MyGeneSet = present_genes)
  melanocyte_subset <- AddModuleScore(object = melanocyte_subset, features = gene_sets, ctrl = 100, name = "GeneSet_")
  score_col <- grep("GeneSet_", colnames(melanocyte_subset@meta.data), value = TRUE)
  colnames(melanocyte_subset@meta.data)[colnames(melanocyte_subset@meta.data) == score_col] <- "Score"
  melanocyte_subset$Module_Score_Z <- scale(melanocyte_subset$Score)[,1]
 #提取元数据
  metadata <- melanocyte_subset@meta.data
#小提琴图
  p_violin <- ggviolin(metadata, 
                       x = "orig.ident", y = "Module_Score_Z", 
                       fill = "orig.ident",
                       add = c("boxplot", "mean_sd"),
                       palette = "npg",
                       xlab = "", ylab = "Gene Set Score") +
    stat_compare_means(method = "t.test", label.y = max(metadata$Module_Score_Z) * 1.1)
  print(p_violin)
#密度图
  median_values <- metadata %>%
    group_by(orig.ident) %>%
    summarise(median_z = median(Module_Score_Z, na.rm = TRUE))
  p_density <- ggplot(metadata, aes(x = Module_Score_Z, fill = orig.ident)) +
    geom_density(alpha = 0.5, color = NA) +
    geom_vline(data = median_values, aes(xintercept = median_z, color = orig.ident), linetype = "dashed", linewidth = 0.8, show.legend = FALSE) +
    scale_fill_brewer(palette = "Set2", name = "Sample") +
    scale_color_brewer(palette = "Set2", guide = "none") +
    labs(x = "Gene Set Score Z-Score", y = "Density", title = "Gene Set Score Distribution with Medians") +
    theme_classic() +
    theme(legend.position = "right", panel.grid.major = element_line(color = "grey90", linewidth = 0.2), plot.title = element_text(hjust = 0.5))
  print(p_density)
  
