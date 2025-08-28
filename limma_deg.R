#' 差异表达分析
#'
#' @title 基于limma进行两组间差异表达分析
#' @param od 结果输出路径
#' @param DEG_exp 表达谱，基因在行样本在列
#' @param DEG_pdata 样本分组文件，第一列为样本，第二列为分组
#' @param controlLabel 对照组的分类标签
#' @param caseLabel 实验组的分类标签
#' @param DEG_FC 等于log(差异倍数阈值)，默认值1
#' @param DEG_P 差异显著性检验P值
#' @param pvalue 是否使用未校正的P，默认NULL是不使用；如果不是NULL，则使用未校正的P
#' @param saveplot 是否生成结果图片，包括火山图和热图，默认为FALSE不生成
#' @param color_fun 热图中样本分组的颜色
#' @return list
#' @examples
#' res <- limma_deg(od = out_dir, DEG_exp = x_exp, DEG_pdata = pdata,
#'                  controlLabel = "normal", caseLabel = "tumor",
#'                  DEG_FC = 1, DEG_P = 0.05, color_fun = color_fun1)
#' @author CY

limma_deg <- function(DEG_exp, DEG_pdata, controlLabel, caseLabel, 
                      DEG_FC = 1, DEG_P = 0.05, pvalue = FALSE, 
                      saveplot = FALSE, color_fun = NULL, od = NULL) {
  library(limma)
  library(ggplot2)
  library(pheatmap)

  DEG_exp <- DEG_exp[, DEG_pdata$sample]

  group_factor <- factor(DEG_pdata$group, levels = c(controlLabel, caseLabel))
  design <- model.matrix(~ group_factor)

  fit <- lmFit(DEG_exp, design)
  fit <- eBayes(fit)

  results <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  results$log2FC <- results$logFC
  results$threshold <- if (pvalue) {
    abs(results$log2FC) > DEG_FC & results$P.Value < DEG_P
  } else {
    abs(results$log2FC) > DEG_FC & results$adj.P.Val < DEG_P
  }

  if (saveplot && !is.null(od)) {
    if (!dir.exists(od)) dir.create(od, recursive = TRUE)

    # Volcano plot
    volcano_data <- data.frame(log2FC = results$log2FC,
                                negLogP = -log10(ifelse(pvalue, results$P.Value, results$adj.P.Val)),
                                Significant = results$threshold)

    volcano_plot <- ggplot(volcano_data, aes(x = log2FC, y = negLogP, color = Significant)) +
      geom_point(alpha = 0.7) +
      theme_minimal() +
      scale_color_manual(values = c("grey", "red")) +
      labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10(P-value)")

    ggsave(filename = file.path(od, "volcano_plot.png"), plot = volcano_plot, width = 6, height = 5)

    # Heatmap
    selected_genes <- rownames(results)[results$threshold]
    if (length(selected_genes) > 1) {
      pdf(file = file.path(od, "heatmap.pdf"))
      pheatmap(DEG_exp[selected_genes, ], 
               annotation_col = data.frame(Group = DEG_pdata$group),
               show_rownames = FALSE,
               cluster_cols = TRUE,
               color = if (is.null(color_fun)) colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100) else color_fun)
      dev.off()
    }
  }

  return(results)
}
