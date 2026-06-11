#---
#title: "SDH inhibition snMultiome RNA & ATAC-seq"
#author: "Dakota Nuttall"
#date: "2025-06-17"
#R-version: 4.3.3
#---
###HELPER FUNCTIONS############################################################################

SCTransform.my <- function(srt_obj, assay = 'RNA',...){
  srt.out <- srt_obj |>
    SCTransform(assay = assay,
                method = "glmGamPoi",
                seed.use = 1234,
                return.only.var.genes = T,
                vars.to.regress = c('percent.mt','percent.ribo'
                ),
                vst.flavor = "v2", variable.features.n = 5000,
                ...)
  return(srt.out)
}

Subset_SCT <- function(global_srt,idents, npcs) {
  pre_SCT_subset <- global_srt
  if (!is.null(idents)){
  pre_SCT_subset <- subset(pre_SCT_subset, idents = idents)
  }

  pre_SCT_subset <- SCTransform(pre_SCT_subset, vars.to.regress = c("percent.mt","percent.ribo"),vst.flavor = "v2", method = "glmGamPoi", variable.features.n = 5000)
  pre_SCT_subset <- RunPCA(pre_SCT_subset, npcs = npcs)
 
  post_SCT_integrated <- IntegrateLayers(object = pre_SCT_subset, method = HarmonyIntegration, new.reduction = "harmony_pca", orig.reduction = "pca",)
  VariableFeatures(post_SCT_integrated) <- rownames(post_SCT_integrated@assays[["SCT"]]@scale.data)[!grepl(pattern="^mt-|^Rps|^Rpl|^Mrpl|^Mrps",rownames(post_SCT_integrated@assays[["SCT"]]@scale.data))]
 return(post_SCT_integrated)
}

proportion_plots <- function(object,cond_color,clust_color,name.clust,name.cond,height,width) {
  pt <- table(Idents(object), object$orig.ident)
  pt <- as.data.frame(pt)
  pt$Var1 <- as.character(pt$Var1)
  p1 <- ggplot(pt, aes(x = factor(Var1, levels = levels(object) ), y = Freq, fill = Var2)) +
    theme_bw(base_size = 15) +
    geom_col(position = "fill", width = 0.5,) +
    xlab("Cluster") +
    ylab("Proportion") +
    theme(legend.title = element_blank(),legend.text  = element_text( size = 20),axis.text=element_text(size=16),axis.title=element_text(size=18)) +
    scale_fill_manual(values = cond_color)
  
  pt$Var2 <- as.character(pt$Var2)
  pt$Var2<- factor(pt$Var2, levels=levels(object$orig.ident))
  p2 <- ggplot(pt, aes(x = Var2, y = Freq, fill = factor(Var1, levels = levels(object) ))) +
    theme_bw(base_size = 15) +
    geom_col(position = "fill", width = 0.5) +
    xlab("Condition") +
    ylab("Proportion") +
    theme(legend.title = element_blank(),axis.text=element_text(size=16),legend.text  = element_text( size = 20),axis.title=element_text(size=18)) +
    scale_fill_manual(values = clust_color)+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  print(p1)
  print(p2)
  ggsave(filename = name.clust,plot = p1, height = height, width = width)
  ggsave(filename = name.cond,plot = p2, height = height, width = width)
}

runPCA_no_mito <- function(srt_object,assay){
  non_mito_genes<-rownames(srt_obj@assays$SCT@counts)[!grepl(pattern="^mt-",rownames(srt_obj@assays$SCT@counts))]
  srt_PCA <- RunPCA(srt_object,assay = assay,features = non_mito_genes)
  return(srt_PCA)
}

quick_CT_plot <- function(object, pt.size, reduction, legend, label,font.size, resolution= NULL){
  if (legend == TRUE){
  p <- DimPlot(object, reduction = reduction, pt.size = pt.size, label = F, cols = clust_colors, group.by = resolution)+ theme(legend.position="right", legend.margin = margin(l=-1,r=0)) + labs(title = NULL,)+
    theme(legend.text  = element_text( size = 10), axis.text=element_blank(),axis.ticks=element_blank(), axis.line = element_line(linewidth =  1),axis.title=element_text(size=16))+ xlab("UMAP 1") + ylab("UMAP 2") 
  }
  else{
    p <- DimPlot(object, reduction = reduction, pt.size = pt.size, label = F, cols = clust_colors, group.by = resolution)+ theme(legend.position="right") + labs(title = NULL,)+
      theme(legend.text  = element_text( size = 10), axis.text=element_blank(),axis.ticks=element_blank(), axis.line = element_line(linewidth =  1),axis.title=element_text(size=16))+ xlab("UMAP 1") + ylab("UMAP 2") + NoLegend()
  }
  if (label == TRUE){
  LabelClusters(p, id = "ident", size = font.size)}
  else{
    p
  }
}

quick_cond_plot <- function(object, pt.size, reduction,legend, group.by = "orig.ident", labels = c("Saline Sham","Malonate Sham","Saline MI","Malonate MI"), cols= treatment_colors) {
  if (legend == TRUE){
  DimPlot(object,reduction = reduction, group.by = group.by, pt.size = pt.size, )+ theme(legend.position="right",axis.text=element_blank(),axis.ticks=element_blank(), axis.line = element_line(linewidth =  1),axis.title=element_text(size=16)) + labs(title = NULL)+
    theme(legend.text  = element_text( size = 10)) +
    scale_color_manual(labels=labels, values = alpha(cols,0.5)
    )+ xlab("UMAP 1") + ylab("UMAP 2")}
  else{
    DimPlot(object,reduction = reduction, group.by = group.by, pt.size = pt.size, )+ theme(legend.position="bottom",axis.text=element_blank(),axis.ticks=element_blank(), axis.line = element_line(linewidth =  1),axis.title=element_text(size=16)) + labs(title = NULL)+
      theme(legend.text  = element_text( size = 14)) +
      scale_color_manual(labels=labels, values = alpha(cols,0.5)
      )+ xlab("UMAP 1") + ylab("UMAP 2")+ NoLegend()
  }
}

Comparison_DEG <- function(seurat_obj,target_cluster,comparison_clusters){
  DEG_list <- list()
  n <- length(comparison_clusters)
  for (i in 1:n){
    markers <- FindMarkers(seurat_obj, ident.1 = target_cluster, ident.2 = comparison_clusters[i],only.pos = F, min.pct = 0.05, logfc.threshold = .25, test.use = "MAST", recorrect_umi = F, assay = "SCT")
    markers$BH_p_val_adj <- p.adjust(markers$p_val, method='BH')
    markers$gene <- rownames(markers)
    markers <- markers[!grepl("mt-",markers$gene),]
    markers <- subset(markers, BH_p_val_adj < 0.05)
    markers$comparison <- paste(target_cluster," vs. ",comparison_clusters[i],sep = "")
    DEG_list[[i]] <- markers
  }
  names(DEG_list) <- comparison_clusters
  return(DEG_list)
}

Comparison_cond_DEG <- function(seurat_obj,target_condition,comparison_conditions){
  DEG_list <- list()
  n <- length(comparison_conditions)
  for (i in 1:n){
    markers <- FindMarkers(seurat_obj, ident.1 = target_condition, ident.2 = comparison_conditions[i],group.by = "orig.ident",only.pos = T, min.pct = 0.05, logfc.threshold = .25, test.use = "MAST", recorrect_umi = F, assay = "SCT")
    markers$BH_p_val_adj <- p.adjust(markers$p_val, method='BH')
    markers$gene <- rownames(markers)
    markers <- markers[!grepl("mt-",markers$gene),]
    markers <- subset(markers, BH_p_val_adj < 0.05)
    markers$comparison <- paste(target_condition," vs. ",comparison_conditions[i],sep = "")
    DEG_list[[i]] <- markers
  }
  names(DEG_list) <- comparison_conditions
  return(DEG_list)
}

cond_Venn <- function(DEG_list,comparison_conditions, colors, file_name,title= NULL){
  gene_vecs <- vector("list", length(DEG_list))
  for(i in 1:length(DEG_list)){
    df <- DEG_list[[i]]
    genes <- df$gene
    gene_vecs[[i]] <-genes
  }
  names(gene_vecs) <- comparison_conditions
  
  # Create Plot
  venn.diagram(
    x = gene_vecs,
    category.names = comparison_conditions,
    filename = file_name,
    output = TRUE,
    imagetype = "tiff",
    fontfamily = "Arial",
    fill = colors,
    cat.fontfamily = "Arial",
    main = title,
    units = "in",
    height = 4,
    width = 4,
    cat.default.pos = "outer",
    cat.dist = -.0125
  )
  
  #Extract vector of intersecting genes
  intersecting_genes <- Reduce(intersect,gene_vecs)
  return(intersecting_genes)
}


Comparison_peaks <- function(seurat_obj,target_cluster,comparison_clusters,annotated_peaks_df ){
  DEG_list <- list()
  n <- length(comparison_clusters)
  for (i in 1:n){
    markers <- FindMarkers(seurat_obj, ident.1 = target_cluster, ident.2 = comparison_clusters[i],only.pos = F, min.pct = 0.05, logfc.threshold = .25, test.use = "wilcox", recorrect_umi = T, assay = "peakunion")
    markers$BH_p_val_adj <- p.adjust(markers$p_val, method='BH')
    markers$peak <- rownames(markers)
    markers <- subset(markers, BH_p_val_adj < 0.05)
    markers <- left_join(markers,annotated_peaks_df[, c("peak", "annotation", "SYMBOL")],by = "peak")
    markers$comparison <- paste(target_cluster," vs. ",comparison_clusters[i],sep = "")
    DEG_list[[i]] <- markers
    
  }
  names(DEG_list) <- comparison_clusters
  return(DEG_list)
}

Comparison_motif <- function(seurat_obj,target_cluster,comparison_clusters){
  DEG_list <- list()
  n <- length(comparison_clusters)
  for (i in 1:n){
    markers <- FindMarkers(seurat_obj, ident.1 = target_cluster, ident.2 = comparison_clusters[i],only.pos = F, recorrect_umi = T,
                           mean.fxn = rowMeans, assay = "chromvar", logfc.threshold = .25, test.use = "wilcox")
    markers$BH_p_val_adj <- p.adjust(markers$p_val, method='BH')
    markers$motif.ID <- rownames(markers)
    markers$motif.names <-  ConvertMotifID(seurat_obj,assay = "peakunion",id =row.names(markers))
    markers <- subset(markers, BH_p_val_adj < 0.05)
    markers$comparison <- paste(target_cluster," vs. ",comparison_clusters[i],sep = "")
    DEG_list[[i]] <- markers
  }
  names(DEG_list) <- comparison_clusters
  return(DEG_list)
}

Comparison_violin <- function(seurat_obj,target_cluster,comparison_clusters,feat,comparison_list, x.font.size=18){
  # Retrieve the identities column name and values safely
  idents <- Idents(seurat_obj)
  p_values <- c()
  label <- c()
  # Loop through each dataframe in the list
  for (i in seq_along(comparison_list)) {
    df <- comparison_list[[i]]
    
    # Find the adjusted p-value for the specific feature
    if (feat %in% df[[7]]) {
      p_values[i] <- df$BH_p_val_adj[df[[7]] == feat]
      label <- c(label, paste0("p < ", signif(df$BH_p_val_adj[df[[7]] == feat], 3)))
    } else {
      p_values[i] <- NA  # Assign NA if feat not found in that dataframe
      label <- c(label, "N.S.")
    }
    
  }
  
  # Create data frame for ggpubr annotation
  stat_df <- data.frame(
    group1 = comparison_clusters,
    group2 = rep(target_cluster,length(comparison_clusters)),
    p.adj = p_values,
    label = c(label
    )
  )
  # Make violin plot using Seurat’s VlnPlot
  
  # Fetch expression and cluster data
  plot_df <- FetchData(seurat_obj, vars = feat)
  plot_df$cluster <- as.factor(idents)
  
  
  # Build violin plot from scratch
  p <- ggplot(plot_df, aes(x = cluster, y = .data[[feat]])) +
    geom_violin(trim = FALSE, scale = "width",aes(fill = cluster)) +
    geom_jitter(width = 0.15,size = .5,alpha = 0.6,color = "black" 
    ) +
    stat_pvalue_manual(data = stat_df,label = "label",  y.position = c(1.25 * max(FetchData(seurat_obj, vars = feat)[, 1]),1.35*max(FetchData(seurat_obj, vars = feat)[, 1])), tip.length = 0.03,step.increase = 0.05,size = 10,bracket.size = 1.5) +
    theme_classic() + labs(title = feat, y = "Expression Level", x = NULL) +
    scale_fill_manual(values = clust_colors) + theme_classic() + ylim(0, 1.5 * max(FetchData(seurat_obj, vars = feat))) + scale_fill_manual(values = clust_colors) + NoLegend() + xlab(NULL) + theme(text=element_text(size=20),axis.text.x = element_text(size=x.font.size),axis.text.y = element_text(size=18),axis.title=element_text(size=18),plot.margin = margin(l = 0,unit = "in")) + labs(title = paste0(feat), y = "Expression Level", x = NULL)
  p
  return(p)
}

Comparison_violin_motif <- function(seurat_obj,target_cluster,comparison_clusters,feat,comparison_list,title, x.font.size=18){
  # Retrieve the identities column name and values safely
  idents <- Idents(seurat_obj)
  p_values <- c()
  label <- c()
  # Loop through each dataframe in the list
for (i in seq_along(comparison_list)) {
    df <- comparison_list[[i]]
    
    # Find the adjusted p-value for the specific feat
    if (feat %in% df$motif.ID) {
      p_values[i] <- df$BH_p_val_adj[df$motif.ID == feat]
      label <- c(label, paste0("p < ", signif(df$BH_p_val_adj[df$motif.ID == feat], 3)))
    } else {
      p_values[i] <- NA  # Assign NA if feat not found in that dataframe
      label <- c(label, "N.S.")
    }
    
  }
  
  # Create data frame for ggpubr annotation
  stat_df <- data.frame(
    group1 = comparison_clusters,
    group2 = rep(target_cluster,length(comparison_clusters)),
    p.adj = p_values,
    label = c(label
    )
  )
  # Make violin plot using Seurat’s VlnPlot
  
  # Fetch expression and cluster data
  plot_df <- FetchData(seurat_obj, vars = feat)
  plot_df$cluster <- as.factor(idents)
  
  
  # Build violin plot from scratch
  p <- ggplot(plot_df, aes(x = cluster, y = .data[[feat]])) +
    geom_violin(trim = FALSE, scale = "width",aes(fill = cluster)) +
    geom_jitter(width = 0.15,size = .5,alpha = 0.6,color = "black" 
    ) +
    stat_pvalue_manual(data = stat_df,label = "label",  y.position = c(1.15 * max(FetchData(seurat_obj, vars = feat)[, 1]),1.30*max(FetchData(seurat_obj, vars = feat)[, 1])), tip.length = 0.03,step.increase = 0.05,size = 10,bracket.size = 1.5) +
    theme_classic() + labs(title = feat, y = "Expression Level", x = NULL) +
    scale_fill_manual(values = clust_colors) + theme_classic() + ylim(NA, 1.55 * max(FetchData(seurat_obj, vars = feat))) + scale_fill_manual(values = clust_colors) + NoLegend() + xlab(NULL) + theme(text=element_text(size=20),axis.text.x = element_text(size=x.font.size),axis.text.y = element_text(size=18),axis.title=element_text(size=18),plot.margin = margin(l = 0,unit = "in")) + labs(title = paste0(title), y = "Expression Level", x = NULL)
  p
  return(p)

  
}

Comparison_violin_peaks <- function(seurat_obj,target_cluster,comparison_clusters,feat,pval_df, x.font.size=18){
  # Retrieve the identities column name and values safely
  idents <- Idents(seurat_obj)
 
  
  # Make violin plot using Seurat’s VlnPlot
  
  # Fetch expression and cluster data
  plot_df <- FetchData(seurat_obj, vars = feat)
  plot_df$cluster <- as.factor(idents)
  
  
  # Build violin plot from scratch
  p <- ggplot(plot_df, aes(x = cluster, y = .data[[feat]])) +
    geom_violin(trim = FALSE, scale = "width",aes(fill = cluster)) +
    geom_jitter(width = 0.15,size = .5,alpha = 0.6,color = "black" 
    ) +
    stat_pvalue_manual(size = 8,
      pval_df,
      label = "label",
      y.position = c(1.25 * max(FetchData(seurat_obj, vars = feat)[, 1]),1.35*max(FetchData(seurat_obj, vars = feat)[, 1])), tip.length = 0.03,step.increase = 0.05,bracket.size = 1.25
    ) +
    theme_classic() + labs(title = feat, y = "Expression Level", x = NULL) +
    scale_fill_manual(values = clust_colors) + theme_classic() + ylim(0, 1.5 * max(FetchData(seurat_obj, vars = feat))) + scale_fill_manual(values = clust_colors) + NoLegend() + xlab(NULL) + theme(text=element_text(size=20),axis.text.x = element_text(size=x.font.size),axis.text.y = element_text(size=18),axis.title=element_text(size=18),plot.margin = margin(l = 0,unit = "in")) + labs(title = paste0(feat), y = "Expression Level", x = NULL)
  p
  return(p)
}

split_into_two_lines <- function(text) {
  words <- unlist(strsplit(text, " "))
  mid <- ceiling(length(words) / 2)
  
  first_part <- paste(words[1:mid], collapse = " ")
  second_part <- paste(words[(mid+1):length(words)], collapse = " ")
  
  return(paste(first_part, second_part, sep = "\n"))
}

split_into_two_lines <- function(text) {
  words <- unlist(strsplit(text, " "))
  total_chars <- nchar(text)
  
  current_chars <- 0
  split_point <- 0
  
  # Find the point where the first part is longer than the second part character-wise
  for (i in 1:length(words)) {
    current_chars <- current_chars + nchar(words[i])
    if (current_chars >= total_chars / 2.5) {
      split_point <- i
      break
    }
  }
  
  first_part <- paste(words[1:split_point], collapse = " ")
  second_part <- paste(words[(split_point + 1):length(words)], collapse = " ")
  
  return(paste(first_part, second_part, sep = "\n"))
}

identify_low_RNA_clusters <- function(srt_object= seurat_object){
  RNA_count <- srt_object@meta.data$nCount_RNA
  clusters_ids <- levels(srt_object@active.ident)
  low_RNA_cluster <- c()
  for (i in clusters_ids){
    cluster_RNA_median <- median(RNA_count[srt_object@active.ident==i])
    if(cluster_RNA_median < 1000){
      low_RNA_cluster <- c(low_RNA_cluster,i)
    }
  }
  return(low_RNA_cluster)
}


prop_heatmap <- function(seurat_obj, filename){
  # Assume `seurat_obj` is your Seurat object
  # Extract metadata
  metadata <- seurat_obj@meta.data
  
  # Ensure relevant metadata columns exist
  # Replace 'orig.ident' and 'Cell_type' with the actual column names in your metadata
  if (!all(c("orig.ident", "Cell_state") %in% colnames(metadata))) {
    stop("Ensure metadata contains 'orig.ident' and 'Cell_state' columns.")
  }
  
  # Tabulate cell counts per treatment and cluster
  cell_counts <- metadata %>%
    dplyr::count(Cell_state, orig.ident) %>%
    pivot_wider(names_from = orig.ident, values_from = n, values_fill = 0)
  
  # Convert to a matrix for heatmap
  heatmap_matrix <- as.matrix(cell_counts[, -1])
  rownames(heatmap_matrix) <- cell_counts$Cell_state
  
  # Swap rows and columns
  heatmap_matrix_swapped <- t(heatmap_matrix)
  
  # Calculate proportions (percentage) by row
  heatmap_matrix_percent <- apply(heatmap_matrix_swapped, 2, function(x) (x / sum(x)) * 100)
  
  # Create the heatmap
  p <- pheatmap(
    mat = heatmap_matrix_percent,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("palegoldenrod", "navy"))(100),
    border_color = "black",
    angle_col = 45,            # Rotate x-axis labels to appear right-side up
    fontsize_row = 12,        # Adjust font size for row labels
    fontsize_col = 12,        # Adjust font size for column labels
    cellwidth = 30,           # Set a fixed cell width for square cells
    cellheight = 30,  # Set a fixed cell height for square cells
    breaks = seq(0, 100, 1),
    legend_breaks = seq(0, 100, 20),  # Set legend to range from 0 to 100%
    legend_labels = paste0(seq(0, 100, 20), "%"), # Format legend labels as percentages
    display_numbers = F,   # Show cell values
    number_format = "%.1f%%", # Format displayed numbers as percentages
    labels_row = rownames(heatmap_matrix_percent), # Ensure y-axis labels are displayed
    labels_col = colnames(heatmap_matrix_percent),  # Ensure x-axis labels are displayed
    title = F, legend = T, show_rownames = T
  )
  ggsave(filename,plot = p,height = 4, width = 8)
  return(p)
}

up_down_volcano <- function(DEG.df,selectLab,italicLab,up_label,down_label,up_x, down_x,volcano_name, cols = c('blue','orange','grey')){
  keyvals <- ifelse(
    DEG.df$avg_log2FC < 0 & DEG.df$BH_p_val_adj < 0.05, cols[1],
    ifelse(DEG.df$avg_log2FC > 0& DEG.df$BH_p_val_adj < 0.05, cols[2], 
           'black'))
  keyvals[is.na(keyvals)] <- cols[3]
  names(keyvals)[keyvals == cols[2]] <- 'Upregulated'
  names(keyvals)[keyvals == cols[3]] <- 'Not Significant'
  names(keyvals)[keyvals == cols[1]] <- 'Downregulated'
  p<- EnhancedVolcano(DEG.df,title = NULL,subtitle = NULL, x = "avg_log2FC", y= "p_val_adj", lab = italicLab, FCcutoff = 0, drawConnectors = T,colAlpha = 1,pCutoff = 0.05, lengthConnectors = 4.0,
                      arrowheads = F,
                      max.overlaps = Inf,
                      maxoverlapsConnectors = Inf,
                      parseLabels = T,
                      boxedLabels = F,
                      legendLabSize = 16,
                      legendIconSize = 5.0,
                      labSize = 5,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      cutoffLineType = 'blank',
                      colCustom = keyvals,
                      vline = 0,
                      vlineType = "solid",
                      selectLab = selectLab, widthConnectors = 1.0,
                      colConnectors = 'black',
                      xlab = expression(paste("Mean log"[2]," fold change")),
                      ylab = expression(paste("-log"[10]," adjusted ",italic("P")," value")))  + NoLegend() + coord_flip() +
    annotate("text", x = up_x, y = 0, label = up_label,  size = 5, hjust = 0)+
    annotate("text", x = down_x, y = 0, label = down_label,  size = 5,hjust = 0)
  
  ggsave(volcano_name,plot = p, width = 6, height = 5)
  p
}


Preprocess <- function(srt_obj, assay = 'RNA', ...) {

  selgenes<-rownames(srt_obj@assays$RNA$counts)[grepl(pattern="^mt-",rownames(srt_obj@assays$RNA$counts))]
  
  srt.out <- srt_obj |>
    NormalizeData(verbose = F) |>
    PercentageFeatureSet(features = selgenes, col.name = 'percent.mt', assay = assay)
  srt.out@meta.data[, 'percent.mt'][is.nan(srt.out@meta.data[, 'percent.mt'])] <- 0
  
  selgenes<-rownames(srt_obj@assays$RNA$counts)[grepl(pattern="^Rps|^Rpl|^Mrpl|^Mrps",rownames(srt_obj@assays$RNA$counts))]
  
  srt.out <- srt.out |>
    NormalizeData(verbose = F) |>
    PercentageFeatureSet(features = selgenes, col.name = 'percent.ribo', assay = assay)
  srt.out@meta.data[, 'percent.ribo'][is.nan(srt.out@meta.data[, 'percent.ribo'])] <- 0
return(srt.out)
  }


FindDimNumber <- function(srt_obj, var.toal = 0.95, reduction = 'pca'){
  if(is.null(srt_obj@reductions[[reduction]])){
    cat("Reduction", reduction, "not found!")
    return(NULL)
  }
  tmp.var <- (srt_obj@reductions[[reduction]]@stdev)^2
  var.cut <- var.toal*sum(tmp.var)
  dimNum = 0
  var.sum = 0
  while(var.sum < var.cut){
    dimNum = dimNum + 1
    var.sum <- var.sum + tmp.var[dimNum]
  }
  return(dimNum)}

DotPlot_SCT <- function(seurat_obj, gene_vector, pt.size=10) {
  
  # Ensure SCT assay is active
  DefaultAssay(seurat_obj) <- "SCT"
  
  # Check that genes exist
  valid_genes <- gene_vector[gene_vector %in% rownames(seurat_obj)]
  missing_genes <- setdiff(gene_vector, valid_genes)
  
  if (length(missing_genes) > 0) {
    message("❗ Missing genes removed: ",
            paste(missing_genes, collapse = ", "))
  }
  
  # Build plot
  p <- DotPlot(seurat_obj, features = valid_genes, assay = "SCT") +
    scale_color_gradient(low = "white", high = "blue") +
    labs(title = NULL,
         y = "Clusters",
         x = "Genes",
         color = "Avg. Scaled Exp",
         size = "% Expressing") +
    theme_bw(base_size = 13)+
    theme(text=element_text(size=pt.size), axis.text.x = element_text(size=pt.size+2))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Swap axes
  p <- p + coord_flip()
  
  return(p)
}

# a function that returns the position of n-th largest
maxn <- function(n) function(x) order(x, decreasing = TRUE)[n]

Scrublet_merged <- function(seurat, conditions){
  for (i in conditions){
    condition_seurat <- subset(seurat, orig.ident == i)
    write10xCounts(path = paste(i,"_counts",sep = ""),x=condition_seurat@assays$RNA$counts,  overwrite = T)
  }
}

plot_peaks <- function(srt,nearest_genes,peaks_of_interest,comp_list,target_clust,compare_clust,abbrev,x.font=18){
  index <- c(1:length(peaks_of_interest))
  DefaultAssay(srt) <- 'peakunion'
  for (i in index){
    #visualize differentially accessible peaks that have nearby significantly upregulated gene
    
    regions_highlight <- StringToGRanges(peaks_of_interest[i])
    regions_highlight$color <- "skyblue4"
    
    p<-CoveragePlot(
      srt,
      region = peaks_of_interest[i],
      region.highlight = regions_highlight,
      extend.upstream = 2000,
      extend.downstream = 2000,  
    )
    
    p & scale_fill_manual(values =  clust_colors) & theme( strip.text =element_text(size=16),axis.title=element_text(size=14),legend.text  = element_text( size = 12))
    ggsave(paste(abbrev,"_",nearest_genes[i],"_peak.tif",sep = ""),height = 5,width = 6)
    
    #Plot expression of peaks of interest
    Comparison_violin(seurat_obj = srt,target_cluster = target_clust,comparison_clusters = compare_clust,feat  = peaks_of_interest[i],comparison_list = comp_list,x.font.size = x.font)
    ggsave(paste(abbrev,"_",nearest_genes[i],"_peak_violin.tif", sep = ""),height = 5,width = 5)
  }
}

plot_DEGs <- function(srt,gene_list,comp_list,target_clust,compare_clust,x.font=18,feature.pt.size=2){
  index <- c(1:length(gene_list))
  for(i in index){
    if (gene_list[i] %in% Features(srt)){
  Comparison_violin(seurat_obj = srt,target_cluster = target_clust,comparison_clusters = compare_clust,feat = gene_list[i],comparison_list = comp_list,x.font.size = x.font)
  ggsave(paste(target_clust,gene_list[i],"vln_DEG.tif",sep ="_"),height = 5, width = 5)
  
  p <- FeaturePlot(
    object = srt,
    features = gene_list[i],
    pt.size = feature.pt.size,
    reduction = "WNN_umap"
  ) +
    theme(text=element_text(size=24),legend.text  = element_text( size = 16), axis.text=element_blank(),axis.ticks=element_blank(), axis.line = element_line(linewidth =  1))+ xlab(NULL) + ylab(NULL) +theme(axis.text=element_blank(),axis.title=element_text(size=20),legend.text  = element_text( size = 16),plot.title = element_text(face = "italic",size = 40))
  p
  ggsave(paste(target_clust,gene_list[i],"umap_DEG.tif",sep ="_"),height = 5, width = 5)
  }else{print(paste("No data on ",gene_list[i],sep=""))}
} 
}
