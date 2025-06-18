---
title: "SDH inhibition snMultiome RNA & ATAC-seq"
author: "Dakota Nuttall"
date: "2025-06-17"
R-version: 4.3.3
---
###HELPER FUNCTIONS############################################################################

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
    theme(legend.title = element_blank(),legend.text  = element_text(face="bold", size = 14),axis.text=element_text(size=16, face = "bold"),axis.title=element_text(size=18,face="bold")) +
    scale_fill_manual(values = cond_color)
  
  pt$Var2 <- as.character(pt$Var2)
  pt$Var2<- factor(pt$Var2, levels=levels(object$orig.ident))
  p2 <- ggplot(pt, aes(x = Var2, y = Freq, fill = factor(Var1, levels = levels(object) ))) +
    theme_bw(base_size = 15) +
    geom_col(position = "fill", width = 0.5) +
    xlab("Condition") +
    ylab("Proportion") +
    theme(legend.title = element_blank(),axis.text=element_text(size=16, face = "bold"),legend.text  = element_text(face="bold", size = 14),axis.title=element_text(size=18,face="bold")) +
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
    theme(legend.text  = element_text(face="bold", size = 10), axis.text=element_blank(),axis.ticks=element_blank(), axis.line = element_line(linewidth =  1),axis.title=element_text(size=16,face="bold"))+ xlab("UMAP 1") + ylab("UMAP 2") 
  }
  else{
    p <- DimPlot(object, reduction = reduction, pt.size = pt.size, label = F, cols = clust_colors)+ theme(legend.position="right") + labs(title = NULL,)+
      theme(legend.text  = element_text(face="bold", size = 10), axis.text=element_blank(),axis.ticks=element_blank(), axis.line = element_line(linewidth =  1),axis.title=element_text(size=16,face="bold"))+ xlab("UMAP 1") + ylab("UMAP 2") + NoLegend()
  }
  if (label == TRUE){
  LabelClusters(p, id = "ident",  fontface = "bold", size = font.size)}
  else{
    p
  }
}

quick_cond_plot <- function(object, pt.size, reduction,legend, group.by = "orig.ident", labels = c("Saline Sham","Malonate Sham","Saline MI","Malonate MI"), cols= treatment_colors) {
  if (legend == TRUE){
  DimPlot(object,reduction = reduction, group.by = group.by, pt.size = pt.size, )+ theme(legend.position="right",axis.text=element_blank(),axis.ticks=element_blank(), axis.line = element_line(linewidth =  1),axis.title=element_text(size=16,face="bold")) + labs(title = NULL)+
    theme(legend.text  = element_text(face="bold", size = 10)) +
    scale_color_manual(labels=labels, values = alpha(cols,0.5)
    )+ xlab("UMAP 1") + ylab("UMAP 2")}
  else{
    DimPlot(object,reduction = reduction, group.by = group.by, pt.size = pt.size, )+ theme(legend.position="bottom",axis.text=element_blank(),axis.ticks=element_blank(), axis.line = element_line(linewidth =  1),axis.title=element_text(size=16,face="bold")) + labs(title = NULL)+
      theme(legend.text  = element_text(face="bold", size = 14)) +
      scale_color_manual(labels=labels, values = alpha(cols,0.5)
      )+ xlab("UMAP 1") + ylab("UMAP 2")+ NoLegend()
  }
}


quick_feature_plots <- function(object,feature, min.cutoff,reduction,pt.size, umap_name,vln_sub_name,vln_cond_name, group.by = "orig.ident",cond_cols= treatment_colors ){
  p1<- FeaturePlot(
    object = object,
    features = feature,
    min.cutoff = min.cutoff,
    pt.size = pt.size,
    reduction = reduction
  ) +
    theme(text=element_text(size=24),legend.text  = element_text(face="bold", size = 16), axis.text=element_blank(),axis.ticks=element_blank(), axis.line = element_line(linewidth =  1))+ xlab(NULL) + ylab(NULL) +theme(axis.text=element_blank(),axis.title=element_text(size=20,face="bold"),legend.text  = element_text(face="bold", size = 16))
  
  p2 <- VlnPlot(object, feature,pt.size = 0.1, cols = clust_colors) + NoLegend()+ xlab(NULL) + theme(text=element_text(size=20),axis.text=element_text(size=18,face="bold"),axis.title=element_text(size=18,face="bold"))
  p3 <- VlnPlot(object, feature,group.by = group.by,pt.size = 0.1, cols = cond_cols) + NoLegend()+ xlab(NULL) + theme(text=element_text(size=20),axis.text=element_text(size=18,face="bold"),axis.title=element_text(size=18,face="bold"))
 
  ggsave(filename = umap_name,plot = p1, height = 6, width = 6)

  ggsave(filename = vln_sub_name,plot = p2, height = 5, width = 5)

  ggsave(filename = vln_cond_name,plot = p3, height = 6, width = 5)
  return(list(p1,p2,p3))
}

quick_motif_plots <- function(object,motif, min.cutoff,reduction,pt.size,title.name){
  p1<- FeaturePlot(
    object = object,
    features = motif,
    min.cutoff = min.cutoff,
    pt.size = pt.size,
    reduction = reduction,
     cols = c("lightgrey", "darkred")
  ) + labs(title = title.name) + xlab("UMAP 1") + ylab("UMAP 2")+theme(axis.text=element_text(size=18),axis.title=element_text(size=18,face="bold"),legend.text  = element_text(face="bold", size = 12))
  
  p2 <- VlnPlot(object, motif,pt.size = 0.1, cols = clust_colors,assay = "chromvar") + NoLegend()+ xlab(NULL) + theme(text=element_text(size=14),axis.text=element_text(size=18,face="bold"),axis.title=element_text(size=18,face="bold"))+ labs(title = title.name)
  p3 <- VlnPlot(object, motif,group.by = "orig.ident",pt.size = 0.1, cols = treatment_colors,assay = "chromvar") + NoLegend()+ xlab(NULL) + theme(text=element_text(size=18),axis.text=element_text(size=18,face="bold"),axis.title=element_text(size=18,face="bold"))+ labs(title = title.name)
  print(p1)
  print(p2)
  print(p3)
}



quick_umaps <- function(object,reduction, pt.size,labels, label.size){
  p1 <- DimPlot(object, reduction = reduction, pt.size = pt.size, label = F, cols = clust_colors)+ theme(legend.position="bottom",axis.text=element_text(size=16),axis.title=element_text(size=18)) + labs(title = NULL)+
    theme(legend.text  = element_text(face="bold", size = 16))+xlab("UMAP 1")+ylab("UMAP 2") 
  p1<- LabelClusters(p1, id = "ident",  fontface = "bold", size = label.size)
  p2<-DimPlot(object,reduction = reduction, group.by = "orig.ident", pt.size = pt.size, )+ theme(legend.position="bottom",axis.text=element_text(size=16),axis.title=element_text(size=18)) + labs(title = NULL)+
    theme(legend.text  = element_text(face="bold", size = 18)) +
    scale_color_manual(labels=labels, values = alpha(treatment_colors,0.5)
    )+xlab("UMAP 1")+ylab("UMAP 2")
  print(p1)
  print(p2)
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

identify_high_mito_clusters <- function(srt_object= seurat_object){
  mt_percent <- srt_object@meta.data$percent.mt
  clusters_ids <- levels(srt_object@active.ident)
  high_mito_cluster <- c()
  for (i in clusters_ids){
    cluster_mt_percent_mean <- mean(mt_percent[srt_object@active.ident==i])
    if(cluster_mt_percent_mean > 4.5){
      high_mito_cluster <- c(high_mito_cluster,i)
    }
  }
  return(high_mito_cluster)
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

Random_IDs <- function(seurat_obj){
  
  # Randomly sample 100 cell IDs from each cluster
  set.seed(1234) # Setting a seed for reproducibility
  cluster_levels <- levels(seurat_obj@active.ident)
  sampled_ids <-  c()
  for (i in cluster_levels) {
    cluster_ids <- Cells(subset(seurat_obj, idents = i))
    if (length(cluster_ids) >= 100) {
      sampled_ids <- c(sampled_ids,sample(cluster_ids, 100)) # Sample 100 IDs if there are enough
    } else {
      sampled_ids <- c(sampled_ids,cluster_ids) # Return all IDs if fewer than 100
    }
  }
  
  # Display the final vector of sampled IDs
  length(sampled_ids)
  return(sampled_ids)
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
    color = colorRampPalette(c("palegoldenrod", "navy"))(50),
    border_color = "black",
    angle_col = 45,            # Rotate x-axis labels to appear right-side up
    fontsize_row = 12,        # Adjust font size for row labels
    fontsize_col = 12,        # Adjust font size for column labels
    cellwidth = 30,           # Set a fixed cell width for square cells
    cellheight = 30,          # Set a fixed cell height for square cells
    legend_breaks = seq(0, 100, 20),  # Set legend to range from 0 to 100%
    legend_labels = paste0(seq(0, 100, 20), "%"), # Format legend labels as percentages
    display_numbers = F,   # Show cell values
    number_format = "%.1f%%", # Format displayed numbers as percentages
    labels_row = rownames(heatmap_matrix_percent), # Ensure y-axis labels are displayed
    labels_col = colnames(heatmap_matrix_percent),  # Ensure x-axis labels are displayed
    title = F, legend = F, show_rownames = F
  )
  ggsave(filename,plot = p,height = 4, width = 8)
  return(p)
}

up_down_volcano <- function(DEG.df,selectLab,italicLab,up_label,down_label,up_x, down_x,volcano_name){
  keyvals <- ifelse(
    DEG.df$avg_log2FC < 0 & DEG.df$p_val_adj < 0.05, 'blue',
    ifelse(DEG.df$avg_log2FC > 0& DEG.df$p_val_adj < 0.05, 'orange', 
           'black'))
  keyvals[is.na(keyvals)] <- 'grey'
  names(keyvals)[keyvals == 'orange'] <- 'Upregulated'
  names(keyvals)[keyvals == 'grey'] <- 'Not Significant'
  names(keyvals)[keyvals == 'blue'] <- 'Downregulated'
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
    annotate("text", x = up_x, y = 0, label = up_label,  size = 5, color = "orange", hjust = 0)+
    annotate("text", x = down_x, y = 0, label = down_label,  size = 5, color = "blue",hjust = 0)
  
  ggsave(volcano_name,plot = p, width = 6, height = 5)
  p
  

####Functions borrowed from Xiao Li and Yi Zhao

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

GetDoublet <- function(srt_obj, doublet_rate, dimN.var.toal){
  ## Scrublet (run via reticulate)
  mtx <- srt_obj@assays$RNA@counts
  mtx <- t(mtx)
  scrub_model <- scr$Scrublet(mtx, expected_doublet_rate = doublet_rate)
  rst <- scrub_model$scrub_doublets(min_gene_variability_pctl = dimN.var.toal*100,
                                    n_prin_comps = 30L,
                                    min_counts = 2, min_cells = 3)
  rst[[2]] <- scrub_model$call_doublets(threshold = 0.25) ## adjusted based on histogram
  sc_doublets <- Cells(srt_obj)[rst[[2]]]
  sc_singlets <- Cells(srt_obj)[!rst[[2]]]
  srt_obj$Scrublet_doublet <- 'Singlet'
  srt_obj$Scrublet_doublet[rst[[2]]] <- 'Doublet'
  Scrublet <- rst[[1]]
  names(Scrublet) <- Cells(srt_obj)
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


}

