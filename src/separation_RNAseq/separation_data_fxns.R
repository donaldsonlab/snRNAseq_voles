library(tidyverse)
library(ggplot2)
library(ggfortify)
library(ggrepel)
library(GGally)
library(wesanderson)
library(viridis)
library(readxl)
library(car)
library(gridExtra)
library(ez)
library(ggpubr)
library(DESeq2)
library(EnhancedVolcano)
library(RRHO2)
library(ggpubr)
library(UpSetR)
library(scales)
library(stringi)
library(SuperExactTest)
library(ComplexHeatmap)
library(colorRamp2)
library(ggridges)
library(gprofiler2)


get_norm_cts <- function(mod_inputs, want_metrics) {
  
  dds <- DESeqDataSetFromMatrix(countData = mod_inputs, colData = want_metrics, design = ~Cohort) # + ParentCode + PairingAgeCat + Cohort
  dds <- DESeq(dds)
  
  norm_cts <- counts(dds, normalized = TRUE) %>% as.data.frame()
  
  return(norm_cts)
}

transpose_norm_cts <- function(norm_cts) {
  #transform normalized counts
  mod_inputs_animal <- t(norm_cts) %>% as.data.frame() #take out scale in middle %>% scale() 
  mod_inputs_animal$Animal <- readr::parse_number(rownames(mod_inputs_animal))
  mod_inputs_meta <- merge(metrics, mod_inputs_animal, on = "Animal")
  
  return(mod_inputs_meta)
}

zscore_norm_cts <- function(mod_inputs_meta, all_genes) {

  norm_new <- mod_inputs_meta %>% select(c("Animal", "Cohort", "Separation", all_genes))
  norm_new[all_genes] <- scale(norm_new[all_genes])
  
  return(norm_new)
}

avg_zscored_cts <- function(norm_new) {

  norm_new_summ <- norm_new %>% group_by(Cohort) %>% dplyr::summarise_if(is.numeric, mean)
  
  norm_new_summ_mtx <- norm_new_summ[,-1] %>% as.matrix()
  rownames(norm_new_summ_mtx) <- norm_new_summ$Cohort
  
  return(norm_new_summ_mtx)
}

make_heatmap <- function(norm_new_summ_mtx, mod, save = FALSE) {
  genes_in_mod <- mod_genes$Gene[mod_genes$Module==mod]
  genes_in_mod <- intersect(genes_in_mod, rownames(want_inputs))
  
  mod_mtx <- norm_new_summ_mtx[,genes_in_mod]
  set.seed(123)
  heat_mod <- Heatmap(t(mod_mtx), column_order = c("SS48P", "SS4I", "OS48P", "OS4I"))
  order <- row_order(heat_mod)
  
  #get out dendrogram info
  og_dend <- row_dend(heat_mod)
  
  #customize colors
  col_fun = colorRamp2(c(min(mod_mtx), 0, max(mod_mtx)), c("#053061","white", "#67001F"))
  fname = paste0("output/heatmap_clusters_newzscore_excl2599_mod_", mod, ".pdf")
  
  if (save == TRUE) {
    pdf(fname, width = 6, height = 20)
    draw(Heatmap(t(mod_mtx), name = "foo", cluster_rows = og_dend, column_order = c("SS48P", "SS4I", "OS48P", "OS4I"), col = col_fun, split = 3)) #split = 8, 
    dev.off()
  }
  
  else {
    draw(Heatmap(t(mod_mtx), name = "foo", cluster_rows = og_dend, column_order = c("SS48P", "SS4I", "OS48P", "OS4I"), col = col_fun, split = 3),
         heatmap_width = unit(6, "in"), heatmap_height = unit(20, "in"))
  }
  
}

make_ridgelines <- function(norm_new_summ_mtx, norm_new, norm_cts, metrics, mod, save = FALSE) {
  #get module genes
  genes_in_mod <- mod_genes$Gene[mod_genes$Module==mod]
  genes_in_mod <- intersect(genes_in_mod, rownames(want_inputs))
  
  #make heatmap to get row order
  mod_mtx <- norm_new_summ_mtx[,genes_in_mod]
  set.seed(123)
  heat_mod <- Heatmap(t(mod_mtx), column_order = c("SS48P", "SS4I", "OS48P", "OS4I"))
  order <- row_order(heat_mod)
  
  norm_new <- norm_new[,c("Animal", "Cohort", genes_in_mod)]
  norm_cts$gene <- rownames(norm_cts)
  
  to_plot <- norm_new %>% pivot_longer(!c(Animal, Cohort), names_to = "gene", values_to = "count")
  to_plot$base_order <- 1:nrow(to_plot)
  to_plot <- to_plot %>% merge(metrics, on= "Animal") %>% merge(norm_cts, on = "gene") %>% arrange(factor(base_order, levels = order))
  ordered_genes <- unique(to_plot$gene)
  
  p <- to_plot %>%
    filter(Cohort %in% c("SS48P", "SS4I")) %>%
    mutate(gene_ord = factor(gene, levels = rev(ordered_genes))) %>%
    ggplot(aes(x = count, y = gene_ord, fill = Cohort)) +
    stat_density_ridges(
      aes(x = count, fill = Cohort, color = Cohort), 
      alpha = 0.6, scale = 0.9, size = 0.9,
      quantile_lines = TRUE, quantiles = 2,
      rel_min_height = 0.001
    ) +
    scale_fill_manual(values = c("SS4I" = "lightgreen", 
                                 "SS48P" = "darkgreen")) +
    scale_color_manual(values = c("SS4I" = "lightgreen", 
                                  "SS48P" = "darkgreen")) +
    xlim(c(-5, 5)) +
    # guides(fill = guide_legend(reverse = TRUE)) +
    # scale_x_continuous(expand = c(0,0)) +
    theme_classic()
  
  print(p)
  
  if (save == TRUE) {
    fname = paste0("output/ridgeline_SS_excl2599_mod", mod, ".pdf")
    ggsave(fname, p, width = 6, height = 20)
  }
  
  
  p <- to_plot %>%
    filter(Cohort %in% c("OS48P", "OS4I")) %>%
    mutate(gene_ord = factor(gene, levels = rev(ordered_genes))) %>%
    ggplot(aes(x = count, y = gene_ord, fill = Cohort)) +
    stat_density_ridges(
      aes(x = count, fill = Cohort, color = Cohort), 
      alpha = 0.6, scale = 0.9, size = 0.9,
      quantile_lines = TRUE, quantiles = 2,
      rel_min_height = 0.001
    ) +
    scale_fill_manual(values = c("OS4I" = "lightpink", 
                                 "OS48P" = "deeppink3")) +
    scale_color_manual(values = c("OS4I" = "lightpink", 
                                  "OS48P" = "deeppink3")) +
    xlim(c(-5, 5)) +
    # guides(fill = guide_legend(reverse = TRUE)) +
    # scale_x_continuous(expand = c(0,0)) +
    theme_classic()
  
  print(p)
  
  if (save == TRUE) {
    fname = paste0("output/ridgeline_OS_excl2599_mod", mod, ".pdf")
    ggsave(fname, p, width = 6, height = 20)
  }
  
}

make_clusters <- function(norm_new_summ_mtx, mod, save = FALSE) {
  set.seed(123)
  heat_mod <- Heatmap(t(norm_new_summ_mtx), column_order = c("SS48P", "SS4I", "OS48P", "OS4I"))
  draw(heat_mod)
  #get out row order from heatmap
  order <- row_order(heat_mod)
  #get out dendrogram info
  og_dend <- row_dend(heat_mod)
}

get_gene_change_scores <- function(mod_genes, metrics, norm_new, n_mods = 23) {
  all_gene_change <- data.frame()
  if (n_mods == "skip") {
    genes_mod <- mod_genes
    want_anis <- metrics$Animal[metrics$Cohort %in% c("SS48P", "SS4I", "OS48P", "OS4I")] %>% unfactor()
    want_idxs <- match(want_anis, metrics$Animal)
    
    want_inputs <- inputs[,want_idxs]
    want_metrics <- metrics %>% filter(Animal %in% want_anis)
    
    thismod_genes <- intersect(genes_mod, rownames(want_inputs))
    
    genes_mod_inputs <- genes_mod
    
    norm_summ <- norm_new %>% select(c("Animal", "Cohort", genes_mod_inputs)) %>% group_by(Cohort) %>% dplyr::summarise_if(is.numeric, median)
    rownames(norm_summ) <- norm_summ$Cohort
    
    # norm_summ_sep <- norm_new %>% select(c("Animal", "Separation", genes_mod_inputs)) %>% group_by(Separation) %>% dplyr::summarise_if(is.numeric, median)
    
    # norm_summ_mtx <- as.matrix(norm_summ[,-1])
    
    # norm_zscore <- scale(norm_summ_mtx)
    # rownames(norm_zscore) <- rownames(norm_summ)
    # norm_zscore <- t(norm_zscore)
    
    norm_summ_t <- t(norm_summ)
    norm_summ_t <- norm_summ_t[-1,]
    norm_summ_t <- as.data.frame(norm_summ_t)
    norm_summ_t2 <- lapply(norm_summ_t, as.numeric) %>% as.data.frame()
    rownames(norm_summ_t2) <- rownames(norm_summ_t)
    norm_summ_t2$OS_change_rel <- (norm_summ_t2$OS4I-norm_summ_t2$OS48P) #/(norm_summ_t2$OS4I+norm_summ_t2$OS48P)
    norm_summ_t2$abs_OS_change_rel <- abs(norm_summ_t2$OS_change_rel)
    norm_summ_t2$SS_change_rel <- (norm_summ_t2$SS4I-norm_summ_t2$SS48P) #/(norm_summ_t2$SS4I+norm_summ_t2$SS48P)
    norm_summ_t2$abs_SS_change_rel <- abs(norm_summ_t2$SS_change_rel)

    #try instead to separate by up and down regulated genes
    norm_summ_t2$OS_ID_change <- if_else(norm_summ_t2$OS_change_rel > 0, "positive", "negative")
    norm_summ_t2$SS_ID_change <- if_else(norm_summ_t2$SS_change_rel > 0, "positive", "negative")
    
    norm_summ_t2$Module <- "random_genes"
    all_gene_change <- rbind(all_gene_change, norm_summ_t2)
    
    
  }
  else {
    for (mod in 1:n_mods) {
      genes_mod <- mod_genes$Gene[mod_genes$Module==mod]
      
      want_anis <- metrics$Animal[metrics$Cohort %in% c("SS48P", "SS4I", "OS48P", "OS4I")] %>% unfactor()
      want_idxs <- match(want_anis, metrics$Animal)
      
      want_inputs <- inputs[,want_idxs]
      want_metrics <- metrics %>% filter(Animal %in% want_anis)
      
      thismod_genes <- intersect(genes_mod, rownames(want_inputs))
      
      genes_mod_inputs <- intersect(genes_mod, all_genes)
      
      norm_summ <- norm_new %>% select(c("Animal", "Cohort", genes_mod_inputs)) %>% group_by(Cohort) %>% dplyr::summarise_if(is.numeric, median)
      rownames(norm_summ) <- norm_summ$Cohort
      
      norm_summ_mtx <- as.matrix(norm_summ[,-1])
      
      norm_zscore <- scale(norm_summ_mtx)
      rownames(norm_zscore) <- rownames(norm_summ)
      norm_zscore <- t(norm_zscore)
      
      norm_summ_t <- t(norm_summ)
      norm_summ_t <- norm_summ_t[-1,]
      norm_summ_t <- as.data.frame(norm_summ_t)
      norm_summ_t2 <- lapply(norm_summ_t, as.numeric) %>% as.data.frame()
      rownames(norm_summ_t2) <- rownames(norm_summ_t)
      norm_summ_t2$OS_change_rel <- (norm_summ_t2$OS4I-norm_summ_t2$OS48P) #/(norm_summ_t2$OS4I+norm_summ_t2$OS48P)
      norm_summ_t2$abs_OS_change_rel <- abs(norm_summ_t2$OS_change_rel)
      norm_summ_t2$SS_change_rel <- (norm_summ_t2$SS4I-norm_summ_t2$SS48P) #/(norm_summ_t2$SS4I+norm_summ_t2$SS48P)
      norm_summ_t2$abs_SS_change_rel <- abs(norm_summ_t2$SS_change_rel)
      
      
      #try instead to separate by up and down regulated genes
      norm_summ_t2$OS_ID_change <- if_else(norm_summ_t2$OS_change_rel > 0, "positive", "negative")
      norm_summ_t2$SS_ID_change <- if_else(norm_summ_t2$SS_change_rel > 0, "positive", "negative")
      
      norm_summ_t2$Module <- mod
      all_gene_change <- rbind(all_gene_change, norm_summ_t2)
  }
 
    
  }
  return(all_gene_change)
  
}

ridgeline_from_hm_clusts <- function(norm_new_summ_mtx, all_gene_change, mod, save = FALSE) {
  #get module genes
  genes_in_mod <- mod_genes$Gene[mod_genes$Module==mod]
  genes_in_mod <- intersect(genes_in_mod, rownames(want_inputs))
  
  mat <- norm_new_summ_mtx[,genes_in_mod]
  mat <- t(mat)
  
  heat_mod = Heatmap(mat, column_order = c("SS48P", "SS4I", "OS48P", "OS4I"))
  order <- row_order(heat_mod)
  og_dend <- row_dend(heat_mod)
  
  col_fun = colorRamp2(c(min(mat), 0, max(mat)), c("#053061","white", "#67001F"))
  
  HM <- Heatmap(mat, name = "foo", cluster_rows = og_dend, column_order = c("SS48P", "SS4I", "OS48P", "OS4I"), col = col_fun, split = 3) 
  HM <- draw(HM)
  
  r.dend <- row_dend(HM)  #Extract row dendrogram
  rcl.list <- row_order(HM)  #Extract clusters (output is a list)
  
  lapply(rcl.list, function(x) length(x))  #check/confirm size clusters
  
  # loop to extract genes for each cluster.
  for (i in 1:length(row_order(HM))){
    if (i == 1) {
      clu <- t(t(row.names(mat[row_order(HM)[[i]],])))
      out <- cbind(clu, paste("cluster", i, sep=""))
      colnames(out) <- c("gene", "Cluster")
    } else {
      clu <- t(t(row.names(mat[row_order(HM)[[i]],])))
      clu <- cbind(clu, paste("cluster", i, sep=""))
      out <- rbind(out, clu)
    }
  }
  
  out <- out %>% as.data.frame()
  
  plots <- list()
  for (c in unique(out$Cluster)) {
    c_genes <- out$gene[out$Cluster == c]
    c_df <- all_gene_change[c_genes,]
    c_df$gene <- rownames(c_df)
    
    c_df2 <- c_df %>% pivot_longer(c(OS_change_rel, SS_change_rel), names_to = "SSOS", values_to = "rel_change")
    
    p <- c_df2 %>% ggplot(aes(x = rel_change, y = 0, fill = SSOS, color = SSOS, group = SSOS, alpha = 0.6)) + 
      stat_density_ridges(quantile_lines = TRUE, quantiles = 2) + 
      xlim(-3, 3) +
      ylim(0, 3) +
      geom_vline(xintercept = 0) + 
      ggtitle(paste0("Module-", mod, " ", c, " n_genes = ", nrow(c_df2)/2, " % genes = ", round(nrow(c_df)/nrow(out), 3))) +
      theme_classic()
    plots[[c]] <- p
  }
  # p2 <- do.call(grid.arrange, plots)
  p2 <- do.call(arrangeGrob, plots)
  # print(p2)
  
  if (save == TRUE) {
    fname <- paste0("output/hm_clust_ridges_mod", mod, ".pdf")
    ggsave(fname, p2, width = 8, height = 10)
  }

  
}

get_hm_clusters <- function(norm_new_summ_mtx, all_gene_change, mod) {
  hm_clust_df <- data.frame()
  for (mod in 1:23) {
    #get module genes
    genes_in_mod <- mod_genes$Gene[mod_genes$Module==mod]
    genes_in_mod <- intersect(genes_in_mod, rownames(want_inputs))
    
    mat <- norm_new_summ_mtx[,genes_in_mod]
    mat <- t(mat)
    
    heat_mod = Heatmap(mat, column_order = c("SS48P", "SS4I", "OS48P", "OS4I"))
    order <- row_order(heat_mod)
    og_dend <- row_dend(heat_mod)
    
    col_fun = colorRamp2(c(min(mat), 0, max(mat)), c("#053061","white", "#67001F"))
    
    HM <- Heatmap(mat, name = "foo", cluster_rows = og_dend, column_order = c("SS48P", "SS4I", "OS48P", "OS4I"), col = col_fun, split = 3) 
    # HM <- draw(HM)
    
    r.dend <- row_dend(HM)  #Extract row dendrogram
    rcl.list <- row_order(HM)  #Extract clusters (output is a list)
    
    lapply(rcl.list, function(x) length(x))  #check/confirm size clusters
    
    # loop to extract genes for each cluster.
    for (i in 1:length(row_order(HM))){
      if (i == 1) {
        clu <- t(t(row.names(mat[row_order(HM)[[i]],])))
        out <- cbind(clu, paste("cluster", i, sep=""))
        colnames(out) <- c("gene", "Cluster")
      } else {
        clu <- t(t(row.names(mat[row_order(HM)[[i]],])))
        clu <- cbind(clu, paste("cluster", i, sep=""))
        out <- rbind(out, clu)
      }
    }
    
    out <- out %>% as.data.frame()
    out$Module <- mod
    hm_clust_df <- rbind(hm_clust_df, out)
  }
  
  return(hm_clust_df)
}

run_clust_go <- function(hm_clust_df, mod) {
  all_go_df <- data.frame()
  for (clust in unique(hm_clust_df$Cluster)) {
    #for mod6 up/up genes
    glist <- hm_clust_df$gene[hm_clust_df$Module == mod & hm_clust_df$Cluster == clust]

    gostres <- gost(query = glist, 
                       organism = "mochrogaster", 
                       domain_scope = "known",
                       evcodes = TRUE,
                       correction_method = "fdr")
    gos_df <- gostres$result %>% filter(source == "GO:BP")
    if (nrow(gos_df) > 0) {
      gos_df$Module <- mod
      gos_df$Cluster <- clust
      all_go_df <- rbind(all_go_df, gos_df)
    }

  }
  return(all_go_df)
}
