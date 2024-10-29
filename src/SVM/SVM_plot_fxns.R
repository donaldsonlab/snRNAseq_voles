library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(forcats)
library(ggpubr)
library(glmmTMB)
library(DHARMa)
library(gridExtra)

process_df_withself <- function(df) {
  df$group <- if_else(df$Var1 == df$y_pred, "self", 
                             if_else((paste(df$Var1, df$y_pred, sep = "x") %in% unique(metadata$pair)), "partner",
                                     if_else((paste(df$y_pred, df$Var1, sep = "x") %in% unique(metadata$pair)), "partner", "other")))
  
  
  summary_classes <- df %>% group_by(Cluster, Var1) %>% dplyr::summarize(sum_classifications = sum(Freq))
  
  
  
  df2 <- merge(df, summary_classes, on = c("Cluster", "Var1"))
  
  df2$norm_classes <- df2$Freq/df2$sum_classifications
  
  
  df2$animal <- df2$Var1
  df2 <- merge(df2, metadata, on = "animal")
  
  return(df2)
}

process_df_excludeself <- function(df) {
  df$group <- if_else(df$Var1 == df$y_pred, "self", 
                                if_else((paste(df$Var1, df$y_pred, sep = "x") %in% unique(metadata$pair)), "partner",
                                        if_else((paste(df$y_pred, df$Var1, sep = "x") %in% unique(metadata$pair)), "partner", "other")))
  
  df <- df %>% filter(group != "self")
  
  summary_classes <- df %>% group_by(Cluster, Var1) %>% dplyr::summarize(sum_classifications = sum(Freq))
  
  
  
  df2 <- merge(df, summary_classes, on = c("Cluster", "Var1"))
  
  df2$norm_classes <- df2$Freq/df2$sum_classifications
  
  
  df2$animal <- df2$Var1
  df2 <- merge(df2, metadata, on = "animal")
  
  return(df2)
}


plot_ridgeline <- function(df, grouping = FALSE, save = FALSE) {
  #define colors and clusters
  f.levels <- c('Drd1Pdyn', 'Drd1PdynOprm1', 'Drd1Penk', 'Drd2Penk', 'Drd2NoPenk', 
                'GABAergicNeurons', 'Dlx2ImmatureNeurons', 'SstNpyInterneurons', 
                'PvalbInterneurons', 'CholinergicInterneurons', 'MatureOligos', 
                'ImmatureOligos', 'Astrocytes', 'Microglia')
  
  clust_cols <- c("Drd1Pdyn" = "#304880",
                  "Drd1PdynOprm1" = "#4898D0",
                  "Drd1Penk" = "#88C0F0",
                  "Drd2Penk" = "#9060A0",
                  "Drd2NoPenk" = "#C078C0",
                  "GABAergicNeurons" = "#389078",
                  "Dlx2ImmatureNeurons" = "#70A830",
                  "SstNpyInterneurons" = "#98D048",
                  "PvalbInterneurons" = "#60D0A0",
                  "CholinergicInterneurons" = "#80E8C0",
                  "MatureOligos" = "#783028",
                  "ImmatureOligos" = "#B82820",
                  "Astrocytes" = "#D04058",
                  "Microglia" = "#F87878")
  
  if (grouping != FALSE) {

    # grouping <- enquo(grouping)
    df <- df[df[[grouping]]==type,]
    # df <- df %>% filter(!!grouping == type)
  }
  
  p <- df %>%
    mutate(ClustFct = fct_rev(factor(Cluster, levels = f.levels))) %>%
    ggplot(aes(y = ClustFct)) +
    stat_density_ridges(
      aes(x = norm_classes, fill = ClustFct, color = group), 
      alpha = 0.8, from = 0, to = 1, scale = 1.5, size = 1.5,
      quantile_lines = TRUE, quantiles = 2,
      # # if you want points on ridgeline plots
      # jittered_points = TRUE,
      # position = position_points_jitter(width = 0, height = 0),
      # point_shape = '|', point_size = 5, point_alpha = 1, alpha = 0.7
      ) +
    geom_vline(xintercept = 0.026, color = "red", linetype = "dashed", size = 1.5) + #for chance line
    scale_color_manual(values = c("self"="black",
                                  "partner"="gray50",
                                  "other"="gainsboro")) +
    scale_fill_manual(values = clust_cols) + 
    guides(fill = guide_legend(reverse = TRUE)) +
    scale_x_continuous(expand = c(0,0)) +
    theme_classic()
  
  print(p)
  
  if (save == TRUE) {
    if ("self" %in% df$group) {
      ggsave("ridgeline_withself_chanceline_points.pdf", p, width = 8, height = 8, units = "in", device = "pdf", bg = "white")
    }
    else {
      ggsave("ridgeline_excludeself_chanceline_points.pdf", p, width = 8, height = 8, units = "in", device = "pdf", bg = "white")
    }
  }
  
}


widen_df_withself <- function(df) {
  
  with.self.wide <- df %>% pivot_wider(id_cols = c("Var1", "Cluster"), names_from = "group", values_from="Freq")
  
  with.self.wide$other <- lapply(with.self.wide$other, sum) %>% as.numeric()
  with.self.wide$partner <- as.numeric(with.self.wide$partner)
  with.self.wide$self <- as.numeric(with.self.wide$self)
  
  with.self.wide$self.pct <- with.self.wide$self/(with.self.wide$self + with.self.wide$other + with.self.wide$partner)
  with.self.wide$partner.pct <- with.self.wide$partner/(with.self.wide$partner + with.self.wide$other)
  
  with.self.wide <- with.self.wide %>% rename("animal" = "Var1")
  
  with_summ <- with.self.wide %>% group_by(Cluster) %>% summarise(mean_self = mean(na.omit(self.pct)), mean_partner = mean(na.omit(partner.pct)))
  
  with.self.wide <- with.self.wide %>% merge(data, on = "animal")
  
  return(with.self.wide)
  
}

widen_df_excludeself <- function(df) {
  
  exclude.self.wide <- df %>% pivot_wider(id_cols = c("Var1", "Cluster"), names_from = "group", values_from="Freq")
  
  exclude.self.wide$other <- lapply(exclude.self.wide$other, sum) %>% as.numeric()
  exclude.self.wide$partner <- as.numeric(exclude.self.wide$partner)
  
  exclude.self.wide$partner.pct <- exclude.self.wide$partner/(exclude.self.wide$partner + exclude.self.wide$other)
  
  exclude.self.wide <- exclude.self.wide %>% rename("animal" = "Var1") %>% merge(data, on = "animal")
  
  return(exclude.self.wide)
}


plot_dotcloud <- function(df) {

  p <- ggplot(df, aes(x = group_type, y = norm_classes*100, color = group_type, fill = group_type, alpha = 0.8)) +
    # geom_violin() +
    # geom_point(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.9)) +
    geom_jitter(size = 0.5, alpha = 1, width = 0.2) +
    scale_color_manual(values = c("FM_other" = "gray69",
                                  "FM_partner"="coral",
                                  "FF_other"="gray69",
                                  "FF_partner"="slateblue",
                                  "MM_other"="gray69",
                                  "MM_partner" = "#1B9E77")) +
    geom_abline(intercept = (1/37)*100, slope = 0, color = "red", linetype = "dashed") +
    xlab("Group") +
    ylab("% of Classifications") +
    theme_classic()
  
  p2 <- p + facet_wrap(~ Cluster, nrow = 5)

  
  print(p2)

  
}

plot_SVM_corrs <- function(df, save = FALSE) {
  all_plots <- list()

  for (clust in f.levels) {
    clust_df <- df %>% filter(Cluster == clust) %>% merge(metadata, on = "animal")
    clust_wide <- clust_df %>% pivot_wider(id_cols = c(pair, pair_type), names_from = color, values_from = partner.pct)
    
    cor <- rcorr(clust_wide$O, clust_wide$B, type = "spearman")
    lab <- paste("R = ", round(cor$r[2], digits = 4), " p = ", round(cor$P[2], digits = 6), sep = "")
    
    
    p <- ggplot(clust_wide, aes(x = O*100, y = B*100)) +
      geom_smooth(method = "lm", color = "black") +
      geom_point(aes(color = pair_type)) +
      ggtitle(paste(clust, lab)) +
      xlab("Partner 1 % Partner Classifications") +
      ylab("Partner 2 % Partner Classifications") +
      scale_color_manual(values = c("FM" = "coral",
                                    "FF"="slateblue",
                                    "MM"="#1B9E77")) +
      coord_fixed(xlim = c(0, 100), ylim = c(0, 100)) +
      theme_classic() +
      theme(plot.title = element_text(size=8), axis.title = element_text(size = 8))
    
    print(p)
    
    all_plots[[clust]] <- p
    
  }
  
  all_plt <- grid.arrange(grobs = all_plots, nrow = 5)
  
  if (save == TRUE) {
    setwd("output/")
    ggsave("all_partner_corrs_200cells.pdf", all_plt, width = 12, height = 12, units = "in")
  }

}

GLM_SVM <- function(df) {
  svm_stats <- data.frame()
  df$animal <- as.character(df$animal)
  
  for (clust in f.levels) {
    clust_df <- df %>% filter(Cluster == clust)
    
    fit <- glmmTMB(norm_classes~group, data = clust_df, family = ordbeta(link = "logit")) #, ziformula=~group
    # print(summary(fit))

    ## simulate and plot residuals if you want
    # simres <- simulateResiduals(fit)
    # plot(simres, title = clust)
    
    # diagnose(fit)
    
    EMM <- emmeans(fit, ~ group)
    coef <- contrast(EMM, "pairwise")
    coef2 <- summary(coef, adjust = "fdr")
    # print(coef2)
    mini.df <- coef2 %>% as.data.frame()
    mini.df$Cluster <- clust

    svm_stats <- rbind(svm_stats, mini.df)
    
  }

  return(svm_stats)
}

GLM_SVM_pairtype <- function(df) {
  
  all_contrasts <- data.frame(contrast = character(), 
                              estimate = numeric(), 
                              SE = numeric(), 
                              df = numeric(), 
                              t.ratio = numeric(), 
                              p.value = numeric(), 
                              Module = character())
  
  for (clust in f.levels) {
    clust_df <- df %>% filter(Cluster == clust)
    # form <- formula(paste0(mod, "~pair_type*group"))
    fit <- glmmTMB(norm_classes~pair_type*group, data = clust_df, family = ordbeta(link = "logit"))
    # print(summary(fit))
    
    # simres <- simulateResiduals(fit)
    # plot(simres, title = clust)
    
    EMM <- emmeans(fit, ~ pair_type*group)
    # print(summary(pairwise_comparisons))
    coef <- contrast(EMM, "pairwise")[c(3, 8, 12, 13, 14, 15)]
    coef2 <- summary(coef, adjust = "fdr")
    # print(coef2)
    mini.df <- coef2 %>% as.data.frame()
    mini.df$Cluster <- clust
    
    all_contrasts <- rbind(all_contrasts, mini.df)
    
  }
  
  return(all_contrasts)
}