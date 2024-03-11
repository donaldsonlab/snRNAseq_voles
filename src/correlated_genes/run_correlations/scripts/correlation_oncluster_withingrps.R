library(dplyr)
library(tidyr)
library(Hmisc)

args = commandArgs(trailingOnly=TRUE)
clust = args[1]
gene = args[2]
outdir = args[3]

#if gene is one of those weird genes that starts with a number
if (grepl("^[[:digit:]]+", gene) == TRUE) {
  gene <- paste0("X", gene)
}

#import data of average expression per cluster per animal
avg_expr <- read.csv("../indir/avg_SCT_cts_per_cluster.csv")
metadata <- read.csv("../indir/seq_beh_metadata.csv")
metadata <- metadata %>% filter(pair != "4918x4967")

#filter dataframe to be only one cluster
# clust <- "ImmatureOligos"

clust_df <- avg_expr %>% filter(cluster == clust)
clust_df <- clust_df %>% merge(metadata, on = "animal") %>% filter(animal != "4918")


#pick gene
# gene <- "Lhfpl3"
gene_clust_df <- clust_df %>% select("animal", gene)


#find true pairs
pairs <- unique(clust_df$pair)
pairs <- pairs[pairs != "4918x4967"]

#find all combos of not-true pairs
all_animals <- metadata$animal
all_animals <- all_animals[all_animals != c("4918", "4967")]
all_animals <- as.character(all_animals)

#to find true correlation between partners
gene_cts_df <- clust_df %>% pivot_wider(id_cols = pair, 
                                       names_from = color, 
                                       values_from = gene)


cor <- rcorr(gene_cts_df$O, gene_cts_df$B, type = "spearman")
Rho_true <- cor$r[2]
p_val_true <- cor$P[2]


shuffled_df <- data.frame(Cluster = character(), 
                          gene = character(), 
                          iteration = character(), 
                          Rho = numeric(), 
                          p_val = numeric(),
                          mean_expr = numeric())

beh_cohorts <- unique(metadata$beh_cohort)
t1 <- Sys.time()
for (i in 1:1000) {

  ani1_list <- list()
  ani2_list <- list()
  for (coh in beh_cohorts) {
    coh_anis <- metadata$animal[metadata$beh_cohort == coh]
    coh_set_one <- sample(coh_anis, (length(coh_anis)/2))
    coh_set_two <- coh_anis[!coh_anis %in% coh_set_one] %>% sample()
    ani1_list <- ani1_list %>% append(coh_set_one)
    ani2_list <- ani2_list %>% append(coh_set_two)
  }
  # set_one <- sample(all_animals, 19)
  # set_two <- all_animals[!all_animals %in% set_one] %>% sample()
  
  iter_pairs <- data.frame(ani1 = unlist(ani1_list), ani2 = unlist(ani2_list))
  
  iter_pairs$pair <- paste0(iter_pairs$ani1, "x", iter_pairs$ani2)
  long_df <- iter_pairs %>% pivot_longer(names_to = "animal_num", cols = c("ani1", "ani2"))
  long_df <- long_df %>% rename(animal = value)
  long_df <- long_df %>% merge(gene_clust_df, on = "animal")
  
  
  gene_cts_df <- long_df %>% pivot_wider(id_cols = pair, 
                                         names_from = animal_num, 
                                         values_from = gene)
  cor <- rcorr(gene_cts_df$ani1, gene_cts_df$ani2, type = "spearman")
  
  mini_shuffled_df <- data.frame(Cluster = clust, 
                                 gene = gene, 
                                 iteration = i, 
                                 Rho = cor$r[2], 
                                 p_val = cor$P[2],
                                 mean_expr = mean(long_df[[gene]])) 
  shuffled_df <- rbind(shuffled_df, mini_shuffled_df)
  
  
}
t2 <- Sys.time()
print(t2-t1)
#get summary stats
#confidence interval
# Calculate the mean and standard error
l.model <- lm(Rho ~ 1, shuffled_df)

# Calculate the confidence interval
confint(l.model, level=0.95)


#manually calculate confidence interval
alpha = 0.05
degrees.freedom <- 1000-1
sterr = sd(shuffled_df$Rho)/sqrt(nrow(shuffled_df))
t_score <- qt(p=alpha/2, df=degrees.freedom, lower.tail=F)
margin_error <- t_score*sterr

summary_stats <- summarise(shuffled_df, 
                           n_samps = n(),
                           mean = mean(Rho), 
                           stdev = sd(Rho), 
                           sterr = sd(Rho)/sqrt(n()),
                           Rho_true = Rho_true,
			   pval_true = p_val_true)
summary_stats$mean_expr <- mean(long_df[[gene]])
summary_stats$gene <- gene
summary_stats$cluster <- clust
print(summary_stats)

# hist(shuffled_df$Rho, breaks=50)
# abline(v = Rho_true, col = 'red', lwd = 2, lty = 'dashed')
# abline(v = (mean(shuffled_df$Rho) + sd(shuffled_df$Rho)*2), col = "blue", lwd=2, lty = "dashed")

fname <- paste0(outdir, clust, "_", gene, "_correlations.txt")
write.csv(summary_stats, file=fname)

