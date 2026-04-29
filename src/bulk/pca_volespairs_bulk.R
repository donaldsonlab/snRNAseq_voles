setwd("/Users/allenma/voles/oct25seq/diffexp/")
#install.packages("renv")
#BiocManager::install("clusterProfiler")

#renv::init(force = TRUE) 
###############################################################################
# 00. SETUP
###############################################################################
###############################################################################

library("RColorBrewer")
library(tidyverse)
library(ggplot2)
library(DESeq2)
library("pheatmap")
library("vsn")
library(ComplexHeatmap)
library(UpSetR)
library(FSA)
library(timereg)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rstatix)




###############################################################################
# 01. INPUT FILES AND CONSTANTS
###############################################################################

outdir <- "/Users/allenma/voles/oct25seq/diffexp/"
indir <- "/Users/allenma/voles/tableoftabletables/"
gtfvole = paste0(indir, "Microtus_ochrogaster.MicOch1.0.115.gtf")
prematingfile <- paste0(indir,"premattingcage.csv")


###############################################################################
# 02. LOAD ANNOTATION AND DESEQ OBJECT
###############################################################################



prematingdf <- read.csv(prematingfile)

# Read the GTF file; adjust the file name as needed
gtf <- read.csv(gtfvole, sep = "\t", header = FALSE, comment.char = "#")

# Add column names as per GTF spec
colnames(gtf) <-  c("chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

gtf <- gtf %>%
  mutate(
    gene_id = str_extract(attributes, "(?<=gene_id )[^; ]+"),
    gene_name = str_extract(attributes, "(?<=gene_name )[^; ]+")
  )

head(gtf)

switchnamesdf <- 
  as.data.frame(gtf) |>
  dplyr::select(gene_id, gene_name) |>
  base::unique()

table(gtf$gene_name)

load(paste0(indir,"deseq2.dds.RData"))
ls()

colData(dds)


###############################################################################
# 03. LOAD AND PREPARE SAMPLE METADATA
###############################################################################

#dds on allenma@login.rc.colorado.edu:/pl/active/DonaldsonLab/Mary_Allen/voles2025/nfcoresrnaseq10092025/star_salmon/salmon.merged.gene_lengths.tsv .
genelengths = read.csv(paste0(indir,"salmon.merged.gene_lengths.tsv"), sep="\t")

ss = read.csv(paste0(indir,"samplesheet.csv"))

ss$fastq_1 <- basename(ss$fastq_1)
ss$fastq_2 <- basename(ss$fastq_2)




coldata <- read.csv(paste0(indir,"donaldson_vole_meta2025.csv"))
dim(coldata)
coldata <- coldata[match(colData(dds)$sample, coldata$samplename), ]
coldata_F <- coldata %>% filter(Sex=="F") %>% dplyr::select(samplename,Pair)
coldata_M <- coldata %>% filter(Sex=="M") %>% dplyr::select(samplename,Pair)
coldata_M <- coldata_M %>% filter(Pair %in% coldata_F$Pair) %>% dplyr::select(samplename,Pair)
mergepairs <- merge(coldata_F,coldata_M, by="Pair")
mergepairs$pairname <- paste0(mergepairs$samplename.x,"x",mergepairs$samplename.y)
coldata <- merge(coldata, mergepairs, by="Pair", all = TRUE)
coldata <- coldata |> 
  mutate(still_together = if_else(!grepl('sep', Condition), Pair, "sep"))


coldata <- coldata[match(colData(dds)$sample, coldata$samplename), ]

colData(dds) <- cbind(colData(dds), coldata)
colData(dds)$samplename==colData(dds)$sample

colData(dds)$Condition <- factor(colData(dds)$Condition, levels=c("2wk_pair","4wk_sep","6wk_pair"))
colData(dds)$Sex <- factor(colData(dds)$Sex)
colData(dds)$Harvest.Date <- factor(colData(dds)$Harvest.Date)
colData(dds)$Cohort <- factor(colData(dds)$Cohort)
colData(dds)$Pair <- factor(colData(dds)$Pair)

coldata <- as.data.frame(colData(dds))

coldata_simple <- coldata %>% dplyr::select(sample, Pair, Sex, Cohort, Condition)
#library_name	title	library_strategy	organism	tissue	treatment	cohort_batch	molecule	single_or_paired-end	instrument model	description	fastq1 fastq2																															

#library_name vole number
#title vole number plus Condition (paring type)
#treatment Condition
#cohort_batch
#seq_batch
#gender
#Pair
#fastq1
#fastq2

ss <- merge(ss, coldata_simple, on="sample")

#write.csv(ss,paste0(outdir,"DataGEO.csv"))
###############################################################################
# 04. INITIAL DESEQ2 MODEL
###############################################################################




#colData(dds)$new_variable <- new_variable_data

design(dds)

design(dds) <- ~ Sex+Condition

dds <- DESeq(dds)
###############################################################################
# 05. QC BEFORE FILTERING
###############################################################################

len_mat <- as.matrix(genelengths[ , -(1:2)])
rownames(len_mat) <- genelengths$gene_id
#len_mat <- len_mat[rownames(counts), ]

countsdf <- as.matrix(counts(dds, normalize=FALSE))

stopifnot( all(colnames(countsdf) == colnames(len_mat)) )



# length-normalized rates per sample
rate <- countsdf / len_mat          # element-wise division

# avoid division by zero if any length==0
rate[!is.finite(rate)] <- 0

# scale each sample to 1e6
tpm <- t( t(rate) * 1e6 / colSums(rate) )

#write.csv(tpm, "tpm_all_samples.csv")
#How do samples look?

ntd <- normTransform(dds)

meanSdPlot(assay(ntd))

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Condition","Sex")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

vsd <- vst(dds, blind=FALSE)

sampleDists <- dist(t(assay(vsd)))


sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


plotPCA(vsd, intgroup=c("Condition", "Sex"))

#three samples look like crap (same three as the fastq qc said)
# 8178, 8184, 8583
pca <- prcomp(t(assay(vsd)))
distances <- sqrt(rowSums((pca$x[,1:2] - colMeans(pca$x[,1:2]))^2))
outliers <- which(distances > mean(distances) + 2 * sd(distances))


###############################################################################
# 06. REMOVE OUTLIERS AND 6-WEEK SAMPLES
###############################################################################


outlier_samples <- c("v8178", "v8184", "v8583")

tpm_nooutliers <- as.data.frame(tpm) %>% dplyr::select(-v8178, -v8184, -v8583)

#write.csv(tpm_nooutliers, "tpm_nooutlier_samples.csv")

remove6week  <- coldata %>% filter(Condition=="6wk_pair")

the6weeknames <- remove6week$samplename

# Subset to keep all other samples
dds <- dds[, !(colnames(dds) %in% outlier_samples)]
dds <- dds[, !(colnames(dds) %in% the6weeknames)]

design(dds) <- ~ Sex+Condition

colData(dds)$Condition <- factor(colData(dds)$Condition, levels=c("2wk_pair","4wk_sep"))

dds <- DESeq(dds)

#How do the samples look now that 6 week and bad quablity samples are removed?
coldata <- as.data.frame(colData(dds))
Mvoles <- coldata %>% filter(Sex=="M")
Mvoles_names <- unique(Mvoles$sample)
###############################################################################
# 07. QC AFTER FILTERING
###############################################################################


res <- results(dds, contrast = c("Condition", "4wk_sep", "2wk_pair"))
res <- as.data.frame(res)

res$gene_id <- rownames(res)

res_sig <- res %>% arrange(padj)  %>% filter(padj<0.1)

res_sig$gene_id <- rownames(res_sig)
res_sig <- merge(res_sig, switchnamesdf, by="gene_id")

res_sig <- res_sig %>% arrange(padj)

plotCounts(dds, gene="ENSMOCG00000000071", intgroup = "Condition")

write.csv(res_sig, paste0(outdir,"Deseq2_genes_sig_4wk_sep_vs_2wk_paired.csv"))

ntd <- normTransform(dds)

meanSdPlot(assay(ntd))


select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Condition","Sex")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

vsd <- vst(dds, blind=FALSE)

###############################################################################
# 08. PCA FOR CUSTOM PAIR DISTANCE ANALYSIS
###############################################################################


sampleDists <- dist(t(assay(vsd)))


sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotPCA(vsd, intgroup=c("Condition", "Sex"))
plotPCA(vsd, intgroup=c("Condition"))
plotPCA(vsd, intgroup=c("Sex"))
plotPCA(vsd, intgroup=c("Cohort"))
plotPCA(vsd, intgroup=c("Harvest.Date"))
plotPCA(vsd, intgroup=c("still_together"))


p <- plotPCA(vsd, intgroup=c("Condition"))  # returns a ggplot object [web:18]

#ggsave(
#  filename = paste0(outdir,"plotPCA_colorCondition.svg"),
#  plot     = p,
#  width    = 6,   # inches
#  height   = 5,
#  dpi      = 300  # ignored for pure vector but fine to set
#)

p <- plotPCA(vsd, intgroup=c("Pair"))  # returns a ggplot object [web:18]

#ggsave(
#  filename = paste0(outdir,"plotPCA_colorPair.svg"),
#  plot     = p,
#  width    = 6,   # inches
#  height   = 5,
#  dpi      = 300  # ignored for pure vector but fine to set
#)


#get all the gene counts not scaled
normcounts = counts(dds, normalize=TRUE)

ncdf_vole <- normcounts

#How variable are the genes 
row_stats <- t(apply(normcounts, 1, function(x) {
  # Calculate the min, max, 10th and 90th percentiles, 25th, 50th, and 75th percentiles,
  # count of 0s, and the standard deviation
  c(
    min = min(x),
    max = max(x),
    mean = mean(x),
    median = median(x),
    `10th_percentile` = quantile(x, probs = 0.1),
    `90th_percentile` = quantile(x, probs = 0.9),
    `25th_percentile` = quantile(x, probs = 0.25),
    `50th_percentile` = quantile(x, probs = 0.50),
    `75th_percentile` = quantile(x, probs = 0.75),
    count_zeros = sum(x == 0),
    std_dev = sd(x)  # Add standard deviation
  )
}))

row_stats <- as.data.frame(row_stats)
row_stats$percent_sample_0 <- row_stats$count_zeros/ncol(normcounts)
row_stats$CV <- row_stats$std_dev/row_stats$mean

row_stats$gene_id <- rownames(row_stats)
dim(row_stats)
row_stats <- merge(row_stats, switchnamesdf)
dim(row_stats)
#write.csv(row_stats, paste0(outdir,"vole_variation_per_gene.csv"))

#transpose for downstram anaysis
normcountsT <- as.data.frame(t(normcounts))
colnames(normcountsT) <- rownames(normcounts)
rownames(normcountsT) <- colnames(normcounts)

#write.csv(as.data.frame(normcounts), "normcounts_Deseq2_qualitysamples_remove6weeks.csv")

# Remamke the PCA with my own code so I can make 2 week and 4 week different shapes
#Using the most variable genes, like the plotpca does

ntop <- 500
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t( assay(vsd)[select, ] )

pca<-prcomp(mat)
pca <- as.data.frame(pca$x)

pca_with_meta <- pca

metasubset <- coldata %>% dplyr::select(samplename, Pair, Condition)
pca_with_meta$samplename <- rownames(pca_with_meta)

pca_with_meta<- merge(pca_with_meta, metasubset, by="samplename")


ggplot(pca_with_meta, aes(x=PC1, y=PC2, color=Pair, size=10))+geom_point(aes(shape = Condition))+theme_classic()


p<- ggplot(pca_with_meta,
           aes(x = PC1, y = PC2, color = Pair)) +
  geom_point(aes(shape = Condition), size = 3, stroke = 1.1) +
  scale_shape_manual(values = c(17, 21)) +  
  labs(
    x = "PC1 (66% of variance)",
    y = "PC2 (4% of variance)"
  ) +
  theme_classic()

p

#ggsave(paste0(outdir,"pcatop500_shapeandcolor.svg"), plot = p, width = 6, height = 6, dpi = 300)



###############################################################################
# 09. REAL VS SCRAMBLED PAIR DISTANCES
###############################################################################


#distance between every animal and every other animal

distance_matrix <- dist(pca[, c("PC1")])
distance_matrix_full <- as.matrix(distance_matrix)
colnames(distance_matrix_full) <-rownames(pca)
rownames(distance_matrix_full) <-rownames(pca)


#scramble the samples

coldata <- as.data.frame(colData(dds))

stdf <- coldata %>% dplyr::select(pairname, Condition) %>% unique()
stdf <- stdf %>%
  mutate(sample1 = str_extract(pairname, "^[^x]+"),
         sample2 = str_extract(pairname, "(?<=x).+$"),
         dist_value = usedist::dist_get(distance_matrix_full, sample1, sample2))

stdf$which <- stdf$Condition
stdf$type <- stdf$Condition

stdftypes <- stdf

stdf_foroutput <- stdf

pc1_for_output <- pca %>% dplyr::select(PC1, PC2)
pc1_for_output$sample1 <- rownames(pc1_for_output)
pc1_for_output$sample2 <- rownames(pc1_for_output)
pc1_for_output$ani1_pc1 <- pc1_for_output$PC1
pc1_for_output$ani2_pc1 <- pc1_for_output$PC1
pc1_for_output$ani1_pc2 <- pc1_for_output$PC2
pc1_for_output$ani2_pc2 <- pc1_for_output$PC2

pc1_for_output_ani1 <- pc1_for_output %>% dplyr::select(sample1, ani1_pc1, ani1_pc2)
pc1_for_output_ani2 <- pc1_for_output %>% dplyr::select(sample2, ani2_pc1, ani2_pc2)

stdf_foroutput <- merge(stdf_foroutput, pc1_for_output_ani1, by="sample1")
stdf_foroutput <- merge(stdf_foroutput, pc1_for_output_ani2, by="sample2")
  

#write.csv(stdf_foroutput, paste0(outdir,"pc_distances_pairs.csv"))
scramblen = 10

for (i in seq(scramblen)){
  halfanimals <- coldata %>% filter(Sex=="F") %>% dplyr::select(Animal, Condition) 
  halfanimals$Animal <- sample(halfanimals$Animal)
  halfanimals$Condition <- sample(halfanimals$Condition)
  otherhalf = coldata %>% filter(! Animal %in% halfanimals$Animal)%>% dplyr::select(Animal)
  otherhalf$Animal <- sample(otherhalf$Animal)
  scramblepairs = cbind(halfanimals, otherhalf)

  colnames(scramblepairs) <- c("Animal", "Condition","Partner")
  scramblepairs$pairname <- paste0("v",scramblepairs$Animal,"xv", scramblepairs$Partner)
  scramblepairs <- scramblepairs %>% dplyr::select(pairname, Condition)


  stdfscrmble <- scramblepairs %>%
    mutate(sample1 = str_extract(pairname, "^[^x]+"),
         sample2 = str_extract(pairname, "(?<=x).+$"),
         dist_value = usedist::dist_get(distance_matrix_full, sample1, sample2))

  stdfscrmble$which = paste0("scramble", i)
  stdfscrmble$type = "scrambled_partmers"
  
  
  stdftypes <- rbind(stdftypes, stdfscrmble)

}

ggplot(stdftypes, aes(x = Condition, y = dist_value, color=type)) +
  geom_violin(trim = FALSE, position = position_dodge(width = 0.8)) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    size = 2,
    alpha = 0.6
  ) + geom_boxplot(width=0.1)+
  theme_classic()



ggplot(stdftypes, aes(x = Condition, y = dist_value, color=type)) +
  geom_violin(trim = FALSE, position = position_dodge(width = 0.8)) +
  geom_boxplot(width=0.1)+
  theme_classic()


dodge <- position_dodge(width = 0.8)


p<-ggplot(stdftypes,
          aes(x = type, y = dist_value, color = type)) +
  geom_violin(trim = FALSE, position = dodge) +
  geom_boxplot(
    aes(x = type, color = type),
    width = 0.2,
    position = dodge,
    outlier.shape = NA
  ) +
  geom_point(
    position = dodge,
    size = 2,
    alpha = 0.6
  )+
  theme_classic()

p

#ggsave(paste0(outdir,"distance_between_pairs_pcatop500.svg"), plot = p, width = 6, height = 4, dpi = 300)



table(stdftypes$pairname)

###############################################################################
# 10. PREMATING CAGE PAIR COMPARISONS
###############################################################################


prematingdf$samplename <- paste0("v",prematingdf$Animal)

prematingdf <- prematingdf %>% filter(samplename %in% coldata$samplename)

premating_pairs_df <- prematingdf %>%
  group_by(Cage.prepairing) %>%
  summarise(
    pairs = list(
      if (dplyr::n() >= 2) combn(samplename, 2, simplify = FALSE) else list()
    ),
    .groups = "drop"
  ) %>%
  unnest_longer(pairs) %>%
  filter(!map_lgl(pairs, is.null)) %>%
  mutate(
    sample1 = map_chr(pairs, 1),
    sample2 = map_chr(pairs, 2)
  ) %>%
  dplyr::select(Cage.prepairing, sample1, sample2)

premating_pairs_df$pairname <- paste0(premating_pairs_df$sample1, "x", premating_pairs_df$sample2)

premating_pairs_df_match_for_stdf <- premating_pairs_df %>% dplyr::select(pairname, sample1, sample2)

premating_pairs_df_match_for_stdf$Condition <- "precagemate"


premating_pairs_df_match_for_stdf <- premating_pairs_df_match_for_stdf %>%
  mutate(dist_value = usedist::dist_get(distance_matrix_full, sample1, sample2))

premating_pairs_df_match_for_stdf$which <- "precagepartners"
premating_pairs_df_match_for_stdf$type <- "precage_partners"

#stdftypes <- rbind(stdftypes, premating_pairs_df_match_for_stdf)

premating_pairs_df_match_for_stdf_gendered <-
  premating_pairs_df_match_for_stdf %>%
  mutate(
    type = ifelse(sample1 %in% Mvoles_names,
                  "precage_partners_M",
                  "precage_partners_F")
  )

stdftypes_not_gendered <- rbind(stdftypes, premating_pairs_df_match_for_stdf)

#first test they are normalaly distributed
shapiro.test(
  stdftypes_not_gendered$dist_value[stdftypes_not_gendered$type == "2wk_pair"]
)

shapiro.test(
  stdftypes_not_gendered$dist_value[stdftypes_not_gendered$type == "4wk_sep"]
)
shapiro.test(
  stdftypes_not_gendered$dist_value[stdftypes_not_gendered$type == "scrambled_partmers"]
)

shapiro.test(
  stdftypes_not_gendered$dist_value[stdftypes_not_gendered$type == "precage_partners"]
)



kruskal.test(dist_value ~ type, data = stdftypes_not_gendered)


stdftypes_not_gendered <- stdftypes_not_gendered %>%
  mutate(
    type2 = ifelse(
      type == "scrambled_partmers",  "scrambled_partmers",
      ifelse(
        type == "precage_partners", "precage_partners",
        "matepaired")))



rstatix::pairwise_wilcox_test(stdftypes_not_gendered, 
  dist_value ~ type2,
  p.adjust.method = "holm",
  detailed = TRUE
)


rstatix::pairwise_wilcox_test(stdftypes_not_gendered, 
                              dist_value ~ type,
                              p.adjust.method = "holm",
                              detailed = TRUE
)



stdftypes <- rbind(stdftypes, premating_pairs_df_match_for_stdf_gendered)

p<-ggplot(stdftypes,
          aes(x = type, y = dist_value, color = type)) +
  geom_violin(trim = FALSE, position = dodge) +
  geom_boxplot(
    aes(x = type, color = type),
    width = 0.2,
    position = dodge,
    outlier.shape = NA
  ) +
  geom_point(
    position = dodge,
    size = 2,
    alpha = 0.6
  )+
  theme_classic()+theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p

#ggsave(paste0(outdir,"distance_between_pairs_pcatop500_wprecage.svg"), plot = p, width = 6, height = 4, dpi = 300)


#closest precage
precagepartners <- stdftypes %>% filter(type=="precage_partners_M" | type=="precage_partners_F" )

pairtype2week <- coldata %>% filter(Condition=="2wk_pair")
pairtype2week <- c(unique(pairtype2week$samplename.x),unique(pairtype2week$samplename.y))
pairtype4weekseq <- coldata %>% filter(Condition=="4wk_sep")
pairtype4weekseq <- c(unique(pairtype4weekseq$samplename.x), unique(pairtype4weekseq$samplename.y))


precagepartners$Condition2week <- precagepartners$sample1 %in% pairtype2week
precagepartners$Mvole <- (precagepartners$sample1 %in% Mvoles_names) |
  (precagepartners$sample2 %in% Mvoles_names)

ggplot(precagepartners, aes(x=Condition2week, y=dist_value))+geom_violin()+geom_boxplot(width=0.1)+geom_jitter()
ggplot(precagepartners, aes(x=Mvole, y=dist_value))+geom_violin()+geom_boxplot(width=0.1)+geom_jitter()


 



###############################################################################
# 11. LOAD MODULE GENES
###############################################################################



mod_genes <- read.csv(paste0(indir,"new_clusts_hotspot-gene-modules.csv"))


whichmod = 6

onemodgenes <- mod_genes %>% filter(Module==whichmod)
dim(onemodgenes)

res_sig_modgenes <- res_sig %>% filter(gene_name %in% onemodgenes)
#none of the genes that change sig with pairtype are in the module 6

geneshalf <- switchnamesdf %>% filter(gene_name %in% onemodgenes$Gene)
genesotherhalf <- switchnamesdf %>% filter(gene_id %in% onemodgenes$Gene)
mod_switchnamesdf <- rbind(geneshalf, genesotherhalf)
dim(mod_switchnamesdf)

#are mod_genes different then other genes

modonly <- mod_genes %>% dplyr::select(Gene, Module)

row_stats <- merge(
  row_stats,
  modonly,
  by.x = "gene_name",
  by.y = "Gene",
  all.x = TRUE
)
row_stats$Module <- as.factor(row_stats$Module )

row_stats$log10mean <- log10(row_stats$mean)
row_stats$log10median <- log10(row_stats$median)
row_stats$trimmedfc <- row_stats$`90th_percentile.90%`/row_stats$`10th_percentile.10%`

table(row_stats$Module)

ggplot(row_stats, aes(x=Module, y=std_dev))+geom_violin()+coord_cartesian(ylim = c(0, 300))
ggplot(row_stats, aes(x=Module, y=CV))+geom_violin()

p <- ggplot(row_stats, aes(x=log10median,color=Module, y=CV))+geom_point()+facet_wrap(~ Module)+coord_cartesian(ylim = c(0, 3))

#ggsave(
#  filename = paste0(outdir,"CV_and_median_exp.svg"),
#  plot     = p,
#  width    = 6,   # inches
#  height   = 5,
#  dpi      = 300  # ignored for pure vector but fine to set
#)


ggplot(row_stats, aes(x=log10median,color=Module, y=trimmedfc))+geom_point()+facet_wrap(~ Module)+coord_cartesian(ylim = c(0, 10))

ggplot(row_stats, aes(x=CV,color=Module, y=trimmedfc))+geom_point()+facet_wrap(~ Module)+coord_cartesian(ylim = c(0, 15), xlim=c(0,3))


#ggsave(
#  filename = paste0(outdir,"trimmedfc_and_median_exp.svg"),
#  plot     = p,
#  width    = 6,   # inches
#  height   = 5,
#  dpi      = 300  # ignored for pure vector but fine to set
#)

#ggplot(row_stats, aes(x=Module, y=log10mean))+geom_violin()
#ggplot(row_stats, aes(x=Module, y=log10median))+geom_violin()
#ggplot(row_stats, aes(x=Module, y=trimmedfc))+geom_violin()+
#  coord_cartesian(ylim = c(0, 10))


vsddf <- as.data.frame(assay(vsd)) 

vsddf <- vsddf %>% filter(rownames(vsddf) %in% mod_switchnamesdf$gene_id)

norm_counts_modt <- normcountsT %>% dplyr::select(mod_switchnamesdf$gene_id)



modgenenames <-colnames(norm_counts_modt)

norm_counts_modt$Animal <- rownames(norm_counts_modt)
  

n_shuffles = 100

pairsetnames <- c("realpairs",seq(n_shuffles))

#collect real pairs
realpairs <- as.data.frame(unique(coldata$pairname))
colnames(realpairs) <- c("pair")
realpairs <- realpairs %>%
  separate(pair, into = c("Ani1", "Ani2"), sep = "x", remove = FALSE)
realpairs$Ani1 <- paste(realpairs$Ani1 , sep="")
realpairs$Ani2 <- paste( realpairs$Ani2 , sep="")

all_ani2 <- unique(realpairs$Ani2)

# build all combinations of Ani1 with every Ani2
all_combos <- expand.grid(
  Ani1 = realpairs$Ani1,
  Ani2 = all_ani2,
  stringsAsFactors = FALSE
)

all_combos$pair <- paste0(all_combos$Ani1, "x", all_combos$Ani2)

# remove the true pairs
fake_pairs <- all_combos %>% filter(!pair %in% realpairs$pair)

###############################################################################
# 12. GENE-LEVEL HELPER FUNCTIONS
###############################################################################

shufflepairsfunction <- function(realpairs){
  shuffle <- realpairs
  shuffle$Ani2 <- sample(shuffle$Ani2)
  shuffle$pair <- paste0(shuffle$Ani1, "x", shuffle$Ani2)
  shuffle <- shuffle %>% arrange(Ani1) 
  shuffle
}

pairlists <- c(
  list(realpairs),
  replicate(n_shuffles, shufflepairsfunction(realpairs), simplify = FALSE)
)


###############################################################################
# 13. RUN GENE-LEVEL CORRELATION ANALYSES
###############################################################################



calc_a_corr <- function(agenename, apair_long){
  cols_to_select <- c("Animal", "pairname", "Variable", agenename)
  apair_long_onegene <- apair_long %>% dplyr::select(all_of(cols_to_select))
  # Set the genename being used to a common name
  colnames(apair_long_onegene)[colnames(apair_long_onegene) == agenename] <- "expression"
  # Pivot to one row per pair, wide-style with Ani1 and Ani2 expression in two columns
  apair_wide <- apair_long_onegene %>%
    pivot_wider(
      id_cols = pairname,
      names_from = Variable,
      values_from = expression
    )
  # Calculate Spearman correlation between Ani1 and Ani2 for this gene across all pairs
  genecor = cor(apair_wide$Ani1, apair_wide$Ani2, method = "spearman", use = "pairwise.complete.obs")
  return(genecor)
}

corrs_for_genes <- function(pairset, norm_counts_modt, gene_names){
  # pairset: dataframe with Ani1 and Ani2 columns
  # norm_counts_modt: expression matrix/dataframe, rows are Animal, columns are genes
  # gene_names: vector of gene names to assess
  
  apairfinal <- pairset %>% base::unique()
  apairfinal <- as.data.frame(apairfinal)
  apairfinal$pairname <- paste(apairfinal$Ani1, "x", apairfinal$Ani2, sep = "")
  # Long format: each row is one animal from a pair
  apair_long <- apairfinal %>%
    pivot_longer(
      cols = c(Ani1, Ani2),
      names_to = "Variable",
      values_to = "Animal"
    )
  # Add expression values: make sure norm_counts_modt rownames are in a column called 'Animal'
  
  apair_long <- merge(apair_long, norm_counts_modt, by = "Animal", all.x = TRUE)
  
  # Calculate correlations for each gene
  gene_corrs <- sapply(gene_names, function(g) calc_a_corr(g, apair_long))
  return(as.data.frame(gene_corrs))
}

#calucating correlations
all_corrs <- list()  # Initialize an empty list

for (i in seq_along(pairlists)){
  pairset = pairlists[[i]]
  pairset_name = pairsetnames[i]
  corrs <- corrs_for_genes(pairset, norm_counts_modt, modgenenames)
  corrs <- as.data.frame(corrs)
  colnames(corrs) <- pairset_name   # Label the column for this pairset
  all_corrs[[i]] <- corrs  
  # Store in list
}


# Combine into one big table by columns
corrs_combined <- do.call(cbind, all_corrs)

corrs_long <- corrs_combined %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "pairset",
    values_to = "correlation"
  )

corrs_long <- corrs_long %>%
  mutate(
    typepairs = if_else(pairset == "realpairs", "realpairs", "randompairs")
  )

corrs_long$abs_corr <- abs(corrs_long$correlation)

# Ensure n is atomic and ungrouped
corrs_long_with_n <- corrs_long %>%
  group_by(gene, typepairs) %>%
  mutate(n = n()) %>%
  ungroup()

corrs_long_with_n <- merge(corrs_long_with_n, switchnamesdf, by.x="gene", by.y="gene_id")


# Data for violins/boxplots
multi_group <- corrs_long_with_n %>% filter(n >= 2)
# Data for singletons (points)
single_group <- corrs_long_with_n %>% filter(n == 1)

gene_order <- single_group %>%
  arrange(correlation) %>%   # Or arrange(desc(correlation)) for reverse order
  pull(gene_name)


# Set gene factor levels 
multi_group$gene <- factor(multi_group$gene_name, levels = gene_order)
single_group$gene <- factor(single_group$gene_name, levels = gene_order)

p<-ggplot() +
  geom_violin(
    data = multi_group,
    aes(
      y = gene, 
      x = correlation, 
      color = typepairs
    ),
    fill = NA
  ) +
  geom_boxplot(
    data = multi_group,
    aes(
      y = gene, 
      x = correlation, 
      color = typepairs
    ),
    width = 0.2, outlier.shape = NA
  ) +
  geom_point(
    data = single_group,
    aes(
      y = gene, 
      x = correlation, 
      color = typepairs
    ),
    size = 2
  ) +
  theme_minimal()

p

#function to find out the number of genes above a certain cuttoff of spearmans
percent_above_and_below <- function(corrs_long,cutoff){
  
  pabove <- corrs_long %>%
    group_by(pairset) %>%
    summarise(percent_above = mean(correlation > cutoff, na.rm = TRUE) * 100) %>%
    arrange(percent_above)
  
  pbelow<- corrs_long %>%
    group_by(pairset) %>%
    summarise(percent_below = mean(correlation < (-1*cutoff), na.rm = TRUE) * 100) %>%
    arrange(percent_below)
  
  pabove <- pabove %>%
    mutate(
      typepairs = if_else(pairset == "realpairs", "realpairs", "randompairs")
    )
  
  pbelow <- pbelow %>%
    mutate(
      typepairs = if_else(pairset == "realpairs", "realpairs", "randompairs")
    )
  
  real_value_above <- pabove %>%
    filter(typepairs == "realpairs") %>%
    pull(percent_above)
  
  real_value_below<- pbelow %>%
    filter(typepairs == "realpairs") %>%
    pull(percent_below)
  
  random_summary_above <- pabove %>%
    filter(typepairs == "randompairs") %>%
    summarise(
      mean_randompairs_above = mean(percent_above),
      median_randompairs_above = median(percent_above),
      sd_randompairs_above = sd(percent_above)
    )
  
  random_summary_below <- pbelow %>%
    filter(typepairs == "randompairs") %>%
    summarise(
      mean_randompairs_below = mean(percent_below),
      median_randompairs_below = median(percent_below),
      sd_randompairs_below = sd(percent_below)
    )
  random_summary_above <- as.data.frame(random_summary_above)
  random_summary_below <- as.data.frame(random_summary_below)
  
  
  random_summary<- cbind(random_summary_above, random_summary_below)
  
  random_summary$real_value_above <- real_value_above
  random_summary$real_value_below <- real_value_below
  random_summary$cuttoff <- cutoff
  return(random_summary)
}

#calcuating genes above a cuttoff of spearmans
# List of cutoffs
cutoffs <- seq(0.1, 0.6, by = 0.05)

# Apply function to each cutoff and combine results
random_summary <- map_dfr(cutoffs, ~percent_above_and_below(corrs_long, .x))

#write.csv(random_summary, file=paste(outdir,"mod_",whichmod,"_bulk_percent_genes_above_corr.csv",sep=""))


ggplot(random_summary, aes(
  x = cuttoff, 
  y = mean_randompairs_above
)) +
  geom_point() +
  geom_errorbar(aes(
    ymin = mean_randompairs_above - sd_randompairs_above,
    ymax = mean_randompairs_above + sd_randompairs_above
  ), width = 0.01) +
  geom_point(aes(y = real_value_above), color = "blue", size = 3) +
  ylim(0, 60)+
  xlim(0, 0.6)+
  ylab("percent of genes above spearman correlation of X")+
  xlab("spearman correlation")

ggplot(random_summary, aes(
  x = cuttoff, 
  y = mean_randompairs_above
)) +
  geom_ribbon(aes(
    ymin = mean_randompairs_above - sd_randompairs_above,
    ymax = mean_randompairs_above + sd_randompairs_above
  ),
  fill = "grey70", alpha = 0.4) +
  geom_point() +
  geom_errorbar(aes(
    ymin = mean_randompairs_above - sd_randompairs_above,
    ymax = mean_randompairs_above + sd_randompairs_above
  ), width = 0.01) +
  geom_point(aes(y = real_value_above), color = "blue", size = 3) +
  ylim(0, 100) +
  xlim(0, 0.6) +
  ylab("percent of genes above spearman correlation of X") +
  xlab("spearman correlation") +
  theme_classic()

ggsave(paste(outdir,"mod_",whichmod,"_bulk_percent_genes_above_corr.svg",sep=""), width = 8, height = 6, dpi = 300)

###############################################################################
# 16. HEATMAP
###############################################################################


#Heatmaps of mod 6
### binary


Mvoles <- coldata %>% filter(Sex=="M")
Fvoles <- coldata %>% filter(Sex=="F")

scaled_data <- norm_counts_modt %>% dplyr::select(-Animal)
scaled_data <- scale(scaled_data)
stilldf <- coldata %>% dplyr::select(still_together) %>% filter(rownames(coldata) %in% rownames(scaled_data))
Pairdf <-  coldata %>% dplyr::select(pairname) %>% filter(rownames(coldata) %in% rownames(scaled_data))
Conditiondf <- coldata %>% dplyr::select(Condition) %>% filter(rownames(coldata) %in% rownames(scaled_data))
Genderdf <- coldata %>% dplyr::select(Sex) %>% filter(rownames(coldata) %in% rownames(scaled_data))
  
ht0 <- Heatmap(as.matrix(stilldf), name="still",width = unit(1, "cm"))
ht1 <- Heatmap(
  as.matrix(scaled_data),
  name = "counts",
  clustering_distance_rows = "binary",   # distance
  clustering_method_rows  = "centroid"     # linkage method
)
ht2 <- Heatmap(as.matrix(Pairdf), name="Pair",width = unit(1, "cm"))
ht3 <- Heatmap(as.matrix(Conditiondf), name="Condition",width = unit(1, "cm"))
ht4 <- Heatmap(as.matrix(Genderdf), name="Gender",width = unit(1, "cm"))
ht_list <- ht1+ht2+ht3+ht4
ht_list = draw(ht_list, main_heatmap = "counts")

pdf(paste0(outdir,"Bothgendersheatmapbianary.pdf"), width = 8, height = 6)  # width/height in inches
draw(ht_list)                                  # or draw(ht_list, merge_legend = TRUE, ...)
dev.off()

scaled_data <- norm_counts_modt %>% dplyr::select(-Animal) %>% filter(rownames(coldata) %in% rownames(Mvoles))
scaled_data <- scale(scaled_data)
stilldf <- coldata %>% dplyr::select(still_together) %>% filter(rownames(coldata) %in% rownames(scaled_data))
Pairdf <-  coldata %>% dplyr::select(pairname) %>% filter(rownames(coldata) %in% rownames(scaled_data))
Conditiondf <- coldata %>% dplyr::select(Condition) %>% filter(rownames(coldata) %in% rownames(scaled_data))
Genderdf <- coldata %>% dplyr::select(Sex) %>% filter(rownames(coldata) %in% rownames(scaled_data))

ht0 <- Heatmap(as.matrix(stilldf), name="still",width = unit(1, "cm"))
ht1 <- Heatmap(
  as.matrix(scaled_data),
  name = "counts",
  clustering_distance_rows = "binary",   # distance
  clustering_method_rows  = "centroid"     # linkage method
)
ht2 <- Heatmap(as.matrix(Pairdf), name="Pair",width = unit(1, "cm"))
ht3 <- Heatmap(as.matrix(Conditiondf), name="Condition",width = unit(1, "cm"))
ht4 <- Heatmap(as.matrix(Genderdf), name="Gender",width = unit(1, "cm"))
ht_list <- ht1+ht2+ht3+ht4
ht_list = draw(ht_list, main_heatmap = "counts")



pdf(paste0(outdir,"Maleheatmapbinary.pdf"), width = 8, height = 6)  # width/height in inches
draw(ht_list)                                  # or draw(ht_list, merge_legend = TRUE, ...)
dev.off()

scaled_data <- norm_counts_modt %>% dplyr::select(-Animal) %>% filter(rownames(coldata) %in% rownames(Fvoles))
scaled_data <- scale(scaled_data)
stilldf <- coldata %>% dplyr::select(still_together) %>% filter(rownames(coldata) %in% rownames(scaled_data))
Pairdf <-  coldata %>% dplyr::select(pairname) %>% filter(rownames(coldata) %in% rownames(scaled_data))
Conditiondf <- coldata %>% dplyr::select(Condition) %>% filter(rownames(coldata) %in% rownames(scaled_data))
Genderdf <- coldata %>% dplyr::select(Sex) %>% filter(rownames(coldata) %in% rownames(scaled_data))

ht0 <- Heatmap(as.matrix(stilldf), name="still",width = unit(1, "cm"))
ht1 <- Heatmap(
  as.matrix(scaled_data),
  name = "counts",
  clustering_distance_rows = "binary",   # distance
  clustering_method_rows  = "centroid"     # linkage method
)
ht2 <- Heatmap(as.matrix(Pairdf), name="Pair",width = unit(1, "cm"))
ht3 <- Heatmap(as.matrix(Conditiondf), name="Condition",width = unit(1, "cm"))
ht4 <- Heatmap(as.matrix(Genderdf), name="Gender",width = unit(1, "cm"))
ht_list <- ht1+ht2+ht3+ht4
ht_list = draw(ht_list, main_heatmap = "counts")

pdf(paste0(outdir,"Femaleheatmapbinary.pdf"), width = 8, height = 6)  # width/height in inches
draw(ht_list)                                  # or draw(ht_list, merge_legend = TRUE, ...)
dev.off()

### euclidean


scaled_data <- norm_counts_modt %>% dplyr::select(-Animal)
scaled_data <- scale(scaled_data)
stilldf <- coldata %>% dplyr::select(still_together) %>% filter(rownames(coldata) %in% rownames(scaled_data))
Pairdf <-  coldata %>% dplyr::select(pairname) %>% filter(rownames(coldata) %in% rownames(scaled_data))
Conditiondf <- coldata %>% dplyr::select(Condition) %>% filter(rownames(coldata) %in% rownames(scaled_data))
Genderdf <- coldata %>% dplyr::select(Sex) %>% filter(rownames(coldata) %in% rownames(scaled_data))

ht0 <- Heatmap(as.matrix(stilldf), name="still",width = unit(1, "cm"))
ht1 <- Heatmap(
  as.matrix(scaled_data),
  name = "counts",
  clustering_distance_rows = "euclidean",   # distance
  clustering_method_rows  = "centroid"     # linkage method
)
ht2 <- Heatmap(as.matrix(Pairdf), name="Pair",width = unit(1, "cm"))
ht3 <- Heatmap(as.matrix(Conditiondf), name="Condition",width = unit(1, "cm"))
ht4 <- Heatmap(as.matrix(Genderdf), name="Gender",width = unit(1, "cm"))
ht_list <- ht1+ht2+ht3+ht4
ht_list = draw(ht_list, main_heatmap = "counts")

pdf(paste0(outdir,"Bothgendersheatmapeuclidean.pdf"), width = 8, height = 6)  # width/height in inches
draw(ht_list)                                  # or draw(ht_list, merge_legend = TRUE, ...)
dev.off()

scaled_data <- norm_counts_modt %>% dplyr::select(-Animal) %>% filter(rownames(coldata) %in% rownames(Mvoles))
scaled_data <- scale(scaled_data)
stilldf <- coldata %>% dplyr::select(still_together) %>% filter(rownames(coldata) %in% rownames(scaled_data))
Pairdf <-  coldata %>% dplyr::select(pairname) %>% filter(rownames(coldata) %in% rownames(scaled_data))
Conditiondf <- coldata %>% dplyr::select(Condition) %>% filter(rownames(coldata) %in% rownames(scaled_data))
Genderdf <- coldata %>% dplyr::select(Sex) %>% filter(rownames(coldata) %in% rownames(scaled_data))

ht0 <- Heatmap(as.matrix(stilldf), name="still",width = unit(1, "cm"))
ht1 <- Heatmap(
  as.matrix(scaled_data),
  name = "counts",
  clustering_distance_rows = "euclidean",   # distance
  clustering_method_rows  = "centroid"     # linkage method
)
ht2 <- Heatmap(as.matrix(Pairdf), name="Pair",width = unit(1, "cm"))
ht3 <- Heatmap(as.matrix(Conditiondf), name="Condition",width = unit(1, "cm"))
ht4 <- Heatmap(as.matrix(Genderdf), name="Gender",width = unit(1, "cm"))
ht_list <- ht1+ht2+ht3+ht4
ht_list = draw(ht_list, main_heatmap = "counts")



pdf(paste0(outdir,"Maleheatmapeuclidean.pdf"), width = 8, height = 6)  # width/height in inches
draw(ht_list)                                  # or draw(ht_list, merge_legend = TRUE, ...)
dev.off()

scaled_data <- norm_counts_modt %>% dplyr::select(-Animal) %>% filter(rownames(coldata) %in% rownames(Fvoles))
scaled_data <- scale(scaled_data)
stilldf <- coldata %>% dplyr::select(still_together) %>% filter(rownames(coldata) %in% rownames(scaled_data))
Pairdf <-  coldata %>% dplyr::select(pairname) %>% filter(rownames(coldata) %in% rownames(scaled_data))
Conditiondf <- coldata %>% dplyr::select(Condition) %>% filter(rownames(coldata) %in% rownames(scaled_data))
Genderdf <- coldata %>% dplyr::select(Sex) %>% filter(rownames(coldata) %in% rownames(scaled_data))

ht0 <- Heatmap(as.matrix(stilldf), name="still",width = unit(1, "cm"))
ht1 <- Heatmap(
  as.matrix(scaled_data),
  name = "counts",
  clustering_distance_rows = "euclidean",   # distance
  clustering_method_rows  = "centroid"     # linkage method
)
ht2 <- Heatmap(as.matrix(Pairdf), name="Pair",width = unit(1, "cm"))
ht3 <- Heatmap(as.matrix(Conditiondf), name="Condition",width = unit(1, "cm"))
ht4 <- Heatmap(as.matrix(Genderdf), name="Gender",width = unit(1, "cm"))
ht_list <- ht1+ht2+ht3+ht4
ht_list = draw(ht_list, main_heatmap = "counts")

pdf(paste0(outdir,"Femaleheatmapeuclidean.pdf"), width = 8, height = 6)  # width/height in inches
draw(ht_list)                                  # or draw(ht_list, merge_legend = TRUE, ...)
dev.off()

###############################################################################
# 17. Genes that correlate
###############################################################################


corrs_pvals_geneexp <- corrs_long %>%
  group_by(gene) %>%
  summarise(
    real_correlation = correlation[typepairs == "realpairs"],
    rand_correlation = list(correlation[typepairs == "randompairs"]),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    p_empirical = (sum(unlist(rand_correlation) >= real_correlation) + 1) /
      (length(unlist(rand_correlation)) + 1)
  ) %>%
  ungroup()

anticorrs_pvals_geneexp <- corrs_long %>%
  group_by(gene) %>%
  summarise(
    real_correlation = correlation[typepairs == "realpairs"],
    rand_correlation = list(correlation[typepairs == "randompairs"]),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    p_empirical = (sum(unlist(rand_correlation) <= real_correlation) + 1) /
      (length(unlist(rand_correlation)) + 1)
  ) %>%
  ungroup()





anticorrs_pvals_geneexp_withnames = merge(anticorrs_pvals_geneexp, switchnamesdf, by.x="gene", by.y="gene_id")
anticorrs_pvals_geneexp_withnames$rand_correlation_str <- vapply(
  anticorrs_pvals_geneexp_withnames$rand_correlation,
  function(x) paste(x, collapse = ";"),
  FUN.VALUE = character(1)
)

anticorrs_pvals_geneexp_withnames <- anticorrs_pvals_geneexp_withnames %>% dplyr::select(-rand_correlation)
anticorrs_pvals_geneexp_withnames$p_empirical <- as.numeric(anticorrs_pvals_geneexp_withnames$p_empirical)
anticorrs_pvals_geneexp_withnames$real_abs <- as.numeric(anticorrs_pvals_geneexp_withnames$real_correlation)
#write.csv(anticorrs_pvals_geneexp_withnames, paste0(outdir,"antricorrs_pvals_geneexp_withnames.csv"))


ggplot(
  anticorrs_pvals_geneexp_withnames,
  aes(x = p_empirical, y = real_correlation)
) +
  geom_point()


corrs_pvals_geneexp_withnames = merge(corrs_pvals_geneexp, switchnamesdf, by.x="gene", by.y="gene_id")
corrs_pvals_geneexp_withnames$rand_correlation_str <- vapply(
  corrs_pvals_geneexp_withnames$rand_correlation,
  function(x) paste(x, collapse = ";"),
  FUN.VALUE = character(1)
)


corrs_pvals_geneexp_withnames <- corrs_pvals_geneexp_withnames %>% dplyr::select(-rand_correlation)
corrs_pvals_geneexp_withnames$p_empirical <- as.numeric(corrs_pvals_geneexp_withnames$p_empirical)
corrs_pvals_geneexp_withnames$real_abs <- as.numeric(corrs_pvals_geneexp_withnames$real_correlation)
corrs_pvals_geneexp_withnames$padj_BH <- p.adjust(
  corrs_pvals_geneexp_withnames$p_empirical,
  method = "BH"
)
#write.csv(corrs_pvals_geneexp_withnames, paste0(outdir,"corrs_pvals_geneexp_withnames.csv"))

ggplot(
  corrs_pvals_geneexp_withnames,
  aes(x = p_empirical, y = real_correlation)
) +
  geom_point()

graphageneAn1Ani2 <- function(agenename){
  apairfinal <- pairlists[[1]] %>% unique()
  apairfinal <- as.data.frame(apairfinal)
  apairfinal$pairname <- paste(apairfinal$Ani1, "x", apairfinal$Ani2, sep = "")
  
  # Long format: each row is one animal from a pair
  apair_long <- apairfinal %>%
    pivot_longer(
      cols = c(Ani1, Ani2),
      names_to = "Variable",
      values_to = "Animal"
    )
  apair_long <- merge(apair_long, norm_counts_modt, by = "Animal", all.x = TRUE)
  behaviordf <- coldata %>% dplyr::select(Free.interaction, Condition)
  behaviordf$Ani2 <- rownames(behaviordf)
  cols_to_select <- c("Animal", "pairname", "Variable", agenename)
  apair_long_onegene <- apair_long %>% dplyr::select(all_of(cols_to_select))
  # Set the genename being used to a common name
  colnames(apair_long_onegene)[colnames(apair_long_onegene) == agenename] <- "expression"
  # Pivot to one row per pair, wide-style with Ani1 and Ani2 expression in two columns
  apair_wide <- apair_long_onegene %>%
    pivot_wider(
      id_cols = pairname,
      names_from = Variable,
      values_from = expression
    )
  apair_wide$Ani1_expression <- log(apair_wide$Ani1)
  apair_wide$Ani2_expression <- log(apair_wide$Ani2)
  apair_wide <- apair_wide %>%
    separate(
      col   = pairname,
      into  = c("Ani1", "Ani2"),
      sep   = "x",
      remove = FALSE   # keep original pairname; set TRUE if you don't need it
    )
  apair_wide <- merge(apair_wide, behaviordf, by="Ani2")
  apair_wide  <- apair_wide %>% arrange(Ani1)
  write.csv(apair_wide, paste0(outdir,agenename,"_countspaired.csv"))
  gn <- switchnamesdf %>% 
    filter(gene_id == agenename) %>% 
    distinct(gene_name) %>%        # keeps unique rows of gene_name
    pull(gene_name)
    p <- ggplot(apair_wide, aes(x=Ani1_expression, y=Ani2_expression, color=Free.interaction,shape = Condition))+
      geom_point()+scale_shape_manual(values = c(
        "4wk_sep" = 16,   # filled circle
        "2wk_pair" = 17    # filled triangle
      ))+ggtitle(paste(agenename, gn))
  
}


###############################################################################
# 18. SINGLE GENES
###############################################################################

#agenename = "ENSMOCG00000005133"
#agenename = "ENSMOCG00000007634"
#agenename = "ENSMOCG00000020806"
#agenename = "ENSMOCG00000020339"
#agenename = "ENSMOCG00000012214"
agenename = "ENSMOCG00000021921"
#agenename = "ENSMOCG00000016899"
#agenename = "ENSMOCG00000020207"
#agenename = "ENSMOCG00000012103"
#agenename = "ENSMOCG00000018972"
#agenename = "ENSMOCG00000022537"

apw = graphageneAn1Ani2(agenename)
apw


listofgenesvector <- c("ENSMOCG00000000769", "ENSMOCG00000002279", "ENSMOCG00000022537", "ENSMOCG00000006851")

for (agenename in listofgenesvector){apw = graphageneAn1Ani2(agenename)
ggsave(
  filename = paste0(outdir,"Female_Male correlation", agenename, ".svg"),
  plot     = apw,
  width    = 5,
  height   = 5,
  dpi      = 300
)}

###############################################################################
# 19. orthologs so I can use go terms
###############################################################################


ensembl <- useEnsembl(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl",
  mirror  = "useast"      # or "uswest", "asia"
)

# Use current Ensembl
ensembl <- useEnsembl(biomart = "ensembl")  


# Prairie vole gene dataset
moc <- useEnsembl(biomart = "ensembl", dataset = "mochrogaster_gene_ensembl")  

# Human gene dataset
hs  <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")    


orth_moc_hs <- getBM(
  attributes = c("ensembl_gene_id",                     # vole ENSMOCG...
                 "external_gene_name",                  # vole gene symbol (if any)
                 "hsapiens_homolog_ensembl_gene",       # human ENSG...
                 "hsapiens_homolog_associated_gene_name",
                 "hsapiens_homolog_orthology_type"),    # one2one, one2many, etc.
  mart = moc
)



orth_1to1 <- subset(orth_moc_hs,
                    hsapiens_homolog_orthology_type == "ortholog_one2one")


res_human <- res %>%
  inner_join(orth_1to1,
             by = c("gene_id" = "ensembl_gene_id"))

# Foreground = significant DE genes mapped to human
fg_hs <- res_human %>%
  filter(!is.na(padj), padj < 0.1) %>%
  pull(hsapiens_homolog_associated_gene_name) %>%
  unique()


fg_hs_down_4week <- res_human %>%
  filter(!is.na(padj), padj < 0.1, log2FoldChange < 0) %>%
  pull(hsapiens_homolog_associated_gene_name) %>%
  unique()

fg_hs_up_4week <- res_human %>%
  filter(!is.na(padj), padj < 0.1, log2FoldChange > 0) %>%
  pull(hsapiens_homolog_associated_gene_name) %>%
  unique()

# Background = all tested genes with a human ortholog
bg_hs <- res_human %>%
  filter(!is.na(hsapiens_homolog_associated_gene_name)) %>%
  pull(hsapiens_homolog_associated_gene_name) %>%
  unique()

#gocat = "BP"
#gocat = "MF"
gocat = "CC"


# If your IDs are ENSEMBL:
ego_bp <- enrichGO(
  gene          = fg_hs,
  universe      = bg_hs,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = gocat,             # "BP", "MF", or "CC"
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE           
)

head(as.data.frame(ego_bp), 20)

p <- dotplot(ego_bp, showCategory = 20)
p

ggsave(
  filename = paste0(outdir,"plot_all_go_",gocat,".png"),
  plot     = p,
  width    = 6,   # inches
  height   = 5,
  dpi      = 300  
)


# If your IDs are ENSEMBL:
ego_bp_down_4week <- enrichGO(
  gene          = fg_hs_down_4week,
  universe      = bg_hs,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = gocat,             # "BP", "MF", or "CC"
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE           
)

# Look at results
head(as.data.frame(ego_bp_down_4week), 20)

p <- dotplot(ego_bp_down_4week, showCategory = 20)
p

ggsave(
  filename = paste0(outdir,"plot_down_4weeks_go_",gocat,".png"),
  plot     = p,
  width    = 6,   # inches
  height   = 5,
  dpi      = 300  # ignored for pure vector but fine to set
)

ego_bp_up_4week <- enrichGO(
  gene          = fg_hs_up_4week,
  universe      = bg_hs,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = gocat,             # "BP", "MF", or "CC"
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE    
)

# Look at results
head(as.data.frame(ego_bp_up_4week), 20)

p <- dotplot(ego_bp_up_4week, showCategory = 20)
p

ggsave(
  filename = paste0(outdir,"plot_up_4weeks_go_",gocat,".png"),
  plot     = p,
  width    = 6,   # inches
  height   = 5,
  dpi      = 300  
)


genes_GO_up <- ego_bp_up_4week %>%
  as.data.frame() %>%                     
  filter(ID == "GO:0048562") %>%
  pull(geneID) %>%
  strsplit("/") %>%                       
  unlist()

genes_GO_down <- ego_bp_down_4week %>%
  as.data.frame() %>%                      
  filter(ID == "GO:0008038") %>%
  pull(geneID) %>%
  strsplit("/") %>%                        
  unlist()


