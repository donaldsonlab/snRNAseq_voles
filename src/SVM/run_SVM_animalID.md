SVM_animalID
================
Liza Brusman
2024-03-08

``` r
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(caTools)
library(e1071)
library(caret)
library(ComplexHeatmap)
library(forcats)
library(ggridges)
```

read in seurat object

``` r
SCT_norm <- readRDS("../seurat_clustering/output/SCT_norm.rds")
```

import metadata

``` r
meta <- read.csv("../../docs/seq_beh_metadata.csv")

meta.pair <- pivot_wider(meta, id_cols = c(pair, pair_type, SS_OS), names_from = color, values_from = c(LT_phuddle, ST_phuddle))
meta.pair <- meta.pair[meta.pair$pair != '4918x4967',]
meta.pair <- meta.pair %>% arrange(pair)

rownames(meta.pair) <- meta.pair$pair
```

see if SVM can identify individual animals

``` r
setwd("output/")

#cluster names
clusters <- c("Drd1Pdyn", "Drd1PdynOprm1", "Drd1Penk", "Drd2Penk", "Drd2NoPenk", "GABAergicNeurons", "Dlx2ImmatureNeurons", "SstNpyInterneurons", "PvalbInterneurons", "CholinergicInterneurons", "MatureOligos", "ImmatureOligos", "Astrocytes", "Microglia")

clust.pred.df.ani <- data.frame(Cluster = character(), Misclass = numeric())
diag.list.ani <- list()
cm.list.ani <- list()
Anis <- unique(SCT_norm$Ani)
Anis <- Anis[Anis != "4918"]


for (clust in clusters) {

  print(clust)
  print(Sys.time())
  #subset cluster
  Idents(SCT_norm) <- "new_clusts"
  SCT_mini <- subset(x = SCT_norm, idents = clust)
  
  #downsample per animal
  Idents(SCT_mini) <- "Ani"

  #remove animal 4918 (no partner in dataset)
  SCT_mini <- subset(SCT_mini, idents = Anis)
  SCT_mini <- subset(x = SCT_mini, downsample = 200)
  
  #turn counts into matrix - top 3000 genes used for clustering
  genes3000 <- GetAssayData(SCT_mini, assay = "SCT", slot = "scale.data") %>% as.matrix()
  genes3000 <- t(genes3000)
  metadata <- SCT_mini@meta.data
  
  gc()
  
    
  dat = data.frame(genes3000, Ani = as.factor(metadata$Ani))
  print(nrow(dat))

  #split data into train and test
  set.seed(123)
  split = sample.split(dat$Ani, SplitRatio = 0.75)
    
  training_set = subset(dat, split == TRUE)
  test_set = subset(dat, split == FALSE)
  
  classifier = svm(formula = Ani ~ .,
                 data = training_set,
                 type = 'C-classification',
                 kernel = 'radial',
                 probability = TRUE,
                 cost = 10)
  
  # Predicting the Test set results
  y_pred = predict(classifier, newdata = test_set, decision.values = TRUE, probability = TRUE)
  
  # Making the Confusion Matrix
  cm = table(test_set[, 3001], y_pred)
  
  print(cm)
  
  # Missclassification Rate
  misclass <- 1-sum(diag(cm))/sum(cm)
  print(misclass)
  
  mini.df <- data.frame(Cluster = clust, Misclass = misclass)
  clust.pred.df.ani <- rbind(clust.pred.df.ani, mini.df)
  
  #add info about confusion matrix diagonals to list
  diag.list.ani[[clust]] <- diag(cm)
  
  #add actual confusion matrices to list
  cm.list.ani[[clust]] <- cm


  }


clust.pred.df.ani$Accuracy = 1 - clust.pred.df.ani$Misclass
```

rearranging data into one df

``` r
all.classifications <- data.frame(Cluster = character(), Var1 = character(), y_pred = character(), Freq = numeric())

for (i in names(cm.list.ani)) {
  test <- data.frame(cm.list.ani[[i]])
  mini.df.all <- test
  mini.df.all$Cluster = i
  
  #this just binds all the classifications (animal-animal)
  all.classifications <- rbind(all.classifications, mini.df.all)

}
```

save SVM output

``` r
setwd("ouput/")
write.csv(all.classifications, "ani_pairwise_classifications_withself.csv")
```

train on 37 animals, test on held-out animal

``` r
setwd("output/")
clust.pred.ani.exclude <- data.frame(Cluster = character(), animal = character(), prediction = character(), Misclass = numeric())
diag.list.ani <- list()
cm.list.ani <- list()
Anis <- unique(SCT_norm$Ani)
Anis <- Anis[Anis != "4918"]

clusters <- c('Drd1Pdyn', 'Drd1PdynOprm1', 'Drd1Penk', 'Drd2Penk', 'Drd2NoPenk', 
                                      'GABAergicNeurons', 'Dlx2ImmatureNeurons', 'SstNpyInterneurons', 
                                      'PvalbInterneurons', 'CholinergicInterneurons', 'MatureOligos', 
                                      'ImmatureOligos', 'Astrocytes', 'Microglia')

indiv.classifications <- data.frame(animal = character(), Cluster = character(), Var1 = character(), y_pred = character(), Freq = numeric())

for (clust in clusters) {

  print(clust)
  print(Sys.time())
  #subset cluster
  Idents(SCT_norm) <- "new_clusts"
  SCT_mini <- subset(x = SCT_norm, idents = clust)
  
  #downsample per animal
  Idents(SCT_mini) <- "Ani"
  SCT_mini <- subset(x = SCT_mini, idents = Anis)
  SCT_mini <- subset(x = SCT_mini, downsample = 200)
  
  #turn counts into matrix - top 3000 genes used for clustering
  genes3000 <- GetAssayData(SCT_mini, assay = "SCT", slot = "scale.data") %>% as.matrix()
  genes3000 <- t(genes3000)
  metadata <- SCT_mini@meta.data
  
  gc()
  
  dat = data.frame(genes3000, Ani = as.factor(metadata$Ani))
  print(nrow(dat))
  
  for (ani in Anis) {
    print(ani)
    ani_pair <- meta$pair[meta$animal==ani]
    print(ani_pair)
    ani_partner <- strsplit(ani_pair, split = "x")[[1]]
    ani_partner <- ani_partner[ani_partner != ani] 
    print(ani_partner)
    
    
    ani_dat <- dat %>% filter(Ani == ani)
    other_dat <- dat %>% filter(Ani != ani)
    
    training_set <- other_dat
    test_set <- ani_dat
    
    classifier = svm(formula = Ani ~ .,
                 data = training_set,
                 type = 'C-classification',
                 kernel = 'radial',
                 probability = TRUE,
                 cost = 10)
  
    # Predicting the Test set results
    y_pred = predict(classifier, newdata = test_set, decision.values = TRUE, probability = TRUE)
    
    # Making the Confusion Matrix
    cm = table(test_set[, 3001], y_pred)
    
    print(cm)
    
    mini.df <- data.frame(cm)
    mini.df$animal <- ani
    mini.df$Cluster <- clust
    
    
    indiv.classifications <- rbind(indiv.classifications, mini.df)

  }
}
```

save SVM output

``` r
setwd("output/")
write.csv(indiv.classifications, "all_pairwise_classifications_excludeself_200cells.csv")
```
