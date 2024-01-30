#all analysis as one script. commented out saving plots

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("svMisc"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("ivmte"))


#import every matrix file and assign its own name
samples <- c("NAc_4905_F",
             "NAc_4975_M",
             "NAc_4918_F",
             # "NAc_4967_M",
             "NAc_4893_M",
             "NAc_4910_M",
             "NAc_4921_F",
             "NAc_4909_M",
             "NAc_4968_F",
             "NAc_4931_M",
             "NAc_4907_F",
             "NAc_5228_F",
             "NAc_5021_F",
             "NAc_5204_M",
             "NAc_4940_F",
             "NAc_4916_M",
             "NAc_4958_F",
             "NAc_4901_F",
             "NAc_4960_F",
             "NAc_4928_M",
             "NAc_5225_F",
             "NAc_5121_M",
             "NAc_4917_M",
             "NAc_4932_M",
             "NAc_5209_M",
             "NAc_4974_M",
             "NAc_4896_F",
             "NAc_4894_M",
             "NAc_5023_F",
             "NAc_4963_F",
             "NAc_5227_F",
             "NAc_5122_M",
             "NAc_4898_F",
             "NAc_4920_F",
             "NAc_4976_M",
             "NAc_5124_M",
             "NAc_4908_M",
             "NAc_4970_M",
             "NAc_4947_F",
             "NAc_923_M"
)

metadata <- read.csv("/scratch/alpine/libr8020/r_analysis/seurat_metadata.csv")

mito_genes <- c('ND1', 'ND2', 'COX1', 'COX2', 'ATP8', 'ATP6', 'COX3', 'ND3', 'ND4L', 'ND4', 'ND5', 'ND6', 'CYTB')

obj.list <- list()

for (s in samples) {{
  nam <- paste(s, "data", sep="_")
  print(nam)
  path = paste("/scratch/alpine/libr8020/r_analysis/seq_data/", s, "/filtered_feature_bc_matrix/", sep="")
  data <- Read10X(data.dir = path)
  assign(nam, data)
  
  #create seurat object from matrix
  seurobj <- CreateSeuratObject(counts = data,min.cells = 3,min.features = 200)
  
  #remove raw data from environment to conserve space
  rm(data)
  
  seurobj <- PercentageFeatureSet(seurobj, features = mito_genes, col.name = "percent_mito")
  
  seurobj  <- subset(x = seurobj, subset  =  nFeature_RNA > 200 & percent_mito < 5)
  
  #pull out animal number from name of sample and add to seurat obj
  seurobj$Ani <- sapply(strsplit(s, "_"), "[", 2)
  #pull out sex of animal from name of sample and add to seurat obj
  seurobj$Sex <- sapply(strsplit(s, "_"), "[", 3)
  
  #add more metadata from meta df
  mini.meta <- metadata[which(metadata$sample_id==s),]
  
  seurobj$Pair <- mini.meta$pair[1]
  
  seurobj$PairType <- mini.meta$pair_type[1]
  
  seurobj$SSOS <- mini.meta$SS_OS[1]
  
  seurobj$Age <- mini.meta$age[1]
  
  seurobj$Preg <- mini.meta$preg[1]
  
  seurobj$SeqCohort <- mini.meta$seq_cohort[1]
  
  seurobj$BehCohort <- mini.meta$beh_cohort[1]
  
  seurobj$SeqConc <- mini.meta$seq_conc[1]
  
  

  #SCTransform normalization
  seurobj <- SCTransform(seurobj)
  
  seurobj<- FindVariableFeatures(seurobj, selection.method = "vst", nfeatures = 2000)
  
  #make variable feature plot
  top20 <- head(VariableFeatures(seurobj), 20)

  #plot variable features with labels
  plot1 <- VariableFeaturePlot(seurobj)
  plot2 <- LabelPoints(plot=plot1, points=top20, repel=TRUE)
  print(plot2)
  #save variable feature plot
  setwd("/scratch/alpine/libr8020/r_analysis/output/")
  fname.varfeat <- paste(s, "var_feat.png", sep="_")
  ggsave(fname.varfeat, plot2, device="png", bg="white")
  
  fname.rds <- paste(s, ".rds", sep="")
  fname.h5 <- paste(s, ".h5Seurat", sep="")

  seurobj <- ScaleData(seurobj,verbose = FALSE)
  seurobj <- RunPCA(seurobj,npcs = 30,verbose = FALSE)
  
  seurobj <- FindNeighbors(seurobj, reduction = "pca", dims=1:17)
  seurobj <- FindClusters(seurobj, resolution=0.2)
  seurobj <- RunUMAP(seurobj, reduction = "pca", dims = 1:17)

  
  elb.plt <- ElbowPlot(object = seurobj, ndims = 30) + ggtitle(s)

  #save elbow plot
  setwd("/scratch/alpine/libr8020/r_analysis/output/")
  fname.elb <- paste(s, "elbow.png", sep="_")
  ggsave(fname.elb, elb.plt, device="png", bg="white")


  setwd("/scratch/alpine/libr8020/r_analysis/output/")
  saveRDS(seurobj, fname.rds)

  umap.plot <- DimPlot(object=seurobj, reduction = "umap") + ggtitle(s)
  umap.vln <- VlnPlot(seurobj, features = c("DRD1A", "DRD2", "Pdyn", "Olig2", "Aspa", "Pdgfra"), pt.size=0)
  setwd("/scratch/alpine/libr8020/r_analysis/output/")
  fname.umap <- paste(s, "umap.pdf", sep="_")
  ggsave(fname.umap, umap.plot, device="pdf", bg="white")
  fname.vln <- paste(s, "umap_vln.pdf", sep="_")
  ggsave(fname.vln, umap.vln, device="pdf", bg="white")

  markers <- FindAllMarkers(seurobj, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
markers <- markers %>% group_by(cluster) %>% arrange(desc(avg_log2FC))

  setwd("/scratch/alpine/libr8020/r_analysis/output")
  fname <- paste(s, "markers.csv", sep="_")
  write.csv(markers, fname)

  obj.list[[s]] <- seurobj
  
  rm(seurobj)
  
}}

features <- SelectIntegrationFeatures(object.list = obj.list)

obj.list <- PrepSCTIntegration(obj.list, anchor.features = features)

all_samps.anchors  <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features, normalization.method = "SCT", reduction = "rpca", dims = 1:30) 

gc()

all_samps <- IntegrateData(anchorset=all_samps.anchors, normalization.method="SCT", dims=1:30)
all_samps <- RunPCA(all_samps,npcs = 30)


all_samps <- FindNeighbors(all_samps, reduction = "pca", dims = 1:17)
all_samps <- FindClusters(all_samps, resolution = 0.2)
all_samps <- RunUMAP(all_samps, reduction = "pca", dims = 1:17)

#SAVE ALL SAMPS
setwd("/scratch/alpine/libr8020/r_analysis/output/")
saveRDS(all_samps, "all_samples_integrated.rds")

setwd("/scratch/alpine/libr8020/r_analysis/output/")
umap.plot <- DimPlot(object = all_samps, reduction = "umap")

umap.vln <- VlnPlot(all_samps, features=c("DRD1A", "DRD2", "Pdyn", "Olig2", "Aspa", "Pdgfra"), pt.size = 0, combine = TRUE)

setwd("/scratch/alpine/libr8020/r_analysis/output/")
ggsave("all_samps_umap.pdf", umap.plot, device="pdf", bg="white")
ggsave("all_samps_umap_vln.pdf", umap.vln, device="pdf", bg="white")

markers <- FindAllMarkers(all_samps, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

markers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)

#write to csv
setwd("/scratch/alpine/libr8020/r_analysis/output/")
fname <- paste("all_samps_markers.csv")
write.csv(markers, fname)


