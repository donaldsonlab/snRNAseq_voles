# snRNAseq_voles
snRNA-seq analysis of nucleus accumbens tissue from prairie voles

All raw sequencing data is on GEO with the accession number GSE255620. It is embargoed until Aug 1, 2025. \

The Seurat object with the snRNA-seq data (named SCT_norm) is on Figshare. Other important files (metadata, and tissue-level RNA-seq data) are in "docs/" directory \

Link for Seurat object (private link - should change to public link upon publication): https://figshare.com/s/13648f94bab3f8cda994 \
Metadata is in "docs/" directory (and on Figshare). File is "seq_beh_metadata.csv" \
RNA-seq data is in "docs/" directory. Files are "Input_RNAseq_metrics.xlsx" and "Merged_all_inputs.txt" \
CSV containing mean(log(counts)) per gene per animal per cluster (for correlated genes analysis which was included in first submission, but not resubmission) is on Figshare (private link): https://figshare.com/s/647301ed42293deed2f5 \

All analysis was run as R or Jupyter notebooks. For R notebooks, the .Rmd file and the .md file (with figure outputs) are both included.

## To run this analysis
1. Download sequencing data from GEO.
2. Run code in src/seurat_integration to create a separate Seurat object per animal and then integrate all animals into one Seurat object. I did this on the computing cluster (RC at CU Boulder).
3. Run code in src/seurat_clustering.
   - First, run seurat_batch_correction_filtering.Rmd to run batch correction and filter out animal/cells that should be excluded.
   - Then, run UMAP_clustering.Rmd to create UMAP and analyze cell type proportions.
   - UMAP_clustering.Rmd also creates the final Seurat object used for all downstream analysis. This file (SCT_norm) is also on Figshare as noted above, so this part of the analysis does not need to be run every time. Note: UMAP is stochastic, so UMAP generated may vary slightly from the UMAP presented in the paper.
4. Run code in src/hotspot.
   - First, run the Jupyter notebook that runs Hotspot.
   - Then, run hotspot_forpaper.Rmd to analyze the output from Hotspot.
   - To analyze behavior data (partner preference and free interaction) with the Hotspot output (e.g. to examine correlations), run behavior_hotspot.Rmd.
5. Run code in src/SVM.
   - First, run run_SVM_animalID.Rmd to run the SVM.
   - Then, run SVM_plots_forpaper.Rmd to visualize the SVM outputs.
6. Run code in separation_RNAseq.
   - Run separation_data_analysis.Rmd to create all plots for separation data.

**Note:** Analysis using src/correlated_genes was not included in the resubmission to Science. This code finds which genes are correlated between partners in each cluster, runs GO analysis, and plots GO terms in a heatmap.



