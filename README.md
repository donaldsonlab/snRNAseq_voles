# snRNAseq_voles
snRNA-seq analysis of nucleus accumbens tissue from prairie voles

All raw sequencing data is on GEO with the accession number GSE255620. It is embargoed until Aug 1, 2025. 

All analysis was run as R or Jupyter notebooks. For R notebooks, the .Rmd file and the .md file (with figure outputs) are both included. 

## Files on Figshare
- SCT_norm.rds : The Seurat object for all snRNA-seq analyses. All animals are combined in this object. Private link (should change to public link upon publication): https://figshare.com/s/13648f94bab3f8cda994
- Animal metadata (seq_beh_metadata.csv) : metadata from snRNA-seq experiment. Also on Github in "docs/"
- Mean SCT-normalized counts per gene, per animal, per cluster (avg_SCT_cts_per_cluster.csv) : CSV containing mean(log(counts)) per gene per animal per cluster (for correlated genes analysis which was included in first submission, but not resubmission). Private link: https://figshare.com/s/647301ed42293deed2f5

## The "docs/" directory
- [Input_RNAseq_metrics.xlsx](docs/Input_RNAseq_metrics.xlsx): metadata for animals in separation experiment (data from [Sadino, et al.](https://doi.org/10.7554/eLife.80517)).
- [Merged_all_inputs.txt](docs/Merged_all_inputs.txt): RNA-seq gene counts for each animal in separation experiment (data from [Sadino, et al.](https://doi.org/10.7554/eLife.80517)).
- [PPTMetrics_coh1234_updated.csv](docs/PPTMetrics_coh1234_updated.csv): Partner preference test behavior data for animals in snRNA-seq experiment. Note: not all animals in this file were sequenced.
- [ani_mod_scores_allcells_lognorm_counts.csv](docs/ani_mod_scores_allcells_lognorm_counts.csv): mean(log(counts)) values from Hotspot analysis for each animal. This file is generated in run_hotspot.ipynb but provided here to avoid having to re-run Hotspot.
- [free_int_beh.xlsx](docs/free_int_beh.xlsx): Behavior data from free interaction test on snRNA-seq experiment animals.
- [new_clusts_hotspot-gene-modules.csv](docs/new_clusts_hotspot-gene-modules.csv): Gene module membership data generated in run_hotspot.ipynb but provided here to avoid re-running Hotspot.
- [seq_beh_metadata.csv](docs/seq_beh_metadata.csv): metadata from snRNA-seq experiment

## To run this analysis
1. Download sequencing data from GEO.
2. Run code in [src/seurat_integration](src/seurat_integration) to create a separate Seurat object per animal and then integrate all animals into one Seurat object. I did this on the computing cluster (RC at CU Boulder).
3. Run code in [src/seurat_clustering](src/seurat_clustering).
   - First, run [seurat_batch_correction_filtering.Rmd](src/seurat_clustering/seurat_batch_correction_filtering.Rmd) to run batch correction and filter out animal/cells that should be excluded.
   - Then, run [UMAP_clustering.Rmd](src/seurat_clustering/UMAP_clustering.Rmd) to create UMAP and analyze cell type proportions.
   - UMAP_clustering.Rmd also creates the final Seurat object used for all downstream analysis. This file (SCT_norm) is also on Figshare as noted above, so this part of the analysis does not need to be run every time. Note: UMAP is stochastic, so UMAP generated may vary slightly from the UMAP presented in the paper.
4. Run code in [src/hotspot](src/hotspot).
   - First, run the Jupyter notebook [run_hotspot.ipynb](src/hotspot/run_hotspot.ipynb) that runs Hotspot.
   - Then, run [hotspot_forpaper.Rmd](src/hotspot/hotspot_forpaper.Rmd) to analyze the output from Hotspot.
   - To analyze behavior data (partner preference and free interaction) with the Hotspot output (e.g. to examine correlations), run [behavior_hotspot.Rmd](src/hotspot/behavior_hotspot.Rmd).
5. Run code in [src/SVM](src/SVM).
   - First, run [run_SVM_animalID.Rmd](src/SVM/run_SVM_animalID.Rmd) to run the SVM.
   - Then, run [SVM_plots_forpaper.Rmd](src/SVM/SVM_plots_forpaper.Rmd) to visualize the SVM outputs.
6. Run code in [src/separation_RNAseq](src/separation_RNAseq).
   - Run [separation_data_analysis.Rmd](src/separation_RNAseq/separation_data_analysis.Rmd) to create all plots for separation data.

**Note:** Analysis using [src/correlated_genes](src/correlated_genes) was not included in the resubmission to Science. This code finds which genes are correlated between partners in each cluster, runs GO analysis, and plots GO terms in a heatmap.

## renv.lock
This is the R environment with the package versions used for this analysis. An R environment can be created from the [renv.lock](renv.lock) file.




