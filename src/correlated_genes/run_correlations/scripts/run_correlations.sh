outdir=/projects/libr8020/run_R_scripts/outdir/
while IFS=, read -r clust gene
do
sleep 20

sbatch --export=gene=$gene,clust=$clust,outdir=$outdir run_correlations.sbatch
done < ../indir/correlation_genes_perclust_minavgcts.csv
