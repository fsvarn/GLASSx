#!/bin/sh
#SBATCH --cpus-per-task=12 
#SBATCH --time 72:00:00
#SBATCH --mem=72g
#SBATCH --output /projects/verhaak-lab/GLASS-III/logs/slurm/CIBERSORTx_HiRes.out
#SBATCH --error /projects/verhaak-lab/GLASS-III/logs/slurm/CIBERSORTx_HiRes.err

#sbatch cibersortx_hires_run.sh --job-name cibersortx_hires
#srun -q batch -N 1 --cpus-per-task=12 --mem 36G -t 8:00:00 --pty bash 

module load singularity

singularity exec \
-B /projects/verhaak-lab/GLASS-III/data/subset/cibersortx_hires_toil/:/src/data \
-B /projects/verhaak-lab/GLASS-III/results/cibersortx/hires/toil/:/src/outdir \
/projects/verhaak-lab/varnf/docker/CIBERSORTx/hires_latest.sif "/src/CIBERSORTxHiRes" \
--username Frederick.S.Varn.Jr.GR@dartmouth.edu \
--token  ef745da3ebf38e86c4b115ea37df2362 \
--mixture TcgaTargetGtex_GBMLGGBRAIN_rsem_gene_symbol_tpm_antilog.txt \
--sigmatrix scgp_sig.txt \
--label 'toil' \
--rmbatchSmode TRUE \
--refsample 10x_scgp_cibersortx_ref_06092020.txt \
--groundtruth klemm_purified_avg_ground_truth_cpm_subset.txt \
--subsetgenes common_genes.txt \
--threads 12 \
--cluster  true