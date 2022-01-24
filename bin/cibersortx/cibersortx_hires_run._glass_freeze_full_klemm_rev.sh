#!/bin/sh
#SBATCH --cpus-per-task=12 
#SBATCH --time 8:00:00
#SBATCH --mem=36g
#SBATCH --output /projects/verhaak-lab/GLASS-III/logs/slurm/CIBERSORTx_HiRes.out
#SBATCH --error /projects/verhaak-lab/GLASS-III/logs/slurm/CIBERSORTx_HiRes.err

#sbatch cibersortx_hires_run._glass_freeze_full_klemm_rev.sh --job-name cibersortx_hires
#srun -q batch -N 1 --cpus-per-task=12 --mem 36G -t 8:00:00 --pty bash 

module load singularity

singularity exec \
-B /projects/verhaak-lab/GLASS-III/data/cibersortx/glass_20210923/:/src/data \
-B /projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-freeze_validation_rev/:/src/outdir \
/projects/verhaak-lab/varnf/docker/CIBERSORTx/hires_latest.sif "/src/CIBERSORTxHiRes" \
--username Frederick.S.Varn.Jr.GR@dartmouth.edu \
--token  ef745da3ebf38e86c4b115ea37df2362 \
--mixture gene_tpm_matrix_all_samples.tsv \
--sigmatrix scgp_sig.txt \
--label 'GLASS' \
--rmbatchSmode TRUE \
--refsample 10x_scgp_cibersortx_ref_06092020.txt \
--groundtruth klemm_purified_avg_ground_truth_cpm.txt \
--subsetgenes common_genes.txt \
--threads 12 \
--cluster  true