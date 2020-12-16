#!/bin/sh
#SBATCH --cpus-per-task=12 
#SBATCH --time 8:00:00
#SBATCH --mem=36g
#SBATCH --output /projects/verhaak-lab/GLASS-III/logs/slurm/CIBERSORTx_HiRes.out
#SBATCH --error /projects/verhaak-lab/GLASS-III/logs/slurm/CIBERSORTx_HiRes.err

#sbatch cibersortx_hires_run.sh --job-name cibersortx_hires
#srun -q batch -N 1 --cpus-per-task=12 --mem 36G -t 8:00:00 --pty bash 

module load singularity

singularity exec \
-B /projects/verhaak-lab/GLASS-III/data/subset/cibersortx_fractions_glass/:/src/data \
-B /projects/verhaak-lab/GLASS-III/results/cibersortx/fractions/GLASS/:/src/outdir \
/projects/verhaak-lab/varnf/docker/CIBERSORTx/fractions_latest.sif "/src/CIBERSORTxFractions" \
--username Frederick.S.Varn.Jr.GR@dartmouth.edu \
--token  ef745da3ebf38e86c4b115ea37df2362 \
--mixture gene_tpm_matrix_all_samples.tsv \
--sigmatrix neftel_cibersortx_sig_08032020.txt \
--perm 100 \
--label 'GLASS' \
--rmbatchBmode TRUE \
--refsample neftel_cibersortx_cell_type_sourceGEP_08032020.txt \
--verbose TRUE
