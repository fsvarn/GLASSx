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
-B /projects/verhaak-lab/GLASS-III/data/subset/cibersortx_group_tcga/:/src/data \
-B /projects/verhaak-lab/GLASS-III/results/cibersortx/group/tcga/:/src/outdir \
/projects/verhaak-lab/varnf/docker/CIBERSORTx/gep_latest.sif Rscript /src/R_modules/CIBERSORTxGEP.R  \
--username Frederick.S.Varn.Jr.GR@dartmouth.edu \
--token  ef745da3ebf38e86c4b115ea37df2362 \
--mixture mes_tcga_mixture.txt \
--sigmatrix scgp_sig.txt \
--label 'TCGA_MES' \
--rmbatchSmode TRUE \
--refsample 10x_scgp_cibersortx_ref_06092020.txt \
--threads 12 

singularity exec \
-B /projects/verhaak-lab/GLASS-III/data/subset/cibersortx_group_tcga/:/src/data \
-B /projects/verhaak-lab/GLASS-III/results/cibersortx/group/tcga/OTH:/src/outdir \
/projects/verhaak-lab/varnf/docker/CIBERSORTx/gep_latest.sif Rscript /src/R_modules/CIBERSORTxGEP.R  \
--username Frederick.S.Varn.Jr.GR@dartmouth.edu \
--token  ef745da3ebf38e86c4b115ea37df2362 \
--mixture other_tcga_mixture.txt \
--sigmatrix scgp_sig.txt \
--label 'TCGA_OTH' \
--rmbatchSmode TRUE \
--refsample 10x_scgp_cibersortx_ref_06092020.txt \
--threads 12 
