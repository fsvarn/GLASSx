#!/bin/bash
#SBATCH --export=all
#SBATCH --time=48:00:00
#SBATCH --job-name="ssGSEA_MSigDB"
#SBATCH --cpus-per-task=8
#SBATCH --output=ssGSEA_MSigDB.out
#SBATCH --mem=48g

cd /projects/verhaak-lab/GLASS-III/R/expression/collab/
Rscript msigdb_idhwt_glass.r