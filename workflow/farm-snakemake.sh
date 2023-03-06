#!/usr/bin/env bash
# run snakemake on a slurm cluster

#SBATCH -J main_admixcov_sims
#SBATCH -A gmcoopgrp
#SBATCH -p high2
#SBATCH -t 10-00:00
#SBATCH --output main_snakemake.log
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 1G
#SBATCH --mail-type END
#SBATCH --mail-user acpsimon@ucdavis.edu

module load miniconda3

snakemake --profile farm-profile -n