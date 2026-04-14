#!/bin/bash -l
#SBATCH -D /quobyte/jrigrp/cstark/JackGWAS/
#SBATCH --job-name=DH_phg_Jack_GWAS
#SBATCH --array=1-11
#SBATCH -o /quobyte/jrigrp/cstark/JackGWAS/logs/%x_%j_%A_%a.sbatch.out
#SBATCH -e /quobyte/jrigrp/cstark/JackGWAS/logs/%x_%j_%A_%a.sbatch.err
#SBATCH -t 24:00:00
#SBATCH --partition=high
#SBATCH --mem=120000
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=crstark@ucdavis.edu

module load gemma

trait=$(sed -n "${SLURM_ARRAY_TASK_ID}p" trait_names.txt)

gemma -bfile DH_phg_Jack_bed \
  -p DH_phg_Jack_phenotypes_BLUPs_gemma.pheno \
  -n ${SLURM_ARRAY_TASK_ID} \
  -gk 1 \
  -o DH_phg_Jack_bed_${trait}

#loop over traits with -n ### gotta run this in sbatch. otherwise takes too long
gemma -bfile DH_phg_Jack_bed \
  -p DH_phg_Jack_phenotypes_BLUPs_gemma.pheno \
  -k output/DH_phg_Jack_bed_${trait}.cXX.txt \
  -lmm 4 \
  -n ${SLURM_ARRAY_TASK_ID} \
  -o DH_phg_Jack_gwas_trait_${trait}