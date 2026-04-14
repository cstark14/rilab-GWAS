### Commands used to filter and LD prune geno set
module load htslib
module load plink2
module load gemma

bgzip -c /quobyte/jrigrp/jri/DH_phg.vcf > /quobyte/jrigrp/cstark/DH_phg.vcf.gz ### i guess forrest gzipped it after I had done it. so this line isn't needed anymore
tabix -C /quobyte/jrigrp/cstark/DH_phg_Jack.vcf.gz

bcftools view /quotbyte/jrigrp/cstark/DH_phg_Jack.vcf.gz \
  -S phenoParents.txt \
  --min-af 0.05:minor \        # MAF > 1%
  --min-ac 1 \
  --exclude 'F_MISSING > 0.05' \      #<5% missingness
  -Oz -o DH_phg_Jack_OnlyPheno_rmMono_MAF_0.05_Missing_0.05.vcf.gz
  
plink2 --vcf DH_phg_Jack_OnlyPheno_rmMono_MAF_0.05_Missing_0.05.vcf.gz \
  --make-pgen \
  --out DH_phg_Jack \
  --set-missing-var-ids '@:#' \
  --new-id-max-allele-len 1000 \
  --allow-extra-chr
  
# LD pruning first... probably not needed
plink2 --pfile DH_phg_Jack \
  --indep-pairwise 50 10 0.2 \
  --out DH_phg_Jack_LD_50_10_0_2 \
  --allow-extra-chr

# Compute PCA
plink2 --pfile DH_phg_Jack \
  --extract DH_phg_Jack_LD_50_10_0_2.prune.in \
  --pca 20 \
  --out DH_phg_Jack_pca \
  --allow-extra-chr

## Convert to bed out for GEMMA  
plink2 --pfile DH_phg_Jack --make-bed --out DH_phg_Jack_bed --allow-extra-chr

#### see phenosToBLUPs for BLUP conversion

#### see gemmaGWAS_sbatch script
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


#### gapit instead. requires numeric
plink2 --vcf DH_phg_Jack_OnlyPheno_rmMono_MAF_0.05_Missing_0.05.vcf.gz \
  --export A \
  --out DH_phg_Jack_filtered_numeric \
  --allow-extra-chr
