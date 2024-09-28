#!/bin/bash

#SBATCH --account=eon

#SBATCH --time=00:12:00

#SBATCH --nodes=1

#SBATCH --ntasks-per-node=5

#SBATCH --cpus-per-task=6

#SBATCH --mem-per-cpu=128G

#SBATCH --partition=teton-massmem

#SBATCH --job-name=cellrangerRun

#SBATCH -o outs_file

#SBATCH -e error_file

#SBATCH --mail-type=END

#SBATCH --mail-user nboldon@uwyo.edu

cd /project/eon/nboldon/MCS2023/cellrangerRun

module load cell-ranger-atac > module-load.out 2> module-load.err

cellranger-atac count --id=C301_cr \
--reference=/apps/u/opt/linux/cell-ranger-atac/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
--fastqs=/project/eon/nboldon/MCS2023/cellrangerRun/C301 \
 --sample=C301 --localcores=6 --localmem=128 > cellranger.out 2> cellranger.err
