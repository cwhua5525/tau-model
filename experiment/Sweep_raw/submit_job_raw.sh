#!/bin/bash
#SBATCH --job-name=run8054final
#SBATCH --output=8054final_output.txt
#SBATCH --error=8054final_error.txt
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --account=STAT8054
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hua00076@umn.edu

# 移動到你的工作資料夾
cd ~/8054final/experiment/Sweep_raw

# 載入R模組
module purge
module load R/4.4.2-openblas-rocky8


# 執行你的R script
Rscript task_raw.R
