#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 5
#$ -l h=neotera
echo $NSLOTS

#Variables:
venv_path="/home/andreza.cunha/.conda/envs/te_sorter_venv"
input_te_file="/Storage/data2/andreza.cunha/results/flTE/GreenFoxtail.flTE.fa"
out_dir="/Storage/data2/andreza.cunha/results/TEsorter/Sviridis"

#Activate virtual environment in conda:
module load miniconda3
conda activate "${venv_path}"

#Run TEsorter using REXdb:
TEsorter "${input_te_file}" -db rexdb-plant -p 5 -pre "${out_dir}/GreenFoxtail"




