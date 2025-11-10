#!/usr/bin/env bash
#$ -pe smp 10
#$ -jc long
#$ -cwd
#$ -p -100
#qsub -l gpu=0 -mods m_mem_free 16G

source /cluster/ddu/pibrahim001/miniconda3/bin/activate
conda activate FMO_JCPE
 
# Holo-FMOPhore
fmophore -dir ./ -d 5 -qm DFTB -t 1 -c 10 -lib -align

# Apo-scan-FMOPhore
# fmophore -prot prot/*.pdb -ligs ligands/ -d 5 -qm DFTB -t 1 -c 10 -lib -align

conda deactivate
conda deactivate