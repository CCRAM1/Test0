#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000M

python PYTHON/correlations.py -f Genus_subsys3.tsv -l | awk -F '\t' '$4<0.05' | awk -F '\t' '$3>0.5' > GenusSubsys3.correlations.tsv

