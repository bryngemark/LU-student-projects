#!/bin/sh
#SBATCH -t 03:00:00 
#SBATCH -A lu2025-2-27
#SBATCH -J PhotoNeutrons
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --mail-user=da1583sa-s@student.lu.se
#SBATCH --mail-type=END

cd OutputFiles
denv fire ../singleNconfig.py 5000 $1

#Here is how you run your bash file: 
#sbatch JobLauncher.sh 
