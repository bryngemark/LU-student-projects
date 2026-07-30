#!/bin/bash
#SBATCH -t 100:05:00

#SBATCH -J es

#Uncomment this if you start running into OOM issues
##SBATCH --mem-per-cpu=10600 #Default is 5300 MB

#I prefer to have both the output and the error logs combined, to make it easier to spot runtime errors
#SBATCH -o slurm/slurm_%j.out
#SBATCH -e slurm/slurm_%j.out

#To stop a failed job from requeueing (probably best to leave this unchanged)
#SBATCH --no-requeue

#Uncomment this if you want your config file to run on a node disk directly
#Faster I/O, but take care to copy your output files before the job ends
##cp $SLURM_SUBMIT_DIR/SophieFile.py $SNIC_TMP/config.py

just -f $pathtojustfile/justfile fire config.py $runnumber $eventsperjob
