#!/bin/bash
#SBATCH --job-name=JOBNAME
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1

### THIS IS A BLANK SCRIPT "SKELETON" TO GET YOU STARTED IF YOU WANT TO MAKE YOUR OWN SLURM SCRIPT...

## JOB LOG HEADER
  perl -E 'say"="x80'; echo "JOB STARTED: `date`"; echo "NODE: `hostname`"; echo "SCRIPT ${0}:"; cat $0; perl -E 'say"="x80'


## SOFTWARE REQUIREMENTS


## PRIMARY VARIABLES


## COMMAND(s) TO RUN


#JOB LOG FOOTER
perl -E 'say"="x80'; echo "JOB COMPLETED: `date`"; perl -E 'say"="x80'
