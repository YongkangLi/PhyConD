#!/bin/bash

 #SBATCH --partition=scavenger
 #SBATCH --cpus-per-task=45
 #SBATCH --mem-per-cpu=1G
 #SBATCH --output=out/slurm.out
 #SBATCH --error=err/slurm.err

 module load Java/11.0.8

 java -Xss1024m -jar smc.jar -s 32 -p 1024 -m 32 -l ${SLURM_ARRAY_TASK_ID} 1> out/"${SLURM_ARRAY_TASK_ID}.out" 2> err/"${SLURM_ARRAY_TASK_ID}.err"
