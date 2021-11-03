#!/bin/bash
#SBATCH --job-name=get_optimal_strat_livestock_disease_simn
#SBATCH --time=05:00:00
#SBATCH --mem=16g
#SBATCH --export=ALL
#SBATCH --partition=ntd

# Submit a job array
#SBATCH --array=1
#SBATCH --ntasks=1

# ${SLURM_SUBMIT_DIR} points to the path where this script was submitted from
cd ${SLURM_SUBMIT_DIR}

module load julia/1.6

# ARGS list
# ARGS[1] job_ID
# ARGS[2] ConfigFn: Location and pathogen configuration
# ARGS[3] BatchID_offset: Value that BatchID is offset by
# ARGS[4] scen_offset: Value that scen_itr is offset by
julia ComputeOptimalRiskThreshold_VaccBehav_Script.jl ${SLURM_ARRAY_TASK_ID} CumbriaFMDconfig 5000 300

echo "Finished"
