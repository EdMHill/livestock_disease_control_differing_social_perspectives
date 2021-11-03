#!/bin/bash
#SBATCH --job-name=livestock_spatial_sellke_run
#SBATCH --time=01:00:00
#SBATCH --mem=16g
#SBATCH --export=ALL
#SBATCH --partition=ntd

# Submit a job array
#SBATCH --array=1-1000%25
#SBATCH --ntasks=1

# ${SLURM_SUBMIT_DIR} points to the path where this script was submitted from
cd ${SLURM_SUBMIT_DIR}

module load julia/1.6

# ARGS list
# ARGS[1] job_ID
# ARGS[2] job_ID_foundation_value: Set value that job_ID will be incremented upon
# ARGS[3] BatchID_offset: Value that BatchID is offset by
# ARGS[4] RNGseed: To be used to initialise the random number generator
# ARGS[5] ConfigFn: Location and pathogen configuration
# ARGS[6] nReps: Number of replicates requested
# ARGS[7] precalc_dist_and_kernel_array_flag: Set if should precalculate distance and kernel arrays
julia run_generic_livestock_disease_spatial_model_sellke_vers.jl ${SLURM_ARRAY_TASK_ID} 0 5000 99 CumbriaSellkeFMDconfig 10 true

echo "Outbreak simulation runs finished"

cd /home/edhill/FEED/src/SpatialSimns/OptimBehaviour_GenericModel

# args/ARGS list
# args[1] job_ID
# args[2] job_ID_foundation_value: Set value that job_ID will be incremented upon
# args[3] BatchID_offset: Value that BatchID is offset by
# args[4] ConfigFn: Location and pathogen configuration
# args[5] nReps: Number of replicates requested
julia ComputeVaccButNoInfecStats.jl ${SLURM_ARRAY_TASK_ID} 0 5000 CumbriaSellkeFMDconfig 10

echo "Finished"
