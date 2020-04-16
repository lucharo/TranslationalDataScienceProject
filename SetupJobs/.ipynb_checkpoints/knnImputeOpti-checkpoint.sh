#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=32:mem=124gb
#PBS -N knnImputeOptiPerParam
#PBS -J 1:50

cd /rdsgpfs/general/user/lc5415/home/BEES_TDS/Scripts
module load anaconda3/personal
source activate TDS
seed=$PBS_ARRAY_INDEX

time Rscript knnOptimizationSingle.R $seed
conda deactivate