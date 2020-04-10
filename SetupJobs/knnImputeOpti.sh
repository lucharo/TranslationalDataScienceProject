#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=20:mem=240gb
#PBS -N knnImputeOpti

module load anaconda3/personal
source activate TDS

cd /rdsgpfs/general/user/lc5415/home/BEES_TDS/Scripts

seed=1

time Rscript knnOptimizationSingle.R $seed
conda deactivate