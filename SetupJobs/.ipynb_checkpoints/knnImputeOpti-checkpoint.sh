#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=8:mem=96gb
#PBS -N knnImputeOpti
#PBS -J 1-20
module load anaconda3/personal

cd /rdsgpfs/general/user/lc5415/home/BEES_TDS/Scripts

seed=$PBS_ARRAY_INDEX

time Rscript knnOptimizationSingle.R $seed
