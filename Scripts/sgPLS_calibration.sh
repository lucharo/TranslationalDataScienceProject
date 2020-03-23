#PBS -l walltime=30:00:00
#PBS -l select=1:ncpus=1:mem=10gb
#PBS -N sevnodes1core
#PBS -J 1-100

module load anaconda3/personal
source activate BEES_env

cd $PBS_O_WORKDIR

ichunk=$PBS_ARRAY_INDEX
nchunks=100
seed = $PBS_ARRAY_INDEX

Rscript sgPLS_calibration.R $seed

conda deactivate
