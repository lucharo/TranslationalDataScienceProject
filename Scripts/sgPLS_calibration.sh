#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=8:mem=20gb
#PBS -N sgPLS_calibration
#PBS -J 1-100

module load anaconda3/personal
module load gsl
source activate BEES_env

cd $PBS_O_WORKDIR

seed=$PBS_ARRAY_INDEX
echo Job number $seed

time Rscript sgPLS_calibration.R $seed

conda deactivate
