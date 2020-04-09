#PBS -l walltime=30:00:00
#PBS -l select=1:ncpus=8:mem=20gb
#PBS -N sgPLS_calibration
#PBS -J 1-2

module load anaconda3/personal
source activate BEES_env

cd $PBS_O_WORKDIR
echo ../$PBS_O_WORKDIR
seed=$PBS_ARRAY_INDEX

time Rscript sgPLS_calibration.R $seed

conda deactivate
