#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=1:mem=15gb

module load anaconda3/personal
source activate BEES_env

cd $PBS_O_WORKDIR

Rscript PolygenicRiskScore.R

conda deactivate
