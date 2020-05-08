#PBS -l walltime=18:00:00
#PBS -l select=1:ncpus=4:mem=40gb
#PBS -N MLAnalysis
#PBS -J 1-27
module load anaconda3/personal
source activate TDS

combi=$PBS_ARRAY_INDEX
cd /rdsgpfs/general/user/lc5415/home/BEES_TDS/Scripts

time Rscript Analysis.R $combi
conda deactivate
