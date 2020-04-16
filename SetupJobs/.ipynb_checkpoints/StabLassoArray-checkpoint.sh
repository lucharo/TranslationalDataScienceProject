#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=8:mem=8gb
#PBS -N StabLassoArrayALL
#PBS -J 1-100
module load anaconda3/personal
source activate TDS

cd /rdsgpfs/general/user/lc5415/home/BEES_TDS/Scripts

seed=$PBS_ARRAY_INDEX

time Rscript StabLassoWCovarsNoBio.R $seed
time Rscript StabLassoWCovars.R $seed
conda deactivate