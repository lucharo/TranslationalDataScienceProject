#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=1:mem=96gb

module load anaconda3/personal

cd /rdsgpfs/general/user/lc5415/home/BEES_TDS/Scripts

Rscript EDA.R
