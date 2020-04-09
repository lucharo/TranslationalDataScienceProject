#PBS -l walltime=6:00:00
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -N EDA

module load anaconda3/personal

cd /rdsgpfs/general/user/lc5415/home/BEES_TDS/Scripts

time Rscript EDA.R
