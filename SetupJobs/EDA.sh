#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=1:mem=96gb
#PBS -N EDAtry1
module load anaconda3/personal

cd /rdsgpfs/general/user/lc5415/home/BEES_TDS/Scripts

time Rscript EDA.R
