#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=24gb
#PBS -N knnImpute

module load anaconda3/personal

cd /rdsgpfs/general/user/lc5415/home/BEES_TDS/Scripts

Rscript knnImputation.R
