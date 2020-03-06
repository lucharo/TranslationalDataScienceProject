#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=20:mem=240gb
#PBS -N knnImpute

module load anaconda3/personal

cd /rdsgpfs/general/user/lc5415/home/BEES_TDS/Scripts

Rscript knnImputation.R
