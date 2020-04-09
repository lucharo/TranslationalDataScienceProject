#PBS -l walltime=4:00:00
#PBS -l select=1:ncpus=1:mem=16gb
#PBS -N PenalisedReg
module load anaconda3/personal

cd /rdsgpfs/general/user/lc5415/home/BEES_TDS/Scripts

time Rscript LASSO_ENet.R
