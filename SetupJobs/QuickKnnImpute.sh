#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=16:mem=64gb
#PBS -N QuickKNN

module load anaconda3/personal

cd /rdsgpfs/general/user/lc5415/home/BEES_TDS/Scripts

Rscript QuickknnImpute.R
