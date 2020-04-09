#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=32:mem=96gb
#PBS -N StabLasso
module load anaconda3/personal

cd /rdsgpfs/general/user/lc5415/home/BEES_TDS/Scripts

time Rscript StabLasso.R


