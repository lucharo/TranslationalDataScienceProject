#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=96gb
#PBS -N StabLassoArray
#PBS -J 1-2
module load anaconda3/personal

cd /rdsgpfs/general/user/lc5415/home/BEES_TDS/Scripts

seed=$PBS_ARRAY_INDEX

time Rscript StabLasso.R $seed
