#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=6:mem=24gb
#PBS -N UnivAnalysis
module load anaconda3/personal

cd /rdsgpfs/general/user/lc5415/home/BEES_TDS/Scripts

time Rscript univariate_analysis.R
