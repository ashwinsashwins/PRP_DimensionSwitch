#!/bin/bash

# ---  Note -----------------------------
## #SBATCH commands must come before the program
## comment out: use ## 
## submit job to slurm:  sbatch run_BI_DIM.sh
## job ID: -J DimensionSwitchPRP

# --- Start of slurm commands -----------
## use Badre lab account
##SBATCH --account=carney-dbadre-condo

# Request time
#SBATCH -t 32:00:00

# Amount of memory
#SBATCH --mem=250G

# Number of nodes
#SBATCH -n 1 
#SBATCH -c 20

## Take advantage of all the threads (linear algebra)
## $SLURM_CPUS_ON_NODE returns actual number of cores on node
## rather than $SLURM_JOB_CPUS_PER_NODE, which returns what --cpus-per-task asks for
##export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# Email yourself
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ashwin_srinivasan@brown.edu

# Name your job
#SBATCH -J DimensionSwitchPRP

# Output file
#SBATCH -o DimensionSwitchPRP%j.txt
#SBATCH -e DimensionSwitchPRP%j.txt

# Job array
## 200,204,205,301,304,305,307,308,309,310,312,313,315,316,317,318,321,322,323,324,325,326,327,328,330,331,333,334,335,337,338,340,341,342,343,344,345,346,347,348,350,351,352,353,354,356,357,358,360,362,363,364,365,366,367,368,369,370,371,372,373,375,376
## 400,405,501,504,505,507,508,509,510,512,513,515,516,517,518,521,522,523,524,525,528,530,531,533,534,535,537,538,540,541,542,543,544,545,546,548,550,551,552,553,554,556,557,558,560,562,563,564,565,566,567,568,569,570,571,572,573,575
#SBATCH --array=219,220,221,222,223,224


##601,602,603,604,605,606,607,608,609,610,611,612,613,614,615,616,617,618,619,620,621,622,623,624

##201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224 
##


#----Start of script --------------------
module load r/4.5.1-iikl
##echo " #$SLURM_NTASKS cores detected on `hostname`"
##srun -n $SLURM_NTASKS > hostlist.txt
##R --no-save < DimensionSwitchPRP_Decode.R
##R --no-save < MERGINGDimensionSwitchPRP_Decode.R
R --no-save < PRPDimensionSwitch_Decode_CROSS.R
##R --no-save < HAL2017_ERROR12_Aggregate_DBi.R
