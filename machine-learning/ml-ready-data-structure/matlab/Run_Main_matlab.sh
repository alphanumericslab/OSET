#!/bin/bash
#SBATCH --job-name=ML_Ready_Feature_Extraction_MATLAB
#SBATCH --output=ML_Ready_MATLAB.out
#SBATCH --error=ML_Ready_MATLAB.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=170:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=65G
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=bmemousavi@gmail.com

# Run the Python script
export PYTHONUNBUFFERED=1 
python Generate_MD5_Checksum.py

# Set the path to your MATLAB executable
MATLAB_EXECUTABLE="/usr/local/MATLAB/R2022b/bin/matlab"

# Set the path to the folder containing functions
ADDITIONAL_PATH1="/labs/samenilab/team/somayyeh_mousavi/codes/OSET/matlab/tools/ecg/"
ADDITIONAL_PATH2="/labs/samenilab/team/somayyeh_mousavi/codes/OSET/matlab/tools/generic/"

# Run MATLAB script
$MATLAB_EXECUTABLE -nodisplay -r "addpath('$ADDITIONAL_PATH1'); addpath('$ADDITIONAL_PATH2'); addpath(pwd); main_matlab; exit;"

echo "Job finished!"
