#!/bin/bash
#SBATCH --job-name=ML_Ready_Feature_Extraction_Python
#SBATCH --output=ML_Ready_Python.out
#SBATCH --error=ML_Ready_Python.err
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
python main_python.py

