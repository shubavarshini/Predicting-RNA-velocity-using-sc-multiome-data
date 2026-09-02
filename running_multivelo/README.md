MultiVelo, a mechanistic model of gene expression that extends the popular RNA velocity framework by incorporating epigenomic data.
Citation: Li, C., Virgilio, M.C., Collins, K.L. & Welch J.D. Multi-omic single-cell velocity models epigenome–transcriptome interactions and improves cell fate prediction. Nat Biotechnol 41, 387-398 (2023).

In this project Multivelo v0.1.3 was established and run to analyse RNA velocity using RNA expression, ATAC chromatin accessibility and LIAM's batch corrected cell clusters.

File description:

1. multivelo-0.1.3_environment.yaml - This YAML file defines the Conda environment required to reproducibly run Multivelo v0.17. Generated using Miniforge3 (Conda version 24.11.3) on a Linux kernel (version 5.14).
2. multivelo_working_env_bash_patching_code.sh - Upon the Conda environment creation, this bash script patches some of the python package code required for the environment to be excecuted smoothly for analysis.
  
   ```bash
   conda env create -f multivelo-0.1.3_environment.yaml
   chmod +x multivelo_working_env_bash_patching_code.sh
   ./multivelo_working_env_bash_patching_code.sh
   #Provide path to conda env folder and name of the environment 
   ```
4. MultiVelo_reorganising_filtering_normalisation_5K_19Jun2026.py - This python script organises, filters and normalises the provided RNA expressions as LOOM files, ATAC data and incorporating LIAM corrected cell cluster information. Once the RNA expressions and ATAC data are prepared the Multivelo model is run to predic RNA Velocities.
5. run_slurm_Multivelo_reorganising_filtering_normalisation_5K_19June2026.sh - SLURM-friendly wrapper for running MultiVelo and the steps defined in the python script: MultiVelo_reorganising_filtering_normalisation_5K_19Jun2026.py. Execute using sbatch. 

 NOTE: sbatch parameters need to updated as per your cluster requirements or other cluster management and job scheduling systems. The paths in the scripts are user defined that needs to be updated.
