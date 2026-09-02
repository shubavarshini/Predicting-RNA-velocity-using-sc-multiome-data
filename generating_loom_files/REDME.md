Velocyto estimates RNA velocities of single cells by distinguishing unspliced and spliced mRNAs in standard single-cell RNA sequencing protocols. 

In this project Velocyto python package is used to run the velocity analysis outputting a LOOM file from scRNA-seq aligned to reference BAM files and GTF annotation file.

File description:
1. velocyto-0.17_environment.yaml - This YAML file defines the Conda environment required to reproducibly run Velocyto v0.17. Generated using Miniforge3 (Conda version 24.11.3) on a Linux kernel (version 5.14).
2. velocyto.sh - SLURM submission script for running the Velocyto package for individual BAM files. 
3. merge_looms_with_union_and_combine.py - For using the individual LOOM files downstream, they need to be merged. This is a customized python script to merge multiple .loom files into one using loompy, compatible with loompy v3.x APIs observed in the user's environment.

NOTE: sbatch parameters need to updated as per your cluster requirements or other cluster management and job scheduling systems. The paths in the scripts are user defined that needs to be updated.
