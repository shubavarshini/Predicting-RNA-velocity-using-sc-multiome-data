#!/bin/bash
#SBATCH --job-name=multivelo_19June2026
#SBATCH --output=multivelo_19June2026.%j.out
# #SBATCH --partition=highmem         # change this to your cluster's large-memory partition
#SBATCH --cpus-per-task=32         # adjust
#SBATCH --mem=500G                 # adjust to 200-500G depending on dataset
#SBATCH --time=48:00:00
# ## Optional GPU if you plan to use GPU-accelerated steps (multivelo itself is CPU-bound):
# ##SBATCH --gres=gpu:1
# #SBATCH --mail-type=END,FAIL
# #SBATCH --mail-user=you@example.com

# Load modules (adjust to your cluster)
source ~/work/miniforge3/etc/profile.d/conda.sh
# activate your conda env
conda activate multivelo-0.1.3

# Set threads to match requested CPUs
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Path adjustments
PY_SCRIPT=/data/cephfs-1/home/users/shal11_c/work/Documents/multivelo_velocyto_Dec27/MultiVelo_shuba_slurm_ready_reorganising_filtering_normalisation_5K_19Jun2026.py
RNA_LOOM=/data/cephfs-1/home/users/shal11_c/work/Documents/multivelo_velocyto_Dec27/merged.loom
ATAC_DIR=/data/cephfs-1/home/users/shal11_c/work/Documents/multivelo_velocyto_Dec27/filtered_feature_matrix_manual_made
PEAK_ANNO=/data/cephfs-1/work/projects/schuelke-cubi-muscle-dev/BtE_P07_P08_analyses/MULTIOME/outputs/reanalyzed_samples_combined_jamm_peaks/outs/atac_peak_annotation.tsv 
LINKAGE=/data/cephfs-1/work/projects/schuelke-cubi-muscle-dev/BtE_P07_P08_analyses/MULTIOME/outputs/reanalyzed_samples_combined_jamm_peaks/outs/analysis/feature_linkage/feature_linkage.bedpe
LIAM_EMB=/data/cephfs-1/work/projects/schuelke-cubi-muscle-dev/BtE_P07_P08_analyses/MULTIOME/outputs/adata_objects/adata_merged_clusters_ATAC_on_jamm_peaks.h5ad
OUT=/data/cephfs-1/home/users/shal11_c/work/Documents/multivelo_velocyto_Dec27/final_version2_multivelo_analysis_19June2026/multivelo_result_hvg_5K_reorganising_filtering_normalisation_19June2026.h5ad

# Print Header to SLURM Log
cat << "EOF"
==============================================================================
PROJECT: MultiVelo v0.1.3 Workflow and Analysis
==============================================================================
Title:       Predicting RNA-velocity using multiome data generated from iPSCs established for human muscle development.
Description: SLURM-friendly wrapper for running MultiVelo v0.1.3 recover_dynamics_chrom 
             on multiome with alignment of all the objects, i.e., RNA, ATAC and LIAM.

Author(s):   Dr. Shuba Varshini Alampalli
Developer:   Implemented in CUBI's SODAR-Kiosc by Federico Marotta
SODAR Data Steward: Dr. Shuba Varshini Alampalli and Dr. Miha Milek
Institute:   Core Unit Bioinformatics (CUBI), Berlin Institute of Health, Germany AND Charité Universitätsmedizin Berlin, Germany 
Contact:     cubi-helpdesk@bih-charite.de

Date:        2026-07-22
Version:     1.0
Environment: conda activate multivelo-0.1.3
==============================================================================
EOF

srun python $PY_SCRIPT \
  --rna_loom $RNA_LOOM \
  --atac_dir $ATAC_DIR \
  --liam_emb $LIAM_EMB \
  --peak_annot $PEAK_ANNO \
  --linkage $LINKAGE \
  --out $OUT \
  --n_jobs $SLURM_CPUS_PER_TASK
