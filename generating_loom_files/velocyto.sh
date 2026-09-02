#!/bin/bash

# ==============================================================================
# Author      : Dr. Shuba Varshini Alampalli
# Date        : Feb 2026
# Description : SLURM submission script for running the Velocyto to generate LOOM files from BAM.
# ==============================================================================

#SBATCH --export=ALL
#SBATCH --job-name=velocyto
#SBATCH --mem=200G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
# #SBATCH -N 1
#SBATCH --output="%x_%j.log"
#SBATCH --time="5-00:00:00"
# #BATCH --partition=medium	# change this to your cluster's large-memory partition

source ~/work/miniforge3/etc/profile.d/conda.sh
conda activate velocyto-0.17

GTF='/data/cephfs-1/work/projects/schuelke-cubi-muscle-dev/BtE_P07_P08_analyses/shared_files/ref_genome/hg38/Homo_sapiens.GRCh38.111.gtf'

BAMS=(
  /data/cephfs-1/work/projects/schuelke-cubi-muscle-dev/BtE_P07_P08_analyses/MULTIOME/outputs/D12_REP1_run1/outs/gex_possorted_bam.bam
  /data/cephfs-1/work/projects/schuelke-cubi-muscle-dev/BtE_P07_P08_analyses/MULTIOME/outputs/D12_REP2_run1/outs/gex_possorted_bam.bam
  /data/cephfs-1/work/projects/schuelke-cubi-muscle-dev/BtE_P07_P08_analyses/MULTIOME/outputs/D20_REP1_run1/outs/gex_possorted_bam.bam
  /data/cephfs-1/work/projects/schuelke-cubi-muscle-dev/BtE_P07_P08_analyses/MULTIOME/outputs/D20_REP2_run1/outs/gex_possorted_bam.bam
  /data/cephfs-1/work/projects/schuelke-cubi-muscle-dev/BtE_P07_P08_analyses/MULTIOME/outputs/D22-15_REP1_run1/outs/gex_possorted_bam.bam
  /data/cephfs-1/work/projects/schuelke-cubi-muscle-dev/BtE_P07_P08_analyses/MULTIOME/outputs/D22-15_REP2_run1/outs/gex_possorted_bam.bam
  /data/cephfs-1/work/projects/schuelke-cubi-muscle-dev/BtE_P07_P08_analyses/MULTIOME/outputs/D8_REP1_run1/outs/gex_possorted_bam.bam
  /data/cephfs-1/work/projects/schuelke-cubi-muscle-dev/BtE_P07_P08_analyses/MULTIOME/outputs/D8_REP2_run1/outs/gex_possorted_bam.bam
)

for BAM in "${BAMS[@]}"; do
  # sample name, e.g. D22-15_REP1_run1
  SAMPLE=$(basename "$(dirname "$(dirname "$BAM")")")

  # barcode path derived from BAM location
  BARCODES="$(dirname "$BAM")/filtered_feature_bc_matrix/barcodes.tsv.gz"

  OUTDIR="~/work/Documents/multivelo_velocyto_Dec27/${SAMPLE}_velocyto"
  
  # create output directory if it does not exist
  mkdir -p "${OUTDIR}"

  echo "Running velocyto for ${SAMPLE}"
  echo "  BAM: ${BAM}"
  echo "  Barcodes: ${BARCODES}"
  
  rsync -av "${BAM}" "${OUTDIR}/"
  BAM_copy="~/work/Documents/multivelo_velocyto_Dec27/${SAMPLE}_velocyto/gex_possorted_bam.bam"
  
  samtools sort -t CB  -m 10G -@ 20 -O BAM -o ${OUTDIR}/cellsorted_gex_possorted_bam.bam ${BAM}

  samtools index -@ 20 $BAM_copy

  velocyto run \
    --samtools-memory 160000 \
    -@ 20 \
    -b "${BARCODES}" \
    -o "${OUTDIR}" \
    "${BAM_copy}" \
    "${GTF}"

done


