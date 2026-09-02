#!/usr/bin/env python3
# coding: utf-8

# ==============================================================================
# PROJECT: MultiVelo v0.1.3 Workflow and Analysis
# ==============================================================================
# Title:       Predicting RNA-velocity using multiome data generated from iPSCs established for human muscle development.
# Description: SLURM-friendly wrapper for running MultiVelo v0.1.3 recover_dynamics_chrom 
#              on multiome with alignment of all the objects, i.e., RNA, ATAC and LIAM.
#
# Author(s):   Dr. Shuba Varshini Alampalli
# Developer:   Implemented in CUBI's SODAR-Kiosc by Federico Marotta
#SODAR Data Steward: Dr. Shuba Varshini Alampalli and Dr. Miha Milek
# Institute:   Core Unit Bioinformatics (CUBI), Berlin Institute of Health, Germany AND Charité Universitätsmedizin Berlin, Germany 
# Contact:     cubi-helpdesk@bih-charite.de
#
# Date:        2026-07-22
# Version:     1.0
#
# Environment: conda activate multivelo-0.1.3
# ==============================================================================

"""
SLURM-friendly wrapper for running MultiVelo recover_dynamics_chrom on multiome.
Updated with "Master Alignment" logic to ensure cell/gene synchronization and 
proper transfer of corrected neighbor graphs from Liam's embedding.

Usage example:
  python MultiVelo_fixed.py --rna_loom /path/merged.loom \
    --atac_dir /path/filtered_feature_matrix_manual_made \
    --liam_emb /path/adata_merged_clusters_ATAC.h5ad \
    --peak_annot /path/atac_peak_annotation.tsv \
    --linkage /path/feature_linkage.bedpe \
    --out multivelo_result.h5ad \
    --n_jobs 16
"""

import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import anndata as ad
import scipy.sparse as sp
import matplotlib.pyplot as plt
import matplotlib as mpl

# Settings
scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo')
np.set_printoptions(suppress=True)
pd.set_option('display.max_columns', 100)

# --- Helper Functions ---

def ensure_sparse(adata):
    """Ensures X and layers are sparse matrices to save memory."""
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    for k, mat in adata.layers.items():
        if not sp.issparse(mat):
            adata.layers[k] = sp.csr_matrix(mat)

def strip_prefix_keep_after_colon_str(s):
    """Cleans RNA cell names (e.g., 'files:AGCT-1' -> 'AGCT-1')"""
    s = str(s)
    return s.split(":")[-1] if ":" in s else s

def convert_minus1_before_suffix_to_x(name):
    """Cleans ATAC cell names (e.g., 'AGCT-1' -> 'AGCTx')"""
    s = str(name)
    idx = s.rfind("-1_")
    if idx != -1:
        return s[:idx] + "x_" + s[idx + 3 :]
    if s.endswith("-1"):
        return s[:-2] + "x"
    return s

def save_and_show_pdf(filename, fig=None):
    """Rasterizes complex vectors, saves as high-res PDF, and displays in Jupyter."""
    if fig is None:
        fig = plt.gcf()
        
    # Iterate through ALL axes (crucial for subplots) to rasterize dots/arrows
    # Text, legends, and axis lines remain fully editable vector objects!
    for ax in fig.axes:
        for c in ax.collections:
            c.set_rasterized(True)
        for p in ax.patches:
            p.set_rasterized(True)
            
    # Save with 300 DPI to ensure the rasterized parts are publication-quality
    plt.savefig(filename, bbox_inches='tight', dpi=300)
    #plt.show() # This renders in Jupyter AND clears the figure memory
    print(f"  -> Saved & Displayed: {filename}\n")

# --- Main Pipeline ---

def main(args):
    # 1. Set Threading Environment
    if args.n_jobs:
        os.environ['OMP_NUM_THREADS'] = str(args.n_jobs)
        os.environ['MKL_NUM_THREADS'] = str(args.n_jobs)

    # ------------------------------------------------------------------
    # DATA LOADING & INITIAL CLEANING
    # ------------------------------------------------------------------
    print(f"\n--- Loading Data ---")
    
    # LOAD RNA
    print(f"Loading RNA: {args.rna_loom}")
    adata_rna = scv.read(args.rna_loom, cache=True)
    adata_rna.obs_names = [strip_prefix_keep_after_colon_str(x) for x in adata_rna.obs_names]
    adata_rna.var_names_make_unique()
    #ensure_sparse(adata_rna)

    # LOAD ATAC
    print(f"Loading ATAC: {args.atac_dir}")
    adata_atac = sc.read_10x_mtx(args.atac_dir, var_names='gene_symbols', cache=True, gex_only=False)
    adata_atac = adata_atac[:, adata_atac.var['feature_types'] == "Peaks"].copy()
    
    # Clean ATAC names *before* aggregation to match RNA style if needed, 
    # but usually aggregation relies on the raw peak names. 
    # We clean obs_names (cells) here.
    adata_atac.obs_names = [convert_minus1_before_suffix_to_x(x) for x in adata_atac.obs_names]
    #ensure_sparse(adata_atac)

    # AGGREGATE PEAKS TO GENES (MultiVelo specific)
    print("Aggregating ATAC peaks to gene scores (this may take time)...")
    adata_atac = mv.aggregate_peaks_10x(
        adata_atac, 
        args.peak_annot, 
        args.linkage,
        verbose=True
    )
    # Result: adata_atac.var_names are now GENES, matching adata_rna.var_names

    # LOAD LIAM EMBEDDING
    print(f"Loading Liam Embedding: {args.liam_emb}")
    adata_liam_emb = ad.read_h5ad(args.liam_emb)
    adata_liam_emb.obs_names = [convert_minus1_before_suffix_to_x(x) for x in adata_liam_emb.obs_names]
    adata_liam_emb.obs_names_make_unique()


    # ------------------------------------------------------------------
    # STEP 1: THE GRAND ALIGNMENT (Intersection)
    # ------------------------------------------------------------------
    print(f"\n--- Step 1: Aligning Cells and Genes ---")
    
    # Intersect Cells (3-way)
    common_cells = sorted(list(
        set(adata_rna.obs_names) & 
        set(adata_atac.obs_names) & 
        set(adata_liam_emb.obs_names)
    ))

    # Intersect Genes (RNA vs ATAC)
    common_genes = sorted(list(
        set(adata_rna.var_names) & 
        set(adata_atac.var_names)
    ))

    print(f"  > Common Cells: {len(common_cells)}")
    print(f"  > Common Genes: {len(common_genes)}")

    if len(common_cells) == 0:
        raise ValueError("Error: No common cells found. Check your name cleaning functions.")

    # Apply Subset (Sorts them identically)
    adata_rna = adata_rna[common_cells, common_genes].copy()
    adata_atac = adata_atac[common_cells, common_genes].copy()
    adata_liam_emb = adata_liam_emb[common_cells].copy() # Liam usually has different vars

    # ------------------------------------------------------------------
    # STEP 2: TRANSFER METADATA & EMBEDDINGS
    # ------------------------------------------------------------------
    print(f"\n--- Step 2: Transferring Metadata & Embeddings ---")
    
    cols_to_transfer = ['batch', 'sample', 'day', 'leiden_merged']
    for col in cols_to_transfer:
        if col in adata_liam_emb.obs.columns:
            adata_rna.obs[col] = adata_liam_emb.obs[col]

    # Transfer UMAP / Embeddings
    if 'X_umap' in adata_liam_emb.obsm:
        adata_rna.obsm['X_umap'] = adata_liam_emb.obsm['X_umap'].copy()
    
    # Try 'embedding' (integrated) or fallback to PCA
    if 'embedding' in adata_liam_emb.obsm:
        adata_rna.obsm['X_liam_raw'] = adata_liam_emb.obsm['embedding'].copy()
        print("  > LIAM embedding stored as 'X_liam_raw'")
    elif 'X_pca' in adata_liam_emb.obsm:
        adata_rna.obsm['X_liam_raw'] = adata_liam_emb.obsm['X_pca'].copy()
        print("  > LIAM PCA stored as 'X_liam_raw'")
    else:
        raise ValueError("No LIAM embedding found in adata_liam_emb.obsm")

    print("  > Metadata and embeddings transferred.")


    # ------------------------------------------------------------------
    # STEP 3: FILTERING (Lock-Step)
    # ------------------------------------------------------------------
    print(f"\n--- Step 3: Filtering & Normalization ---")
    
    # Filter RNA Cells
    sc.pp.filter_cells(adata_rna, min_counts=500)
    
    # Sync Deletions to ATAC & Liam
    valid_cells = adata_rna.obs_names
    adata_atac = adata_atac[valid_cells].copy()
    adata_liam_emb = adata_liam_emb[valid_cells].copy()
    print(f"  > Cells remaining: {len(valid_cells)}")

    # Filter Genes
    sc.pp.filter_genes(adata_rna, min_cells=20)
    
    # Save raw counts for Seurat v3
    adata_rna.layers['counts'] = adata_rna.X.copy()
    
    # Normalize RNA
    scv.pp.filter_and_normalize(adata_rna, min_shared_counts=10, log=True)
    
    # HVG Selection
    print("  > Selecting Top 5000 HVGs (Seurat v3)...")
    sc.pp.highly_variable_genes(
        adata_rna, 
        n_top_genes=5000, 
        flavor="seurat_v3", 
        #layer='counts', 
        subset=True
        #batch_key='batch'
    )
    
    # Sync ATAC Genes to match HVGs
    final_genes = adata_rna.var_names
    adata_atac = adata_atac[:, final_genes].copy()
    
    # TF-IDF Normalize ATAC (Standard for MultiVelo)
    mv.tfidf_norm(adata_atac)

    # ------------------------------------------------------------------
    # STEP 5: MOMENTS & DYNAMICS
    # ------------------------------------------------------------------
    print(f"\n--- Step 5: Calculating Moments & Running Dynamics ---")
    
    # 1. PCA (needed for moments calculation dimensions)
    sc.pp.pca(adata_rna, n_comps=30)
    # 2. calculated neighbors on raw RNA → This finds cells that are similar in the current raw data space (mostly within the same batch).
    sc.pp.neighbors(adata_rna, n_neighbors=30, n_pcs=30)
    # 3. Moments
    # calculated moments → This smooths the data using those "safe" local neighbors, preserving the clean splicing kinetics within each batch.
    scv.pp.moments(adata_rna, n_pcs=args.n_pcs, n_neighbors=30)
    
    # Smooth ATAC using the same graph
    mv.knn_smooth_chrom(adata_atac, conn=adata_rna.obsp['connectivities'], n_neighbors=30)

    print("\n Ready for MultiVelo velocity calculation!")
    print(f"   - RNA shape: {adata_rna.shape}")
    print(f"   - ATAC shape: {adata_atac.shape}")
    print(f"   - Using LIAM batch-corrected neighborhoods")

    # ==============================================================================
    # STEP 7: scVelo of the RNA seq
    # ==============================================================================
    print("\n--- Step 7: scVelo on the RNA and Recover Dynamics ---")
    adata_rna_scv = adata_rna.copy()
    scv.tl.recover_dynamics(adata_rna_scv,n_jobs=32)

    adata_rna_scv.write("/data/cephfs-1/home/users/shal11_c/work/Documents/multivelo_velocyto_Dec27/multivelo_results_Jan12_2026/scvelo_result_hvg_5K_recoverdynamics_19June2026.h5ad")


    # Velocity
    print("  > Running Velocity (Stochastic)...")
    out_dir = "/data/cephfs-1/home/users/shal11_c/work/Documents/multivelo_velocyto_Dec27/final_version2_multivelo_analysis_19June2026/scVelo_fig"
    scv.tl.velocity(adata_rna_scv, mode="stochastic", min_r2=0.02)
    scv.tl.velocity_graph(adata_rna_scv, n_jobs=32)

    # Latent time
    scv.tl.latent_time(adata_rna_scv)

    adata_rna_scv.write("/data/cephfs-1/home/users/shal11_c/work/Documents/multivelo_velocyto_Dec27/final_version2_multivelo_analysis_19June2026/scVelo_fig/scvelo_result_hvg_5K_stochastic_latentTime_19Jun2026.h5ad")

    scv.pl.velocity_embedding_stream(adata_rna_scv, basis='umap', legend_loc="right", color='leiden_merged')
    save_and_show_pdf(f"{out_dir}/Fig_scv_5K_Differentiation_Trajectories_rna_5K.pdf")

    scv.pl.velocity_embedding_grid(adata_rna_scv, basis='umap', arrow_length=3, arrow_size=2,color='leiden_merged')
    save_and_show_pdf(f"{out_dir}/Fig_scv_5K_Differentiation_Trajectories_grid_rna_5K.pdf")
    
    #scv.pl.scatter(adata_rna_scv, color="latent_time", color_map="gnuplot")
    #save_and_show_pdf(f"{out_dir}/Fig_scv_5K_latent_time_rna.pdf")

    print("\n--- Step 8: scVelo on the RNA for MultiVelo ---")
    scv.tl.recover_dynamics(adata_rna,n_jobs=32)
    scv.tl.velocity(adata_rna, mode="dynamical", min_r2=0.02)
    scv.tl.velocity_graph(adata_rna, n_jobs=32)

    adata_rna.write("/data/cephfs-1/home/users/shal11_c/work/Documents/multivelo_velocyto_Dec27/final_version2_multivelo_analysis_19June2026/scVelo_fig/RNA_scVelo_forMultivelo_5K_19June2026.h5ad")

    print("scVelo Done.")
    
    # --- Fix for Missing Colors ---
    # Check if Liam has the colors defined
    if 'leiden_merged_colors' in adata_liam_emb.uns:
        print("Transferring colors from Liam's embedding...")
        adata_rna.uns['leiden_merged_colors'] = adata_liam_emb.uns['leiden_merged_colors'].copy()
    else:
        # Fallback: If Liam doesn't have them either, generate new ones
        print("Liam's object is missing colors too. Generating new palette...")
        # Force categorical dtype just in case
        if 'leiden_merged' in adata_rna.obs:
            adata_rna.obs['leiden_merged'] = adata_rna.obs['leiden_merged'].astype('category')
            # Use default scanpy palette
            sc.pl.palettes.default_20
            # Determine number of categories
            n_cats = len(adata_rna.obs['leiden_merged'].cat.categories)
            # Assign colors
            adata_rna.uns['leiden_merged_colors'] = sc.pl.palettes.default_20[:n_cats]

    # ==============================================================================
    # STEP 9: Running multi-omic dynamical model
    # ==============================================================================
    print("\n--- Step 9: multi-omic dynamical model MultiVelo ---")
    #MultiVelo incorporates chromatin accessibility information into RNA velocity and achieves better lineage predictions.
    adata_result = mv.recover_dynamics_chrom(
        adata_rna,
        adata_atac,
        max_iter=5,
        init_mode='invert',
        parallel=True,
        save_plot=False,
        rna_only=False,
        fit=True,
        n_anchors=500,
        n_pcs=30,
        extra_color_key='leiden_merged',
        n_jobs=32
        )

    print("Writing result to:", args.out)
    adata_result.write(args.out)
    print("MultiVelo Done.")


if __name__ == '__main__':
    p = argparse.ArgumentParser(description="MultiVelo Wrapper with Master Alignment")
    
    # Input Files
    p.add_argument('--rna_loom', required=True, help="Path to RNA Loom file")
    p.add_argument('--atac_dir', required=True, help="Path to ATAC 10x directory")
    p.add_argument('--liam_emb', required=True, help="Path to .h5ad containing Liam's corrected embedding/graph")
    
    # Paths for Aggregate Peaks (Made arguments for flexibility, defaults provided)
    p.add_argument('--peak_annot', default='/data/cephfs-1/work/projects/schuelke-cubi-muscle-dev/BtE_P07_P08_analyses/MULTIOME/outputs/reanalyzed_samples_combined_jamm_peaks/outs/atac_peak_annotation.tsv',
                   help="Path to atac_peak_annotation.tsv")
    p.add_argument('--linkage', default='/data/cephfs-1/work/projects/schuelke-cubi-muscle-dev/BtE_P07_P08_analyses/MULTIOME/outputs/reanalyzed_samples_combined_jamm_peaks/outs/analysis/feature_linkage/feature_linkage.bedpe',
                   help="Path to feature_linkage.bedpe")
    
    # Output
    p.add_argument('--out', required=True, help="Output .h5ad filename")

    # Parameters
    p.add_argument('--n_pcs', type=int, default=30)
    p.add_argument('--n_neighbors', type=int, default=50) # Used if Liam graph is missing
    p.add_argument('--n_anchors', type=int, default=300)
    p.add_argument('--max_iter', type=int, default=5)
    p.add_argument('--init_mode', type=str, default='invert')
    p.add_argument('--n_jobs', type=int, default=16)
    
    args = p.parse_args()
    main(args)
