#!/usr/bin/env python3
# ==============================================================================
# Author      : Dr. Shuba Varshini Alampalli
# Date        : Feb 2026
# Description : Custom python script to merge the loom files generated for multiple scRNA and chnging barcodes to include replicate metadata
# ==============================================================================
"""
merge_looms_with_union_and_combine.py

Robust merging of multiple .loom files into one using loompy, compatible with loompy v3.x APIs observed in the user's environment.

Approach:
 - Compute union of genes and union of layers across inputs. 
 - For each input, create temporary on-disk numpy.memmap files (one per layer) sized (n_genes_union, ncols_for_this_input). Fill those memmaps by streaming input data in column chunks and mapping input rows -> union rows.
 - Call loompy.create(tmp_loom, layers=memmap_dict, row_attrs=..., col_attrs=...) to create a proper .loom temp file that has the '' (empty-string) main layer and other named layers. This avoids API mismatches with loompy.new and allows writing the main matrix cleanly.
 - After creating all temp looms, call loompy.combine(temp_paths, output_path).
 - Clean up memmap backing files and temp looms (unless --keep-temp).

Notes and trade-offs:
 - This approach uses on-disk memmaps for each layer per-temp loom. That requires temporary disk space roughly equal to the sum of the temp loom sizes (which can be large). It avoids large RAM usage.
 - It preserves layer dtypes where possible, writes missing layers as zeros, and writes column attributes (ca) and row attribute with gene names (ra).
 - Barcodes written are ORIGINAL_BARCODE + SEP + SAMPLE_ID (per my request).

Usage:
    python merge_looms_with_union_and_combine.py -o merged.loom file1.loom file2.loom
    python merge_looms_with_union_and_combine.py --map mapping.csv -o merged.loom #Prefered usage

Requires:
    loompy (v3.x), numpy

"""
import argparse
import csv
import os
import re
import sys
import tempfile
import shutil
from collections import OrderedDict
from typing import List, Dict, Tuple

import numpy as np
import loompy
import scipy.sparse as sp

PREFERRED_GENE_KEYS = ["Gene", "gene", "GeneName", "gene_name", "Name", "name"]
PREFERRED_BARCODE_KEYS = [
    "CellID", "CellIDs", "cell_id", "barcode", "Barcode",
    "Cell", "CellName", "CellNames", "barcodes",
]


def sanitize_sample_id(s: str) -> str:
    return re.sub(r"[^0-9A-Za-z\-_]", "_", s)


def read_mapping_csv(path: str) -> Dict[str, str]:
    mapping = {}
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        for r in reader:
            if len(r) < 2:
                continue
            infile = r[0].strip()
            sid = r[1].strip()
            if infile:
                mapping[infile] = sanitize_sample_id(sid or os.path.splitext(os.path.basename(infile))[0])
    return mapping


def find_gene_key(loom: loompy.LoomConnection) -> str:
    for k in PREFERRED_GENE_KEYS:
        if k in loom.ra:
            return k
    if len(loom.ra.keys()) > 0:
        return list(loom.ra.keys())[0]
    return None


def find_barcode_key(loom: loompy.LoomConnection) -> str:
    for k in PREFERRED_BARCODE_KEYS:
        if k in loom.ca:
            return k
    return None


def gather_union_genes(input_paths: List[str], user_gene_key: str = None) -> Tuple[List[str], Dict[str, str]]:
    union = OrderedDict()
    gene_key_used = {}
    for p in input_paths:
        with loompy.connect(p, "r") as ds:
            gk = user_gene_key or find_gene_key(ds)
            if gk is None:
                raise RuntimeError(f"Cannot determine gene row attribute for file {p}")
            gene_key_used[p] = gk
            genes = np.asarray(ds.ra[gk], dtype=str)
            for g in genes:
                if g not in union:
                    union[g] = True
    return list(union.keys()), gene_key_used


def build_col_attr_placeholders(all_col_keys: List[str], ncols: int) -> Dict[str, np.ndarray]:
    placeholders = {}
    for k in all_col_keys:
        arr = np.empty(ncols, dtype=object)
        arr[:] = ""
        placeholders[k] = arr
    return placeholders


def make_memmap(path: str, shape: Tuple[int, int], dtype):
    # Create memmap file on disk with zeros
    mm = np.memmap(path, dtype=dtype, mode="w+", shape=shape)
    # Ensure zero-init
    mm[:] = 0
    # Flush to disk
    mm.flush()
    return mm


def close_and_delete_memmap(mm: np.memmap, path: str):
    try:
        mm.flush()
        del mm
    except Exception:
        pass
    try:
        os.remove(path)
    except Exception:
        pass


def main():
    parser = argparse.ArgumentParser(description="Merge loom files by union of genes; append sampleID to barcodes (ORIGINAL + SEP + SAMPLE_ID) using memmap-backed temp creation")
    parser.add_argument("inputs", nargs="*", help="Input loom files (ignored if --map is provided)")
    parser.add_argument("-o", "--output", required=True, help="Output loom file path")
    parser.add_argument("--map", help="CSV mapping file: input_file,sample_id (use this instead of listing inputs).")
    parser.add_argument("--sep", default="_", help="Separator appended between original barcode and sample ID (default: '_')")
    parser.add_argument("--barcode-key", default="", help="Force a specific column-attribute key to be used as barcode (auto-detect otherwise)")
    parser.add_argument("--gene-key", default="", help="Force a specific row-attribute key to use as gene names (auto-detect otherwise)")
    parser.add_argument("--keep-temp", action="store_true", help="Keep temporary files (do not delete)")
    parser.add_argument("--tmpdir", default=None, help="Directory to create temporary files in (default: system temp dir)")
    parser.add_argument("--chunk", type=int, default=2000, help="Column chunk size when streaming (default: 2000)")
    args = parser.parse_args()

    # Build mapping
    if args.map:
        mapping = read_mapping_csv(args.map)
        if not mapping:
            raise RuntimeError(f"Mapping CSV {args.map} contains no valid rows.")
        if args.inputs:
            inputs = args.inputs
            sample_ids = []
            for p in inputs:
                if p not in mapping:
                    sample_ids.append(sanitize_sample_id(os.path.splitext(os.path.basename(p))[0]))
                else:
                    sample_ids.append(mapping[p])
        else:
            inputs = list(mapping.keys())
            sample_ids = [mapping[k] for k in inputs]
    else:
        if not args.inputs:
            parser.error("Either supply input looms on the command line or provide --map mapping.csv")
        inputs = args.inputs
        sample_ids = [sanitize_sample_id(os.path.splitext(os.path.basename(p))[0]) for p in inputs]

    for p in inputs:
        if not os.path.exists(p):
            raise FileNotFoundError(f"Input file not found: {p}")

    print(f"Found {len(inputs)} input files.")

    # union genes
    print("Scanning files to build union of genes...")
    union_genes, gene_key_used = gather_union_genes(inputs, user_gene_key=(args.gene_key or None))
    n_genes_union = len(union_genes)
    union_genes_arr = np.asarray(union_genes, dtype=str)
    print(f"Union genes count: {n_genes_union}")

    # union layers & col keys
    layers_union = set()
    col_keys_union = set()
    ncols_list = []
    per_file_ncols = {}
    for p in inputs:
        with loompy.connect(p, "r") as ds:
            layers_union.update(list(ds.layers.keys()))
            col_keys_union.update(list(ds.ca.keys()))
            ncols_list.append(ds.shape[1])
            per_file_ncols[p] = ds.shape[1]
    layers_union = list(layers_union)
    col_keys_union = list(col_keys_union)
    total_cols = sum(ncols_list)
    print(f"Union layers: {layers_union}")
    print(f"Union column-attribute keys: {col_keys_union}")
    print(f"Total columns in final file will be: {total_cols}")

    tmpdir = args.tmpdir or tempfile.mkdtemp(prefix="loom_merge_")
    created_tmpdir_by_us = args.tmpdir is None
    print(f"Temporary files and memmaps will be created in: {tmpdir}")

    temp_paths = []
    union_index = {g: i for i, g in enumerate(union_genes)}

    CHUNK = max(1, int(args.chunk))
    barcode_key_forcing = args.barcode_key if args.barcode_key else None
    row_key_name = args.gene_key or "Gene"

    # For each input, create memmaps for each layer then create temp loom via loompy.create
    for idx, p in enumerate(inputs):
        sample = sample_ids[idx]
        basename = os.path.splitext(os.path.basename(p))[0]
        tmpname = os.path.join(tmpdir, f"tmp_{idx}_{basename}.loom")
        temp_paths.append(tmpname)
        print(f"Preparing temp memmaps for '{p}' -> will create '{tmpname}' with sampleID='{sample}'")

        with loompy.connect(p, "r") as ds_in:
            ncols = ds_in.shape[1]

            # prepare column attributes for this temp file
            col_attrs = build_col_attr_placeholders(col_keys_union, ncols)

            if barcode_key_forcing:
                barcode_key = barcode_key_forcing
            else:
                candidate = find_barcode_key(ds_in)
                barcode_key = candidate if candidate else "CellID"

            if barcode_key in ds_in.ca:
                original_barcodes = np.asarray(ds_in.ca[barcode_key], dtype=str)
            else:
                original_barcodes = np.array([f"cell{i}" for i in range(ncols)], dtype=str)

            # new barcode format: ORIGINAL + SEP + SAMPLE_ID
            new_barcodes = np.array([f"{b}{args.sep}{sample}" for b in original_barcodes], dtype=str)
            col_attrs[barcode_key] = new_barcodes

            # copy other col attrs
            for k, v in ds_in.ca.items():
                if k == barcode_key:
                    continue
                col_attrs[k] = np.asarray(v, dtype=object)

            # Create memmap files for each layer. Keep track for cleanup.
            memmap_paths = {}
            memmaps = {}
            try:
                for layer in layers_union:
                    # dtype: try to use input layer dtype if present, otherwise float32
                    if layer in ds_in.layers:
                        try:
                            dtype = ds_in.layers[layer].dtype
                        except Exception:
                            dtype = np.float32
                    else:
                        dtype = np.float32

                    # create a memmap file for this layer
                    safe_layer_name = "main" if layer == "" else layer.replace("/", "_")
                    mm_path = os.path.join(tmpdir, f"mm_{idx}_{safe_layer_name}.dat")
                    mm = make_memmap(mm_path, (n_genes_union, ncols), dtype=dtype)
                    memmap_paths[layer] = mm_path
                    memmaps[layer] = mm

                # Fill memmaps by streaming input data in column chunks
                for layer in layers_union:
                    mm = memmaps[layer]
                    is_present = (layer in ds_in.layers)
                    if is_present:
                        in_genes = np.asarray(ds_in.ra[gene_key_used[p]], dtype=str)
                        idx_map = np.array([union_index[g] for g in in_genes], dtype=np.int64)
                        for c0 in range(0, ncols, CHUNK):
                            c1 = min(ncols, c0 + CHUNK)
                            if layer == "" or layer == getattr(ds_in, "main_layer", None):
                                in_block = ds_in[:, c0:c1]
                            else:
                                in_block = ds_in.layers[layer][:, c0:c1]

                            # convert sparse to dense if needed
                            if sp.issparse(in_block):
                                in_block = in_block.toarray()
                            else:
                                in_block = np.asarray(in_block)

                            # place rows into union positions directly into memmap
                            mm[idx_map, c0:c1] = in_block
                    else:
                        # leave zeros (already zero-initialized)
                        continue

                    # flush periodically
                    mm.flush()

                # Build layers dict expected by loompy.create: keys should be layer names ('' allowed)
                layers_for_create = {}
                for layer, mm_path in memmap_paths.items():
                    # open memmap read-only view for create (numpy.memmap works as array)
                    arr = np.memmap(mm_path, mode="r", dtype=np.dtype(memmaps[layer].dtype), shape=(n_genes_union, ncols))
                    layers_for_create[layer] = arr

                # Now call loompy.create with layers dict, row_attrs and col_attrs
                row_attrs = {row_key_name: union_genes_arr}
                # loompy.create expects row_attrs and col_attrs as dict of numpy arrays
                col_attrs_for_create = {k: np.asarray(v) for k, v in col_attrs.items()}
                print(f"Calling loompy.create for temp loom: {tmpname} (this writes the memmaps into a loom file)")
                loompy.create(tmpname, layers_for_create, row_attrs, col_attrs_for_create, file_attrs={"shape": (n_genes_union, ncols)})
                print(f"Created temp loom: {tmpname}")
            finally:
                # close and remove memmaps
                for layer, mm in list(memmaps.items()):
                    path = memmap_paths.get(layer)
                    try:
                        close_and_delete_memmap(mm, path)
                    except Exception:
                        pass

    # Combine temp looms into final merged loom
    print("Combining temporary looms into final output using loompy.combine...")
    try:
        loompy.combine(temp_paths, args.output)
    except Exception as e:
        # fallback manual concatenation
        print(f"loompy.combine failed with: {e}. Falling back to manual concatenation.", file=sys.stderr)
        col_start = 0
        final_file_attrs = {"shape": (n_genes_union, total_cols)}
        with loompy.new(args.output, file_attrs=final_file_attrs) as ds_out:
            ds_out.ra[row_key_name] = union_genes_arr
            placeholders = build_col_attr_placeholders(col_keys_union, total_cols)
            for k_c, v_c in placeholders.items():
                ds_out.ca[k_c] = v_c

            for tp in temp_paths:
                with loompy.connect(tp, "r") as ds_temp:
                    ncols = ds_temp.shape[1]
                    col_end = col_start + ncols
                    ds_out[:, col_start:col_end] = ds_temp[:, :]
                    for layer in layers_union:
                        is_main_layer = (layer == "" or layer == ds_temp.main_layer)
                        if is_main_layer:
                            continue
                        ds_out.layers[layer][:, col_start:col_end] = ds_temp.layers[layer][:, :]
                    for k in col_keys_union:
                        if k in ds_temp.ca:
                            ds_out.ca[k][col_start:col_end] = ds_temp.ca[k]
                    col_start = col_end
        print("Fallback manual concatenation finished.")

    print(f"Output written to: {args.output}")

    # cleanup temp looms
    if not args.keep_temp:
        print("Removing temporary files...")
        for tp in temp_paths:
            try:
                os.remove(tp)
            except Exception:
                pass
        if created_tmpdir_by_us:
            try:
                shutil.rmtree(tmpdir, ignore_errors=True)
            except Exception:
                pass
        print("Temporary files removed.")
    else:
        print(f"Kept temporary files in: {tmpdir}")

    print("Done.")


if __name__ == "__main__":
    main()
