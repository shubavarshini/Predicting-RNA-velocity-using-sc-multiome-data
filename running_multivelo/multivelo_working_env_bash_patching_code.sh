#!/bin/bash

echo "================================================="
echo "   MultiVelo Environment Patcher   "
echo "================================================="

# 1. Prompt the user for their base Conda environments path
echo "We need the path to where your Conda environments are stored."
read -p "Enter your envs directory (e.g., ~/work/miniforge3/envs): " ENVS_PATH

# Expand the tilde (~) to the actual home directory path if used
ENVS_PATH="${ENVS_PATH/#\~/$HOME}"

# 2. Prompt the user for the environment name
echo ""
echo "We also need the exact name of the Conda environment."
read -p "Enter the environment name (e.g., multivelo_working_env_patched): " ENV_NAME

# 3. Construct the full paths
SITE_PACKAGES_DIR="${ENVS_PATH}/${ENV_NAME}/lib/python3.12/site-packages/"
PYTHON_EXEC="${ENVS_PATH}/${ENV_NAME}/bin/python"

# 4. Verify the path actually exists before trying to patch
echo ""
if [ ! -d "$SITE_PACKAGES_DIR" ]; then
    echo "Error: Could not find the directory at:"
    echo "$SITE_PACKAGES_DIR"
    echo "Please check your path and environment name, then try again."
    exit 1
fi

if [ ! -f "$PYTHON_EXEC" ]; then
    echo "Error: Could not find Python executable at:"
    echo "$PYTHON_EXEC"
    echo "Is the environment fully built?"
    exit 1
fi

echo "Found site-packages! Applying manual patches to ${ENV_NAME}..."

# 5. Export the path so the Python snippet can read it safely
export PATCH_BASE_DIR="$SITE_PACKAGES_DIR"

# 6. Run the Python patcher using the specific environment's Python
"$PYTHON_EXEC" -c '
import os

# Grab the dynamic path passed from Bash
base_dir = os.environ.get("PATCH_BASE_DIR")

def apply_patch(file_path, old_code, new_code):
    full_path = os.path.join(base_dir, file_path)
    if os.path.exists(full_path):
        with open(full_path, "r") as f:
            content = f.read()
        
        # Only write if the old code is actually in the file
        if old_code in content:
            with open(full_path, "w") as f:
                f.write(content.replace(old_code, new_code))
            print(f"  [PATCHED] {file_path}")
        else:
            print(f"  [SKIPPED] {file_path} (Code not found or already patched)")
    else:
        print(f"  [ERROR] File not found: {file_path}")

# 1. Patch settings.py
apply_patch("scvelo/settings.py", 
    "warnings.filterwarnings(\"ignore\", category=cbook.mplDeprecation)", 
    "import matplotlib\nwarnings.filterwarnings(\"ignore\", category=matplotlib.MatplotlibDeprecationWarning)")

# 2. Patch neighbors.py
apply_patch("scvelo/preprocessing/neighbors.py", 
    "neighbors.compute_neighbors(write_knn_indices=True, **kwargs)", 
    "neighbors.compute_neighbors(**kwargs)")

# 3. Patch utils.py (Categories)
apply_patch("scvelo/plotting/utils.py", 
    "obs_vals.cat.categories = obs_vals.cat.categories.astype(str)", 
    "new_categories = obs_vals.cat.categories.astype(str)\n    adata.obs[value_to_plot] = pd.Categorical(obs_vals, categories=new_categories)")

# 4. Patch auxiliary.py (Missing import)
apply_patch("multivelo/auxiliary.py", 
    "import sys", 
    "import sys\nfrom .dynamical_chrom_func import top_n_sparse")

# 5. Patch _parallelize.py (NumPy arrays)
apply_patch("scvelo/core/_parallelize.py", 
    "res = np.array(res) if as_array else res", 
    "res = np.array(res, dtype=object) if as_array else res")
'

echo "All done! Your environment is ready to use."
